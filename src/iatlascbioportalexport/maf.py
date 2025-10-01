import argparse
import logging
from multiprocessing import Pool
import os
import subprocess
from typing import Dict

import pandas as pd

import utils

syn = utils.synapse_login()

def read_and_merge_maf_files(input_folder_synid: str) -> pd.DataFrame:
    """Read in and merge MAF files from a specified folder

    Args:
        folder: Synapse id of folder containing MAF files

    Return:
        pd.DataFrame: Merged maf of all mafs in input folder
    """
    entities = syn.getChildren(input_folder_synid)
    # Filter for files ending in .maf
    dfs = []
    for item in entities:
        if item["name"].endswith(".maf"):
            df = pd.read_csv(syn.get(item["id"]).path, sep="\t", comment="#")
            dfs.append(df)

    if not dfs:
        print(f"No MAF files found in  {input_folder_synid:}")
        return None
    else:
        merged_df = pd.concat(dfs, ignore_index=True)
        return merged_df


def split_into_chunks(
    dataset_name: str,
    input_df: pd.DataFrame,
    datahub_tools_path: str,
    max_rows: int = 40000,
) -> int:
    """Splits the maf file into designated number of equal (if possible) chunks.
       Maf files don't need to be randomnly shuffled prior because the
       variants are independent of one another.

       Genome nexus cannot handle maf file sizes greater than ~55 MB

    Args:
        dataset_name (str): Name of the dataset
        input_df (pd.DataFrame): Input maf data to be split
        datahub_tools_path (str): Path to the datahub tools repo
        max_rows (int): maximum rows to split each chunk into

    Returns:
        int: Number of maf chunks
    """
    dataset_dir = os.path.join(
        f"{datahub_tools_path}/add-clinical-header/", dataset_name
    )
    n_maf_chunks = (len(input_df) + max_rows - 1) // max_rows  # ceiling division
    for i in range(n_maf_chunks):
        maf_chunk = input_df[i * max_rows : (i + 1) * max_rows]
        maf_chunk.to_csv(
            f"{dataset_dir}/data_mutations_{i + 1}.txt",
            sep="\t",
            index=False,
            float_format="%.12g",
        )
    return n_maf_chunks


def run_genome_nexus(
    dataset_name: str, n_maf_chunks: int, n_workers: int, datahub_tools_path: str
) -> None:
    """Runs genome nexus annotator on each of the split maf chunks. Logging
        is saved for each chunk.
        
        This will parallelize the workflow if n_workers is defined or > 1,
        otherwise it will run genome nexus on each of the chunk(s) serially

    Args:
        dataset_name (str): Name of the dataset
        n_maf_chunks (int): Number of chunks the maf was split into
        n_workers (int): Number of worker processes to use to
            run the genome nexus annotator in parallel.
        datahub_tools_path (str): Path to the datahub tools repo
    """
    dataset_dir = os.path.join(
        f"{datahub_tools_path}/add-clinical-header/", dataset_name
    )

    if n_workers and n_workers > 1:
        args_list = [(i, dataset_dir) for i in range(1, n_maf_chunks + 1)]
        with Pool(processes=n_workers) as pool:
            pool.starmap(run_genome_nexus_for_one_chunk, args_list)
    else:
        for i in range(1, n_maf_chunks + 1):
            run_genome_nexus_for_one_chunk(i, dataset_dir)


def run_genome_nexus_for_one_chunk(i: int, dataset_dir: str) -> None:
    """Worker script for the run_genome_nexus main script. Runs
        genome nexus on one maf chunk

    Args:
        i (int): maf chunk number
        dataset_dir (str): the directory to pull the maf chunks from
    """
    filename = f"data_mutations_{i}.txt"
    cmd = f"""
    docker run --rm -it -e GENOMENEXUS_BASE=https://grch38.genomenexus.org \
        -v {dataset_dir}:/data genomenexus/gn-annotation-pipeline:master java \
        -jar annotationPipeline.jar \
        --filename /data/{filename} \
        --output-filename /data/data_mutations_annotated_{i}.txt \
        --error-report-location /data/data_mutations_error_report_{i}.txt \
        --isoform-override mskcc
    """
    with open(f"{dataset_dir}/genome_nexus_log_chunk_{i}.txt", "w") as f:
        subprocess.run(
            cmd,
            shell=True,
            executable="/bin/bash",
            check=True,
            stdout=f,
            stderr=subprocess.STDOUT,
        )


def concatenate_mafs(
    dataset_name: str, n_maf_chunks: int, datahub_tools_path: str
) -> Dict[str, pd.DataFrame]:
    """Concatenates the annotated mafs together
        and the error reports of the failed annotations together

    Args:
        dataset_name (str): Name of dataset
        n_maf_chunks (int): Number of maf chunks
        datahub_tools_path (str): Path to the datahub tools repo

    Returns:
        Dict[str, pd.DataFrame]: A dictionary with the following keys:
            "annotated_maf" : the annotated maf dataset concatenated together
            "error_maf" : error report outputted by genome nexus of failed annotations
    """
    dataset_dir = os.path.join(
        f"{datahub_tools_path}/add-clinical-header/", dataset_name
    )
    annotated_mafs = []
    error_mafs = []
    for i in range(1, n_maf_chunks + 1):
        annotated_mafs.append(
            pd.read_csv(
                f"{dataset_dir}/data_mutations_annotated_{i}.txt", sep="\t", comment="#"
            )
        )
        error_mafs.append(
            pd.read_csv(f"{dataset_dir}/data_mutations_error_report_{i}.txt", sep="\t")
        )

    annotated_mafs_all = pd.concat(annotated_mafs)
    error_mafs_all = pd.concat(error_mafs)

    annotated_mafs_all_processed = postprocessing(input_df=annotated_mafs_all)
    annotated_mafs_all_processed.to_csv(
        f"{dataset_dir}/data_mutations.txt", sep="\t", index=False, float_format="%.12g"
    )

    error_mafs_all.to_csv(
        f"{dataset_dir}/data_mutations_error_report.txt",
        sep="\t",
        index=False,
        float_format="%.12g",
    )
    return {"annotated_maf": annotated_mafs_all, "error_maf": error_mafs_all}


def postprocessing(input_df: pd.DataFrame) -> pd.DataFrame:
    """Additional postprocessing steps we have to do for the
       datasets. This can be removing variants that
       cannot be annotated, etc.

    Args:
        input_df (pd.DataFrame): input annotated maf file

    Returns:
        pd.DataFrame: output maf file postprocessed
    """
    annotated_maf = input_df.copy()
    # remove the failed variants on MT chromosome
    annotated_maf = annotated_maf[annotated_maf.Chromosome != "chrM"]
    return annotated_maf


def generate_meta_files(dataset_name: str, datahub_tools_path: str) -> None:
    """Generates the meta* files for the given dataset

    Args:
        dataset_name (str): name of the iatlas dataset
        datahub_tools_path (str): Path to the datahub tools repo
    """
    dataset_dir = os.path.join(
        f"{datahub_tools_path}/add-clinical-header/", dataset_name
    )

    cmd = f"""
    python3 {datahub_tools_path}/generate-meta-files/generate_meta_files.py \
        -d {dataset_dir} \
        -s iatlas_{dataset_name} \
        -m {datahub_tools_path}/generate-meta-files/datatypes.txt
    """
    # Run in shell to allow sourcing
    subprocess.run(cmd, shell=True, executable="/bin/bash")


def validate_export_files(
    input_df: pd.DataFrame, output_df: pd.DataFrame, **kwargs
) -> None:
    """Validates the export files, doing basic checks

        Validation # 1: Check that number of rows are equal
        Validation # 2: Check that there are no duplicates
        Validation # 3: Check that the output's Tumor_Sample_Barcode values
         exist in the input df

    Args:
        input_df (pd.DataFrame): input maf data
        output_df (pd.DataFrame): output annotated maf data
    """
    logger = kwargs.get("logger", logging.getLogger(__name__))
    if len(input_df) != len(output_df):
        logger.error(
            f"Output rows {len(output_df)} are not equal to input rows {len(input_df)}."
        )
    # no dups
    if len(output_df[output_df.duplicated()]) > 0:
        logger.error("There are duplicates in the output.")
    # check that the Tumor_Sample_Barcode exists in original maf
    if set(list(output_df.Tumor_Sample_Barcode.unique())) != set(
        list(input_df.Tumor_Sample_Barcode.unique())
    ):
        logger.error(
            "The Tumor_Sample_Barcode values are not equal in the output compared to input."
        )


def validate_that_allele_freq_are_not_na(
    input_df: pd.DataFrame,
    **kwargs,
) -> None:
    """Validation: Checks that there are no NAs in the
        maf allele frequency columns: t_ref_count, t_alt_count
        as these are required to calculate the allele frequency(AF):
            AF = t_alt_count / (t_alt_count + t_ref_count)

    Args:
        input_df (pd.DataFrame): input dataframe with allele freq columns
    """
    logger = kwargs.get("logger", logging.getLogger(__name__))
    # check that allele _freq are present
    allele_freq_cols = ["t_ref_count", "t_alt_count"]
    if set(allele_freq_cols) <= set(input_df.columns):
        if input_df[allele_freq_cols].isna().any().any():
            logger.error(
                f"There are NAs in the allele frequency columns: {allele_freq_cols}"
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset",
        type=str,
        help="Name of dataset to run processing for",
    )
    parser.add_argument(
        "--input_folder_synid",
        type=str,
        help="Synapse id for the input folder containing all of the mafs to be merged",
    )
    parser.add_argument(
        "--max_rows",
        type=int,
        help="Max rows per maf chunk used to split maf file. Default: 40000",
        default=40000,
    )
    parser.add_argument(
        "--datahub_tools_path",
        type=str,
        help="Path to datahub-study-curation-tools repo",
    )
    parser.add_argument(
        "--n_workers",
        type=int,
        default=None,
        help="Number of worker processes to run genome nexus in parallel. Optional. Defaults to None.",
    )
    parser.add_argument(
        "--clear_workspace",
        action="store_true",
        default=False,
        help="Whether to clear local directory of files or not",
    )
    
    args = parser.parse_args()
    if args.clear_workspace:
        utils.clear_workspace(dir_path=f"{args.datahub_tools_path}/add-clinical-header")
    dataset_flagger = utils.ErrorFlagHandler()
    dataset_logger = utils.create_logger(
        dataset_name=args.dataset,
        datahub_tools_path=args.datahub_tools_path,
        log_file_name="iatlas_maf_validation_log.txt",
        flagger=dataset_flagger
    )
    maf_df = read_and_merge_maf_files(input_folder_synid=args.input_folder_synid)
    n_maf_chunks = split_into_chunks(
        dataset_name=args.dataset,
        input_df=maf_df,
        datahub_tools_path=args.datahub_tools_path,
        max_rows=args.max_rows,
    )
    run_genome_nexus(
        dataset_name=args.dataset,
        n_maf_chunks=n_maf_chunks,
        n_workers=args.n_workers,
        datahub_tools_path=args.datahub_tools_path,
    )
    mafs = concatenate_mafs(
        dataset_name=args.dataset,
        n_maf_chunks=n_maf_chunks,
        datahub_tools_path=args.datahub_tools_path,
    )
    validate_export_files(
        input_df=maf_df, output_df=mafs["annotated_maf"], logger=dataset_logger
    )
    validate_that_allele_freq_are_not_na(mafs["annotated_maf"], logger=dataset_logger)
    generate_meta_files(
        dataset_name=args.dataset, datahub_tools_path=args.datahub_tools_path
    )
    if dataset_flagger.had_error:
        dataset_logger.error("FAILED: Validation of study failed")


if __name__ == "__main__":
    main()
