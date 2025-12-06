import argparse
import os
import shutil
from typing import Dict

from iatlascbioportalexport import utils

syn = utils.synapse_login()


def parse_metadata_kv(text: str) -> Dict[str, str]:
    """
    Parse a simple 'key: value' metadata file into a dict.
    This matches the format you showed (one key per line).
    Lines starting with '#' or blank lines are ignored.

    Args:
        text (str): input text of read in metadata file

    Returns:
        Dict[str, str]: Output metadata dictionary of form:
            'key: value',
            'key: value'
            ...
    """
    meta = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        meta[k.strip()] = v.strip()
    return meta


def write_metadata_kv(meta: Dict[str, str]) -> str:
    """
    Write dict back to 'key: value' lines.
    Preserves no ordering guarantees (unless Python 3.7+ insertion order
    and you don't reorder keys).

    Args:
        meta (Dict[str, str]): Input metadata dictionary of form:
            'key: value',
            'key: value'
            ...

    Returns:
        str: the metadata as text
    """
    return "\n".join(f"{k}: {v}" for k, v in meta.items()) + "\n"


def download_and_patch_files(
    syn: "synapseclient.Synapse",
    data_synid: str,
    meta_synid: str,
    out_dir: str,
    out_data_filename: str,
    out_meta_filename: str,
    cancer_study_identifier: str,
    data_filename_in_meta: str,
) -> None:
    """Downloads the listed file and associated metadata file
        from Synapse, does some modifying of the fields and saves as
        standardized names

    Args:
        syn (synapseclient.Synapse): Synapse client connection
        data_synid (str): Synapse id of input data
        meta_synid (str): Synapse id of metadata file
        out_dir (str): Output directory
        out_data_filename (str): Output data filename
        out_meta_filename (str): output metadata filename
        cancer_study_identifier (str): value to replace with in the
            cancer_study_identifier key in the metadata file
        data_filename_in_meta (str): value to replace with in the
            data_filename key in the metadata file
    """
    # Download both entities
    data_ent = syn.get(data_synid, downloadLocation=out_dir)
    meta_ent = syn.get(meta_synid, downloadLocation=out_dir)

    downloaded_data_path = data_ent.path
    downloaded_meta_path = meta_ent.path

    # Read + patch metadata
    with open(downloaded_meta_path, "r", encoding="utf-8") as f:
        meta_text = f.read()

    meta = parse_metadata_kv(meta_text)

    # update these fields
    meta["cancer_study_identifier"] = cancer_study_identifier
    meta["data_filename"] = data_filename_in_meta

    # Save final data file
    final_data_path = os.path.join(out_dir, out_data_filename)
    shutil.copyfile(downloaded_data_path, final_data_path)

    # Save final metadata file
    final_meta_path = os.path.join(out_dir, out_meta_filename)
    with open(final_meta_path, "w", encoding="utf-8") as f:
        f.write(write_metadata_kv(meta))

    print(f"Wrote:")
    print(f"  data: {final_data_path}")
    print(f"  meta: {final_meta_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Download and patch two (data, metadata) pairs from Synapse."
    )
    parser.add_argument("--dataset", required=True, help="Dataset name")

    # Pair 1: generic assay
    parser.add_argument("--generic-assay-data-synid", required=True)
    parser.add_argument("--generic-assay-metadata-synid", required=True)

    # Pair 2: gene expression
    parser.add_argument("--gene-expression-data-synid", required=True)
    parser.add_argument("--gene-expression-metadata-synid", required=True)

    parser.add_argument(
        "--datahub_tools_path",
        type=str,
        help="Path to datahub-study-curation-tools repo",
    )

    args = parser.parse_args()

    out_dir = utils.get_local_dataset_output_folder_path(
        args.dataset, args.datahub_tools_path
    )

    # generic assay files
    download_and_patch_files(
        syn=syn,
        data_synid=args.generic_assay_data_synid,
        meta_synid=args.generic_assay_metadata_synid,
        out_dir=out_dir,
        out_data_filename="data_rna_seq_mrna.txt",
        out_meta_filename="meta_rna_seq_mrna.txt",
        cancer_study_identifier=f"iatlas_{args.dataset}",
        data_filename_in_meta="data_rna_seq_mrna.txt",
    )

    # gene expression files
    download_and_patch_files(
        syn=syn,
        data_synid=args.gene_expression_data_synid,
        meta_synid=args.gene_expression_metadata_synid,
        out_dir=out_dir,
        out_data_filename="data_gene_signatures.txt",
        out_meta_filename="meta_gene_signatures.txt",
        cancer_study_identifier=f"iatlas_{args.dataset}",
        data_filename_in_meta="data_gene_signatures.txt",
    )


if __name__ == "__main__":
    main()
