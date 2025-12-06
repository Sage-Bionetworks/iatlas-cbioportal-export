import argparse
import os
import shutil
from typing import Dict

import utils

syn = utils.synapse_login()

def parse_metadata_kv(text: str) -> Dict[str, str]:
    """
    Parse a simple 'key: value' metadata file into a dict.
    This matches the format you showed (one key per line).
    Lines starting with '#' or blank lines are ignored.

    Args:
        text (str): _description_

    Returns:
        Dict[str, str]: _description_
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
        meta (Dict[str, str]): _description_

    Returns:
        str: _description_
    """
    return "\n".join(f"{k}: {v}" for k, v in meta.items()) + "\n"


def download_and_patch_files(
    syn : "synapseclient.Synapse",
    data_synid: str,
    meta_synid: str,
    out_dir: str,
    out_data_filename: str,
    out_meta_filename: str,
    cancer_study_identifier: str,
    data_filename_in_meta: str,
) -> None:
    """_summary_

    Args:
        syn (synapseclient.Synapse): _description_
        data_synid (str): _description_
        meta_synid (str): _description_
        out_dir (str): _description_
        out_data_filename (str): _description_
        out_meta_filename (str): _description_
        cancer_study_identifier (str): _description_
        data_filename_in_meta (str): _description_
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
