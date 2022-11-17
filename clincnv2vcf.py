import argparse
import os.path
from time import strftime

import pandas as pd
from pyfaidx import Fasta


def main(reference_genome: Fasta, input_file: str, sample_id: str) -> None:
    df_in = pd.read_csv(
        input_file, delimiter="\t", skiprows=check_skip_lines(input_file)
    )

    df_out = pd.DataFrame()
    df_out["#CHROM"] = df_in["#chr"]
    df_out["POS"] = df_in["start"]
    df_out["ID"] = "."
    df_out["REF"] = df_in.apply(
        lambda x: get_base(reference_genome, x["#chr"], x["start"]), axis=1
    )
    df_out["ALT"] = df_in["CN_change"].apply(lambda x: f"<{get_DELDUP_by_CN(x)}>")
    df_out["QUAL"] = df_in["loglikelihood"]
    df_out["FILTER"] = "PASS"
    df_out["INFO"] = (
        df_in["CN_change"].apply(lambda x: f"SVTYPE={get_DELDUP_by_CN(x)};")
        + df_in["end"].astype(str)
        + ";SVLEN="
        + (df_in["end"] - df_in["start"]).astype(str)
        + ";"
        + "CIPOS=-500,500;CIEND=-500,500"
    )
    df_out["FORMAT"] = "GT:CN"

    df_out[sample_id] = "./.:" + df_in["CN_change"].astype(str)

    df_out.fillna(".", inplace=True)

    with open(f"{sample_id}_cnvs_clincnv.vcf", "w") as f:
        f.write(VCF_HEAD + df_out.to_csv(sep="\t", index=False))


def get_base(reference_genome: Fasta, chromosome: str, position: int) -> str:
    return reference_genome[chromosome][position - 1 : position].seq


def check_skip_lines(file: str) -> list[int]:
    skip_lines = []
    with open(file) as f:
        for i, line in enumerate(f):
            if line.startswith("##"):
                skip_lines.append(i)
    return skip_lines


def get_DELDUP_by_CN(cn: int) -> str:
    if cn < 2:
        return "DEL"
    elif cn > 2:
        return "DUP"
    else:
        return "."


if __name__ == "__main__":
    VCF_HEAD = f"""##fileformat=VCFv4.3
##fileDate={strftime("%Y%m%d")}
##reference=file:///NGS_Daten_nfs/ngs-bits/data/genomes/GRCh37.fa
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=CNV,Description="Copy number variable region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
"""

    parser = argparse.ArgumentParser(
        description="Convert ClinCNV TSV file(s) to VCF-format."
    )
    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        help="reference genome in FASTA format.",
        type=str.strip,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        nargs="+",
        help="input file(s) in ClinCNV TSV format.",
        type=str.strip,
    )
    parser.add_argument(
        "-s",
        "--sampleid",
        nargs="+",
        help="sample id(s). Will be infered from input filename if not given.",
        type=str.strip,
    )
    args = parser.parse_args()

    if os.path.isfile(args.reference):
        reference_genome = Fasta(args.reference)
    else:
        raise FileNotFoundError(f"Reference genome not found: {args.reference}")

    if args.sampleid and len(args.input) != len(args.sampleid):
        raise ValueError("Number of input files and sample ids do not match.")

    for index, file in enumerate(args.input):
        if os.path.isfile(file):
            if not args.sampleid:
                sample_id = (
                    os.path.basename(file).split(".")[0].removesuffix("_cnvs_clincnv")
                )
            else:
                sample_id = args.sampleid[index]
            main(reference_genome, file, sample_id)
        else:
            raise FileNotFoundError(f"Input file not found: {file}")
