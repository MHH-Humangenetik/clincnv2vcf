import argparse
import os.path
import re
from time import strftime

import pandas as pd
from pyfaidx import Fasta


def main(reference_genome: Fasta, input_file: str, sample_id: str, ci: int) -> None:
    folder = os.path.dirname(input_file)
    skip_lines, metadata = parse_metadata(input_file)
    df_in = pd.read_csv(input_file, delimiter="\t", skiprows=skip_lines)
    if metadata["gender of sample"].lower() == "m":
        male = True
    else:
        male = False

    df_out = pd.DataFrame()
    df_out["#CHROM"] = df_in["#chr"]
    df_out["POS"] = df_in["start"]
    df_out["ID"] = "."
    df_out["REF"] = df_in.apply(
        lambda x: get_base(reference_genome, x["#chr"], x["start"]), axis=1
    )
    df_out["ALT"] = df_in.apply(
        lambda x: get_DELDUP_by_CN(
            x["CN_change"],
            ltgt=True,
            male_x=(male and (x["#chr"] == "X") or (x["#chr"] == "chrX")),
        ),
        axis=1,
    )
    df_out["QUAL"] = df_in[
        "loglikelihood"
    ].abs()  # Use absolute value of loglikelihood because clincnv sometimes reports negative values
    df_out["FILTER"] = df_in.apply(
        lambda x: "PASS"
        if get_DELDUP_by_CN(
            x["CN_change"],
            male_x=(male and (x["#chr"] == "X") or (x["#chr"] == "chrX")),
        )
        != "."
        else f"CN={x['CN_change']}",
        axis=1,
    )
    df_out["INFO"] = (
        df_in.apply(
            lambda x: f"SVTYPE={get_DELDUP_by_CN(x['CN_change'], male_x=(male and ((x['#chr'] == 'X') or (x['#chr'] == 'chrX'))))}",
            axis=1,
        )
        + ";END="
        + df_in["end"].astype(str)
        + ";SVLEN="
        + (df_in["end"] - df_in["start"]).astype(str)
        + f";CIPOS=-{ci},{ci}"
        + f";CIEND=-{ci},{ci}"
        + ";NOREGIONS="
        + df_in["no_of_regions"].astype(str)
        + ";PAF="
        + df_in["potential_AF"].astype(str)
        + ";GENES="
        + df_in["genes"].apply(clean_info)
        + ";QVALUE="
        + df_in["qvalue"].astype(str)
        + ";OLAPAFIMGAG="
        + df_in["overlap af_genomes_imgag"].apply(clean_info)
        + ";CNPATO="
        + df_in["cn_pathogenic"].apply(clean_info)
        + ";CLINVAR="
        + df_in["clinvar_cnvs"].apply(clean_info)
        + ";GENEINFO="
        + df_in["gene_info"].apply(clean_info)
        + ";OWNPATO="
        + df_in["ngsd_pathogenic_cnvs"].apply(clean_info)
    )
    df_out["FORMAT"] = "GT:CN"

    df_out[sample_id] = "./.:" + df_in["CN_change"].astype(str)

    df_out.fillna(".", inplace=True)

    with open(os.path.join(folder, f"{sample_id}_cnvs_clincnv.vcf"), "w") as f:
        VCF_HEAD = f"""##fileformat=VCFv4.3
##fileDate={strftime("%Y%m%d")}
##reference=file://{reference_genome.filename}
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant",Source="ClinCNV",Version={metadata["ClinCNV version"]}>
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record",Source="ClinCNV",Version={metadata["ClinCNV version"]}>
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles",Source="ClinCNV",Version={metadata["ClinCNV version"]}>
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants",Source="ClinCNV",Version={metadata["ClinCNV version"]}>
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants",Source="ClinCNV",Version={metadata["ClinCNV version"]}>
##INFO=<ID=NOREGIONS,Number=1,Type=Integer,Description="Number of regions in the CNV",Source="ClinCNV",Version={metadata["ClinCNV version"]}>
##INFO=<ID=PAF,Number=1,Type=Float,Description="Potential allelic fraction",Source="ClinCNV",Version={metadata["ClinCNV version"]}>
##INFO=<ID=GENES,Number=.,Type=String,Description="List of genes affected by the CNV",Source="ClinCNV",Version={metadata["ClinCNV version"]}>
##INFO=<ID=QVALUE,Number=1,Type=Float,Description="Q-value of the CNV",Source="ClinCNV",Version={metadata["ClinCNV version"]}>
##INFO=<ID=OLAPAFIMGAG,Number=.,Type=String,Description="Overlap with inhouse database of IMGAG and GnomAD",Source="ClinCNV",Version={metadata["ClinCNV version"]}>
##INFO=<ID=CNPATO,Number=.,Type=String,Description="Known pathogenic CNVs",Source="ClinCNV",Version={metadata["ClinCNV version"]}>
##INFO=<ID=CLINVAR,Number=.,Type=String,Description="ClinVar overlapping CNVs",Source="ClinCNV",Version={metadata["ClinCNV version"]}>
##INFO=<ID=GENEINFO,Number=.,Type=String,Description="Gene informations (i.e. stringency, region, etc.)",Source="ClinCNV",Version={metadata["ClinCNV version"]}>
##INFO=<ID=OWNPATO,Number=.,Type=String,Description="Known CNVs in inhouse database",Source="ClinCNV",Version={metadata["ClinCNV version"]}>
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=CNV,Description="Copy number variable region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
"""
        f.write(VCF_HEAD + df_out.to_csv(sep="\t", index=False))


def get_base(reference_genome: Fasta, chromosome: str, position: int) -> str:
    return reference_genome[chromosome][position - 1 : position].seq


def parse_metadata(file: str) -> tuple[list[int], dict[str, str]]:
    skip_lines = []
    metadata = {}
    metadata_pattern = re.compile(r"##(.*?):\s*(.*)")
    with open(file) as f:
        for i, line in enumerate(f):
            if line.startswith("##"):
                skip_lines.append(i)
                metadata_match = re.match(metadata_pattern, line)
                if metadata_match:
                    metadata[metadata_match.group(1)] = metadata_match.group(2)
    return skip_lines, metadata


def get_DELDUP_by_CN(cn: int, ltgt: bool = False, male_x: bool = False) -> str:
    if male_x:
        cn_decision = 1
    else:
        cn_decision = 2
    if cn == cn_decision:
        return "."

    return_value = ""
    if cn < cn_decision:
        return_value = "DEL"
    elif cn > cn_decision:
        return_value = "DUP"

    if ltgt:
        return_value = f"<{return_value}>"

    return return_value


def clean_info(info: str) -> str:
    if pd.isna(info):
        return "."
    else:
        return re.sub(r"[^a-zA-Z0-9\.]", "_", info)


if __name__ == "__main__":

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
    parser.add_argument(
        "-c",
        "--confidenceinterval",
        default="500",
        help="confidence interval around POS/END for imprecise variants. Will be used in in both directions. Defaults to `500` resulting in `CIPOS=-500,500;CIEND=-500,500`.",
        type=int,
    )
    args = parser.parse_args()

    if os.path.isfile(args.reference):
        reference_genome = Fasta(filename=args.reference)
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
            main(reference_genome, file, sample_id, args.confidenceinterval)
        else:
            raise FileNotFoundError(f"Input file not found: {file}")
