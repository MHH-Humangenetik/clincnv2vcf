# CNV2VCF

CNV2VCF is a tool to convert output of [ClinCNV](https://github.com/imgag/ClinCNV) to VCF format.

DEL/DUP interpretation is not (yet) gender aware for the X chromosome.

## Custom info fields

The following custom INFO fields are added to the VCF:
- NOREGIONS: Number of regions in the CNV (no_of_regions)
- PAF: Potential allelic fraction (potential_AF)
- GENES: List of genes affected by the CNV (genes)
- QVALUE: Q-value of the CNV (qvalue)
- OLAPAFIMGAG: Overlap with inhouse database of IMGAG and GnomAD (overlap af_genomes_imgag)
- CNPATO: Known pathogenic CNVs (cn_pathogenic)
- CLINVAR: ClinVar overlapping CNVs (clinvar_cnvs)
- GENEINFO: Gene informations (i.e. stringency, region, etc.) (gene_info)
- IHPATO: Known CNVs in inhouse database (ngsd_pathogenic_cnvs)

## Usage

- `-r REFERENCE, --reference REFERENCE` reference genome in FASTA format.
- `-i INPUT [INPUT ...], --input INPUT [INPUT ...]` input file(s) in ClinCNV TSV format.
- `-s SAMPLEID [SAMPLEID ...], --sampleid SAMPLEID [SAMPLEID ...]` sample id(s). Will be infered from input filename if not given.
- `-c CONFIDENCEINTERVAL, --confidenceinterval CONFIDENCEINTERVAL` confidence interval around POS/END for imprecise variants. Will be used in in both directions. Defaults to `500` resulting in `CIPOS=-500,500;CIEND=-500,500`.

## Attribution

- Uses [pyfaidx](https://pypi.org/project/pyfaidx/) (see Shirley MD, Ma Z, Pedersen BS, Wheelan SJ. 2015. Efficient "pythonic" access to FASTA files using pyfaidx. *PeerJ PrePrints* 3:e970v1 https://doi.org/10.7287/peerj.preprints.970v1)