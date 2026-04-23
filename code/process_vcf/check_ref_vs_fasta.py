import sys
import pysam
from cyvcf2 import VCF

VCF_PATH = "/home/jqian54/sulab/enformer_fine_tuning/data/GTEx_v8_SNP_array_VCF/gtex_genotypes_SNPsOnly.bcf.gz"
FASTA_PATH = "/home/jqian54/sulab/enformer_fine_tuning/data/hg38_nochr.fa"

OUT_TSV = "ref_fasta_mismatches.tsv"
OUT_TXT = "ref_fasta_mismatch_count.txt"

vcf = VCF(VCF_PATH)
fasta = pysam.FastaFile(FASTA_PATH)

mismatch_count = 0

with open(OUT_TSV, "w") as tsv:
    tsv.write("CHROM\tPOS\tVCF_REF\tFASTA_SEQ\n")

    for rec in vcf:
        chrom = rec.CHROM
        pos = rec.POS            # 1-based
        ref = rec.REF

        # FASTA fetch is 0-based, half-open
        fasta_seq = fasta.fetch(chrom, pos - 1, pos - 1 + len(ref))

        if fasta_seq.upper() != ref.upper():
            mismatch_count += 1
            tsv.write(f"{chrom}\t{pos}\t{ref}\t{fasta_seq}\n")

with open(OUT_TXT, "w") as txt:
    txt.write(f"TOTAL_MISMATCHES\t{mismatch_count}\n")

print(f"Done.")
print(f"Mismatch TSV : {OUT_TSV}")
print(f"Mismatch count: {OUT_TXT}")
