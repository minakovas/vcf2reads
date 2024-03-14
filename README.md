# vcf2reads
This is a CLI tool for generating reads in fasta format based on variants from vcf file and reference genome.

It was written in February 2023 as a test assignment for the position of a bioinformatics specialist. It gets a reference genome and a vcf file with genomic variants as input and generates a new fasta file with random fragments of the genome of a given length, where the positions containing genomic variants are replaced with alleles from the vcf file.  

At the moment it supports SNPs and indels and filtering variants by allele frequency (AF). Only python standard library is used.

## Input
File with genomic varians in vcf format. Example:
```
##fileformat=VCFv4.2
##fileDate=20240313
##source=PLINKv1.90
##contig=<ID=1,length=1002>
##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1    sample2
1       1000    rs0     C       .       .       .       PR      GT      ./.     0/0
1       1001    rs10    A       G       .       .       PR      GT      0/0     0/1
```
File with reference genome in fasta format. Example:
```
>chr1
AGCAGGTGGGTAGGGAGGGGTTCGTGAGCTGCAGTGTGGGGGGAGGCTCCTGGCTGTGCC
CAGTCCTCCTGCCCTGGAGGGAAAGGCTCTGGAAAACACAGCCCAGAGTGCTCACCCTAC
>chr2
TTCTTCAGCTTATCCAGCTGCGTCTGGGCATAACGCATCTGGTTTTCAACACTCTTCAGC
TCATCATCCAGCTCCAGCTGCTGTCGTTTAAATTCACTGTATTCTCTCTGATACCTGTGA
```
## Usage
```
vcf2reads.py [-h] --genome GENOME --vcf VCF --read_len READ_LEN --N N --out OUT --sample SAMPLE [--random_seed RANDOM_SEED] [--af_min AF_MIN] [--af_max AF_MAX]
```
Arguments:
* `genome` - path to <genome.fasta>, required.
* `vcf` - path to <file.vcf.vgz>, required.
* `read_len` - length of read to generate, required.
* `N` - number of reads to generate, required.
* `out` - path to output file, required.
* `sample` - name of the sample in vcf, required.
* `random_seed` - optional, default 42.
* `af_min` - lower threshold for filtering by AF, optional. If not specified, filtering by AF is not applied.
* `af_max` - higher threshold for filtering by AF, optional. If not specified, filtering by AF is not applied.

Example:
```
python vcf2reads.py --genome hg38.fasta \
--vcf variants.vcf \
--sample SAMPLE1 \
--read_len 150 \
--N 100000 \
--out results/reads.fasta \
--random_seed 88
```
## Output
Returns fasta file with random fragments of genome with variants from vcf file. Format:
```
>Sample=<sample name>:Chr=<chromosome number>:Pos=<position of variant>:Ref=<reference allele>:Alt=<alternative allele>:GT=<genotype 0 or 1>
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>...
```
Example:
```
>Sample=HG02771:Chr=chrY:Pos=14385140:Ref=C:Alt=T:GT=1
TAAATCTTCACAACAGCCATGAGTTAGAAGATATCATTATAACAACATTG
>Sample=HG02771:Chr=chrY:Pos=16955589:Ref=T:Alt=C:GT=1
tCTCACCAAATTACAAGTTCAACAAACTTTTGCAACTTGCATCGCTTGCA
...
```
