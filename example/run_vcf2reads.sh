#!/usr/bin/bash

# You should specify your own paths
GENOME=/home/alex/Desktop/BTK_task/fasta/chrY_PA.fna
VCF1=/home/alex/Desktop/BTK_task/vcf/gnomad_chrY_10samples_50000entries_sorted.vcf.vgz
VCF2=/home/alex/Desktop/BTK_task/vcf/gnomad.genomes.v3.1.2.hgdp_tgp.chrY_samples_10.vcf.vgz

python vcf2reads.py --genome $GENOME --vcf $VCF2 --sample HG03660 --read_len 150 --N 30000 --out result1.fasta --random_seed 13
python vcf2reads.py --genome $GENOME --vcf $VCF3 --sample HG02771 --read_len 50 --N 100000 --out result2.fasta --random_seed 12 --af_min 0.5 --af_max 0.99
