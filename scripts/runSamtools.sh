#!/bin/bash

REF=$1
BAM=$2
REGION=$3

samtools mpileup -uf $REF $BAM -r $REGION | bcftools view -bvcg - > var.raw.bcf  
bcftools view var.raw.bcf | /home/rimmer/Downloads/samtools-0.1.16/bcftools/vcfutils.pl varFilter -D100 > var.flt.vcf  
