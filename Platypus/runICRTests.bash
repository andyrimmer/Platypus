#!/bin/bash

python Platypus.py callVariants --bamFiles=dataTests/17-ICRTests/BAMs/V009.bam --refFile=/home/rimmer/Genomes/human_g1k_v37_ebv.fa --regions=14:95572021-95572221 --output=ICR_Test_V009.vcf $*
python Platypus.py callVariants --bamFiles=dataTests/17-ICRTests/BAMs/V025.bam --refFile=/home/rimmer/Genomes/human_g1k_v37_ebv.fa --regions=13:32905038-32905238 --output=ICR_Test_V025.vcf $*
python Platypus.py callVariants --bamFiles=dataTests/17-ICRTests/BAMs/V028.bam --refFile=/home/rimmer/Genomes/human_g1k_v37_ebv.fa --regions=17:41246233-41246433 --output=ICR_Test_V028.vcf $*
python Platypus.py callVariants --bamFiles=dataTests/17-ICRTests/BAMs/V046.bam --refFile=/home/rimmer/Genomes/human_g1k_v37_ebv.fa --regions=5:176720830-176721030 --output=ICR_Test_V046.vcf $*
python Platypus.py callVariants --bamFiles=dataTests/17-ICRTests/BAMs/V051.bam --refFile=/home/rimmer/Genomes/human_g1k_v37_ebv.fa --regions=11:32410500-32410700 --output=ICR_Test_V051.vcf $*
python Platypus.py callVariants --bamFiles=dataTests/17-ICRTests/BAMs/V059.bam --refFile=/home/rimmer/Genomes/human_g1k_v37_ebv.fa --regions=10:104263846-104264046 --output=ICR_Test_V059.vcf $*
python Platypus.py callVariants --bamFiles=dataTests/17-ICRTests/BAMs/V069.bam --refFile=/home/rimmer/Genomes/human_g1k_v37_ebv.fa --regions=7:148504682-148504882 --output=ICR_Test_V069.vcf $*
python Platypus.py callVariants --bamFiles=dataTests/17-ICRTests/BAMs/V070.bam --refFile=/home/rimmer/Genomes/human_g1k_v37_ebv.fa --regions=17:41246483-41246683 --output=ICR_Test_V070.vcf $*
python Platypus.py callVariants --bamFiles=dataTests/17-ICRTests/BAMs/V900.bam --refFile=/home/rimmer/Genomes/human_g1k_v37_ebv.fa --regions=6:79424333-79424533 --output=ICR_Test_V900.vcf $*
python Platypus.py callVariants --bamFiles=dataTests/17-ICRTests/BAMs/V901.bam --refFile=/home/rimmer/Genomes/human_g1k_v37_ebv.fa --regions=6:79424333-79424533 --output=ICR_Test_V901.vcf $*
