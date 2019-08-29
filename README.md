# RealSVsim (SV simulator by spiking in real events)

A simulator to generate genome sequences with real or random genomic variants.

Author: Peng Xu

Email: pxu@uabmc.edu

Draft date: August. 29, 2019

## Description

RealSVsim can simulate a genome with a variety of genomic variants, including SNP, INDEL and structural variations (e.g. deletion, insertion, inversion and translocation). It can spike in real variant events the user provides and also generate random variant events.

## System requirements and dependency

Require Python 3 environment

## Installation

```
git clone git@github.com:penguab/RealSVsim.git
```

## Usage

A reference genome sequence (fasta file format) is required as an input. Users can provide their own real variant events (bed file format), or provide the number of total variants to be simulated. Users can also provide both real events and total variant numbers. If the total variant number is smaller than the number of real events, RealSVsim will randomly choose real events to match the total variant number. If the total variant number is larger than the number of real events, RealSVsim will automatically generate random events to match the total variant number.
```
python3 RealSVsim.py -b SNP_INDEL_SV.bed(optional) -g genome.fa

        -----------optional-----------
                -s snp number
                -d deletion number
                -i insertion number
                -u duplication number
                -v inversion number
                -p translocation pair (How many chromosome pairs are formed)
                -n translocation number (How many translocation events occur within each chromosome pair)
```

Example 1. Real SV events spike-in using GIAB HG002 dataset.
```
#download real SV events from GIAB
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
gzip -d HG002_SVs_Tier1_v0.6.vcf.gz

#Convert vcf to bed file format
perl -lane '@a=split //,$F[3];@b=split //,$F[4];next unless $F[6]eq"PASS" and ($#a==0 or $#b==0);next if $F[3]=~/N/ or $F[4]=~/N/;print join "\t","chr".$F[0],$F[1],$F[1]+abs($#a-$#b),"DEL",abs($#a-$#b),(join "",@a[1..$#a]) if $#b==0;print join "\t","chr".$F[0],$F[1]-1,$F[1],"INS",abs($#a-$#b),(join "",@b[1..$#b]) if $#a==0;' HG002_SVs_Tier1_v0.6.vcf  >HG002_SVs_Tier1_v0.6.SV.bed

#simulate genome by providing real SV events
python3 RealSVsim.py -b HG002_SVs_Tier1_v0.6.SV.bed -g genome.fa
```

Example 1. Random SV events simulation.
```
python3 RealSVsim.py -g genome.fa -d 5000 -i 5000 -u 500 -v 200 -p 10 -n 3
```

## News


