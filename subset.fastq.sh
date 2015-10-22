#! /bin/bash -l

#Author: Sebastian DiLorenzo
#Date: 2015-03-03

#Grab first x rows from fastq file.
#Call as sh subset.fastq.sh filename.fastq.gz #lines

zcat $1 | head -n $2 - | gzip - > tiny$1


