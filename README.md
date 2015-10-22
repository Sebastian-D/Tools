# Tools
I thought maybe I should commit some of my more useful tools to git for easier future access (and the good of my fellow man?)

## remove_unpaired.py
For the special case of having a SAM/BAM file with a mix of paired-end and single-end reads where the single-end reads dont have the correct flags for easy SAMTOOLS removal. This file separates the reads into two files.

## subset.fastq.sh
Grab the first nr of lines from a fastq.gz file so you can dryrun your pipeline/tool/method on this subset of data.
