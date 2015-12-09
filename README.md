# Tools
I thought maybe I should commit some of my more useful tools to git for easier future access (and the good of my fellow man?)

#### remove_unpaired.py
For the special case of having a SAM/BAM file with a mix of paired-end and single-end reads where the single-end reads dont have the correct flags for easy SAMTOOLS removal. This file separates the reads into two files.

#### subset.fastq.sh
Grab the first nr of lines from a fastq.gz file so you can dryrun your pipeline/tool/method on this subset of data.

#### Create_Cosmic.sh
Create a COSMIC file from any version of their database such as exists in the GATK resource files. Useful for example when using MuTect.

#### mergeConvertCpGReport.r
Merge Single-end and Paired-end aligned bismark_methylation_extractor Cytosine reports. (*.CpG_report.txt). Optional conversion to [methylKit](https://github.com/al2na/methylKit) input format.
Useful for when you had a low mapping efficiency and want as much aligned content as possible.

#### CpG_report2methylKit.r
Convert a single CpG report (cytosine report) to methylKit input format.

#### CpG_island_coverage_histogram.r
Produce coverage plots of methylation data.

Usage: `Rscript CpG_island_coverage_histograms.r methylRawObject.Rdata Hsapien_CpG_islands_hg19.txt`

The Hsapien_CpG_islands_hg19.txt can be downloaded from here. Any file with same format will work.
How to create the methylrawobject, which must be named myobj:
```R
library(methylKit)
myobj=read(<arguments>)
save(myobj,file="methylRawObject.Rdata")
```


