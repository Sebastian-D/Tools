# NGS pipeline
This is a pipeline for NGS data implemented in bpipe. It does not have normal variant calling but will rather include a module for cancer/normal variant calling in the future. It is closely based on the GATK best practices and runs multiple paired-end lanes in parallel with configuration options for multithreading and memory allocation.

#### Data
This pipeline handles whole-genome sequencing DNA data and exome data.

#### Installation & guide

1. Install Bpipe.
2. Download NGSpipe.pipe.
3. Edit NGSpipe.pipe with your own constants such as indel/dbsnp files and paths to the required programs.
4. Edit the runcommand so the pattern matches your files. At the moment it is L00%_R*.fastq.gz which will match a typical sample pair to several pairs of files in the format:
```
samplename_L001_R1.fastq.gz
samplename_L001_R2.fastq.gz
samplename_L003_R1.fastq.gz
samplename_L003_R2.fastq.gz
```

Execute `$bpipe run NGSpipe.pipe sample*.fastq.gz` to run on all files, sorted by lanes and read number.


