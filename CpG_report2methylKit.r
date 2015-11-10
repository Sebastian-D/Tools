# Description: Convert the cytosine_report of bismark_methylation_extractor to methylKit input format.

# Usage: Rscript CpG_report2methylKit.r CpG_report.txt

# Author: Sebastian DiLorenzo, sebastian.dilorenzo@bils.se

# Date: 2015-11-09

#Input arguments
arg <- commandArgs(trailingOnly = TRUE)

#Load file
file1 <- read.table(arg[1],colClasses = c("character","integer","character","integer","integer","character","character"))
#file1 <- read.table("BRNO1003.merged.dedupped.picard_sort.CpG_report.txt",colClasses = c("character","integer","character","integer","integer","character","character"))

#Set column names
colnames(file1) <- c("chr","pos","strand","count_methylated","count_unmethylated","C-context","trinucleotide-context")

#Get total coverage
file1$tot_coverage <- file1$count_unmethylated + file1$count_methylated

#Remove rows with zero tot_coverage
file1 <- file1[!(file1$tot_coverage == 0),]

#Reassign strand -/+ to R/F
file1$strand[file1$strand=="-"] <- "R"
file1$strand[file1$strand=="+"] <- "F"

#Insert chr infront of all chromosome numbers if it is needed
if (!any(grep("chr",file1$chr[1])))
{
	file1$chr <- sub("^","chr",file1$chr)
}

#Create freqC and freqT
file1$freqC <- (file1$count_methylated / file1$tot_coverage) * 100
file1$freqT <- (file1$count_unmethylated / file1$tot_coverage) * 100

#Drop tot_counts
file1 <- file1[ , -which(names(file1) %in% c("count_methylated","count_unmethylated","C-context","trinucleotide-context"))]

#Sort columns by chromosome and position
file1 <- file1[with(file1,order(chr,pos)), ]

#Create chrBase column
file1$chrBase <- paste(file1$chr,file1$pos,sep=".")

#Reorder columns
file1 <- file1[,c(7,1,2,3,4,5,6)]

#Get a good sample name
tmp_name <- unlist(strsplit(arg[1],"[.]"))
tmp_name2 <- unlist(strsplit(tmp_name[1],"/"))
rm(tmp_name)

write.table(file1,file=paste(tmp_name2[length(tmp_name2)],".methylKit_input.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F)

#cleanup and save image just incase
save.image(paste(tmp_name2[length(tmp_name2)],".methylkitConverted.Rdata",sep=""))

write.table(file1,file="BRNO1003.methylKit_input.txt")
save.image()
