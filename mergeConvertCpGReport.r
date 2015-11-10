# Description: Merge the cytosin_reports of bismark_methylation_extractor from PE and SE aligned reads. Optionally convert input format to methylKit.

# Usage: Rscript merge_CpGreport.r -help SE1.CpG_report.txt SE2.CpG_report.txt PE.CpG_report.txt -methylkit=T/F 

# Author: Sebastian DiLorenzo, sebastian.dilorenzo@bils.se

# Date: 2015-11-02

# This is the input format for methylKit
# chrBase chr base    strand  coverage    freqC   freqT
# chr21.9764539   chr21   9764539 R   12  25.00   75.00
# chr21.9764513   chr21   9764513 R   12  0.00    100.00
# chr21.9820622   chr21   9820622 F   13  0.00    100.00
# chr21.9837545   chr21   9837545 F   11  0.00    100.00
# chr21.9849022   chr21   9849022 F   124 72.58   27.42
# chr21.9853326   chr21   9853326 F   17  70.59   29.41

#Input arguments
arg <- commandArgs(trailingOnly = TRUE)

#Show usage
if (arg[1]=="-help" | arg[1]=="-h"){
	stop("Usage: Rscript mergeConvertCpGReport.r SE1.CpG_report.txt SE2.CpG_report.txt PE.CpG_report.txt -methylkit=<T or F>")
}

cat("Reading files.\n")
#Load the files
file1 <- read.table(arg[1],colClasses = c("character","integer","character","integer","integer","character","character"))
file2 <- read.table(arg[2],colClasses = c("character","integer","character","integer","integer","character","character"))
file3 <- read.table(arg[3],colClasses = c("character","integer","character","integer","integer","character","character"))

#Set column names
colnames(file1) <- c("chr","pos","strand","count_methylated","count_unmethylated","C-context","trinucleotide-context")
colnames(file2) <- c("chr","pos","strand","count_methylated","count_unmethylated","C-context","trinucleotide-context")
colnames(file3) <- c("chr","pos","strand","count_methylated","count_unmethylated","C-context","trinucleotide-context")

cat("Merging files.\n")
#Merge the files
tmp_merged <- merge(file1,file2,by.x=c("chr","pos","strand","C-context","trinucleotide-context"),by.y=c("chr","pos","strand","C-context","trinucleotide-context"),all=TRUE)
merged <- merge(tmp_merged,file3,by.x=c("chr","pos","strand","C-context","trinucleotide-context"),by.y=c("chr","pos","strand","C-context","trinucleotide-context"),all=TRUE)

#Clean
rm(file1)
rm(file2)
rm(file3)
rm(tmp_merged)

cat("Modifying files.\n")
#Set all missing counts to zero
merged[is.na(merged)] <- 0

#Add the counts of different files
merged$tot_count_methylated <- merged$count_methylated.x + merged$count_methylated.y +merged$count_methylated
merged$tot_count_unmethylated <- merged$count_unmethylated.x + merged$count_unmethylated.y + merged$count_unmethylated

#Remove unecessary columns
merged <- merged[ , -which(names(merged) %in% c("count_methylated.x","count_unmethylated.x","count_methylated.y","count_unmethylated.y","count_methylated","count_unmethylated"))]

#Reorder the columns
merged <- merged[,c(1,2,3,6,7,4,5)]

#Sort columns by chromosome and position
merged <- merged[with(merged,order(chr,pos)), ]

#Get a good sample name
tmp_name <- unlist(strsplit(arg[1],"/"))
tmp_name2 <- unlist(strsplit(tmp_name[length(tmp_name)],"[.]"))
rm(tmp_name)

#Output the merged file
write.table(merged,file=paste(tmp_name2[1],".Merged.CpG_report.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F)

#cleanup and save image just incase
save.image(paste(tmp_name2[1],".Merged.Rdata",sep=""))

## If we want to convert the format to methylKit
if (arg[4]=="-methylkit=T" | arg[4]=="-methylkit=TRUE"){
	cat("Starting methylkit conversion.\n")
	merged$tot_coverage <- merged$tot_count_methylated + merged$tot_count_unmethylated

	#Remove rows with zero tot_coverage
	merged <- merged[!(merged$tot_coverage == 0),]

	#Reassign strand -/+ to R/F
	merged$strand[merged$strand=="-"] <- "R"
	merged$strand[merged$strand=="+"] <- "F"

	#Insert chr infront of all chromosome numbers
	#Check if "chr" already there, if not, insert it NEEDS WORK
	if (!any(grep("chr",merged$chr[1])))
	{
		merged$chr <- sub("^","chr",merged$chr)
	}

	#Create freqC and freqT
	merged$freqC <- (merged$tot_count_methylated / merged$tot_coverage) * 100
	merged$freqT <- (merged$tot_count_unmethylated / merged$tot_coverage) * 100

	#Drop tot_counts
	merged <- merged[ , -which(names(merged) %in% c("tot_count_methylated","tot_count_unmethylated","C-context","trinucleotide-context"))]

	#Create chrBase column
	merged$chrBase <- paste(merged$chr,merged$pos,sep=".")

	#Reorder columns
	merged <- merged[,c(7,1,2,3,4,5,6)]

	#Save output to file
	write.table(merged,file=paste(tmp_name2[1],".MergedAndConverted.CpG_report.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F)

	#cleanup and save image just incase
	save.image(paste(tmp_name2[1],".MergedAndConverted.Rdata",sep=""))
} else {
	cat("Not methylkit converted.\n")
}

#write.table(merged,file="wip.txt",quote=F,sep="\t",row.names=F,col.names=F)

#save.image("wip.Rdata")


