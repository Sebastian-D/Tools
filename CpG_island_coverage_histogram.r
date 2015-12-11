#Description: Produce CpG island coverage histogram. (made for CpG islands from the start, thus all the specific naming)
#
#Usage: Rscript CpG_island_coverage_histograms.r Input.txt Hsapien_CpG_islands_hg19.txt
#
#Author: Sebastian DiLorenzo
#
#Date: 19/11/2015
#
#Input.txt should be a tab separated file with columns chromosome, position, coverage
#example:
#chr1    3000827 8
#chr1    3001007 5
#chr1    3001008 3
#
#Requires: Rcpp

library(Rcpp)

#Input arguments
arg <- commandArgs(trailingOnly = TRUE)

#Read sample coverage data
sample <- read.table(arg[1],colClasses=c("character","integer","integer"),header=F)
colnames(sample) <- c("chr","pos","coverage")

#read CpG island regions
CpG_islands <- read.table(arg[2],colClasses=c("character","integer","integer"),header=T)
#colnames(CpG_islands) <- c("chr","start","end")

#Safeguard wd so that plots are created in correct location
wd <- getwd()

#Get name
#name <- "test.txt"
name <- strsplit(arg[1],"/")
name <- name[[1]][length(name[[1]])]
cat("Starting work on ",name,"\n")

#Get region name
region_name <- strsplit(arg[2],"/")
region_name <- region_name[[1]][length(region_name[[1]])]

#Create Rcpp function
cppFunction('NumericVector CpG_index(NumericVector POS, NumericVector START, NumericVector END) {
        //Try to write a Rcpp function which checks if POS is between START and END
        // should return an index of positions to pick from the coverage of the methylRawObject

        //Create result vector
        NumericVector CpGi_index(POS.size());

        //Iterator over intervals (cpg islands)
        int k = 0;
        for( int i = 0; i < POS.size(); i += 1){
                if(POS[i] >= START[k] && POS[i] <= END[k]){
                        CpGi_index[i] = 1;

                }
                else if(POS[i] > END[k]){
                        k++;
                }
                else{
                        CpGi_index[i] = 0;
                }
        }
        return CpGi_index;

}')


#Only look at chromosomes that have CpG islands according to the CpG_islands list. Remove all others from sample
sample <- sample[which(sample$chr %in% CpG_islands$chr),]

#Only look at chromosomes that exist in this sample. (this is why we reload the list of CpG_islands for each sample)
CpG_islands <- CpG_islands[which(CpG_islands$chr %in% sample$chr),]

#unique(CpG_islands$chr)
#unique(sample$chr)

#order islands and sample the same
CpG_islands <- CpG_islands[order(CpG_islands$chr,CpG_islands$start),]
sample <- sample[order(sample$chr,sample$pos),]

sample_chrom_ends <- NULL
j = 1
#Find largest end for each chromosome
for (k in unique(sample$chr))
{
        sample_chrom_ends[j] <- max(sample$pos[which(sample$chr == k)])
        j = j + 1
}

#add all previous values
sample_chrom_ends <- cumsum(as.numeric(sample_chrom_ends))
#Adjust vector to lengthen to continuous list
sample_chrom_ends <- c(0,sample_chrom_ends)
sample_chrom_ends <- sample_chrom_ends[-length(sample_chrom_ends)]

j = 1
#chain coordinates together to continuous list
for (k in unique(CpG_islands$chr))
{
        CpG_islands$start[which(CpG_islands$chr == k)] <- CpG_islands$start[which(CpG_islands$chr == k)] + sample_chrom_ends[j]
        CpG_islands$end[which(CpG_islands$chr == k)] <- CpG_islands$end[which(CpG_islands$chr == k)] + sample_chrom_ends[j]
        sample$pos[which(sample$chr == k)] <- sample$pos[which(sample$chr == k)] + sample_chrom_ends[j]
        j = j + 1
}

#check that we have a continuous list.
#if true then it is continuous
#all(sample$pos == cummax(sample$pos))

#find local maxima
#local <- c(sample$pos[1],sample$pos) - c(sample$pos,sample$pos[length(sample$pos)])
#If positives exist then the list is not continuous. this also gives the positions of said positive values
#positives <- which(local > 0)

#execute rcpp code to get index of cpgi positions for this sample
idx <- CpG_index(sample$pos,CpG_islands$start,CpG_islands$end)

#subset
cpgi <- sample$coverage[as.logical(idx)]
COV <- subset(cpgi, cpgi <= 200)

#Move to correct wd, just in case
setwd(wd)

#Produce COVERAGE histogram with fixed Y axis
cat("Plotting", name, "\n")

png(paste(name,'.Coverage.png',sep=''),width=1680,height=1050)
        hist(COV ,main=paste("Coverage density of ",name," in regions from ", region_name,sep=""),xlab="Coverage",ylab="Density",
                col="olivedrab3",labels=TRUE,freq=FALSE,breaks=200,xlim=c(1,60),ylim=c(0,0.10),las=1)
axis(side=1,at=seq(0,60,5))

dev.off()





 