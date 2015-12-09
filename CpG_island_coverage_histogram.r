#Description: Companion script for methylKit. Produce histograms from a methylRaw object.
#				Works with any intervals formatted as Hsapien_CpG_islands_hg19.txt is.
#
#Usage: Rscript CpG_island_coverage_histograms.r methylRawObject.Rdata Hsapien_CpG_islands_hg19.txt
#
#Author: Sebastian DiLorenzo
#
#Date: 19/11/2015
#
#Constraints: You will need to have installed the R package methylKit and generated a methylRAWobject (see ?read from methylKit)
#How to create the methylrawobject which must be named myobj:
#library(methylKit)
#myobj=read(<arguments>)
#save(myobj,file="methylRawObject.Rdata")

library(methylKit)
library(Rcpp)

#Input arguments
arg <- commandArgs(trailingOnly = TRUE)

#Load mehtylraw object
load(arg[1])

#read CpG island regions
CpG_islands <- read.table(arg[2],colClasses=c("character","integer","integer",header=T)
#colnames(CpG_islands) <- c("chr","start","end")

#Backup the original
original_CpG_island <- CpG_islands
#Revert
#CpG_islands <- original_CpG_island

#Safeguard wd so that plots are created in correct location
wd <- getwd()

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


#forloop over samples
for (i in 1:length(myobj))
{
        #Get a sample
        sample <- getData(myobj[[i]])
        name = myobj[[i]]@sample.id

        #Send message
        cat("Starting work on ",name,"\n")

        #Revert to original CpG_islands
        CpG_islands <- original_CpG_island

        #Convert factor to character
        sample$chr <- as.character(sample$chr)

        #Only look at chromosomes that have CpG islands according to the CpG_islands list. Remove all others from sample
        sample <- sample[which(sample$chr %in% CpG_islands$chr),]

        #Only look at chromosomes that exist in this sample. (this is why we reload the list of CpG_islands for each sample)
        CpG_islands <- CpG_islands[which(CpG_islands$chr %in% sample$chr),]

        #unique(CpG_islands$chr)
        #unique(sample$chr)

        #order islands and sample the same
        CpG_islands <- CpG_islands[order(CpG_islands$chr,CpG_islands$start),]
        sample <- sample[order(sample$chr,sample$start),]

        sample_chrom_ends <- NULL
        #CpG_chrom_ends <- NULL
        j = 1
        #Find largest end for each chromosome
        for (k in unique(sample$chr))
        {
                #CpG_chrom_ends[j] <- max(CpG_islands$end[which(CpG_islands$chr == i)])
                sample_chrom_ends[j] <- max(sample$end[which(sample$chr == k)])
                j = j + 1
                #cat(j,"\n")
        }

        #backups
        #backup_cpg_ends <- CpG_chrom_ends
        #backup_sample_ends <- sample_chrom_ends
        #backup_CpG_islands <- CpG_islands
        #backup_sample <- sample

        #revert
        #sample <- backup_sample
        #CpG_islands <- backup_CpG_islands


        #add all previous values
        sample_chrom_ends <- cumsum(as.numeric(sample_chrom_ends))
        #Adjust vector to lengthen to continuous list
        sample_chrom_ends <- c(0,sample_chrom_ends)
        sample_chrom_ends <- sample_chrom_ends[-length(sample_chrom_ends)]

        #CpG_chrom_ends <- cumsum(as.numeric(CpG_chrom_ends))
        #CpG_chrom_ends <- c(0,CpG_chrom_ends)
        #CpG_chrom_ends <- CpG_chrom_ends[-length(CpG_chrom_ends)]

        j = 1
        #chain coordinates together to continuous list
        for (k in unique(CpG_islands$chr))
        {
                CpG_islands$start[which(CpG_islands$chr == k)] <- CpG_islands$start[which(CpG_islands$chr == k)] + sample_chrom_ends[j]
                CpG_islands$end[which(CpG_islands$chr == k)] <- CpG_islands$end[which(CpG_islands$chr == k)] + sample_chrom_ends[j]
                sample$start[which(sample$chr == k)] <- sample$start[which(sample$chr == k)] + sample_chrom_ends[j]
                sample$end[which(sample$chr == k)] <- sample$end[which(sample$chr == k)] + sample_chrom_ends[j]
                j = j + 1
                #cat(i,"\n")
        }

        #check that we have a continuous list.
        #if true then it is continuous
        #all(sample$start == cummax(sample$start))

        #find local maxima
        #local <- c(sample$start[1],sample$start) - c(sample$start,sample$start[length(sample$start)])
        #If positives exist then the list is not continuous. this also gives the positions of said positive values
        #positives <- which(local > 0)

        #execute rcpp code to get index of cpgi positions for this sample
        idx <- CpG_index(sample$start,CpG_islands$start,CpG_islands$end)

	#subset
        cpgi <- sample$coverage[as.logical(idx)]
	COV <- subset(cpgi, cpgi <= 200)

        #Move to correct wd, just in case
        setwd(wd)

        #Produce COVERAGE histogram with fixed Y axis
        cat("Plotting", name, "\n")

	png(paste(name,'.Coverage.png',sep=''),width=1680,height=1050)
                hist(COV ,main=paste("Coverage density of ",name," in CpG islands",sep=""),xlab="Coverage",ylab="Density",
                        col="olivedrab3",labels=TRUE,freq=FALSE,breaks=200,xlim=c(1,60),ylim=c(0,0.10),las=1)
	axis(side=1,at=seq(0,60,5))

	dev.off()
}




 