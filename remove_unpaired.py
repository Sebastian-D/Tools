#!/usr/bin/env python

import sys

#Author: Sebastian DiLorenzo
#Date: 2015-22-10

#What this script is for:
#If you for some reason have a SAM or BAM file with broken mate pairs (one mate missing) that is queryname sorted
#and the samtools flags are broken for said reads as well. Note that if the flag is not broken you can use samtools to the same effect as this script. 
# This script will remove the single reads to one file (singles.sam) and the paired reads to another file (paired.sam).
#For me this occurred in my bismark aligned paired end data, perhaps due to deduplicate_bismark but that is not confirmed.

#Usage BAM file:
#samtools view -h input.bam | ./remove_unpaired.py

#Usage SAM file:
#./remove_unpaired.py input.sam

#Take each line from the input file (which has been piped in as sam format using samtools) and
#remove lines that are unpaired.

#open two files for writing
output = open("paired.sam","w")
singles = open("singles.sam","w")

#Grab and write the header (@)
while True:
	line1 = sys.stdin.readline()
	if line1.startswith("H"):
		break
	else:
		output.write(line1)

#Check two lines against eachother, if they match write them as a pair and jump to the next two lines. If they dont write
#the first line as a single and iterate only one step looking at the next pair
while True:
	#read in second line first (we already have line1 from reading the header)
	line2 = sys.stdin.readline()
	#If line2 doesnt exist we are at EOF
	if not line2:
		break
	#Are the two querynames the same?
	if line1.split("\t")[0] == line2.split("\t")[0]:
		output.write(line1)
		output.write(line2)
		line1 = sys.stdin.readline()
	else:
		singles.write(line1)
		line1 = line2

#Close the files
singles.close()
output.close()
