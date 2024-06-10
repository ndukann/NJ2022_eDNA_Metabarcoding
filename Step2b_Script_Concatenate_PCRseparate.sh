#!/bin/bash
#Written by Joran Vanhollebeke & Isolde Cornelis
#Date: 

clear -x

#Input is asked
echo "Enter your input folder: (vb: /home/genomics/sderycke/Minion_run/basecalled_folder)"
read -e INPUTFOLDER
echo

cd $INPUTFOLDER
echo $INPUTFOLDER
echo

#Input is asked
echo "Enter your RUNNAME: (vb: MiFish_UE-S)"
read -e RUN_NAME
echo

#RUN_NAME=MiFish_UE-S
echo $RUN_NAME

# Define the variable GZFILES
echo find ${INPUTFOLDER}/${RUN_NAME}

#possible that you have to change this from fastq to gz depending if the files are already unzipped
GZFILES=$(find $INPUTFOLDER/${RUN_NAME}*/processed-reads/*/trimmed-R*/ -type f -name "*.gz")

echo Contains >> log.txt
echo $GZFILES >> log.txt

#Unzip all files
gunzip $GZFILES

echo $PWD
FILES=$(ls ${INPUTFOLDER}/${RUN_NAME}1/processed-reads/sense/trimmed-R1/*.fastq)
echo $FILES

mkdir -p $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/sense/trimmed-R1/
mkdir -p $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/antisense/trimmed-R1/
mkdir -p $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/sense/trimmed-R2/
mkdir -p $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/antisense/trimmed-R2/

#Copy and zip the files of all three PCR replicates into one concatenated folder

cp $INPUTFOLDER/${RUN_NAME}*/processed-reads/sense/trimmed-R1/*_S*.R1.fastq $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/sense/trimmed-R1
cp $INPUTFOLDER/${RUN_NAME}*/processed-reads/sense/trimmed-R2/*_S*.R2.fastq $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/sense/trimmed-R2
cp $INPUTFOLDER/${RUN_NAME}*/processed-reads/antisense/trimmed-R1/*_S*.R1.fastq $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/antisense/trimmed-R1
cp $INPUTFOLDER/${RUN_NAME}*/processed-reads/antisense/trimmed-R2/*_S*.R2.fastq $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/antisense/trimmed-R2

gzip $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/*/trimmed-R*/*.R*.fastq

echo "finished"
