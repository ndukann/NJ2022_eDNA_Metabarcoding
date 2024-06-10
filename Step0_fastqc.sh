#!/bin/bash
################################################################################
# Author	: Joran Vanhollebeke 
# Usage		: 


# conda activate multiqc
# /home/genomics/jvanhollebeke/01_Basic_detailed_scripts/Step0_fastqc.sh /home/genomics/smaes/genohub####/Raw_data/ /home/genomics/smaes/genohub#####/fastqc/ 4
# Info		: Don't forget to check the extentions of your file (fq.gz or fastq.gz) and adapt the script appropriately  
################################################################################
# VALIDATE INPUT
################################################################################
function usage(){
	errorString="Running this FastQC script requires 3 parameters:\n
		1. Path of the folder with fastq.gz files.\n
		2. Path of the output folder.\n
		3. Number of threads to use.\n";
	echo -e ${errorString};
	exit 1;
}
if [ "$#" -ne 3 ]; then
	usage
fi
################################################################################
# INPUT FOLDER CONTAINING .fastq.gz files
################################################################################
inputFolder=$1;
# Remove trailing slash if this is last char
len=${#inputFolder};
lastPos=$(expr $len - 1);
lastChar=${inputFolder:$lastPos:1};
if [[ $lastChar == '/' ]]; then
	inputFolder=${inputFolder:0:$lastPos};
fi
################################################################################
# CREATE FOLDER (IF NOT EXISTS)
################################################################################
outFolder=$2;
mkdir -p ${outFolder}
################################################################################
# RUN FASTQC # $1 = path fastq files, $2 = output folder, $3 = threads
################################################################################
# Concatenate filenames
inputFiles="";
for i in $( ls ${inputFolder}/ | grep fastq.gz); do
	inputFile="${inputFolder}/$i";
    # Add filename
    inputFiles="$inputFiles $inputFile";
done
echo "### ${inputFiles} ###";
# Compose command
fastqcCommand="fastqc --extract -t $3 -o ${outFolder} ${inputFiles}";
# Show command
echo -e "$fastqcCommand";
# Execute
output=$(eval $fastqcCommand);
# Show output
echo -e "$output\n";
#################################################################################

# Automatically execute a multiqc on the made quality files
echo "cd $outFolder";
echo "multiqc .";
