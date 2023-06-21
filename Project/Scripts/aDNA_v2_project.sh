#!/bin/bash
 
############## NOTICE ##############
# Currently under development      #
# Patrick Cichowicz                #
#                                  #
# aDNA analysis automation script  #
####################################

 # Source Colors for text for better visibilty
source /Users/Patrick/aDNA/temp/colors.txt

 ########################################################################
 # Usage function - how to use script                                   #
 ########################################################################

# Still in development #
# plan to indicate whether SE/PE data, etc
set -e

# Usage() {
#    # Display Help
#    echo
#    echo -e "Usage" 
  
#    echo -e "\t$0 [options] -n <Project Name> -s <single-end or paired-end reads>" 
#    echo
#    echo -e " ${Blue}----- Required ----- ${Color_Off}"
#    echo
#    echo -e "\t -n |\t\t\t ${Green}Name of the project ${Color_Off}";
#    echo -e "\t -s |\t\t\t ${Green}Single-end or paired-end reads for bwa samse/sampe alignment ${Color_Off}"
#    echo -e "\t\t <SE || se || single || single-end>\t\t ${Green}Single-end option ${Color_Off}"
#    echo -e "\t\t <PE || pe || paired || paired-end>\t\t ${Green}Paired-end option ${Color_Off}"
# }
# if [ "$#" -eq 0 ]
# then
#     Usage
#     exit
# fi

# Name=""

# while getopts ":hn:" option; do
#   case $option in 
#     h) #Display help
#       Usage
#       exit;;
    
#     n) Name=$OPTARG;;
    
#     \?) #Incorrect/no option
#       echo -e "\tPlease use -h for help or provide the arguments needed"
#       Usage
#       exit;;
#   esac
# done

  ########################################################################
  # Functions 
  ########################################################################

# Here are functions that will be used frequently during the script and to avoid
# writing the same command. Also idea to print a text argument that can be easily
# grep'd 

# --------------------------------------------- IN DEVELOPMENT -------------------------------------------- #

# general function that exits after printing its text argument

err() {
  echo "$1...exiting";
  exit 1; # any non-0 exit code signals an error
}
# function to check return code of programs and exits with standard message if code is non-zero;
# otherwise displays completiong message and date.
#   arg 1 is the return code (usually $?)
#   arg2 is text describing what command ran
check_res() {
  if [ "$1" == "0" ]; then
    echo "..Done $2 `date`";
  else
    err "$2 returned non-0 exit code $1";
  fi
}
# --------------------------------------------- IN DEVELOPMENT -------------------------------------------- #

# Echo function that displays text in red indicating missing file or error
echo_r() {
  echo -e "${Red}$1${Color_Off}"
}

# Echo function that displays text in green indicating file exists or successful command
echo_g() {
  echo -e "${Green}$1${Color_Off}"
}

# Echo the time and date of command in yellow text
time_y() {
  echo -e "${BIYellow}`Date` ${BIBlue}$1${Color_Off}"
}

# Echo the text in blue
echo_b() {
  echo -e "${BIBlue}$1${Color_Off}"
}
  ########################################################################
  # Reference genome file and path location
  # For other reference genomes, create new variable with path of new reference file
  ########################################################################

reference="/Users/Patrick/aDNA/REFERENCE/UCSC/Homo/genome.fa"
refpath="/Users/Patrick/aDNA/REFERENCE/UCSC/Homo"  

time_y "Starting aDNA pipeline script automation"
sleep 1

# Get user to enter project name inorder to create dir along with subdirs for output files
# and log files.
read -p "$(echo "Enter Project name:" | pr -to30)" PROJECT
sleep 1

 # Making required directories specified from the user.
time_y "Making required input/output and log directories ..."

mkdir -p /Users/Patrick/aDNA/Project/${PROJECT}/{fastp,bwa,bam,final,logs,statistics}
find $PROJECT -print | pr -to30

# Checking for required files and directories before the pipeline starts
time_y "Checking for reference genome directory and files ..."
sleep 1
# list of bwa reference index files created
index_files=("amb" "ann" "bwt" "pac" "sa")
missing_index=()

# Check if the ref path/directory is correct/exists, along with all the necessary index files
# Also checks if the index files themselves are not empty (0), using the -s condition
if [ -d "$refpath" ]
then
  echo_g "Reference directory exists" | pr -to30

  for index in ${index_files[@]}; do 
    # Iterate through index of bwa index files required
    name="${reference}.${index}"

    if [[ (! -f "$name") || (! -s "$name") ]]
    then
      echo_r "Index $index file is missing. Please fix" | pr -to30
      missing_index+=("$name")
    else 
      echo_g "Reference genome index file $index file exist." | pr -to30
    fi
  
  done 
  # Prints all missing or empty reference genome index files in order for user to know which
  # files need to be fixed.
  if [ ${#missing_index[@]} -gt 0 ]
  then
    echo_r "${#missing_index[@]} index files missing, Aborting ..." | pr -to30
    exit 1 # exits with 1 error code indicating script is ending
  fi
  
else
  echo_r "Reference directory does not exist. Please confirm path or directory exists. Aborting .."
  exit 1 # exits with 1 error code indicating script is ending

fi | tee /Users/Patrick/aDNA/Project/${PROJECT}/logs/${PROJECT}_index.log # Writes all outcomes into log file

  # Making sure fastq files are available to be processed

  # path to directory where all project fastq files should be downloaded to
fastq_path="/Users/Patrick/aDNA/user_reduction"

read -p "$(echo "Enter name to project ${PROJECT}'s fastq files:" | pr -to30)" fastqs # consider making this better

echo_b "Checking for fastq files ..." | pr -to30
 
  # Check if the given path for fastq files exists
if [ ! -d "${fastq_path}/$fastqs" ]; then
  echo_r "$PROJECT fastq directory not found."
  exit 1
#else
# echo Sample DIR found ...
fi

  # Make a list of the fastq files 
  # Also assuming that paired end reads, check if both R1 and R2 files exist for each sample
fastq_files=($(find ${fastq_path}/${fastqs} -maxdepth 1 -type f -name "*.fastq.gz" | cut -f1-7 -d_ | sort | uniq))

  # list for samples which their R1 or R2 reads are missing
missing_fastq=()

sleep 1

# Process each Fastq file
for fastq_file in ${fastq_files[@]}; do
  filename=$(basename "$fastq_file")
  #sample_name="${filename%_*}"
  sleep 1

  #echo $fastq_file
  #echo $filename
  # Check if both R1 and R2 files exist for the current sample
  
  r1_file="${fastq_file}_S1_L004_R1_001.fastq.gz"
  r2_file="${fastq_file}_S1_L004_R2_001.fastq.gz"
  # r1_file="${fastq_path}/${fastqs}/${sample_name}_S1_L004_R1_001.fastq.gz"
  # r2_file="${fastq_path}/${fastqs}/${sample_name}_S1_L004_R2_001.fastq.gz"

  #echo $r1_file
  #echo $r2_file
  
  if [[ (! -s "$r1_file") || (! -s "$r2_file") ]]; then # check for both R1 and R2 files & if not empty files
    missing_fastq+=("$fastq_file")
    echo_r "Paired end reads for sample $filename missing" | pr -to30
  else
  
  #echo_g "R1 Fastq file for sample '$sample_name' found."
  #echo_g "R2 Fastq file for sample '$sample_name' found."
    echo_g "Paired end reads for sample $filename found" | pr -to30
  fi

done

  # Check to see how many samples are missing/empty and record samples not processed into log file
if [[ ${#missing_fastq[@]} -gt 0 ]]; then
  echo_r "${#missing_fastq[@]} Sample(s) not processed and found in log directory." | pr -to30
  echo "${missing_fastq}" > /Users/Patrick/aDNA/Project/${PROJECT}/logs/${PROJECT}_missing_fastqs.log | pr -to30

else
  echo_g "All paired end reads found."

fi

# echo -e "$Color_Off Number of fastq files to be processed : $(find ${rawfq}$fastqs -name "*fastq.gz" | wc -l)" | pr -to31

  #ADD IN CHECK FOR PAIRED READS

 

echo_b "All directories and files are ready ...  starting aDNA pipeline analysis" | pr -to30
echo

time_y "fastp - Adapter trimming & length/quality control"  

          # --------------------------------------------------- #
          # ---------------------- fastp ---------------------- #
          # --------------------------------------------------- #

# paths for in and output files required for fastp trimming and quality control
# fastq files are stored in the fastq_files variable above
fastp_out="/Users/Patrick/aDNA/Project/$PROJECT/fastp"

# Statistic files, determine number of reads before and after trimming

# fastq files to be trimmed & get the read lengths
# number of reads before and after trimming
tasks=${#fastq_files[@]}
task_count_fastp=0

# # Creating fastq_statistics text file with headers
 echo Sample_Name Total_Reads Trimmed_Reads Read_Length >> 1_fastq_statistics.txt

for fastq in ${fastq_files[@]}; do

  fastp_base=$(basename "$fastq_file")

  r1_file="${fastq}_S1_L004_R1_001.fastq.gz"
  r2_file="${fastq}_S1_L004_R2_001.fastq.gz"
  
  r1_out="${fastp_out}/${fastp_base}_S1_L004_R1_001.trim.fastq.gz"
  r2_out="${fastp_out}/${fastp_base}_S1_L004_R2_001.trim.fastq.gz"

  ((task_count_fastp++))
  basename=${fastq%.fastq.gz}

  echo "Sample $fastp_base - Sample number $task_count_fastp out of $tasks ..."

# Calculate total number of reads by - total number of lines / 4
  echo Calculating number of reads total ...

  #echo $(($(gzcat $fastq | wc -l ) /4))
  #total_reads=$(($(gzcat $fastq | wc -l ) /4))
  echo $fastq $total_reads >> /Users/Patrick/aDNA/Project/${PROJECT}/statistics/statistics_fastp.txt # writes to statistic file


  #echo Performing fastp adapter and length trimming ...
  # length and quality scores >= 30
  echo
  echo "fastp -l 30 -q 30 --in1 $r1_file --in2 $r2_file --out1 $r1_out --out2 $r2_out" # can use paralell command here
  echo
  sleep 2
  echo Calculating number of reads after fastp trimming ...
  #trimmed_reads=$(($(gzcat ${basename}_trim.fastq.gz | wc -l ) /4))

# Read length of trimmed fastq files
  echo Calculating read length of trimmed fastq files ...
  #read_length=$(gzcat ${basename}_trim.fastq.gz | awk "NR%4==2 {sum+=length($1)} END {print sum/ (NR/4)}")
  echo

 # rm *.json *.html
# Write all statistics from fastq files to fastq_statistics.txt file.
 # echo $basename | cut -f1 -d_ $total_reads $trimmed_reads $read_length >> 1_fastq_statistics.txt
done
#
          # ----------------------------------------------------- #
          # ---------------------- bwa aln ---------------------- #
          # ----------------------------------------------------- #

# trimmed fastq files align to reference genome using bwa aln
# samtools for flagstat report / sort by mapped reads only, index and convert to bam files
# mapdamage for fragmentation and 5' 3' end DNA damage patterns

# trimmed fastq files to be aligned to GRCh37 Reference Genome using bwa aln and samse
# task_count_aln=0
# tasks_bwa=$(ls *_trim* | wc -l)
# #

#
# # Creat bwa_statistics file
# echo Sample_Name Mapped_reads Mapped_length >> 2_bwa_statistics.txt
# #
# for trim_fastq in $(ls *_trim.fastq.gz);do
#   #
#    ((task_count_aln++))
#    basename=${trim_fastq%_trim.fastq.gz}
#
#   echo "Sample number $task_count_aln out of $tasks_bwa ..."
#   echo $basename
#   #
#   #
#   # SM=$(echo $trim_fastq | cut -f1 -d_)
#   # CN=CGG
#   # PL=ILLUMINA
#   # LB=$(echo $trim_fastq | cut -f4 -d_)
#   # ID=$(echo $trim_fastq | cut -f1 -d.)
#   #
#   # bwa aln $reference $trim_fastq > ${basename}_aln.sai
#   #
#   # bwa samse -r "@RG\tID:${ID}\tSM:${SM}\tCN:CGG\tPL:ILLUMINA\tLB:${LB}" \
#   # $reference ${basename}_aln.sai $trim_fastq | samtools sort -o ${basename}.bam
#   # # # Discard all unmapped read and keep mapped reads only.
#   #
#   # echo "Discarding all unmapped reads ..."
#   # samtools view -b -F 4 ${basename}.bam  > ${basename}_mapped.bam
#
#   echo "Calculating number of mapped reads ..."
#   mapped_reads=$(samtools view -c ${basename}_mapped.bam)
#
#   echo "Calculating mean read length of mapped reads ..."
#   mapped_length=$(samtools view ${basename}_mapped.bam | awk '{sum+=length($10)} END {print sum/NR}')
#
#   echo $basename $mapped_reads $mapped_length >> 2_bwa_statistics.txt
# done

          # ---------------------------------------------------- #
          # ---------------------- Picard ---------------------- #
          # ---------------------------------------------------- #

# Marking all duplicate reads using the Picard tool

# task_count_md=0
# tasks_md=$(ls *_mapped.bam | wc -l)
#
# # Create Picard MarkDuplicate statisitc file
# echo Sample_Name Duplicate_Reads Optical_Reads Unique_Reads >> 3_picard_statistics.txt
#
# for mapped in $(ls *_mapped.bam); do
#   ((task_count_md++))
#
#   basename=${mapped%_mapped.bam}
#
#   echo "Sample number $task_count_md out of $tasks_md ..."
#   echo "Marking all duplicates in sample $basename ..."
#
#   #java -jar /usr/local/Cellar/picard/picard.jar
#
#   java -jar $PICARD MarkDuplicates \
#   OPTICAL_DUPLICATE_PIXEL_DISTANCE=12000 \
#   I=${mapped} \
#   O=`basename ${mapped} _mapped.bam`_md.bam \
#   REMOVE_DUPLICATES=false \
#   METRICS_FILE=`basename ${mapped} _mapped.bam`_md.bam.metrics \
#   TAGGING_POLICY=All \
#   VALIDATION_STRINGENCY=LENIENT
#
#       #WHY NOT USE THE SAMTOOLS VIEW -f and -F OPTIONS TO GET THIS INFO...... COULD BE A POSSIBILITY
#       #UNABLE TO POTENTIALLY GET THE OPTICAL READS
#
# # For read pairs, need to add col 6 + 7 for duplicates, then col 8 for optical only
# # For unique, most likely need to add up col 2 and 3 for total number of reads
#   #dupread=$(awk 'NR==8 {print $6}' ${basename}_md.bam.metrics)
#   dupread=$(awk 'NR==8, sum = $6 + $7 {print sum}' ${basename}_md.bam.metrics)
#   optreads=$(awk 'NR==8 {print $8}' ${basename}_md.bam.metrics)
#   uniqueread=$(awk 'NR==8 {sum = $2 - $6 ; print sum}' ${basename}_md.bam.metrics)
#
#   echo $basename $dupread $optreads $uniqueread >> 3_picard_statistics.txt
# done

# From here on out, need to use -F 1024 option to remove all the duplicated reads for further analysis


          # ---------------------------------------------------- #
          # ---------------------- mapDamage ------------------- #
          # ---------------------------------------------------- #

# for mapdam in $(ls *rmd.bam); do
#   echo "mapDamage - Sample $mapdam"
#
#   mapDamage -i $mapdam -r $reference -d results
#
# done





# Take statisitc files and merge them into one main, skip the first column so no duplicates
#paste fastq_statistics.txt <(cut -d " " -f2- bwa_statistics.txt) | column -t
#
  # echo "Merging statistic files ..."
  #
  # arr=([1-9]*statistics.txt)
  # file="${arr[0]}"
  #
  # for f in "${arr[@]:1}"; do
  #   paste "$file" <(cut -d " " -f2- $f) > _file.tmp && mv _file.tmp file.tmp
  #   file=file.tmp
  #
  #
  # done











echo -e "$BIYellow $(Date) $BIBlue aDNA pipeline automation completed ${Color_Off}"
echo -e "$BIYellow $(Date) ${Color_Off} Maybe some logs/errors message here, total time taken to run scripts ... etc"

# Picard metrics file readme

# 1 - LIBRARY
# 2 - UNPAIRED_READS_EXAMINED
# 3 - READ_PAIRS_EXAMINED
# 4 - SECONDARY_OR_SUPPLEMENTARY_RDS
# 5 - UNMAPPED_READS
# 6 - UNPAIRED_READ_DUPLICATES
# 7 - READ_PAIR_DUPLICATES
# 8 - READ_PAIR_OPTICAL_DUPLICATES
# 9 - PERCENT_DUPLICATION
# 10 - ESTIMATED_LIBRARY_SIZE


# Statistics File Columns Readme
# 1 - Sample name
# 2 - Total number of reads
# 3 - Total number of reads after fastp trimming
# 4 - Mean read length of trimmed reads (Column 3)
# 5 - Number of reads mapped to reference genome
# 6 - Mean read length of mapped reads only
# 7 - Number of duplicate reads (Picard)
# 8 - Number of unique mapped reads (excluding duplicate reads)
