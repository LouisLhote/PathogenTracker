#!/bin/bash
#this script is a test pipeline for pathogen screening 


###############################################################################################fonction#####################################################################################################
###########################################################################################################################################################################################################

Help()
{
    # Display Help
    echo "Welcome to PathogenTracker, a pipeline to track pathogen in animal sample"
    echo
    echo  "PathogenTracker v0.1"
    echo  "OPTIONS:"
    echo
    echo
    echo "USAGE bash PathogenTrackerv0.1.sh [input-file] [output-directory] [mode]"
    echo
    echo "-h     Print this Help."
    echo "MODE:"
    echo
    echo "-t     Trimming mode, use it if you have R1 and R2 not trimmed file"
    echo "-c     Collapsed mode, use it if you have alreay trimmed and collapsed reads"
    echo "File options:"
    echo
    echo "-1      R1 file, USE ONLY in trimming mode"
    echo "-2      R2 file, USE ONLY in trimming mode"
    echo "-m      Collapsed file USE ONLY in collapsed mode"
    echo
    echo "-o      MANDATORY, careful, is your path is  /my/output/directory/ write it /my/output/directory  " 
}




########################################################################classic mode with adapter removal#############################################################################################

trimming-mode()
{


echo "#################################################################################################################################################################################################"
echo "#################################################################################################################################################################################################"
echo "###################################################################WELCOME TO PathogenTracker####################################################################################################"
echo "#################################################################################################################################################################################################"
echo "#################################################################################################################################################################################################"
echo
echo
echo 


if [ -n "${s1}" ]
then
  echo "The R1 sequence is ${s1}"
else
  echo "######Please enter a R1 file!!!#########"
  exit
fi

if [ -n "${s2}" ]
then
 echo "The R2 sequence is ${s2}"
else
  echo "######Please enter a R2 file!!!#########"
  exit
fi

if [ -n "${out_dir}" ]
then
 echo "The results directory is ${out_dir}"
else
  echo "######Please enter a output directory file!!!#########"
  exit
fi


surname=`basename ${s1}`
name=`echo ${surname} | cut -f1 -d'_' `
echo ${name}




echo "================================================================Trimming step: Filtering with AdapterRemoval==================================================="

#########QF############
###use adapterremoval
####warning : adapterremoval is install in my personnal conda, if you want to use it, please install it : https://adapterremoval.readthedocs.io/en/stable/



mkdir $out_dir/AdapterRemoval_output

AdapterRemoval --file1 ${s1} --file2 ${s2} --basename $out_dir/AdapterRemoval_output/${name} --threads 6 --collapse --minadapteroverlap 1 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minlength 30 --gzip --trimns --trimqualities


echo "================================================================The trimming step is done!==================================================="


cd $out_dir
cd AdapterRemoval_output
mkdir collapsed_reads
mv *collapsed.gz collapsed_reads/
cd collapsed_reads



echo "================================================================Duplicate Removal step: Removing the dupliacte with prinseq++==================================================="
#######Duplicate removal########

mkdir prinseq_out_good
mkdir prinseq_out_bad

eval "$(conda shell.bash hook)"

conda activate prinseq-plus-plus
prinseq++ -fastq ${name}.collapsed.gz -derep 14 -out_good prinseq_out_good/${name}-passed.fq -out_bad prinseq_out_bad/${name}-not_passed.fq -VERBOSE 2 -threads 6


echo "================================================================Duplicate Removal step is done==================================================="
conda deactivate

#########alignement on host genome###########


cd prinseq_out_good
mkdir bam-file
cd bam-file
mkdir aligned
mkdir unaligned
cd ..


echo "================================================================Alignement on host genome step: The reads are aligned with bowtie2==================================================="

echo "========First step, the reads are aligned on a mix of differents host genome======"
###mapped reads
bowtie2 -x /raid_md0/Reference_Genomes/multi/chimera -U ${name}-passed.fq -p 12 | samtools view -Sb  - -F4 >  bam-file/aligned/${name}_aligned.bam
echo "========Second step, the unaligned reads are extract from the alignement."
###unmapped reads
bowtie2 -x /raid_md0/Reference_Genomes/multi/chimera -U ${name}-passed.fq -p 12 | samtools view -Sb  - -f4 >  bam-file/unaligned/${name}_unaligned.bam

####bam to fastq
cd bam-file
cd unaligned
mkdir fastq




samtools bam2fq ${name}_unaligned.bam -@4 |gzip - > fastq/${name}_unaligned.fastq.gz

echo "================================================================Alignement step is done==================================================="

echo "================================================================Taxonomic assignement step: using Kraken2==================================================="

####kraken2 assignement
cd fastq
mkdir kraken-report
mkdir kraken-output
kraken2 --db /raid_md0/metagenomic_resources/kraken_db_custom/ncbint/ ${name}_unaligned.fastq.gz --report-minimizer-data --minimum-hit-groups 3 --output kraken-output/${name}-unmapped-krak-output --report kraken-report/${name}-unmapped-krak-report.txt --threads 20
echo "================================================================Taxonomic assignement step is done==================================================="

echo "This is the end of the pipeline, the final file is the kraken-report file, you can vizualise it with pavian: https://github.com/fbreitwieser/pavian   "
echo
echo
echo "All the kraken-report are here ${output_dir}/AdapterRemoval_output/prinseq_out_good/bam-file/unaligned/fastq/kraken-report/ "
exit
}


###################################################################################collapsed mode,without the adapterremoval step###########################################################################

collapsed-mode()
{
echo "#################################################################################################################################################################################################"
echo "#################################################################################################################################################################################################"
echo "###################################################################WELCOME TO PathogenTracker####################################################################################################"
echo "#################################################################################################################################################################################################"
echo "#################################################################################################################################################################################################"
echo
echo
echo



if [ -n "${merged}" ]
then
  echo "The collapsed sequence is ${merged}"
else
  echo "######Please enter a collapsed file!!!#########"
  exit
fi

if [ -n "${out_dir}" ]
then
 echo "The results directory is ${out_dir}"
else
  echo "######Please enter a output directory file!!!#########"
  exit
fi




surname=`basename ${merged}`
name=`echo ${surname} | cut -f1 -d'.' `
echo ${name}


echo "================================================================Duplicate Removal step: Removing the dupliacte with prinseq++==================================================="
#######Duplicate removal########

cd $out_dir
mkdir $out_dir/prinseq_out_good
mkdir $out_dir/prinseq_out_bad

eval "$(conda shell.bash hook)"

conda activate prinseq-plus-plus
prinseq++ -fastq ${merged} -derep 14 -out_good $out_dir/prinseq_out_good/${name}-passed.fq -out_bad $out_dir/prinseq_out_bad/${name}-not_passed.fq -VERBOSE 2 -threads 6


echo "================================================================Duplicate Removal step is done==================================================="
conda deactivate

#########alignement on host genome###########

cd $out_dir/prinseq_out_good
mkdir bam-file
cd bam-file
mkdir aligned
mkdir unaligned
cd ..

echo "================================================================Alignement on host genome step: The reads are aligned with bowtie2==================================================="

echo "========First step, the reads are aligned on a mix of differents host genome======"
###mapped reads
bowtie2 -x /raid_md0/Reference_Genomes/multi/chimera -U ${name}-passed.fq -p 12 | samtools view -Sb  - -F4 >  bam-file/aligned/${name}_aligned.bam
echo "========Second step, the unaligned reads are extract from the alignement."
###unmapped reads
bowtie2 -x /raid_md0/Reference_Genomes/multi/chimera -U ${name}-passed.fq -p 12 | samtools view -Sb  - -f4 >  bam-file/unaligned/${name}_unaligned.bam

####bam to fastq
cd bam-file
cd unaligned
mkdir fastq




samtools bam2fq ${name}_unaligned.bam -@4 |gzip - > fastq/${name}_unaligned.fastq.gz

echo "================================================================Alignement step is done==================================================="

echo "================================================================Taxonomic assignement step: using Kraken2==================================================="

####kraken2 assignement
cd fastq
mkdir kraken-report
mkdir kraken-output
kraken2 --db /raid_md0/metagenomic_resources/kraken_db_custom/ncbint/ ${name}_unaligned.fastq.gz --report-minimizer-data --minimum-hit-groups 3 --output kraken-output/${name}-unmapped-krak-output --report kraken-report/${name}-unmapped-krak-report.txt --threads 20
echo "================================================================Taxonomic assignement step is done==================================================="

echo "This is the end of the pipeline, the final file is the kraken-report file, you can vizualise it with pavian: https://github.com/fbreitwieser/pavian   "
echo
echo
echo "All the kraken-report are here ${output_dir}/prinseq_out_good/bam-file/unaligned/fastq/kraken-report/ "
exit

}



###########################################################################input#########################################################################################################################
#########################################################################################################################################################################################################
while getopts ":o:1:2:m: h c t" option; do
     case $option in
          h) # display Help
             Help
             exit;;
          1) # r1
             s1=${OPTARG}
             ;;
          2) # r2
             s2=${OPTARG}
             ;;
          m) # merged
             merged=${OPTARG}
             ;;
          o) # out dir
             out_dir=${OPTARG}
             ;;
          t) # trimming mode
             trimming-mode
            ;;
          c) # collapsed mode
             collapsed-mode
            ;;
          \?) # incorrect option
             echo "Error: Invalid option"
            exit;;
     esac
done



