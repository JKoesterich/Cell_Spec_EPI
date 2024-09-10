#!/bin/bash

#SBATCH --partition=mem       # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --cpus-per-task=8           # Cores per task (>1 if multithread tasks)
#SBATCH --mem=64G                   # Real memory (RAM) required (MB)
######SBATCH --time=3-00:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=ATCRCP_PipelineSE_%j_STDOUT.out     # STDOUT output file
#SBATCH --error=ATCRCP_PipelineSE_%j_STDERR.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env


# Add command line inputs, first is the input file name (entire path), second is the output folder (a folder with the name of the file (minus the pathing and the .fastq.gz) will be created inside this folder)
# Since most steps are similar, I am converting this to run the Cut&Run and ChIPseq options as well, command line input 3 will be either ATAC, CUTRUN, or CHIP
# Example code ATCRCP_Generalized_SingleEnded_pipeline.sh /path/to/SEfile /pathTo/directoryThat/willHold/OutputFolder ATAC
#if [[ $3 != "ATAC" && $3 != "CUTRUN" && $3 != "CHIP" ]]
# Commenting out the above for now as the last step (SEARC only allows paired end data)
if [[ $3 != "ATAC" && $3 != "CHIP" ]]
then
	echo "Input 3 is invalid: currently only ATAC and CHIP are accepted"
	exit
fi
cd /scratch/jhk148
source miniconda3/bin/deactivate
source miniconda3/bin/activate miniconda3/envs/macs-py2.7
INFIL=$1
FILFUL=${INFIL##*/}
FIL=${FILFUL::-9}
OFOL="$2"/"$FIL"
mkdir $OFOL
mkdir "$OFOL"/FASTQC
module load java
java -jar /projects/community/Trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar  SE -phred33 -threads 8 $INFIL "$OFOL"/"$FIL"_Trmd.fq.gz ILLUMINACLIP:Reference_Data/NGS/hg19/adapters/TruSeq3-SE.fa:2:30:10 ILLUMINACLIP:Reference_Data/NGS/hg19/NEBNextPrimers1and2.fa:2:30:10 LEADING:15 TRAILING:15 MINLEN:36 CROP:50
module load FastQC
fastqc  "$OFOL"/"$FIL"_Trmd.fq.gz  -o "$OFOL"/FASTQC
module load bowtie2
if [[ $3 == "ATAC" ]]
then
	bowtie2  -x Bowtie2Index/genome -U "$OFOL"/"$FIL"_Trmd.fq.gz -S "$OFOL"/"$FIL"_Aligned.bam -X 2000 -p 8 --no-mixed --no-discordant 2> "$OFOL"/"$FIL"_bowtie_stats.log
elif [[ $3 == "CUTRUN" ]]
then
	bowtie2  -x Bowtie2Index/genome -U "$OFOL"/"$FIL"_Trmd.fq.gz -S "$OFOL"/"$FIL"_Aligned.bam -I 10 -X 2000 -p 8 --dovetail --phred33 --local --very-sensitive-local --no-unal --no-mixed --no-discordant 2> "$OFOL"/"$FIL"_bowtie_stats.log
elif [[ $3 == "CHIP" ]]
then
	bowtie2  -x Bowtie2Index/genome -U "$OFOL"/"$FIL"_Trmd.fq.gz -S "$OFOL"/"$FIL"_Aligned.bam -X 2000 -p 4 --no-mixed --no-discordant 2> "$OFOL"/"$FIL"_bowtie_stats.log
fi
module load java
/projects/community/GATK/gatk-4.1.4.1/gatk SortSam --INPUT="$OFOL"/"$FIL"_Aligned.bam --OUTPUT="$OFOL"/"$FIL"_sorted.bam -SO=coordinate  --TMP_DIR=/mnt/scratch
python PlotMAPQScoresInBam.py --bam "$OFOL"/"$FIL"_sorted.bam --out "$OFOL"/FASTQC --samplename $FIL
module load samtools
samtools view -bq 30 "$OFOL"/"$FIL"_sorted.bam > "$OFOL"/"$FIL"_filt.bam
#module load GATK/4.1.4.1-yc759
/projects/community/GATK/gatk-4.1.4.1/gatk --java-options '-Xmx56G' MarkDuplicates  --REMOVE_DUPLICATES true  --CREATE_INDEX true  --TMP_DIR /mnt/scratch --METRICS_FILE "$OFOL"/"$FIL"_dups.tab --INPUT "$OFOL"/"$FIL"_filt.bam --OUTPUT "$OFOL"/"$FIL"_dupchkd.bam --ASSUME_SORT_ORDER=coordinate  2>> "$OFOL"/"$FIL"_markdups.log
module load bedtools2
bedtools intersect -v -a "$OFOL"/"$FIL"_dupchkd.bam -b Reference_Data/NGS/hg19/hg19-blacklist-kundaje.bed > "$OFOL"/"$FIL"_noblklst.bam
if [[ $3 == "ATAC" ]]
then
	module load java
	/projects/community/GATK/gatk-4.1.4.1/gatk SortSam --INPUT="$OFOL"/"$FIL"_noblklst.bam --OUTPUT="$OFOL"/"$FIL"_shifted.bam -SO=coordinate  --TMP_DIR=/mnt/scratch
	module load samtools
	samtools index "$OFOL"/"$FIL"_shifted.bam
	module load java
	igvtools count "$OFOL"/"$FIL"_shifted.bam "$OFOL"/"$FIL"_noblklst.tdf hg19
	/usr/bin/perl count_dup.pl "$OFOL"/"$FIL"_noblklst.bam "$OFOL"/"$FIL"_NirDups.txt
	macs2 callpeak -t "$OFOL"/"$FIL"_shifted.bam --outdir $OFOL -f BAM -n "$FIL"_macs2 -g hs
elif [[ $3 == "CUTRUN" ]]
then
	samtools index "$OFOL"/"$FIL"_noblklst.bam
	igvtools count --pairs "$OFOL"/"$FIL"_noblklst.bam "$OFOL"/"$FIL"_noblklst.tdf hg19
	/usr/bin/perl count_dup.pl "$OFOL"/"$FIL"_noblklst.bam "$OFOL"/"$FIL"_NirDups.txt
	macs2 callpeak -t "$OFOL"/"$FIL"_shifted.bam --outdir $OFOL -f BAM -n "$FIL"_macs2 -g hs -q 0.05
elif [[ $3 == "CHIP" ]]
then
	igvtools count --pairs "$OFOL"/"$FIL"_noblklst.bam "$OFOL"/"$FIL"_noblklst.tdf hg19
	/usr/bin/perl count_dup.pl "$OFOL"/"$FIL"_noblklst.bam "$OFOL"/"$FIL"_NirDups.txt
	macs2 callpeak -t "$OFOL"/"$FIL"_noblklst.bam --outdir $OFOL -f BAM -n "$FIL"_macs2 -g hs
fi
