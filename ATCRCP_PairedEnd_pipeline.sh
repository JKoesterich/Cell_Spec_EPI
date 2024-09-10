#!/bin/bash

#SBATCH --partition=mem       # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --cpus-per-task=8           # Cores per task (>1 if multithread tasks)
#SBATCH --mem=64G                   # Real memory (RAM) required (MB)
######SBATCH --time=3-00:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=ATCRCP_PipelinePE_%j_STDOUT.out     # STDOUT output file
#SBATCH --error=ATCRCP_PipelinePE_%j_STDERR.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env


# Add command line inputs, first is the R1 input file name (entire path), second is the R2 input file name (entire path), third is the output folder (a folder with the name of the file (minus the pathing and the .fastq.gz) will be created inside this folder)
# Since most steps are similar, I am converting this to run the Cut&Run and ChIPseq options as well, command line input 4 will be either ATAC, CUTRUN, or CHIP
# Example code ATCRCP_Generalized_SingleEnded_pipeline.sh /path/to/PEfile_R1 /path/to/PEfile_R2 /pathTo/directoryThat/willHold/OutputFolder ATAC
if [[ $4 != "ATAC" && $4 != "CUTRUN" && $4 != "CHIP" ]]
then
	echo "Input 4 is invalid: currently only ATAC, CHIP, and CUTRUN are accepted"
	exit
fi
source miniconda3/bin/deactivate
source miniconda3/bin/activate miniconda3/envs/macs-py2.7
INFIL=$1
FILFUL=${INFIL##*/}
FIL1=$(echo $FILFUL | awk '{gsub(".gz","",$0); gsub(".fastq","",$0);print($0)}')
INFIL2=$2
FILFUL2=${INFIL2##*/}
FIL2=$(echo $FILFUL2 | awk '{gsub(".gz","",$0); gsub(".fastq","",$0);print($0)}')
FIL=$(echo $FIL1 | awk '{gsub("_R1","",$0);print($0)}')
OFOL="$3"/"$FIL"
R1O="$OFOL"/"$FIL1"
R2O="$OFOL"/"$FIL2"
mkdir $OFOL
mkdir "$OFOL"/FASTQC
module load java
java -jar /projects/community/Trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 8 $INFIL $INFIL2 "$R1O"_Trmd.fq.gz "$R1O"_Trmd_Unpaired.fq.gz "$R2O"_Trmd.fq.gz "$R2O"_Trmd_Unpaired.fq.gz ILLUMINACLIP:Reference_Data/NGS/hg19/adapters/TruSeq3-PE.fa:2:15:4:4:true ILLUMINACLIP:Reference_Data/NGS/hg19/NEBNextPrimers1and2.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25

module load FastQC
fastqc  "$R1O"_Trmd.fq.gz  -o "$OFOL"/FASTQC
fastqc  "$R2O"_Trmd.fq.gz  -o "$OFOL"/FASTQC
module load bowtie2
if [[ $4 == "ATAC" ]]
then
	bowtie2  -x Bowtie2Index/genome -1 "$R1O"_Trmd.fq.gz -2 "$R2O"_Trmd.fq.gz -S "$OFOL"/"$FIL"_Aligned.bam -X 2000 -p 8 --no-mixed --no-discordant 2> "$OFOL"/"$FIL"_bowtie_stats.log
elif [[ $4 == "CUTRUN" ]]
then
#	bowtie2  -x Bowtie2Index/genome -1 "$R1O"_Trmd.fq.gz -2 "$R2O"_Trmd.fq.gz -S "$OFOL"/"$FIL"_Aligned.bam -I 10 -X 2000 -p 8 --dovetail --phred33 --local --very-sensitive-local --no-unal --no-mixed --no-discordant 2> "$OFOL"/"$FIL"_bowtie_stats.log
	bowtie2  -x Bowtie2Index/genome -1 "$R1O"_Trmd.fq.gz -2 "$R2O"_Trmd.fq.gz -S "$OFOL"/"$FIL"_Aligned.bam -I 10 -X 700 -p 8 --dovetail --phred33 --local --very-sensitive-local --no-unal --no-mixed --no-discordant 2> "$OFOL"/"$FIL"_bowtie_stats.log
elif [[ $4 == "CHIP" ]]
then
	bowtie2  -x Bowtie2Index/genome -1 "$R1O"_Trmd.fq.gz -2 "$R2O"_Trmd.fq.gz -S "$OFOL"/"$FIL"_Aligned.bam -X 2000 -p 4 --no-mixed --no-discordant 2> "$OFOL"/"$FIL"_bowtie_stats.log
fi
#module load java
#module load GATK/4.1.4.1-yc759
/projects/community/GATK/gatk-4.1.4.1/gatk SortSam --INPUT="$OFOL"/"$FIL"_Aligned.bam --OUTPUT="$OFOL"/"$FIL"_sorted.bam -SO=coordinate  --TMP_DIR=/mnt/scratch
python PlotMAPQScoresInBam.py --bam "$OFOL"/"$FIL"_sorted.bam --out "$OFOL"/FASTQC --samplename $FIL
module load samtools
samtools view -bq 30 "$OFOL"/"$FIL"_sorted.bam > "$OFOL"/"$FIL"_filt.bam
#module load GATK/4.1.4.1-yc759
/projects/community/GATK/gatk-4.1.4.1/gatk --java-options '-Xmx56G' MarkDuplicates  --REMOVE_DUPLICATES true  --CREATE_INDEX true  --TMP_DIR /mnt/scratch --METRICS_FILE "$OFOL"/"$FIL"_dups.tab --INPUT "$OFOL"/"$FIL"_filt.bam --OUTPUT "$OFOL"/"$FIL"_dupchkd.bam --ASSUME_SORT_ORDER=coordinate  2>> "$OFOL"/"$FIL"_markdups.log
module load bedtools2
bedtools intersect -v -a "$OFOL"/"$FIL"_dupchkd.bam -b Reference_Data/NGS/hg19/hg19-blacklist-kundaje.bed > "$OFOL"/"$FIL"_noblklst.bam
if [[ $4 == "ATAC" ]]
then
	#module load java
	#module load GATK/4.1.4.1-yc759
	/projects/community/GATK/gatk-4.1.4.1/gatk SortSam --INPUT="$OFOL"/"$FIL"_noblklst.bam --OUTPUT="$OFOL"/"$FIL"_shifted.bam -SO=coordinate  --TMP_DIR=/mnt/scratch
	#module load samtools
	samtools index "$OFOL"/"$FIL"_shifted.bam
	igvtools count "$OFOL"/"$FIL"_shifted.bam "$OFOL"/"$FIL"_noblklst.tdf hg19
	/usr/bin/perl count_dup.pl "$OFOL"/"$FIL"_noblklst.bam "$OFOL"/"$FIL"_NirDups.txt
	macs2 callpeak -t "$OFOL"/"$FIL"_shifted.bam --outdir $OFOL -f BAM -n "$FIL"_macs2 -g hs
elif [[ $4 == "CUTRUN" ]]
then
	#bedtools bamtobed -bedpe -i "$OFOL"/"$FIL"_noblklst.bam > "$OFOL"/"$FIL"_sample.bed
	samtools sort -n "$OFOL"/"$FIL"_noblklst.bam -o "$OFOL"/"$FIL"_noblklst_sorted.bam
	samtools view "$OFOL"/"$FIL"_noblklst_sorted.bam | awk '{a=substr($1,length($1)-1); if(a == ".1" || a == ".2"){$1 = substr($1,0,length($1)-2); gsub(" ","\t",$0); print($0)}else{print($0)}}' > "$OFOL"/"$FIL"_noblklst_sorted_noRsuffix_nohead.bam	
	samtools view -H "$OFOL"/"$FIL"_noblklst_sorted.bam | cat - "$OFOL"/"$FIL"_noblklst_sorted_noRsuffix_nohead.bam > "$OFOL"/"$FIL"_noblklst_sorted_noRsuffix.bam
	rm "$OFOL"/"$FIL"_noblklst_sorted_noRsuffix_nohead.bam
	samtools fixmate "$OFOL"/"$FIL"_noblklst_sorted_noRsuffix.bam "$OFOL"/"$FIL"_noblklst_sorted_noRsuffix_fixmate.bam
	bedtools bamtobed -bedpe -i "$OFOL"/"$FIL"_noblklst_sorted_noRsuffix_fixmate.bam > "$OFOL"/"$FIL"_sample.bed
	awk '$1==$4 && $6-$2 < 1000 {print $0}' "$OFOL"/"$FIL"_sample.bed > "$OFOL"/"$FIL"_clean.bed
	cut -f 1,2,6 "$OFOL"/"$FIL"_clean.bed | sort -k1,1 -k2,2n -k3,3n > "$OFOL"/"$FIL"_fragements.bed
	bedtools genomecov -bg -i "$OFOL"/"$FIL"_fragements.bed -g Reference_Data/NGS/../bedtools/genomes/human.hg19.genome > "$OFOL"/"$FIL"_fragements.bedgraph
	samtools index "$OFOL"/"$FIL"_noblklst.bam
	igvtools count --pairs "$OFOL"/"$FIL"_noblklst.bam "$OFOL"/"$FIL"_noblklst.tdf hg19
	/usr/bin/perl count_dup.pl "$OFOL"/"$FIL"_noblklst.bam "$OFOL"/"$FIL"_NirDups.txt
	macs2 callpeak -t "$OFOL"/"$FIL"_sorted.bam --outdir $OFOL -n "$FIL"_macs2 -g hs --nomodel --keep-dup all 
	bash SEACR-master/SEACR_1.3.sh "$OFOL"/"$FIL"_fragements.bedgraph 0.01 non stringent "$OFOL"/"$FIL"_NonStringent
elif [[ $4 == "CHIP" ]]
then
	igvtools count --pairs "$OFOL"/"$FIL"_noblklst.bam "$OFOL"/"$FIL"_noblklst.tdf hg19
	/usr/bin/perl count_dup.pl "$OFOL"/"$FIL"_noblklst.bam "$OFOL"/"$FIL"_NirDups.txt
	macs2 callpeak -t "$OFOL"/"$FIL"_noblklst.bam --outdir $OFOL -f BAM -n "$FIL"_macs2 -g hs
fi
