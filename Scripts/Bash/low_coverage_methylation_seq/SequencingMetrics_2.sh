# Summarizing sequencing and analysis metrics
# Output: summary files with sample name and metrics as tab deliminated file
# Maike Bensberg (maike.bensberg@liu.se)
# 2023-11-28

# metrics that will be analyzed:
### number of reads
### number of properly aligned read pairs
### % duplicates
### number of all CpGs covered
### number of merged CpGs covered
### CpG, CHG and CHH methylation without applying any conversion cutoff
### average global CpG methylation from merged strands

# based on output from script: LowPassMeth_preprocess_MB.sh 

#######################################################################################################################################

# Stop script at runtime errors
set -e

# Output directory
out_dir="/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/LowPassEMseq/combined_analysis/AJLS_treated_HMA/seq_metrics"
mkdir ${out_dir}

# Data folder raw fastq files
data_fastq="/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/LowPassEMseq/combined_analysis/AJLS_treated_HMA/data"
# Data folder duplication metrics
data_dup="/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/LowPassEMseq/combined_analysis/AJLS_treated_HMA/BWA/Alignments/MarkDup/Metrics"
# Data folder alignment files
folder_alignmentStats="/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/LowPassEMseq/combined_analysis/AJLS_treated_HMA/BWA/Alignments/stats"
# Data folder for all methylation (no conversion cutoff)
data_meth_all="/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/LowPassEMseq/combined_analysis/AJLS_treated_HMA/Methylation/allMethylation"
# Data folder for all methylation (min conversion cutoff 90%)
data_CpG_merged="/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/LowPassEMseq/combined_analysis/AJLS_treated_HMA/Methylation/CpGmethylation"

#######################################################################################################################################

# Start message
echo "Summarizing sequencing metrics for all samples in folders"
echo $data_fastq
echo $data_dup
echo $folder_alignmentStats
echo $data_meth_all
echo $data_CpG_merged
echo ""
date

#######################################################################################################################################

# total number of reads
echo -e "sample \t total_read_pairs" > ${out_dir}/NumberReads.txt
for i in `ls ${data_fastq}/*_L001_R1_001.fastq.gz`
do sample=$(basename $i _L001_R1_001.fastq.gz)
lines=$(zcat ${data_fastq}/${sample}_L001_R1_001.fastq.gz|wc -l)
reads=$(echo "$lines/4"|bc)
echo -e "$sample \t $reads" >> ${out_dir}/NumberReads.txt
done


# percent duplication
echo -e "sample \t percent_duplicate" > ${out_dir}/Duplication.txt
for i in `ls ${data_dup}/*_L001.bam.txt`
do sample=$(basename $i _L001.bam.txt)
dupl=$(awk 'NR==8{print $10}' ${data_dup}/${sample}_L001.bam.txt)
echo -e "$sample \t $dupl" >> ${out_dir}/Duplication.txt
done


# number of properly aligned read pairs
echo -e "sample \t aligned_read_pairs" > ${out_dir}/AlignedReads.txt
for i in `ls ${folder_alignmentStats}/*_L001.sam.txt`
do sample=$(basename $i _L001.sam.txt)
total_reads=$(awk 'NR==12{print $1}' ${folder_alignmentStats}/${sample}_L001.sam.txt)
aligned_pairs=$(echo "$total_reads/2"|bc)
echo -e "$sample \t $aligned_pairs" >> ${out_dir}/AlignedReads.txt
done

#######################################################################################################################################

# number of all CpGs covered
echo -e "sample \t CpGs_covered_all" > ${out_dir}/CpGs_covered_all.txt
for i in `ls ${data_meth_all}/*_L001.bam_Trim5_CpG.bedGraph`
do sample=$(basename $i _L001.bam_Trim5_CpG.bedGraph)
CpGs_all=$(tail -n +2 ${data_meth_all}/${sample}_L001.bam_Trim5_CpG.bedGraph | wc -l)
echo -e "$sample \t $CpGs_all" >> ${out_dir}/CpGs_covered_all.txt
done


# number of merged CpGs covered
echo -e "sample \t mergedCpGs_covered" > ${out_dir}/CpGs_merged_covered.txt
for i in `ls ${data_CpG_merged}/*_L001.bam_Trim5_merged_CpG.bedGraph`
do sample=$(basename $i _L001.bam_Trim5_merged_CpG.bedGraph)
CpGs_merged=$(tail -n +2 ${data_CpG_merged}/${sample}_L001.bam_Trim5_merged_CpG.bedGraph | wc -l)
echo -e "$sample \t $CpGs_merged" >> ${out_dir}/CpGs_merged_covered.txt
done

#######################################################################################################################################

# CpG, CHG and CHH methylation without applying any conversion cutoff
# summary text file of methylation for all chromosomes except from control DNA (pUC19 and lambda DNA)
# output (tab deliminated): sample name   average CpG methylation   average CHG methylation   average CHH methylation
# mitochondrial DNA is excluded

echo -e "sample \t CpGmethylation \t CHGmethylation \t CHHmethylation" > ${out_dir}/Summary_Methylation_cytosines_all.txt
for i in `ls ${data_meth_all}/*_L001.bam_Trim5_CpG.bedGraph`
do 
sample=$(basename $i _L001.bam_Trim5_CpG.bedGraph)
CpG=$(awk '/chr/ && !/chrM/ { total += $4; count += 1 } END { print total/count}' ${data_meth_all}/${sample}_CpG.bedGraph)
CHG=$(awk '/chr/ && !/chrM/ { total += $4; count += 1 } END { print total/count}' ${data_meth_all}/${sample}_CHG.bedGraph)
CHH=$(awk '/chr/ && !/chrM/ { total += $4; count += 1 } END { print total/count}' ${data_meth_all}/${sample}_CHH.bedGraph)
echo -e "$sample \t $CpG \t $CHG \t $CHH" >> ${out_dir}/Summary_Methylation_cytosines_all.txt
done

#######################################################################################################################################

# average global CpG methylation from merged strands
# mitochondrial DNA and control DNA excluded

echo -e "sample \t CpGmethylation" > ${out_dir}/Summary_CpGmethylation_merged.txt
for i in `ls ${data_CpG_merged}/*.bam_Trim5_merged_CpG.bedGraph`
do 
sample=$(basename $i .bam_Trim5_merged_CpG.bedGraph)
CpG=$(awk '/chr/ && !/chrM/ { total += $4; count += 1 } END { print total/count}' ${data_CpG_merged}/${sample}.bam_Trim5_merged_CpG.bedGraph)
echo -e "$sample \t $CpG" >> ${out_dir}/Summary_CpGmethylation_merged.txt
done

#######################################################################################################################################

# Completion message
echo ""
echo "Done"
date
