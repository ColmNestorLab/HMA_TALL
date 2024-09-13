# Summarizing sequencing and analysis metrics of EM-seq data
# Output: summary files with sample name and metrics as tab deliminated file
# Maike Bensberg (maike.bensberg@liu.se)
# 2023-12-14

# metrics that will be analyzed:
### number of reads
### % duplicates
### number of properly aligned read pairs
### number of CpGs (merged) covered by at least 1 read
### number of CpGs (merged) covered by at least 5 reads
### average CpG, CHG and CHH methylation without minimum coverage (strands not merged)
### average CpG methyation for control DNAs (no minimum coverage)
### average CpG methyation for control DNAs (coverage >= 5)
### average CpG methylation without minimum coverage
### average CpG methylation with minimum 5 reads coverage

# based on output from script: EMseq_CpGm_analysis.sh

#######################################################################################################################################

# Stop script at runtime errors
set -e

# Output directory
out_dir="/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/seq_metrics"
mkdir ${out_dir}

# Data folder with raw fastq files
data_fastq="/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/data2"
# Data folder with duplication metrics
data_dup="/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/BWA/Alignments/MarkDup/Metrics"
# Data folder alignment files
folder_alignmentStats="/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/BWA/Alignments/stats"
# Data folder for all methylation (no coverage cutoff and strands not merged)
data_meth_all="/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/Methylation/allMethylation"
# Data folder for CpG methylation (strands merged)
data_CpG_merged="/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/Methylation/CpGmethylation"
# Data folder for CpG methylation (strands merged and coverageat least 5)
data_CpG_merged_cov5="/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/EMseq/LOUCY_SUPT1_HMAtreated_12_2023/Methylation/CpGmethylation_cov5"

#######################################################################################################################################

# Start message
echo "Summarizing sequencing metrics for all samples in the following folders:"
echo $data_fastq
echo $data_dup
echo $folder_alignmentStats
echo $data_meth_all
echo $data_CpG_merged
echo $data_CpG_merged_cov5
echo ""
date

#######################################################################################################################################

# total number of reads
echo -e "sample \t total_read_pairs" > ${out_dir}/NumberReads.txt
for i in `ls ${data_fastq}/*_R1_001.fastq.gz`
do sample=$(basename $i _R1_001.fastq.gz)
lines=$(zcat ${data_fastq}/${sample}_R1_001.fastq.gz|wc -l)
reads=$(echo "$lines/4"|bc)
echo -e "$sample \t $reads" >> ${out_dir}/NumberReads.txt
done


# percent duplication
echo -e "sample \t percent_duplicates" > ${out_dir}/Duplication.txt
for i in `ls ${data_dup}/*.bam.txt`
do sample=$(basename $i .bam.txt)
dupl=$(awk 'NR==8{print $10}' ${data_dup}/${sample}.bam.txt)
echo -e "$sample \t $dupl" >> ${out_dir}/Duplication.txt
done


# number of properly aligned read pairs
echo -e "sample \t aligned_read_pairs" > ${out_dir}/AlignedReads.txt
for i in `ls ${folder_alignmentStats}/*.sam.txt`
do sample=$(basename $i .sam.txt)
total_reads=$(awk 'NR==12{print $1}' ${folder_alignmentStats}/${sample}.sam.txt)
aligned_pairs=$(echo "$total_reads/2"|bc)
echo -e "$sample \t $aligned_pairs" >> ${out_dir}/AlignedReads.txt
done

#######################################################################################################################################

# number of all CpGs covered (strands merged)
echo -e "sample \t CpGs_covered_all" > ${out_dir}/CpGs_covered_all_merged.txt
for i in `ls ${data_CpG_merged}/*.bam_Trim5_merged_CpG.bedGraph`
do sample=$(basename $i .bam_Trim5_merged_CpG.bedGraph)
CpGs_all=$(tail -n +2 ${data_CpG_merged}/${sample}.bam_Trim5_merged_CpG.bedGraph | wc -l)
echo -e "$sample \t $CpGs_all" >> ${out_dir}/CpGs_covered_all_merged.txt
done


# number of CpGs covered with minimum of 5 reads (strands merged)
echo -e "sample \t mergedCpGs_covered" > ${out_dir}/CpGs_merged_covered.txt
for i in `ls ${data_CpG_merged_cov5}/*.bam_Trim5_merged_cov5_CpG.bedGraph`
do sample=$(basename $i .bam_Trim5_merged_cov5_CpG.bedGraph)
CpGs_merged=$(tail -n +2 ${data_CpG_merged_cov5}/${sample}.bam_Trim5_merged_cov5_CpG.bedGraph | wc -l)
echo -e "$sample \t $CpGs_merged" >> ${out_dir}/CpGs_merged_covered.txt
done

#######################################################################################################################################

# CpG, CHG and CHH methylation without applying any conversion or coveragecutoff
# summary text file of methylation for all chromosomes except from control DNA (pUC19 and lambda DNA)
# mitochondrial DNA is excluded
# output (tab deliminated): sample name   average CpG methylation   average CHG methylation   average CHH methylation

echo -e "sample \t CpGmethylation \t CHGmethylation \t CHHmethylation" > ${out_dir}/Summary_Methylation_cytosines_all.txt
for i in `ls ${data_meth_all}/*_CpG.bedGraph`
do 
sample=$(basename $i _CpG.bedGraph)
CpG=$(awk '/chr/ && !/chrM/ { total += $4; count += 1 } END { print total/count}' ${data_meth_all}/${sample}_CpG.bedGraph)
CHG=$(awk '/chr/ && !/chrM/ { total += $4; count += 1 } END { print total/count}' ${data_meth_all}/${sample}_CHG.bedGraph)
CHH=$(awk '/chr/ && !/chrM/ { total += $4; count += 1 } END { print total/count}' ${data_meth_all}/${sample}_CHH.bedGraph)
echo -e "$sample \t $CpG \t $CHG \t $CHH" >> ${out_dir}/Summary_Methylation_cytosines_all.txt
done

#######################################################################################################################################

# average CpG methyation for control DNA (no minimum coverage)
# pUC19 (M77789.2) should methylated at all CpGs and lambda (J02459.1) should be unmethylated at all CpGs

echo -e "sample \t pUC19_methylation \t lambda_methylation" > ${out_dir}/controlDNA_methylation.txt
for i in `ls ${data_CpG_merged}/*.bam_Trim5_merged_CpG.bedGraph`
do 
sample=$(basename $i .bam_Trim5_merged_CpG.bedGraph)
pUC19_meth=$(awk '/M77789.2/ { total += $4; count += 1 } END { print total/count}' ${data_CpG_merged}/${sample}.bam_Trim5_merged_CpG.bedGraph)
lambda_meth=$(awk '/J02459.1/ { total += $4; count += 1 } END { print total/count}' ${data_CpG_merged}/${sample}.bam_Trim5_merged_CpG.bedGraph)
echo -e "$sample \t $pUC19_meth \t $lambda_meth" >> ${out_dir}/controlDNA_methylation.txt
done


# average CpG methyation for control DNA (coverage >= 5)
# pUC19 (M77789.2) should methylated at all CpGs and lambda (J02459.1) should be unmethylated at all CpGs

echo -e "sample \t pUC19_methylation \t lambda_methylation" > ${out_dir}/controlDNA_methylation_cov5.txt
for i in `ls ${data_CpG_merged_cov5}/*.bam_Trim5_merged_cov5_CpG.bedGraph`
do 
sample=$(basename $i .bam_Trim5_merged_cov5_CpG.bedGraph)
pUC19_meth=$(awk '/M77789.2/ { total += $4; count += 1 } END { print total/count}' ${data_CpG_merged_cov5}/${sample}.bam_Trim5_merged_cov5_CpG.bedGraph)
lambda_meth=$(awk '/J02459.1/ { total += $4; count += 1 } END { print total/count}' ${data_CpG_merged_cov5}/${sample}.bam_Trim5_merged_cov5_CpG.bedGraph)
echo -e "$sample \t $pUC19_meth \t $lambda_meth" >> ${out_dir}/controlDNA_methylation_cov5.txt
done

#######################################################################################################################################

# average CpG methylation without minimum coverage
# mitochondrial DNA and control DNA excluded

echo -e "sample \t CpGmethylation" > ${out_dir}/Summary_CpGmethylation_all.txt
for i in `ls ${data_CpG_merged}/*.bam_Trim5_merged_CpG.bedGraph`
do 
sample=$(basename $i .bam_Trim5_merged_CpG.bedGraph)
CpG=$(awk '/chr/ && !/chrM/ { total += $4; count += 1 } END { print total/count}' ${data_CpG_merged}/${sample}.bam_Trim5_merged_CpG.bedGraph)
echo -e "$sample \t $CpG" >> ${out_dir}/Summary_CpGmethylation_all.txt
done


# average CpG methylation with minimum 5 reads coverage
# mitochondrial DNA and control DNA excluded

echo -e "sample \t CpGmethylation" > ${out_dir}/Summary_CpGmethylation_cov5.txt
for i in `ls ${data_CpG_merged_cov5}/*.bam_Trim5_merged_cov5_CpG.bedGraph`
do 
sample=$(basename $i .bam_Trim5_merged_cov5_CpG.bedGraph)
CpG=$(awk '/chr/ && !/chrM/ { total += $4; count += 1 } END { print total/count}' ${data_CpG_merged_cov5}/${sample}.bam_Trim5_merged_cov5_CpG.bedGraph)
echo -e "$sample \t $CpG" >> ${out_dir}/Summary_CpGmethylation_cov5.txt
done

#######################################################################################################################################

# Completion message
echo ""
echo "Done"
date



