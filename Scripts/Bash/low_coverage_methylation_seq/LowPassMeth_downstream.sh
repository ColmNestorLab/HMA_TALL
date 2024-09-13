# Downstream analysis of low-pass EM-seq data
# Maike Bensberg (maike.bensberg@liu.se) and Sandra Hellberg (sandra.hellberg@liu.se) and Mykolas Malevicius
# 2023-09-28

# input files need to be pre-processed (e.g. LowPassMeth_preprocess_MB.sh)
# Pre-processing includes: quality control, trimming, alignment, mark duplicates, indexing, generating Mbias plots and extracting methylation information


# 0. Set paths to input and output and record package versions 
# 1. Summarizing global CpG methylation
# 2. Summarize CpG methylation per chromosome
# 3. Summarize methylation in CGIs and non-CGIs
# 4. Summarize CpG methylation per region (promoter, intergenic, gene body)
# 5. Summarize methylation at imprinted regions
# 6. Check methylation at LINE, SINE and ERV elements
# 7. Check methylation at HERV elements per chromosome
# 8. Check percent of CpGs covered in each region

# use conda environment
# conda activate HMA_emseq

############################################################################################################

# Stop script at runtime errors
set -e

# Start message
echo "Downstream analysis of low-pass EM-seq data"
echo "###########################################"
echo ""
date
echo ""

### Set input and output directories (change needed!)

# Path to where files can be found (bedgraph files with CpG methylation)
methylation_dir=/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/LowPassEMseq/combined_analysis/AJLS_treated_HMA/Methylation/CpGmethylation

# Path to the output folder for summary files
outpath=/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/LowPassEMseq/combined_analysis/AJLS_treated_HMA/Methylation/CpGmethSummaries
mkdir $outpath

# Path to the folder with regions of interest as bed files
# CGIs, genes, HERVs, imprinted genes, LINE elements, promoters
ROI_path=/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/LowPassEMseq/combined_analysis/AJLS_treated_HMA/ROI_beds



### Output data folder and package versions to log

echo "Methylation will be summarized for files from the following folder:"
echo $methylation_dir
echo "The files are expected to be pre-processed (bedgraph files with methylation information for each CpG)"
echo ""
echo "files with regions of interest should be available here:"
echo $ROI_path
echo ""
echo "Package versions:"
echo "BedTools: $(bedtools --version)"
echo ""

########################################################################################################

### 1. Summarizing global CpG methylation

# summary text file of methylation for all chromosomes except from control DNA (pUC19 and lambda DNA)
# mitochondrial DNA is also excluded
# output (tab deliminated): sample name   average CpG methylation

echo "Summarizing global CpG methylation"
echo "##################################"
echo ""

echo -e "sample \t CpGmethylation" >> $outpath/GlobalCpGmethylation.txt

for i in `ls $methylation_dir/*.bam_Trim5_merged_CpG.bedGraph`
do 
sample=$(basename $i .bam_Trim5_merged_CpG.bedGraph)
CpG=$(awk '/chr/ && !/chrM/ { total += $4; count += 1 } END { print total/count}' $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph)
echo -e "$sample \t $CpG" >> $outpath/GlobalCpGmethylation.txt
done

########################################################################################################

### 2. Summarize CpG methylation per chromosome

echo "Summarizing CpG methylation per chromosome"
echo "##########################################"
echo ""

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY" "chrM")

echo -e "sample \t chr1 \t chr2 \t chr3 \t chr4 \t chr5 \t chr6 \t chr7 \t chr8 \t chr9 \t chr10 \t chr11 \t chr12 \t chr13 \t chr14 \t chr15 \t chr16 \t chr17 \t chr18 \t chr19 \t chr20 \t chr21 \t chr22 \t chrX \t chrY \t chrM" >> $outpath/CpGmethylation_perChr.txt

for i in `ls $methylation_dir/*.bam_Trim5_merged_CpG.bedGraph`
do 
sample=$(basename $i .bam_Trim5_merged_CpG.bedGraph)
l=$(echo -e "$sample")
for c in "${chromosomes[@]}"
do
meth=$(awk -v chr="$c" '$1 == chr { total += $4; count++ } END { if (count > 0) print total/count }' $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph)
l=$(echo -e "$l \t $meth")
done
echo -e $l >> $outpath/CpGmethylation_perChr.txt
done

########################################################################################################

### 3. Summarize methylation in CGIs and non-CGIs

echo "Summarizing CpG methylation in CGIs and non CGIs"
echo "################################################"
echo "bed-file of CGIs (hg38) downloaded from UCSC: 2023-09-29"
echo ""

echo -e "sample \t CGI \t non_CGI" >> $outpath/CGIs_CpGmethylation.txt
for i in `ls $methylation_dir/*.bam_Trim5_merged_CpG.bedGraph`
do
sample=$(basename $i .bam_Trim5_merged_CpG.bedGraph)
CGI_meth=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/CGIs_hg38_UCSC_20230929.bed -wa | awk '{ total += $4; count += 1 } END { print total/count}')
nonCGI_meth=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/CGIs_hg38_UCSC_20230929.bed -v | awk '{ total += $4; count += 1 } END { print total/count}')
echo -e "$sample \t $CGI_meth \t $nonCGI_meth" >> $outpath/CGIs_CpGmethylation.txt
done

########################################################################################################

### 4. Summarize CpG methylation per region (promoter, intergenic, gene body)

echo "Summarizing CpG methylation in promoters, genes and intergenic regions"
echo "######################################################################"
echo "bed-file of genes and promoters (hg38) downloaded from UCSC: 2023-09-29"
echo "promoters defined as 1kb upstream and 500 bp downstream of the TSS"
echo "promoters divided into promoter_all, promoter_CGI and promoter_nonCGI"
echo ""

# extract promoter region from genes bed file
# awk 'BEGIN {OFS="\t"} $6 == "-" {$2 = $3-500} {print}' genes_hg38_ucsc_20230929.bed > temp.bed
# awk 'BEGIN {OFS="\t"} $6 == "-" {$3 = $3+1000} {print}' temp.bed > temp2.bed
# awk 'BEGIN {OFS="\t"} $6 == "+" {$3 = $2+500} {print}' temp2.bed > temp3.bed
# awk 'BEGIN {OFS="\t"} $6 == "+" {$2 = $2-1000} {print}' temp3.bed > temp4.bed
# awk 'BEGIN {OFS="\t"} $2 < 0 {$2 = 0} {print}' temp4.bed > temp5.bed
# cat temp5.bed > promoters_1kbup_500bodown_hg38_ucsc_20240115.bed



# CGI promoters defined as one of the following should be true: 20% of the promoter overlaps with a CGI or 50% of the CGI overlaps with the promoter
# bedtools intersect -a promoters_1kbup_500bodown_hg38_ucsc_20240115.bed -b CGIs_hg38_UCSC_20230929.bed -wa -u -f 0.2 -F 0.5 -e > promoter_CGI.bed

echo -e "sample \t promoters_all \t promoters_CGI \t promoters_nonCGI \t genes \t intergenic" >> $outpath/Genes_CpGmethylation.txt
for i in `ls $methylation_dir/*.bam_Trim5_merged_CpG.bedGraph`
do
sample=$(basename $i .bam_Trim5_merged_CpG.bedGraph)
promoters_all_meth=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/promoters_1kbup_500bodown_hg38_ucsc_20240115.bed -wa -u | awk '{ total += $4; count += 1 } END { print total/count}')
promoters_CGI_meth=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/promoter_CGI.bed -wa -u | awk '{ total += $4; count += 1 } END { print total/count}')
promoters_nonCGI_meth=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/promoter_nonCGI.bed -wa -u | awk '{ total += $4; count += 1 } END { print total/count}')
genes_meth=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/genes_min_hg38_ucsc_20230929.bed -wa -u | awk '{ total += $4; count += 1 } END { print total/count}')
inter_meth=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/promoters_1kbup_500bodown_hg38_ucsc_20240115.bed $ROI_path/genes_min_hg38_ucsc_20230929.bed -v | awk '{ total += $4; count += 1 } END { print total/count}')
echo -e "$sample \t $promoters_all_meth \t $promoters_CGI_meth \t $promoters_nonCGI_meth \t $genes_meth \t $inter_meth" >> $outpath/Genes_CpGmethylation.txt
done

########################################################################################################

### 5. Summarize methylation at imprinted regions

echo "Summarizing CpG methylation at imprinted regions"
echo "################################################"
echo "bed-file of imprinted regions summarized from Akbari et al, 2022 eLife: https://doi.org/10.7554/eLife.77898"
echo ""

echo -e "sample \t imprinted \t non_imprinted" >> $outpath/imprinted_regions_CpGmethylation.txt
for i in `ls $methylation_dir/*.bam_Trim5_merged_CpG.bedGraph`
do
sample=$(basename $i .bam_Trim5_merged_CpG.bedGraph)
imprinted_meth=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/imprinted_Akbari2022_published.bed -wa -u | awk '{ total += $4; count += 1 } END { print total/count}')
non_imprinted_meth=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/imprinted_Akbari2022_published.bed -v | awk '{ total += $4; count += 1 } END { print total/count}')
echo -e "$sample \t $imprinted_meth \t $non_imprinted_meth" >> $outpath/imprinted_regions_CpGmethylation.txt
done

########################################################################################################

### 6. Check methylation at LINE, SINE and ERV elements

# Methylation at LINE elements is often used as an indicator for loss of global methylation
# bed file with locations of HERV elements downloaded from UCSC Table Browser by Sandra
# using her file

echo "Summarizing CpG methylation at LINE elements and ERVs"
echo "######################################################"
echo "bed files with locations of ERVs, LINEs and SINEs downloaded from UCSC Table Browser by Sandra"
echo ""


echo -e "sample \t LINEs \t SINEs \t ERVs" >> $outpath/LINEandHERV_CpGmethylation.txt
for i in `ls $methylation_dir/*.bam_Trim5_merged_CpG.bedGraph`
do sample=$(basename $i .bam_Trim5_merged_CpG.bedGraph)
CpG_LINE=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/hg38_rmsk_LINE.bed -u -wa | awk '{ total += $4; count += 1 } END { print total/count}')
CpG_SINE=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/hg38_rmsk_SINE.bed -u -wa | awk '{ total += $4; count += 1 } END { print total/count}')
CpG_ERV=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/hg38_rmsk_ERV.bed -u -wa | awk '{ total += $4; count += 1 } END { print total/count}')
echo -e "$sample \t $CpG_LINE \t $CpG_SINE \t $CpG_ERV" >> $outpath/LINEandHERV_CpGmethylation.txt
done

########################################################################################################

### 7. Check methylation at ERV elements per chromosome

# bed file with locations of HERV elements downloaded from UCSC Table Browser by Sandra
# using her file

echo "Summarizing CpG methylation at HERVs per chromosome"
echo "###################################################"
echo ""

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

echo -e "sample \t chr1 \t chr2 \t chr3 \t chr4 \t chr5 \t chr6 \t chr7 \t chr8 \t chr9 \t chr10 \t chr11 \t chr12 \t chr13 \t chr14 \t chr15 \t chr16 \t chr17 \t chr18 \t chr19 \t chr20 \t chr21 \t chr22 \t chrX \t chrY" >> $outpath/HERV_CpGmethylation_perChr.txt

for i in `ls $methylation_dir/*.bam_Trim5_merged_CpG.bedGraph`
do 
sample=$(basename $i .bam_Trim5_merged_CpG.bedGraph)
l=$(echo -e "$sample")
HERVintersect_file="$outpath/temp_intersection.bed"
bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/hg38_rmsk_ERV.bed -u -wa > "$HERVintersect_file"
for c in "${chromosomes[@]}"
do
meth=$(awk -v chr="$c" '$1 == chr { total += $4; count++ } END { if (count > 0) print total/count }' $HERVintersect_file)
l=$(echo -e "$l \t $meth")
done
echo -e $l >> $outpath/HERV_CpGmethylation_perChr.txt
rm -f "$HERVintersect_file"
done

########################################################################################################

### 8. Check percent of CpGs covered in each region

# bed file with all CpGs in the genome converted from a MethylDackel --cytosine_report file (awk 'BEGIN {OFS="\t"} NR%2==1 {print $1, $2-1, $2+1}' cytosine_report.txt > CpGs.bed

echo "Summarizing Ccoverage for all regions"
echo "###################################################"
echo "percent of CpGs covered by one read or more"
echo ""

total_global=$(cat $ROI_path/CpGs.bed | wc -l)
total_CGI=$(bedtools intersect -a $ROI_path/CpGs.bed -b $ROI_path/CGIs_hg38_UCSC_20230929.bed -u -wa | wc -l)
total_nonCGI=$(bedtools intersect -a $ROI_path/CpGs.bed -b $ROI_path/CGIs_hg38_UCSC_20230929.bed -v -wa | wc -l)
total_promoter_all=$(bedtools intersect -a $ROI_path/CpGs.bed -b $ROI_path/promoters_1kbup_500bodown_hg38_ucsc_20240115.bed -u -wa | wc -l)
total_promoter_CGI=$(bedtools intersect -a $ROI_path/CpGs.bed -b $ROI_path/promoter_CGI.bed -u -wa | wc -l)
total_promoter_nonCGI=$(bedtools intersect -a $ROI_path/CpGs.bed -b $ROI_path/promoter_nonCGI.bed -u -wa | wc -l)
total_genes=$(bedtools intersect -a $ROI_path/CpGs.bed -b $ROI_path/genes_min_hg38_ucsc_20230929.bed -u -wa | wc -l)
total_inter=$(bedtools intersect -a $ROI_path/CpGs.bed -b $ROI_path/genes_min_hg38_ucsc_20230929.bed $ROI_path/promoters_1kbup_500bodown_hg38_ucsc_20240115.bed -v -wa | wc -l)
total_LINE=$(bedtools intersect -a $ROI_path/CpGs.bed -b $ROI_path/hg38_rmsk_LINE.bed -u -wa | wc -l)
total_SINE=$(bedtools intersect -a $ROI_path/CpGs.bed -b $ROI_path/hg38_rmsk_SINE.bed -u -wa | wc -l)
total_ERV=$(bedtools intersect -a $ROI_path/CpGs.bed -b $ROI_path/hg38_rmsk_ERV.bed -u -wa | wc -l)

echo -e "sample \t cov_global \t cov_CGI \t cov_nonCGI \t cov_promoter_all \t cov_promoter_CGI \t cov_promoter_nonCGI \t cov_genes \t cov_intergenic \t cov_LINEs \t cov_SINEs \t cov_ERVs" >> $outpath/coverage.txt
for i in `ls $methylation_dir/*.bam_Trim5_merged_CpG.bedGraph`
do sample=$(basename $i .bam_Trim5_merged_CpG.bedGraph)
# global
CpGs=$(cat $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph | wc -l)
cov_global=$(echo "scale=5;$CpGs/$total_global"|bc -l)
# CGI
CpGs=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/CGIs_hg38_UCSC_20230929.bed -u -wa | wc -l)
cov_CGI=$(echo "scale=5;$CpGs/$total_CGI"|bc -l)
# nonCGI
CpGs=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/CGIs_hg38_UCSC_20230929.bed -v -wa | wc -l)
cov_nonCGI=$(echo "scale=5;$CpGs/$total_nonCGI"|bc -l)
# all promoters
CpGs=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/promoters_1kbup_500bodown_hg38_ucsc_20240115.bed -u -wa | wc -l)
cov_promoter_all=$(echo "scale=5;$CpGs/$total_promoter_all"|bc -l)
# CGI promoters
CpGs=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/promoter_CGI.bed -u -wa | wc -l)
cov_promoter_CGI=$(echo "scale=5;$CpGs/$total_promoter_CGI"|bc -l)
# nonCGI promoters
CpGs=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/promoter_nonCGI.bed -u -wa | wc -l)
cov_promoter_nonCGI=$(echo "scale=5;$CpGs/$total_promoter_nonCGI"|bc -l)
# gene bodies
CpGs=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/genes_min_hg38_ucsc_20230929.bed -u -wa | wc -l)
cov_genes=$(echo "scale=5;$CpGs/$total_genes"|bc -l)
# intergenic
CpGs=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/genes_min_hg38_ucsc_20230929.bed $ROI_path/promoters_1kbup_500bodown_hg38_ucsc_20240115.bed -v -wa | wc -l)
cov_inter=$(echo "scale=5;$CpGs/$total_inter"|bc -l)
# LINEs
CpGs=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/hg38_rmsk_LINE.bed -u -wa | wc -l)
cov_LINE=$(echo "scale=5;$CpGs/$total_LINE"|bc -l)
# SINEs
CpGs=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/hg38_rmsk_SINE.bed -u -wa | wc -l)
cov_SINE=$(echo "scale=5;$CpGs/$total_SINE"|bc -l)
# ERVs
CpGs=$(bedtools intersect -a $methylation_dir/${sample}.bam_Trim5_merged_CpG.bedGraph -b $ROI_path/hg38_rmsk_ERV.bed -u -wa | wc -l)
cov_ERV=$(echo "scale=5;$CpGs/$total_ERV"|bc -l)
echo -e "$sample \t $cov_global \t $cov_CGI \t $cov_nonCGI \t $cov_promoter_all \t $cov_promoter_CGI \t $cov_promoter_nonCGI \t $cov_genes \t $cov_inter \t $cov_LINE \t $cov_SINE \t $cov_ERV" >> $outpath/coverage.txt
done



########################################################################################################

# Completion message
echo ""
echo "Done"
date





