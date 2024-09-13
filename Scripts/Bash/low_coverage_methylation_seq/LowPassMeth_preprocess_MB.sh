# Preprocessing of EM-seq data
# Maike Bensberg (maike.bensberg@liu.se) and Sandra Hellberg (sandra.hellberg@liu.se) and Mykolas Malevicius
# 2023-09-28

# remember to activate conda environment (HMA_emseq)


# 0. Set paths to input and output and record package versions 
# 1. FastQC analysis of all samples
# 2. Quality trimming using cutadapt wrapped in Trim-Galore
# 3. Alignment using BWA-meth
# 4. Get stats for alignment
# 5. Convert alignment files to bam
# 6. Sort aligned files
# 7. Mark duplicates using picard
# 8. Index files
# 9. Create M-bias plots
# 10. Extracting methylation data for each cytosine
# 11. Extracting methylation data per CpG with strands merged
# 12. Extracting coverage per CpG (strands merged)
 
# Requirement: a reference genome must already have been indexed using BWA-MEM or BWA-MEM2
# bwameth.py index $REF

# use conda environment
# conda activate HMA_emseq

############################################################################################################

# Stop script at runtime errors
set -e

# Start message
echo "Precrocessing EM-seq data"
date
echo ""


### Set input and output directories (change needed!)

# Path to where fastq files can be found 
data_dir=/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/LowPassEMseq/combined_analysis/AJLS_treated_HMA/data

# Path to output folder
outpath=/mnt/wwn-0x5000cca28fd30c7d-part1/HMA_TALL/LowPassEMseq/combined_analysis/AJLS_treated_HMA

# Path to reference genome
# Requirement: a reference genome must already have been indexed using BWA-MEM or BWA-MEM2
# bwameth.py index $REF
genome_path=/mnt/wwn-0x5000cca28fd30c7d-part1/OncoABIS/EM-seq/genomes/BWA_Genome


### Output data folder and package versions to log

echo "Input data folder:"
echo $data_dir
echo ""
echo "Package versions:"
echo "FastQC: $(fastqc --version)"
echo "TrimGalore: $(trim_galore --version)"
echo "BWA-meth: $(bwameth.py --version)"
echo "Samtools: $(samtools --version)"
echo "Picard: $(picard MarkDuplicates --version)"
echo ""

############################################################################################################

### 1.  FastQC analysis of all samples

mkdir $outpath/FastQC
list_fastq=($(find "$data_dir" -type f -name '*fastq.gz' -print))

echo "Start FastQC"
echo ""

for i in "${list_fastq[@]}"
do 
Reads=$(basename "$i")
fastqc -o $outpath/FastQC $data_dir/${Reads}
done

############################################################################################################

### 2. Quality trimming using cutadapt wrapped in Trim-Galore

mkdir $outpath/Trimgalore

echo "Start Trimming"
echo ""

for i in `ls $data_dir/*_R1_001.fastq.gz`
do 
Reads=$(basename $i _R1_001.fastq.gz)
trim_galore -j 3 --paired --fastqc -o $outpath/Trimgalore $data_dir/${Reads}_R1_001.fastq.gz $data_dir/${Reads}_R2_001.fastq.gz
done

############################################################################################################

### 3.  Alignment using BWA-meth

mkdir $outpath/BWA
mkdir $outpath/BWA/Alignments
mkdir $outpath/BWA/Alignments/sam

echo "Start Alignment using BWA-meth"
echo ""

for i in `ls $outpath/Trimgalore/*1.fq.gz`
do 
TrimmedReads=$(basename $i _R1_001_val_1.fq.gz)
bwameth.py --reference ${genome_path}/Combined_Genome.fa -t 20 $outpath/Trimgalore/${TrimmedReads}_R1_001_val_1.fq.gz $outpath/Trimgalore/${TrimmedReads}_R2_001_val_2.fq.gz > $outpath/BWA/Alignments/sam/${TrimmedReads}.sam
done

############################################################################################################

### 4.  Get stats for alignment

mkdir  $outpath/BWA/Alignments/stats

echo "Get Stats for alignment"
echo ""

for i in `ls $outpath/BWA/Alignments/sam/*.sam`
do 
AlignmentsSAM=$(basename $i)
samtools flagstat $outpath/BWA/Alignments/sam/${AlignmentsSAM} > $outpath/BWA/Alignments/stats/${AlignmentsSAM}.txt
done

############################################################################################################

### 5.  Convert alignment files to bam

mkdir $outpath/BWA/Alignments/bam

echo "Convert sam to bam"
echo ""

for i in `ls $outpath/BWA/Alignments/sam/*.sam`
do 
AlignmentsSAM=$(basename $i .sam)
samtools view -b -@ 8 $i > $outpath/BWA/Alignments/bam/${AlignmentsSAM}.bam
done

############################################################################################################

# 6.  Sort aligned files

mkdir $outpath/BWA/Alignments/sortedBam

echo "Sort aligned bam files"
echo ""

for i in `ls $outpath/BWA/Alignments/bam/*.bam`
do 
AlignmentsBAM=$(basename $i)
samtools sort -@ 8 $i > $outpath/BWA/Alignments/sortedBam/sorted.${AlignmentsBAM}
done

############################################################################################################

### 7.  Mark duplicates using picard

mkdir $outpath/BWA/Alignments/MarkDup
mkdir $outpath/BWA/Alignments/MarkDup/Metrics

echo "Mark duplicates using Picard"
echo ""

for i in `ls $outpath/BWA/Alignments/sortedBam/*.bam`
do 
SortedReads=$(basename $i)
picard MarkDuplicates \
I=$i \
O=$outpath/BWA/Alignments/MarkDup/markdup.${SortedReads} \
M=$outpath/BWA/Alignments/MarkDup/Metrics/dupmetrics.${SortedReads}.txt \
TAG_DUPLICATE_SET_MEMBERS=true 
done

############################################################################################################

### 8.  Index files

echo "Index files"
echo ""

for i in `ls $outpath/BWA/Alignments/MarkDup/*.bam`
do 
MarkDupBAM=$(basename $i)
samtools index -@ 8 $i
done

############################################################################################################

### 9. Create M-bias plots

mkdir $outpath/Methylation
mkdir $outpath/Methylation/Mbias
mkdir $outpath/Methylation/Mbias/Mbias
mkdir $outpath/Methylation/Mbias/MbiasTrim5

echo "Creating M-bias plots"
echo ""

# M-bias plots without any trimming
for i in `ls $outpath/BWA/Alignments/MarkDup/*.bam`
do 
MarkedDup=$(basename $i)
MethylDackel mbias -@ 20 $genome_path/Combined_Genome.fa $i $outpath/Methylation/Mbias/Mbias/${MarkedDup}_mbias
done

# M-bias plots with trimming 5 bases from each end
for i in `ls $outpath/BWA/Alignments/MarkDup/*.bam`
do 
MarkedDup=$(basename $i)
MethylDackel mbias --nOT 6,5,6,5 --nOB 6,5,6,5 -@ 20 $genome_path/Combined_Genome.fa $i $outpath/Methylation/Mbias/MbiasTrim5/${MarkedDup}_mbiasTrim5
done

############################################################################################################

### 10. Extracting methylation data for each cytosine

# 5 bases trimmed from each end of each read
# information on methylation divided by cytosine context (CpG, CHG, CHH) with seperate files for each


mkdir $outpath/Methylation/allMethylation

# output as bedgraph; one file for each cytosine context (CG, CHG and CHH)

echo "Extracting information on methylation per cytosine"
echo ""

for i in `ls $outpath/BWA/Alignments/MarkDup/*.bam`
do 
MarkedDup=$(basename $i);
MethylDackel extract --nOT 6,5,6,5 --nOB 6,5,6,5 -@ 20 --CHG --CHH -o $outpath/Methylation/allMethylation/${MarkedDup}_Trim5 $genome_path/Combined_Genome.fa $i
done

############################################################################################################

# 11. Extracting methylation data per CpG with strands merged

# 5 bases trimmed from each end of each read
# information merged for each CpG (no strand-specific information)

mkdir $outpath/Methylation/CpGmethylation

echo "Extracting methylation per CpG"
echo ""


for i in `ls $outpath/BWA/Alignments/MarkDup/*.bam`
do 
MarkedDup=$(basename $i);
MethylDackel extract --nOT 6,5,6,5 --nOB 6,5,6,5 -@ 20 --mergeContext -o $outpath/Methylation/CpGmethylation/${MarkedDup}_Trim5_merged $genome_path/Combined_Genome.fa $i
done

############################################################################################################

# 12. Extracting coverage per CpG (strands merged)

mkdir $outpath/Methylation/coverage

echo "Extracting coverage per CpG"
echo ""

for i in `ls $outpath/BWA/Alignments/MarkDup/*.bam`
do 
MarkedDup=$(basename $i);
MethylDackel extract --nOT 6,5,6,5 --nOB 6,5,6,5 -@ 20 --mergeContext --counts -o $outpath/Methylation/coverage/${MarkedDup}_Trim5_merged $genome_path/Combined_Genome.fa $i
done


# 5 bases trimmed from each end of each read
# information merged for each CpG (no strand-specific information)

############################################################################################################

# Completion message
echo "Preprocessing done"
date
echo ""

