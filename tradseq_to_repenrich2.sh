#!/bin/bash
#SBATCH --job-name=re_pipe
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mkmwong@stanford.edu
#SBATCH -p hns,normal

#### README  ##################################################################
#### This script is for directly adapting thhe raw fastq file from         ####
#### Morrison Lab in-house UV-lesion sequence method TRAD-Seq to repeat    ####
#### quantifying tools Repenrich2.                                         ####
#### - This script is only for running on Sherlock Cluster at Stanford.    #### 
#### - Repeenrich2 need to be installed prior to running this script.      ####
#### - Only python version: 2.x is compatible with Repenrich2.             ####
#### - Please fill in the variables below prior to running the script.     ####
#### - Repenrich2 need 64gb of memory to run so don't change the mem or    ####
####   cpu per task setting, but time setting can be adjusted depending on ####
####   numbe of samples/ how big the fastqs are.                           ####
###############################################################################

#### !!!! Please edit the form before running the script !!!! ####
fastqname=('TP53_200J_UVB_1' 'TP53_200J_UVB_2')    #prefix of each fastq file 
suffix=_r    #input suffix
ori_ext=.fastq    #input extention
new_ext=.fastq    #option: .fastq or .fastq.gz, whether cutadapt output should be zipped
file_path=/scratch/groups/ashbym/mandy/tp53_repenrich2    #path of input file
outfile_path=/scratch/groups/ashbym/mandy/tp53_repenrich2    #path to store all output
genome_path=/home/groups/ashbym/mandy/lambda_human/hg19_lambda    #path to genome sequence and bowtie2 index
picard_path=/home/groups/ashbym/mandy    #path to picard executabe
script_path=/scratch/groups/ashbym/mandy    #path where this script is
repenrich2_path=/home/groups/ashbym/mandy/RepEnrich2    #path of repenrich installation
repmskr_path=/home/groups/ashbym/mandy    #path of repeatmasker file required for running repeat masker
repenrich2_setup=/home/groups/ashbym/mandy/RepEnrich2_setup_hg19    #path of setup of repenrich2

#### Step0: Load required modules
echo "Loading required modules..."
module load python/2.7.13 R java
module load biology py-cutadapt/1.18_py27 bowtie2 samtools 
module load bedtools py-biopython

for barcode in "${fastqname[@]}"; do

    echo "Currnetly running sample" $barcode
#### Step1: trim raw reads of adaptor sequences ####
    echo "Trimming adaptor sequences from reads..."
    cutadapt -m 50:55 -a CTGTCTCTTATACACATCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A "CACTGCNNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;max_error_rate=0.15;min_overlap=11"   -o $outfile_path"/"$barcode"_trimmed"$suffix"1"$new_ext -p $outfile_path"/"$barcode"_trimmed"$suffix"2"$new_ext $file_path"/"$barcode$suffix"1"$ori_ext $file_path"/"$barcode$suffix"2"$ori_ext

#### Step2: Align trimmed reads to concatenated human-lambda genome ####
#### At the same time extract UMI for de-duplication of reads ####
#### Aligned output is stored as sorted bam file ####
    echo "Aligning reads..."
    bowtie2 --very-sensitive --sam-no-qname-trunc -q -p 16 -x $genome_path -1 $outfile_path"/"$barcode"_trimmed"$suffix"1"$new_ext -2 $outfile_path"/"$barcode"_trimmed"$suffix"2"$new_ext | $script_path"/extract_umi.sh" | samtools view -bS - | samtools sort - > $outfile_path"/"$barcode".bam"

#### Step3: Index the bam file
    echo "Indexing bam file..."
    samtools index $outfile_path"/"$barcode".bam"

#### Step4: De-duplication 
    echo "Deduplicating reads..."
    java -Xmx32g -jar $picard_path"/picard.jar" MarkDuplicates I=$outfile_path"/"$barcode".bam" O=$outfile_path"/"$barcode"_rmdup.bam" M=$outfile_path"/"$barcode"_rmdup_metrics.txt" REMOVE_DUPLICATES=true BARCODE_TAG=RX

#### Step5: Remove lambda reads and keeo only the human reads ####
    echo "Separating human and lambda reads..."
    samtools index $outfile_path"/"$barcode"_rmdup.bam"
    samtools idxstats $outfile_path"/"$barcode"_rmdup.bam" | cut -f 1 | grep -v lambda | xargs samtools view -b $outfile_path"/"$barcode"_rmdup.bam" > $outfile_path"/"$barcode"_human.bam"
    samtools view -b $outfile_path"/"$barcode"_rmdup.bam" lambda > $outfile_path"/"$barcode"_lambda.bam"

#### Step6: Get read names of only dipyrimidine reads ####
    echo "Extrcting readnames of  dipyrimidine reads..."
    cp $script_path"/build_signal_track2.R" $outfile_path"/build_signal_track2_temp_"$barcode".R" 
    echo "build_signal_track2(\""$outfile_path"/"$barcode"_human.bam\"", "\""$outfile_path"/"$barcode".bed\"" ,tmp_out= "\""$outfile_path"/"$barcode".tab\")" >> $outfile_path"/build_signal_track2_temp_"$barcode".R"
    Rscript $outfile_path"/build_signal_track2_temp_"$barcode".R"

#### Step7: Filter read by read names ####
    echo "Filtering non-dipyrimidine reads..."
    samtools index $outfile_path"/"$barcode"_human.bam"
    cat $outfile_path"/"$barcode".tab" | cut -d \, -f 1 > $outfile_path"/"$barcode"_readname.tab" 
    java -Xmx32g -jar $picard_path"/picard.jar" FilterSamReads I=$outfile_path"/"$barcode"_human.bam" O=$outfile_path"/"$barcode"_dipy.bam" READ_LIST_FILE=$outfile_path"/"$barcode"_readname.tab" FILTER=includeReadList

#### Step8: Subset reads for repenrich2 ####
    echo "Subsetting reads for Repenrich2..."
    python $repenrich2_path"/RepEnrich2_subset.py" $outfile_path"/"$barcode"_dipy.bam" 30 $outfile_path"/"$barcode --pairedend TRUE

#### Step9: Running repenrich2 ####
    echo "Running Repenrich2..."
    python $repenrich2_path"/RepEnrich2.py" $repmskr_path"/hg19_repeatmasker_clean.txt" $outfile_path"/"$barcode $barcode $repenrich2_setup $outfile_path"/"$barcode"_multimap_R1.fastq" --fastqfile2 $outfile_path"/"$barcode"_multimap_R2.fastq" $outfile_path"/"$barcode"_unique.bam" --cpus 16 --pairedend TRUE

done

