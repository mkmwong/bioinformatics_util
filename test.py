### Purpose of this script: this is to use for when there are multiple      ###
### files in a directory requires same kind of processing. This script will ###
### be used to generaate the bash script required to run the job, as well   ###
### to submit the job the slurm scheduler.                                  ###
### option for alignment: STAR(RNASeq), bowtie2(DNASeq), bwa(DNASeq),       ### 

import os
import subprocess

### paths
INDEX_PATH="/home/groups/ashbym/mandy/combined_genome"
ALIGNER_PATH=""
F_PATH="/scratch/groups/ashbym/mandy/TRAD-Seq/"
SAMTOOLS_PATH=""
OUTPUT_PATH="/scratch/groups/ashbym/mandy/test_dir/"
### job details 
SCRIPT_PREFIX="STAR2run"
ALIGNER_TYPE="STAR"
ENDEDNESS="SINGLE"
F_NAME_SUFFIX="_R"
F_NAME_EXT=".fastq.gz"
N_THREAD="16"
CONDA_ENV="base"
OTHER_OPTS=""
### Slurm job schedule settings
JOB_NAME="TEST_BOWTIE2"
TIME="2:00:00"
MEM_PER_CPU="4G"
CPUS_PER_TASK="16"
MAIL_TYPE="ALL"
MAIL_USER="mkmwong@stanford.edu"
PARTITION="hns,normal"
f_names = ["ATGAGGCT"]

def write_header(f, job_name, time, mem_per_cpu, cpus_per_task, mail_type, mail_user):
    f.write("#!/bin/bash\n")
    f.write("#SBATCH --job-name=" + job_name + "\n")
    f.write("#SBATCH --time=" + time + "\n")
    f.write("#SBATCH --mem-per-cpu=" + mem_per_cpu + "\n")
    f.write("#SBATCH --cpus-per-task=" + cpus_per_task + "\n")
    f.write("#SBATCH --mail-type=" + mail_type + "\n")
    f.write("#SBATCH --mail-user=" + mail_user + "\n")
    f.write("#SBATCH -p " + PARTITION + "\n")

   
def write_bowtie2(f, f_name):
    if ENDEDNESS=="PAIRED":
        f.write( ALIGNER_PATH + ALIGNER_TYPE + " -p " + N_THREAD +   " -x " + INDEX_PATH + \
" -1 " + F_PATH + f_name + F_NAME_SUFFIX + "1" + F_NAME_EXT + "  -2 " + F_PATH + f_name + F_NAME_SUFFIX + "2" + F_NAME_EXT +  \
  OTHER_OPTS + " | " + SAMTOOLS_PATH + "samtools view -bS > " + OUTPUT_PATH + f_name + ".bam\n")
    elif ENDEDNESS=="SINGLE":
        f.write( ALIGNER_PATH + ALIGNER_TYPE + " -p " + N_THREAD + " -x " + INDEX_PATH + \
" -U " + F_PATH + f_name + F_NAME_EXT + OTHER_OPTS  + " | " + SAMTOOLS_PATH + "samtools view -bS > " \
 + OUTPUT_PATH + f_name + ".bam\n") 

def write_bwa(f, f_name):
    if ENDEDNESS=="PAIRED":
        f.write( ALIGNER_PATH + ALIGNER_TYPE + " mem  -t " + N_THREAD + " " + INDEX_PATH + ".fa " + \
F_PATH + f_name + F_NAME_SUFFIX + "1" + F_NAME_EXT + " " + F_PATH + f_name + F_NAME_SUFFIX + "2" + F_NAME_EXT + \
OTHER_OPTS + " | " + SAMTOOLS_PATH + "samtools view -bS > " + OUTPUT_PATH + f_name + ".bam\n"  )
    elif ENDEDNESS=="SINGLE":
        f.write( ALIGNER_PATH + ALIGNER_TYPE + " mem  -t " + N_THREAD + " " + INDEX_PATH + ".fa " + \
F_PATH + f_name + F_NAME_EXT +  \
OTHER_OPTS + " | " + SAMTOOLS_PATH + "samtools view -bS > " + OUTPUT_PATH + f_name + ".bam\n"  )

def write_STAR(f, f_name):
    if ENDEDNESS=="PAIRED" and "gz" in F_NAME_EXT:
        f.write( ALIGNER_PATH + ALIGNER_TYPE + " --runThreadN " + N_THREAD + " --genomeDir " + INDEX_PATH + \
" --readFilesCommand gunzip -c --readFilesIn " + F_PATH + f_name + F_NAME_SUFFIX + "1" + F_NAME_EXT + " " + \
F_PATH + f_name + F_NAME_SUFFIX + "2" + F_NAME_EXT + " --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --quantMode GeneCounts " + \
" --outFilterMultimapNmax 1 --outFilterMatchNmin 35 --outFileNamePrefix " + OUTPUT_PATH + f_name + ".bam\n" )
    elif ENDEDNESS=="PAIRED" and "gz" not in F_NAME_EXT:
        f.write( ALIGNER_PATH + ALIGNER_TYPE + " --runThreadN " + N_THREAD + " --genomeDir " + INDEX_PATH + \
" --readFilesIn " + F_PATH + f_name + F_NAME_EXT + " " + \
F_PATH + f_name + F_NAME_SUFFIX + "2" + F_NAME_EXT + " --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --quantMode GeneCounts " + \
" --outFilterMultimapNmax 1 --outFilterMatchNmin 35 --outFileNamePrefix " + OUTPUT_PATH + f_name + ".bam\n" )
    elif ENDEDNESS=="SINGLE" and "gz" in F_NAME_EXT: 
        f.write( ALIGNER_PATH + ALIGNER_TYPE + " --runThreadN " + N_THREAD + " --genomeDir " + INDEX_PATH + \
" --readFilesCommand gunzip -c --readFilesIn " + F_PATH + f_name + F_NAME_SUFFIX + "1" + F_NAME_EXT + \
" --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --quantMode GeneCounts " + 
" --outFilterMultimapNmax 1 --outFilterMatchNmin 35 --outFileNamePrefix " + OUTPUT_PATH + f_name + ".bam\n" )
    elif ENDEDNESS=="SINGLE" and "gz" not in F_NAME_EXT:
        f.write( ALIGNER_PATH + ALIGNER_TYPE + " --runThreadN " + N_THREAD + " --genomeDir " + INDEX_PATH + \
" --readFilesIn " + F_PATH + f_name + F_NAME_EXT + " " + \
" --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --quantMode GeneCounts " + 
" --outFilterMultimapNmax 1 --outFilterMatchNmin 35 --outFileNamePrefix " + OUTPUT_PATH + f_name + ".bam\n" )

for f_name in f_names:
    if not os.path.isdir(OUTPUT_PATH):
        os.makedirs(OUTPUT_PATH)
    script_name = SCRIPT_PREFIX + "_" + f_name + ".sh"
    f = open(script_name, "w") 
    write_header(f, JOB_NAME, TIME, MEM_PER_CPU, CPUS_PER_TASK, MAIL_TYPE, MAIL_USER)
    if CONDA_ENV:
        f.write("conda activate " + CONDA_ENV + "\n")
    if ALIGNER_TYPE == "bowtie2":
        write_bowtie2(f, f_name)
    elif ALIGNER_TYPE == "bwa":
        write_bwa(f, f_name)
    elif ALIGNER_TYPE == "STAR":
        write_STAR(f, f_name)
    elif ALIGNER_TYPE == "hisat2":
        write_hisat2(f, f_name)

    f.close()      
    os.chmod(script_name, 0o755)
    rc = subprocess.call("sbatch "+script_name, shell=True) 



