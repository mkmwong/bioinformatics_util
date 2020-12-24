### Purpose of this script: this is use for running MELT on all bam files   ###
### in a directory. 

import os
import subprocess

### paths
INDEX_PATH="/home/groups/ashbym/mandy/hg38_genome/grch38_1kgmaj.fa"
MELT_PATH="/home/groups/ashbym/mandy/MELTv2.2.0/MELT.jar"
F_PATH="/scratch/groups/ashbym/mandy/test"
OUTPUT_PATH="/scratch/groups/ashbym/mandy/test_dir/"
### job details 
GENE_PATH="/home/groups/ashbym/mandy/MELTv2.2.0/add_bed_files/Hg38/Hg38.genes.bed"
TRANSPO_F="/home/groups/ashbym/mandy/MELTv2.2.0/transposon_file_list.txt"
CONDA_ENV="dnanexus"
JAVA_MEM="Xmx32G"
### Slurm job schedule settings
JOB_NAME="TEST_MELT"
TIME="2:00:00"
MEM_PER_CPU="2G"
CPUS_PER_TASK="16"
MAIL_TYPE="ALL"
MAIL_USER="mkmwong@stanford.edu"
PARTITION="hns,normal"

def write_header(f, job_name, time, mem_per_cpu, cpus_per_task, mail_type, mail_user):
    f.write("#!/bin/bash\n")
    f.write("#SBATCH --job-name=" + job_name + "\n")
    f.write("#SBATCH --time=" + time + "\n")
    f.write("#SBATCH --mem-per-cpu=" + mem_per_cpu + "\n")
    f.write("#SBATCH --cpus-per-task=" + cpus_per_task + "\n")
    f.write("#SBATCH --mail-type=" + mail_type + "\n")
    f.write("#SBATCH --mail-user=" + mail_user + "\n")
    f.write("#SBATCH -p " + PARTITION + "\n")


f_names = [f for f in os.listdir(F_PATH) if f.endswith(".bam")]

for f_name in f_names:
    sub_f = f_name.split(".")[0]
    script_name = "MELT_" + sub_f + ".sh"
    f = open(script_name, "w")
    write_header(f, JOB_NAME, TIME, MEM_PER_CPU, CPUS_PER_TASK, MAIL_TYPE, MAIL_USER)
    # sort the bam file if it isn't sorted
    f.write("module load java\n")
    f.write("conda activate " + CONDA_ENV + "\n")
    f.write("if samtools view -H " + F_PATH + "/" + f_name + " | head -1 | grep 'SO:coordinate';\n")
    f.write("then\n    echo 'BAM file has been sorted. Skip sorting.';\n")
    f.write("    mv " + F_PATH + "/" + f_name + " " + F_PATH + "/" + sub_f + "_sorted.bam;\n" )
    f.write("else\n    echo 'BAM file has not been sorted. Now sorting...'; \n")
    f.write("    samtools sort -o " + F_PATH + "/" + sub_f + "_sorted.bam " + F_PATH + "/" + f_name + ";\nfi \n")
    # index the sorted bam file
    f.write("echo 'Indexing bam file...'\n")
    f.write("samtools index " + F_PATH + "/" + sub_f + "_sorted.bam \n")
    # preprocess for MELT
    f.write("echo 'Preprocessing for MELT...'\n")
    f.write("java -" + JAVA_MEM + " -jar " + MELT_PATH + " Preprocess \\\n" )
    f.write("    -bamfile " + F_PATH + "/" + sub_f + "_sorted.bam" + " \\\n    -h " + INDEX_PATH +"\n")
    # running MELT
    f.write("echo 'Running MELT...'\n")
    f.write("java -" + JAVA_MEM + " -jar " + MELT_PATH + " Single \\\n")
    f.write("    -bamfile " + F_PATH + "/" + sub_f + "_sorted.bam" + " \\\n    -h " + INDEX_PATH +" \\\n")
    f.write("    -n " + GENE_PATH + " \\\n")
    f.write("    -w " + F_PATH + "/" + sub_f  + " \\\n")
    f.write("    -t " + TRANSPO_F + " \\\n    -k -ac")
    f.close()
    os.chmod(script_name, 0o755)
    rc = subprocess.call("sbatch "+script_name, shell=True) 


