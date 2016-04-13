# MSc - A. Echchiki 
# Illumina data analysis

# data repository
ls /home/jroux/archive/MinION/run_Illumina_2015_11_19 # paired-end reads

# aim: quality check
# program: fastqc-0.11.2
# 
# move to datadir
cd /home/jroux/archive/MinION/run_Illumina_2015_11_19
# run fastqc on all .fastq.gz 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -L /bin/bash -J fastqc$a -N "export PATH=/software/bin:$PATH; module add UHTS/Quality_control/fastqc/0.11.2; fastqc -t 2 -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/ /home/jroux/archive/MinION/run_Illumina_2015_11_19/$i"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_fastqc_raw

# aim: adapter removal 
# program: cutadapt-1.8
# move to datadir
cd /home/jroux/archive/MinION/run_Illumina_2015_11_19
# run cutadapt on all .fastq.gz
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -L /bin/bash -J trimraw$a -N "export PATH=/software/bin:$PATH; module add UHTS/Quality_control/cutadapt/1.8; cutadapt -a TCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/$i /home/jroux/archive/MinION/run_Illumina_2015_11_19/$i"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed

