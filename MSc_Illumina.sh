# MSc - A. Echchiki 
# Illumina data analysis

# data repository
ls /home/jroux/archive/MinION/run_Illumina_2015_11_19 # paired-end reads

# aim: quality check
# program: fastqc-0.11.2
# 
# move to datadir
cd /home/jroux/archive/MinION/run_Illumina_2015_11_19
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -L /bin/bash -J fastqc$a -N "export PATH=/software/bin:$PATH; module add UHTS/Quality_control/fastqc/0.11.2; fastqc -t 2 -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_fastqc_raw /home/jroux/archive/MinION/run_Illumina_2015_11_19/$i"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_fastqc_raw

# aim: adapter removal 
# program: cutadapt/1.8
#
# move to datadir
cd /home/jroux/archive/MinION/run_Illumina_2015_11_19
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -L /bin/bash -J trimraw$a -N "export PATH=/software/bin:$PATH; module add UHTS/Quality_control/cutadapt/1.8; cutadapt -a TCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/$i /home/jroux/archive/MinION/run_Illumina_2015_11_19/$i"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt

# aim: adapter removal 
# program: trimmomatic/0.33
#
# move to datadir
cd /home/jroux/archive/MinION/run_Illumina_2015_11_19
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -M 20971520 -L /bin/bash -J trimmomatic$a -N "export PATH=/software/bin:$PATH; module add UHTS/Analysis/trimmomatic/0.33; trimmomatic SE -threads 8 /home/jroux/archive/MinION/run_Illumina_2015_11_19/$i /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/$i ILLUMINACLIP:/scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_adapter/Illumina_Index22.txt:3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:60"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/trimmomatic

# aim: quality check on cutadapt output
# program: fastqc-0.11.2
# 
# move to datadir
cd /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -L /bin/bash -J cutadaptfastqc$a -N "export PATH=/software/bin:$PATH; module add UHTS/Quality_control/fastqc/0.11.2; fastqc -t 2 -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/fastqc /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/$i"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/fastqc

# aim: quality check on trimmomatic output
# program: fastqc-0.11.2
# 
# move to datadir
cd /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/trimmomatic
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -L /bin/bash -J trimmomaticfastqc$a -N "export PATH=/software/bin:$PATH; module add UHTS/Quality_control/fastqc/0.11.2; fastqc -t 2 -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/trimmomatic/fastqc /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/trimmomatic/$i"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/trimmomatic/fastqc

# download .fasta & annotation 

# aim: extract cds 
# program: cufflinks/2.2.1 
# 
# submit job 
bsub -q dee-hugemem -L /bin/bash -J dmel_cds -N "module add UHTS/Assembler/cufflinks/2.2.1; gffread -g /scratch/cluster/monthly/aechchik/MSc/dmel_files/genome/dm6.fa -x dm6_cds.fa /scratch/cluster/monthly/aechchik/MSc/dmel_files/genome/dm6.gtf"
# how many cds
cd /scratch/cluster/monthly/aechchik/MSc/dmel_files/genome; awk '/>/{ print }'  dm6_cds.fa | wc -l # 30318

# aim: index the reference 
# program: /kallisto/0.42.4
# 
# submit job 
bsub -q dee-hugemem -L /bin/bash -J index_dmel -N "module add UHTS/Analysis/kallisto/0.42.4; kallisto index -i dm6_cds.idx dm6_cds.fa"

# aim: pseudoalignment of reads to ref
# program: /kallisto/0.42.4
# 
# submit job 
#
bsub -J kallisto_pseudobamSA4 "kallisto quant -i  Pseud_protegens_S5_cds.idx -o /scratch/cluster/monthly/mls_2015/aechchik/kallisto -b 100 --paired -l 365 -s 0 --pseudobam /scratch/cluster/monthly/mls_2015/aechchik/SA4_qc/SA4t.fastq.gz > SA4t.sam"





