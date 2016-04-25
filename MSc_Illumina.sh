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

# aim: quality check on cutadapt output
# program: fastqc-0.11.2
# 
# move to datadir
cd /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -L /bin/bash -J cutadaptfastqc$a -N "export PATH=/software/bin:$PATH; module add UHTS/Quality_control/fastqc/0.11.2; fastqc -t 2 -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/fastqc /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/$i"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/fastqc

# aim: adapter removal 
# program: trimmomatic/0.33
#
# move to datadir
cd /home/jroux/archive/MinION/run_Illumina_2015_11_19
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -M 20971520 -L /bin/bash -J trimmomatic$a -N "export PATH=/software/bin:$PATH; module add UHTS/Analysis/trimmomatic/0.33; trimmomatic SE -threads 8 /home/jroux/archive/MinION/run_Illumina_2015_11_19/$i /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/$i ILLUMINACLIP:/scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_adapter/Illumina_Index22.txt:3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:60"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/trimmomatic

# aim: quality check on trimmomatic output
# program: fastqc-0.11.2
# 
# move to datadir
cd /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/trimmomatic
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -L /bin/bash -J trimmomaticfastqc$a -N "export PATH=/software/bin:$PATH; module add UHTS/Quality_control/fastqc/0.11.2; fastqc -t 2 -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/trimmomatic/fastqc /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/trimmomatic/$i"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/trimmomatic/fastqc

# aim: get reference genome and annotation
# 
# create dir for ref files 
mkdir /scratch/cluster/monthly/aechchik/MSc/dmel_files; cd /scratch/cluster/monthly/aechchik/MSc/dmel_files
# download Dmel genome
wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-31/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.31.dna.genome.fa.gz
# unzip genome
gunzip Drosophila_melanogaster.BDGP6.31.dna.genome.fa.gz
# download Dmel annotation
wget ftp://ftp.ensembl.org/pub/release-84/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.84.gtf.gz
# unzip annotation
gunzip Drosophila_melanogaster.BDGP6.84.gtf.gz

# aim: extract transcripts 
# program: cufflinks/2.2.1 
# 
# submit job 
bsub -q dee-hugemem -L /bin/bash -J dmel_transcripts -N "module add UHTS/Assembler/cufflinks/2.2.1; gffread -g Drosophila_melanogaster.BDGP6.31.dna.genome.fa -x Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa Drosophila_melanogaster.BDGP6.84.gtf"
# how many transcripts
awk '/>/{ print }' Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa | wc -l # 30353

# aim: indexing by kallisto
# program: /kallisto/0.42.4
# 
# submit job 
bsub -q dee-hugemem -L /bin/bash -J index_dmel -N "module add UHTS/Analysis/kallisto/0.42.4; kallisto index -i Drosophila_melanogaster.BDGP6.31.dna.transcripts.idx Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa"

# aim: quantification by kallisto
# program: /kallisto/0.42.4
# 
# create dir 
mkdir /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_kallisto
# re-name files to match pairs, to be recognized by kallisto
cd /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/ ; for filename in *.gz; do mv "./$filename" "./${filename:25}";  done
# keep R1/R2_num.fastq.gz
mkdir R1; mkdir R2
mv R1*.gz ./R1; mv R2*gz ./R2
cd ./R1
for filename in *.gz; do mv "./$filename" "./${filename:3}";  done #remove info R1
for file in *fastq.gz; do mv $file $(basename $file .fastq.gz)_R1.fastq.gz; done # add info R1 before extension
mv *.gz ../ # move to original directory
cd ../R2
for filename in *.gz; do mv "./$filename" "./${filename:3}";  done # remove info R2
for file in *fastq.gz; do mv $file $(basename $file .fastq.gz)_R2.fastq.gz; done # add info R2 before extension
mv *.gz ../ # move to original directory
ls  /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/ # files are ordered num1_R1.fastq.gz num1_R2.fastq.gz num2_R1.fastq.gz num2_R2.fastq.gz ...
# kallisto will consider them as pairs? hope so

# sumbit job
bsub -J kallisto_quant "module add UHTS/Analysis/kallisto/0.42.4; kallisto quant -i /scratch/cluster/monthly/aechchik/MSc/dmel_files/Drosophila_melanogaster.BDGP6.31.dna.transcripts.idx -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_kallisto -b 100 -l 243 -s 30 /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/*.fastq.gz"
