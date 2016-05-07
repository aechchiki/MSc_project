# MSc - A. Echchiki 
# Illumina data analysis

# data repository
ls /home/jroux/archive/MinION/run_Illumina_2015_11_19 # paired-end reads

# aim: quality check
# software: fastqc-0.11.2
# 
# move to datadir
cd /home/jroux/archive/MinION/run_Illumina_2015_11_19
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -L /bin/bash -J fastqc$a -N "export PATH=/software/bin:$PATH; module add UHTS/Quality_control/fastqc/0.11.2; fastqc -t 2 -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_fastqc_raw /home/jroux/archive/MinION/run_Illumina_2015_11_19/$i"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_fastqc_raw

# aim: adapter removal 
# software: cutadapt/1.8
#
# move to datadir
cd /home/jroux/archive/MinION/run_Illumina_2015_11_19
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -L /bin/bash -J trimraw$a -N "export PATH=/software/bin:$PATH; module add UHTS/Quality_control/cutadapt/1.8; cutadapt -a TCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/$i /home/jroux/archive/MinION/run_Illumina_2015_11_19/$i"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt

# aim: quality check on cutadapt output
# software: fastqc-0.11.2
# 
# move to datadir
cd /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -L /bin/bash -J cutadaptfastqc$a -N "export PATH=/software/bin:$PATH; module add UHTS/Quality_control/fastqc/0.11.2; fastqc -t 2 -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/fastqc /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/$i"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/fastqc

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
# software: cufflinks/2.2.1 
# 
# submit job 
bsub -q dee-hugemem -L /bin/bash -J dmel_transcripts -N "module add UHTS/Assembler/cufflinks/2.2.1; gffread -g Drosophila_melanogaster.BDGP6.31.dna.genome.fa -x Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa Drosophila_melanogaster.BDGP6.84.gtf"
# how many transcripts
awk '/>/{ print }' Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa | wc -l # 30353

# aim: genome indexing by kallisto
# software: kallisto/0.42.4
# 
# submit job 
bsub -q dee-hugemem -L /bin/bash -J index_dmel -N "module add UHTS/Analysis/kallisto/0.42.4; kallisto index -i Drosophila_melanogaster.BDGP6.31.dna.transcripts.idx Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa"

# aim: transcript quantification by kallisto
# software: kallisto/0.42.4
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
#
# sumbit job
bsub -J kallisto_quant "module add UHTS/Analysis/kallisto/0.42.4; kallisto quant -i /scratch/cluster/monthly/aechchik/MSc/dmel_files/Drosophila_melanogaster.BDGP6.31.dna.transcripts.idx -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_kallisto -b 100 /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/*.fastq.gz"

# note: hisat2 works with unzipped files 
gunzip /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/*.fastq.gz

# aim: genome indexing by hisat2
# software: hisat/2.0.2
#
# move to dir
cd /scratch/cluster/monthly/aechchik/MSc/dmel_files/
# submit job 
bsub "module add UHTS/Aligner/hisat/2.0.2; hisat2-build Drosophila_melanogaster.BDGP6.31.dna.genome.fa"

# aim: extract splice sites from gtf for hisat2
# software: hisat/2.0.2
#
# move to dir
cd /scratch/cluster/monthly/aechchik/MSc/dmel_files/
# submit job 
bsub -q dee-hugemem "module add UHTS/Aligner/hisat/2.0.2; extract_splice_sites.py Drosophila_melanogaster.BDGP6.31.dna.genome.fa"

# aim: map reads vs ref (no known gtf as input)
# software: hisat/2.0.2; samtools/1.3
#
# move to dir
cd /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt
# submit job 
bsub -q dee-hugemem -J hisat2_nogtf "module add UHTS/Aligner/hisat/2.0.2; module add UHTS/Analysis/samtools/1.3; hisat2 -x /scratch/cluster/monthly/aechchik/MSc/dmel_files/bt2_index.idx -1 `ls *_R1* | tr '\n' ','` -2 `ls *_R2* | tr '\n' ','` | samtools view -bS > hisat.bam"

# aim: map reads vs ref (known gtf as input)
# software: hisat/2.0.2
#
# move to dir
cd /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt
# submit job 
bsub -q dee-hugemem -J hisat2_gtf "module add UHTS/Aligner/hisat/2.0.2; module add UHTS/Analysis/samtools/1.3; hisat2 -x /scratch/cluster/monthly/aechchik/MSc/dmel_files/bt2_index.idx -1 `ls *_R1* | tr '\n' ','` -2 `ls *_R2* | tr '\n' ','` --known-splicesite-infile /scratch/cluster/monthly/aechchik/MSc/dmel_files/Drosophila_melanogaster.BDGP6.splices.txt  | samtools view -bS  > hisat_known-splicesite.bam"

# aim: genome indexing by star
# software: STAR/2.5.0b
#
# move to dir
cd /scratch/cluster/monthly/aechchik/MSc/dmel_files/
# submit job 
bsub -M 8388608 -q dee-hugemem -J star_idx "module add UHTS/Aligner/STAR/2.5.0b; STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /scratch/cluster/monthly/aechchik/MSc/dmel_files/ --genomeFastaFiles Drosophila_melanogaster.BDGP6.31.dna.genome.fa --sjdbGTFfile Drosophila_melanogaster.BDGP6.84.gtf"

# aim: map reads vs ref
# software: STAR/2.5.0b
#
# submit job 
bsub -q dee-hugemem -J star_map "module add UHTS/Aligner/STAR/2.5.0b; STAR --runThreadN 8 --genomeDir /scratch/cluster/monthly/aechchik/MSc/dmel_files/ --readFilesIn /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/Reads1.fastq /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/Reads2.fastq"
