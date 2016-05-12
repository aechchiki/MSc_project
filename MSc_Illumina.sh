### MSc - A. Echchiki 
### Illumina data analysis


### sequencing data repository
ls /home/jroux/archive/MinION/run_Illumina_2015_11_19 # paired-end reads


### reference genome and annotation: Ensembl v84

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


### raw data quality check

# aim: quality check
# software: fastqc-0.11.2
# 
# move to datadir
cd /home/jroux/archive/MinION/run_Illumina_2015_11_19
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -J fastqc$a "module add UHTS/Quality_control/fastqc/0.11.2; fastqc -t 2 -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_fastqc_raw /home/jroux/archive/MinION/run_Illumina_2015_11_19/$i"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_fastqc_raw


### adapter removal 

# aim: adapter removal 
# software: cutadapt/1.8
#
# move to datadir
cd /home/jroux/archive/MinION/run_Illumina_2015_11_19
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -J trimraw$a "module add UHTS/Quality_control/cutadapt/1.8; cutadapt -a TCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/$i /home/jroux/archive/MinION/run_Illumina_2015_11_19/$i"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt

# aim: adapter removal 
# software: trimmomatic/0.33
#
# move to datadir
cd /home/jroux/archive/MinION/run_Illumina_2015_11_19; 
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -M 20971520 -L /bin/bash -J trimmomatic$a "module add UHTS/Analysis/trimmomatic/0.33; trimmomatic SE -threads 8 /home/jroux/archive/MinION/run_Illumina_2015_11_19/$i /scratch/cluster/monthly/aechchik/MSc/illumina/trimmed/trimmomatic$i ILLUMINACLIP:/scratch/cluster/monthly/aechchik/MSc/illumina/adapter/Illumina_Index22.txt:3:25:6 MINLEN:60"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/trimmomatic


### trimmed data quality check 

# aim: quality check on cutadapt output
# software: fastqc-0.11.2
# 
# move to datadir
cd /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -L /bin/bash -J cutadaptfastqc$a -N "export PATH=/software/bin:$PATH; module add UHTS/Quality_control/fastqc/0.11.2; fastqc -t 2 -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/fastqc /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/$i"; done
# output in /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/fastqc

# aim: quality check on trimmomatic output
# software: fastqc-0.11.2
# 
# move to datadir
cd /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/trimmomatic
# submit job 
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -L /bin/bash -J trimmomaticfastqc$a -N "export PATH=/software/bin:$PATH; module add UHTS/Quality_control/fastqc/0.11.2; fastqc -t 2 -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/trimmomatic/fastqc /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/trimmomatic/$i"; done


### mapping 

## kallisto 

# aim: genome indexing by kallisto
# software: kallisto/0.42.4
# 
# submit job 
bsub -q dee-hugemem -L /bin/bash -J index_dmel -N "module add UHTS/Analysis/kallisto/0.42.4; kallisto index -i Drosophila_melanogaster.BDGP6.31.dna.transcripts.idx Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa"

# aim: transcript quantification by kallisto on cutadapt reads
# software: kallisto/0.42.4
# 
# create dir 
mkdir /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_kallisto
# re-name files to match pairs, to be recognized by kallisto
cd /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/ ; for filename in *.gz; do mv "./$filename" "./${filename:25}";  done # keep R1/R2_num.fastq.gz
mkdir R1; mkdir R2
mv R1*.gz ./R1; mv R2*gz ./R2
cd ./R1; for filename in *.gz; do mv "./$filename" "./${filename:3}";  done #remove info R1
for file in *fastq.gz; do mv $file $(basename $file .fastq.gz)_R1.fastq.gz; done # add info R1 before extension
mv *.gz ../ # move to original directory
cd ../R2; for filename in *.gz; do mv "./$filename" "./${filename:3}";  done # remove info R2
for file in *fastq.gz; do mv $file $(basename $file .fastq.gz)_R2.fastq.gz; done # add info R2 before extension
mv *.gz ../ # move to original directory
ls  /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/ # files are ordered num1_R1.fastq.gz num1_R2.fastq.gz num2_R1.fastq.gz num2_R2.fastq.gz ...
#
# sumbit job
bsub -J kallisto_quant "module add UHTS/Analysis/kallisto/0.42.4; kallisto quant -i /scratch/cluster/monthly/aechchik/MSc/dmel_files/Drosophila_melanogaster.BDGP6.31.dna.transcripts.idx -o /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_kallisto -b 100 /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/*.fastq.gz"

# aim: transcript quantification by kallisto on trimmomatic reads
# software: kallisto/0.42.4
# 
# create dir 
mkdir /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_kallisto
# re-name files to match pairs, to be recognized by kallisto
cd /scratch/cluster/monthly/aechchik/MSc/illumina/trimmed/trimmomatic ; for filename in *.gz; do mv "./$filename" "./${filename:25}";  done # keep R1/R2_num.fastq.gz
mkdir R1; mkdir R2
mv R1*.gz ./R1; mv R2*gz ./R2
cd ./R1; for filename in *.gz; do mv "./$filename" "./${filename:3}";  done #remove info R1
for file in *fastq.gz; do mv $file $(basename $file .fastq.gz)_R1.fastq.gz; done # add info R1 before extension
mv *.gz ../ # move to original directory
cd ../R2; for filename in *.gz; do mv "./$filename" "./${filename:3}";  done # remove info R2
for file in *fastq.gz; do mv $file $(basename $file .fastq.gz)_R2.fastq.gz; done # add info R2 before extension
mv *.gz ../ # move to original directory
ls  /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/ # files are ordered num1_R1.fastq.gz num1_R2.fastq.gz num2_R1.fastq.gz num2_R2.fastq.gz ...
#
# submit job
bsub -J kallisto_quant_trimmo "module add UHTS/Analysis/kallisto/0.42.4; kallisto quant -i /scratch/cluster/monthly/aechchik/MSc/dmel_files/Drosophila_melanogaster.BDGP6.31.dna.transcripts.idx -o /scratch/cluster/monthly/aechchik/MSc/illumina/kallisto/trimmomatic -b 100 /scratch/cluster/monthly/aechchik/MSc/illumina/trimmed/trimmomatic/*.fastq"


## hisat2

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
bsub -q dee-hugemem "module add UHTS/Aligner/hisat/2.0.2; hisat2_extract_splice_sites.py Drosophila_melanogaster.BDGP6.84.gtf > hisat2_splices.txt"
# note: Will then be useful to extract splice sites in mapping mode, using --known-splicesite-infile <path> : "with this mode, you can provide a list of known splice sites, which HISAT2 makes use of to align reads with small anchors. You can create such a list using python hisat2_extract_splice_sites.py genes.gtf > splicesites.txt, where hisat2_extract_splice_sites.py is included in the HISAT2 package, genes.gtf is a gene annotation file, and splicesites.txt is a list of splice sites with which you provide HISAT2 in this mode.
#  Note that it is better to use indexes built using annotated transcripts (such as genome_tran or genome_snp_tran), which works better than using this option. 

# aim: download indexes built using annotated transcripts
# wget 
# 
cd /scratch/cluster/monthly/aechchik/MSc/dmel_files
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/bdgp6_tran.tar.gz # D. melanogaster, Ensembl BDGP6 
bsub -q dee-hugemem -J unzip_hisat2idx "tar -xvzf bdgp6_tran.tar.gz"

# aim: hisat2 mapping 
# software: hisat/2.0.2
# 
# note: hisat2 works with unzipped files 
gunzip /scratch/cluster/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/*.fastq.gz
# submit job 
bsub -q dee-hugemem -J hisat2_gtf_cutadapt "module add UHTS/Aligner/hisat/2.0.2; hisat2 -q -x /scratch/cluster/monthly/aechchik/MSc/dmel_files/bdgp6_tran/genome_tran -1 `ls /scratch/cluster/monthly/aechchik/MSc/illumina/trimmed/cutadapt/reads_t/*_R1* | tr '\n' ','` -2 `ls /scratch/cluster/monthly/aechchik/MSc/illumina/trimmed/cutadapt/reads_t/*_R2* | tr '\n' ','` > /scratch/cluster/monthly/aechchik/MSc/illumina/hisat2/dmel_hisat2.sam"

# note: hisat2 works with unzipped files 
gunzip /scratch/cluster/monthly/aechchik/MSc/illumina/trimmed/trimmomatic/*.fastq.gz
# submit job 
bsub -q dee-hugemem -J hisat2_gtf_trim " module add UHTS/Aligner/hisat/2.0.2; module add UHTS/Analysis/samtools/1.3; hisat2 -q -x /scratch/cluster/monthly/aechchik/MSc/dmel_files/bdgp6_tran/genome_tran -1 `ls /scratch/cluster/monthly/aechchik/MSc/illumina/trimmed/trimmomatic/*_R1* | tr '\n' ','` -2 `ls /scratch/cluster/monthly/aechchik/MSc/illumina/trimmed/trimmomatic/*_R2* | tr '\n' ','` > /scratch/cluster/monthly/aechchik/MSc/illumina/hisat2/trimmomatic/hisat2_trimmomatic.sam"


### star

# aim: genome indexing by star
# software: STAR/2.5.0b
#
# move to dir
cd /scratch/cluster/monthly/aechchik/MSc/dmel_files/
# submit job 
bsub -M 8388608 -q dee-hugemem -J star_idx "module add UHTS/Aligner/STAR/2.5.0b; STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /scratch/cluster/monthly/aechchik/MSc/dmel_files/star_idx --genomeFastaFiles /scratch/cluster/monthly/aechchik/MSc/dmel_files/Drosophila_melanogaster.BDGP6.31.dna.genome.fa --sjdbGTFfile /scratch/cluster/monthly/aechchik/MSc/dmel_files/Drosophila_melanogaster.BDGP6.84.gtf"

#note: mapping of cutadapt reads was not successful 

# aim: map reads vs ref (trimmomatic)
# software: STAR/2.5.0b
#
# submit job
a=0; for i in {001..070}; do echo $i; a=$((a+1)); R1='_R1.fastq'; R2='_R2.fastq'; bsub -M 8388608 -q dee-hugemem -J starmap_$a -N "module add UHTS/Aligner/STAR/2.5.0b; STAR --outFileNamePrefix /scratch/cluster/monthly/aechchik/MSc/illumina/star/trimmomatic/ --runThreadN 8 --genomeDir /scratch/cluster/monthly/aechchik/MSc/dmel_files/star_idx --readFilesIn /scratch/cluster/monthly/aechchik/MSc/illumina/trimmed/trimmomatic/$i$R1 /scratch/cluster/monthly/aechchik/MSc/illumina/trimmed/trimmomatic/$i$R2" ; done



#####


## not sure

# aim: extract transcripts 
# software: cufflinks/2.2.1 
# 
# submit job 
bsub -q dee-hugemem -L /bin/bash -J dmel_transcripts -N "module add UHTS/Assembler/cufflinks/2.2.1; gffread -g Drosophila_melanogaster.BDGP6.31.dna.genome.fa -x Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa Drosophila_melanogaster.BDGP6.84.gtf"
# how many transcripts
awk '/>/{ print }' Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa | wc -l # 30353
