### MSc - A. Echchiki 
### Illumina data analysis


### GET RAW READS + QC ###

# copy from original to new location 
bsub -q dee-hugemem -J copy_reads "cp -r /home/jroux/archive/MinION/run_Illumina_2015_11_19 /scratch/beegfs/monthly/aechchik/MSc/illumina/reads/"
# move to dir 
cd /scratch/beegfs/monthly/aechchik/MSc/illumina/reads/
# rename reads 
for filename in *.gz; do mv "./$filename" "./${filename:34}"; done # remove non-informative filename
mkdir R1; mkdir R2 # create temp subdir
mv R1*.gz ./R1; mv R2*gz ./R2 # separate R1 from R2 reads
cd ./R1 # prepare to rename R1 reads in R1 subdir
for filename in *.gz; do mv "./$filename" "./${filename:3}"; done # remove info R1 before file number
for filename in *fastq.gz; do mv $filename $(basename $filename .fastq.gz)_R1.fastq.gz; done # add info R1 after file number
mv *.gz ../ # move renamed R1 to parent
cd ../R2 # prepare to rename R1 reads in R2 subdir
for filename in *.gz; do mv "./$filename" "./${filename:3}"; done # remove info R2 before file number
for filename in *fastq.gz; do mv $filename $(basename $filename .fastq.gz)_R2.fastq.gz; done # add info R2 after file number
mv *.gz ../ # move renamed R2 to parent
cd ../ # move to parent
rm -r ./R1/ ; rm -r ./R2/ # remove temp subdir

# quality control on raw reads 
# software: fastqc-0.11.2
#
# move to dir
cd /scratch/beegfs/monthly/aechchik/MSc/illumina/reads/
# submit job
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -J fastqc$a "module add UHTS/Quality_control/fastqc/0.11.2; fastqc /scratch/beegfs/monthly/aechchik/MSc/illumina/reads/$i -o /scratch/beegfs/monthly/aechchik/MSc/illumina/reads/fastqc"; done
# parameters:
# -o: specified output directory, directory must exist


### ADAPTER TRIMMING + QC ###

# note: 
# adapter sequnece was inferred from LGTF report, mentioning index CGTACG. From http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences_1000000002694-01.pdf, this corresponds to TruSeq Adapter, Index 22. Final adapter seq: GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTT.

# adapter trimming only
# software: cutadapt-1.8
#
# move to dir
cd /scratch/beegfs/monthly/aechchik/MSc/illumina/reads/
# submit job
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -J cutadapt_trim_$a "module add UHTS/Quality_control/cutadapt/1.8; cutadapt $i -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTT -o /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/unfiltered/$i"; done
# parameters:
# -a:  sequence of an adapter that was ligated to the 3' end; the adapter itself and anything that follows is trimmed. 
# -o: specified output file

# quality control on trimmed reads
# software: fastqc-0.11.2
#
# move to dir
cd /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/unfiltered/
# submit job
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -J fastqc$a "module add UHTS/Quality_control/fastqc/0.11.2; fastqc /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/unfiltered/$i -o /scratch/beegfs/monthly/aechchik/MSc/illumina/reads/cutadapt/unfiltered/fastqc"; done
# parameters:
# -o: specified output directory, directory must exist

# adapter trimming, read filtering
# software: cutadapt-1.8
#
# move to dir
cd /scratch/beegfs/monthly/aechchik/MSc/illumina/reads/
# submit job
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -J cutadapt_trimfilter_$a "module add UHTS/Quality_control/cutadapt/1.8; cutadapt $i -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTT -o /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/filtered/$i -m 36 --too-short-output /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/filtered/discarded/$i.fastq "; done
# parameters:
# -a:  sequence of an adapter that was ligated to the 3' end; the adapter itself and anything that follows is trimmed. 
# -o: specified output file
# -m: define min length for the read to be kept. 36 was calculated: 100 (read length) - 64 (full adapter length)
# --too-short-output: specified output file for discarded reads
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTT
# quality control on trimmed and filtered reads
# software: fastqc-0.11.2
#
# move to dir
cd /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/filtered/
# submit job
a=0; for i in $(ls *.fastq.gz); do echo $i; a=$((a+1)); bsub -q dee-hugemem -J fastqc$a "module add UHTS/Quality_control/fastqc/0.11.2; fastqc /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/filtered/$i -o /scratch/beegfs/monthly/aechchik/MSc/illumina/reads/cutadapt/filtered/fastqc"; done
# parameters:
# -o: specified output directory, directory must exist


### ALIGNMENT ###

# reference genome and annotation download from Ensembl v84
#
cd /scratch/beegfs/monthly/aechchik/MSc/dmel_files/dmel_bdgp6_fa
# download Dmel genome
wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-31/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.31.dna.genome.fa.gz
# unzip genome
gunzip Drosophila_melanogaster.BDGP6.31.dna.genome.fa.gz
#
cd /scratch/beegfs/monthly/aechchik/MSc/dmel_files/dmel_bdgp6_gtf
# download Dmel annotation
wget ftp://ftp.ensembl.org/pub/release-84/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.84.gtf.gz
# unzip annotation
gunzip Drosophila_melanogaster.BDGP6.84.gtf.gz


# hisat2: idx generation
# software: hisat-2.0.2
# 
# index ref genome 
# move to dir
cd /scratch/beegfs/monthly/aechchik/MSc/dmel_files/hisat2_idx/new_idx
# submit job
bsub -J hisat2_idx "module add UHTS/Aligner/hisat/2.0.2; hisat2-build /scratch/beegfs/monthly/aechchik/MSc/dmel_files/dmel_bdgp6_fa/Drosophila_melanogaster.BDGP6.31.dna.genome.fa  /scratch/beegfs/monthly/aechchik/MSc/dmel_files/hisat2_idx/new_idx/hisat2idx" # input: fasta genome; output: define the basename of idx files 

# hisat2: alignment
# software: hisat-2.0.2
# 
# part 1: use prebuilt indexes, built using annotated transcripts
# download prebuilt indexes for Dmel, Ensembl BDGP6 
# note: seems working better than --known-splicesites-infile, built using gtf (man check)
# move to dir 
cd /scratch/beegfs/monthly/aechchik/MSc/dmel_files/hisat2_idx/prebuilt_idx
# wget from ftp: D. melanogaster, Ensembl BDGP6 
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/bdgp6_tran.tar.gz 
# unzip
bsub -q dee-hugemem -J unzip_hisat2idx "tar -xvzf bdgp6_tran.tar.gz"
# submit job
# note: use UNfiltered cutadapt trimmed reads (error message on filtered cutadapt: LSF)
bsub -q dee-hugemem -J hisat2_prebuiltIDX "module add UHTS/Aligner/hisat/2.0.2; hisat2 -q -x /scratch/beegfs/monthly/aechchik/MSc/dmel_files/hisat2_idx/prebuilt_idx/bdgp6_tran/genome_tran -1 `ls /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/unfiltered/*_R1* | tr '\n' ','` -2 `ls /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/unfiltered/*_R2* | tr '\n' ','` > /scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/hisat2/prebuilt_idx/prebuilt_hisat2.sam"
# parameters: 
# -q: quantification mode (alignemnt)
# -x: dir where to read idx 
# -1: list R1 reads comma sep
# -2: list R2 reads comma sep 
#
# part 2: use new indexes, built using hisat2-build function
# submit job 
# note: use UNfiltered cutadapt trimmed reads (error message on filtered cutadapt: LSF)
bsub -q dee-hugemem -J hisat2_denovoIDX "module add UHTS/Aligner/hisat/2.0.2; hisat2 -q -x /scratch/beegfs/monthly/aechchik/MSc/dmel_files/hisat2_idx/denovo_idx/hisat2idx -1 `ls /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/unfiltered/*_R1* | tr '\n' ','` -2 `ls /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/unfiltered/*_R2* | tr '\n' ','` > /scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/hisat2/denovo_idx/denovo_hisat2.sam"
# parameters: 
# -q: quantification mode (alignemnt)
# -x: dir where to read idx 
# -1: list R1 reads comma sep
# -2: list R2 reads comma sep 

# star: idx generation 
# software: STAR/2.5.0b
#
# move to dir
cd /scratch/beegfs/monthly/aechchik/MSc/dmel_files/
# 
# part 1: use info in gtf
bsub -M 8388608 -q dee-hugemem -J star_GTFidx "module add UHTS/Aligner/STAR/2.5.0b; STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /scratch/beegfs/monthly/aechchik/MSc/dmel_files/star_idx/gtf_idx --genomeFastaFiles /scratch/beegfs/monthly/aechchik/MSc/dmel_files/dmel_bdgp6_fa/Drosophila_melanogaster.BDGP6.31.dna.genome.fa --sjdbGTFfile /scratch/beegfs/monthly/aechchik/MSc/dmel_files/dmel_bdgp6_gtf/Drosophila_melanogaster.BDGP6.84.gtf"
# TODO: describe pars
#
# part 2: don't use info in gtf
bsub -M 8388608 -q dee-hugemem -J star_noGTFidx "module add UHTS/Aligner/STAR/2.5.0b; STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /scratch/beegfs/monthly/aechchik/MSc/dmel_files/star_idx/nogtf_idx --genomeFastaFiles /scratch/beegfs/monthly/aechchik/MSc/dmel_files/dmel_bdgp6_fa/Drosophila_melanogaster.BDGP6.31.dna.genome.fa"
# TODO: describe pars

# star: alignment 
# software: STAR/2.5.0b
#
# part 1: use annotations in gtf (refer to idx built on gtf)
# move to dir
# note: use filtered cutadapt trimmed reads (error message on UNfiltered cutadapt: LSF)
cd /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/filtered/
# submit job 
# note: memory!
a=0; for i in {001..070}; do echo $i; a=$((a+1)); R1='_R1.fastq.gz'; R2='_R2.fastq.gz'; bsub -M 8388608 -q dee-hugemem -J star_wGTF$a -N "module add UHTS/Aligner/STAR/2.5.0b; STAR --outFileNamePrefix /scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/star/star_gtf/$i --runThreadN 8 --genomeDir /scratch/beegfs/monthly/aechchik/MSc/dmel_files/star_idx/gtf_idx --readFilesCommand zcat --readFilesIn /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/$i$R1 /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/$i$R2 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40"  ; done
# parameters:
# --outFileNamePrefix: path to output 
# --runThreadN: n. CPU to use
# --genomeDir: dir to idx files
# --readFilesCommand: specify input reads are .gz
# --readFilesIn: path/to/pairedR1 path/to/pairedR2
# --genomeFastaFiles: dir to fasta
# --outFilterScoreMinOverLread: outFilterScoreMin normalized to read length 
# --outFilterMatchNminOverLread: outFilterMatchNmin normalized to read length
# --outFilterMatchNmin: alignment will be output only if the number of matched bases is higher than this value
#
# part 2: do not use annotations in gtf (refer to idx built without gtf)
# move to dir
# note: use filtered cutadapt trimmed reads (error message on UNfiltered cutadapt: LSF)
cd /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/filtered/
# submit job 
# note: memory!
a=0; for i in {001..070}; do echo $i; a=$((a+1)); R1='_R1.fastq.gz'; R2='_R2.fastq.gz'; bsub -M 8388608 -q dee-hugemem -J star_noGTF$a -N "module add UHTS/Aligner/STAR/2.5.0b; STAR --outFileNamePrefix /scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/star/star_nogtf/$i --runThreadN 8 --genomeDir /scratch/beegfs/monthly/aechchik/MSc/dmel_files/star_idx/nogtf_idx --readFilesCommand zcat --readFilesIn /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/$i$R1 /scratch/beegfs/monthly/aechchik/MSc/illumina/cutadapt/$i$R2 --outFilterMatchNmin 90" ; done


# parameters:
# --outFileNamePrefix: path to output 
# --runThreadN: n. CPU to use
# --genomeDir: dir to idx files
# --readFilesCommand: specify input reads are .gz
# --readFilesIn: path/to/pairedR1 path/to/pairedR2
# --genomeFastaFiles: dir to fasta
# --outFilterScoreMinOverLread: outFilterScoreMin normalized to read length 
# --outFilterMatchNminOverLread: outFilterMatchNmin normalized to read length
# --outFilterMatchNmin: alignment will be output only if the number of matched bases is higher than this value

# convert sam to bam 
cd /scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/star/star_gtf/
a=0; for i in $(ls *.sam); do echo $i; a=$(echo $i | cut -d'.' -f1 | cut -d'_' -f3); echo $a; b=$((b + 1)); bsub -q dee-hugemem -J samtobam$b "module add UHTS/Analysis/samtools/1.3; samtools view -Sb $i > $a'.bam'"; done
/scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/star/star_nogtf/
a=0; for i in $(ls *.sam); do echo $i; a=$(echo $i | cut -d'.' -f1 | cut -d'_' -f3); echo $a; b=$((b + 1)); bsub -q dee-hugemem -J samtobam$b "module add UHTS/Analysis/samtools/1.3; samtools view -Sb $i > $a'.bam'"; done

# sort bam files
cd /scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/star/star_gtf/
a=0; for i in $(ls *.bam); do echo $i; a=$((a+1)); bsub -q dee-hugemem -J sortbam$a "module add UHTS/Analysis/samtools/1.3; samtools sort $i -o $i.sorted.bam"; done
/scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/star/star_nogtf/
a=0; for i in $(ls *.bam); do echo $i; a=$((a+1)); bsub -q dee-hugemem -J sortbam$a "module add UHTS/Analysis/samtools/1.3; samtools sort $i -o $i.sorted.bam"; done

# merge sorted bam 
bsub -q dee-hugemem -J mergebam "module add UHTS/Analysis/samtools/1.3; samtools merge merged_sorted.bam *sorted.bam"; done
# TOO BIG FOR CUFFLINKS!


### TRANSCRIPT ASSEMBLY ###

#TODO: complete with parameters

#CUFFLINKS REFGUIDED bam star gtf
cd /scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/star/star_gtf/ ;
a=0; for i in $(ls *.sorted*); do echo $i; a=$((a+1)); bsub  -M 20971520 -q dee-hugemem -J cufflinks_refguided_star_gtf$a "module add UHTS/Assembler/cufflinks/2.2.1; cufflinks -p 8 -o /scratch/beegfs/monthly/aechchik/MSc/illumina/transcript_assembly/cufflinks/cufflinks_refguided/bam_from_star/star_gtf/$a -g /scratch/beegfs/monthly/aechchik/MSc/dmel_files/dmel_bdgp6_gtf/Drosophila_melanogaster.BDGP6.84.gtf /scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/star/star_gtf/$i"; done

#CUFFLINKS DENOVO bam star gtf
cd /scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/star/star_gtf/ ;
a=0; for i in $(ls *.sorted*); do echo $i; a=$((a+1)); bsub -M 20971520 -q dee-hugemem -J cufflinks_denovo_star_gtf$a "module add UHTS/Assembler/cufflinks/2.2.1; cufflinks -p 8 -o /scratch/beegfs/monthly/aechchik/MSc/illumina/transcript_assembly/cufflinks/cufflinks_denovo/bam_from_star/star_gtf/$a /scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/star/star_gtf/$i"; done

#CUFFLINKS REFGUIDED bam star nogtf
cd /scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/star/star_nogtf/ ;
a=0; for i in $(ls *.sorted*); do echo $i; a=$((a+1)); bsub -M 20971520 -q dee-hugemem -J cufflinks_refguided_star_nogtf$a "module add UHTS/Assembler/cufflinks/2.2.1; cufflinks -p 8 -o /scratch/beegfs/monthly/aechchik/MSc/illumina/transcript_assembly/cufflinks/cufflinks_refguided/bam_from_star/star_nogtf/$a -g /scratch/beegfs/monthly/aechchik/MSc/dmel_files/dmel_bdgp6_gtf/Drosophila_melanogaster.BDGP6.84.gtf /scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/star/star_nogtf/$i"; done

#CUFFLINKS DENOVO bam star nogtf
cd /scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/star/star_nogtf/ ;
a=0; for i in $(ls *.sorted*); do echo $i; a=$((a+1)); bsub -M 20971520 -q dee-hugemem -J cufflinks_denovo_star_nogtf$a "module add UHTS/Assembler/cufflinks/2.2.1; cufflinks -p 8 -o /scratch/beegfs/monthly/aechchik/MSc/illumina/transcript_assembly/cufflinks/cufflinks_denovo/bam_from_star/star_nogtf/$a /scratch/beegfs/monthly/aechchik/MSc/illumina/alignment/star/star_nogtf/$i"; done

# TODO: repeat for hisat2 in split mode




# TODO: subsection 
#TODO: repeat on filtered/unfiltered cutadapt

# load Cufflinks 
module add UHTS/Assembler/cufflinks/2.2.1;
module add UHTS/Analysis/kallisto/0.42.4;



## cufflinks & kallisto 

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

# aim: transcript quantification by kallisto on cutadapt reads
# software: kallisto/0.42.4
# 
# sumbit job
bsub -J kallisto_quant "module add UHTS/Analysis/kallisto/0.42.4; kallisto quant -i /scratch/beegfs/monthly/aechchik/MSc/dmel_files/Drosophila_melanogaster.BDGP6.31.dna.transcripts.idx -o /scratch/beegfs/monthly/aechchik/MSc/MSc_Illumina/Ill_kallisto -b 100 /scratch/beegfs/monthly/aechchik/MSc/MSc_Illumina/Ill_trimmed/cutadapt/*.fastq.gz"
