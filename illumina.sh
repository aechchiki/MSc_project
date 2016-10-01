## 0) original data

# original directory
hiseq_od=/home/jroux/archive/MinION/run_Illumina_2015_11_19
ls -l $hiseq_od | wc -l # 142 (71*2 paired-end)

# working directory
hiseq_raw=/scratch/beegfs/monthly/aechchik/MSc/illumina/reads/data/raw
mkdir -p $hiseq_raw/
# unzip and merge to R1 and R2
cat $hiseq_od/*R1*.fastq.gz > $hiseq_raw/R1_all.fastq.gz
cat $hiseq_od/*R2*.fastq.gz > $hiseq_raw/R2_all.fastq.gz

# fastqc on raw
hiseq_raw_qc=/scratch/beegfs/monthly/aechchik/MSc/illumina/reads/data/raw/qc
mkdir -p $hiseq_raw_qc/
# load fastqc
module add UHTS/Quality_control/fastqc/0.11.2
# run qc
fastqc $hiseq_raw/R1_all.fastq.gz -o $hiseq_raw_qc/
fastqc $hiseq_raw/R2_all.fastq.gz -o $hiseq_raw_qc/


## 1) trimming

# trimming directory
hiseq_trim=/scratch/beegfs/monthly/aechchik/MSc/illumina/reads/data/trimmed/
mkdir -p $hiseq_trim/

# get adapters file
cd $hiseq_trim && { curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip ; unzip Trimmomatic-0.33.zip; cd -; }
# keep TruSeq paired-end adapters, remove non-used files
mv $hiseq_trim/Trimmomatic-0.33/adapters/TruSeq3-PE.fa $hiseq_trim
rm -r $hiseq_trim/Trimmomatic*

# load trimmomatic
module add UHTS/Analysis/trimmomatic/0.33
# run trimming
trimmomatic PE -t 8 -phred33 -trimlog $hiseq_trim/TrimLog $hiseq_trim/R1_all.fastq.gz $hiseq_trim/R2_all.fastq.gz $hiseq_trim/R1_paired.fastq.gz $hiseq_trim/R1_unpaired.fastq.gz $hiseq_trim/R2_paired.fastq.gz $hiseq_trim/R2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:33
# description: 
# PE = paired end input 
# -t <number> = threads 
# -phred33 = base quality encoding
# -trimlog <logfile>
# <input1> <input2> <output_paired1> <output_unpaired1> <output_paired2> <output_unpaired2> 
# ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold> t adapter and other illumina-specific sequences from the read
# LEADING:<quality> Remove low quality bases from the beginning of the read
# TRAILING:<quality> Remove low quality bases from the beginning of the read
# SLIDINGWINDOW:<windowSize>:<requiredQuality> Performs a sliding window trimming approach. It starts scanning at the 5‟ end and clips the read once the average quality within the window falls below a threshold
# MINLEN:<length> Removes reads that fall below the specified minimal length

# trimmed reads for processing
zcat $hiseq_trim/R1_paired.fastq.gz | wc -l #1007718592 reads
zcat $hiseq_trim/R2_paired.fastq.gz | wc -l #1007718592 reads

# quality control on trimmed
# load fastqc
module add UHTS/Quality_control/fastqc/0.11.2
# run qc
fastqc $hiseq_trim/R1_paired.gz -o $hiseq_trim/
fastqc $hiseq_trim/R2_paired.gz -o $hiseq_trim/

# 2) get genome/annotation 

# reference folder
dmel=/scratch/beegfs/monthly/aechchik/MSc/dmel_files

# get chromosomes 
mkdir -p $dmel/chromosomes/
cd $dmel/chromosomes/ && { curl -O ftp://ftp.ensemblgenomes.org/pub/metazoa/release-31/fasta/drosophila_melanogaster/dna/*.dna.chromosome*; gzip -d *.gz; cd -; }

# get full genome
mkdir $dmel/genome/
cd $dmel/genome/ && { curl -O ftp://ftp.ensemblgenomes.org/pub/metazoa/release-31/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.31.dna.genome.fa.gz; gzip -d *.gz; cd -; }

# get gtf annotation
mkdir $dmel/gtf/
cd $dmel/gtf/ && { curl -O ftp://ftp.ensembl.org/pub/release-84/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.84.gtf.gz; gzip -d *.gz; cd -}

# extract transcripts from full genome
mkdir $dmel/transcripts/
# load cufflinks
module add UHTS/Assembler/cufflinks/2.2.1
# extract annotated regions
gffread -g $dmel/genome/Drosophila_melanogaster.BDGP6.31.dna.genome.fa -x $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa $dmel/gtf/Drosophila_melanogaster.BDGP6.84.gtf


## 3) alignment

# prepare wdir 
hiseq_aln=/scratch/beegfs/monthly/aechchik/MSc/illumina/alignment
mkdir -p $hiseq_aln

# hisat2

# load hisat2
module add UHTS/Aligner/hisat/2.0.2
# index genome
hisat2-build $dmel/genome/Drosophila_melanogaster.BDGP6.31.dna.genome.fa $dmel/genome/hisat
# build a list of known splice sites
hisat2_extract_splice_sites.py $dmel/gtf/Drosophila_melanogaster.BDGP6.84.gtf > $dmel/gtf/gtf_splicesites.txt
# prepare outdir
mkdir -p $hiseq_aln/hisat/
# alignment 1-pass 
# note: tested with -M 20971520 on LSF
hisat2 -q --phred33 --known-splicesite-infile $dmel/gtf/gtf_splicesites.txt --novel-splicesite-outfile $dmel/gtf/hisat1_new_splicesites.txt --rna-strandness RF --dta-cufflinks -x $dmel/genome/hisat -1 $hiseq_trim/R1_paired.gz -2 $hiseq_trim/R2_paired.gz -S $hiseq_aln/hisat/hisat_1pass.sam
# alignment 2-pass
# note: tested with -M 20971520 on LSF
hisat2 -q --phred33 --known-splicesite-infile $dmel/gtf/gtf_splicesites.txt --novel-splicesite-infile $dmel/gtf/hisat1_new_splicesites.txt --rna-strandness RF --dta-cufflinks -x $dmel/genome/hisat -1 $hiseq_trim/R1_paired.gz -2 $hiseq_trim/R2_paired.gz -S $hiseq_aln/hisat/hisat_2pass.sam


# star 

# load star
module add UHTS/Aligner/STAR/2.5.0b
# index genome 
STAR --runMode genomeGenerate --genomeDir $dmel/genome/star/ --genomeFastaFiles $dmel/genome/Drosophila_melanogaster.BDGP6.31.dna.genome.fa --sjdbGTFfile $dmel/gtf/Drosophila_melanogaster.BDGP6.84.gtf --sjdbOverhang 100
# prepare outdir
mkdir -p $hiseq_aln/star/
# alignment 1-pass
# note: tested with -M 20971520 on LSF
STAR --genomeDir $dmel/genome/star/ --readFilesIn $hiseq_trim/R1_paired.gz $hiseq_trim/R2_paired.gz --outFilterMismatchNoverLmax 0.04 --alignSJDBoverhangMin 1 --outFileNamePrefix $hiseq_aln/star/star
# prepare alignment 2-pass
# note: need to re-index genome with new splice junctions
# note: tested with -M 8388608 on LSF
STAR --runMode genomeGenerate --genomeDir $dmel/genome/star2/ --genomeFastaFiles $dmel/genome/Drosophila_melanogaster.BDGP6.31.dna.genome.fa --sjdbGTFfile $dmel/gtf/Drosophila_melanogaster.BDGP6.84.gtf --sjdbFileChrStartEnd $hiseq_aln/star/star/STARSJ.out.tab --sjdbOverhang 100
# alignment 2-pass
# note: tested with -M 20971520 on LSF
STAR --genomeDir $dmel/genome/star2/ --readFilesIn $hiseq_trim/R1_paired.gz $hiseq_trim/R2_paired.gz --outFilterMismatchNoverLmax 0.04 --alignSJDBoverhangMin 1 --outFileNamePrefix $hiseq_aln/star/star2 --quantMode TranscriptomeSAM


# mapsplice

# load mapsplice 
module add UHTS/Analysis/MapSplice/2.1.5
# load bowtie
# note: use bowtie1: bowtie2 idx is not supported
module add UHTS/Aligner/bowtie/0.12.9
# index genome in chromosome files with bowtie
bowtie-build $dmel/chromosomes/*.dna.chromosome* $dmel/chromosomes/dmel_chr_bowtie
# get rid of spaces in ref IDs
mkdir $dmel/genome/bowtie/
cp $dmel/genome/Drosophila_melanogaster.BDGP6.31.dna.genome.fa $dmel/genome/bowtie/dmel_ref.fa
sed -i 's/ /_/g' $dmel/genome/bowtie/dmel_ref.fa
# prepare outdir
mkdir -p $hiseq_aln/mapsplice/
# alignment 
mapsplice.py -p 8 -c $dmel/genome/bowtie/ -x $dmel/chromosomes/dmel_chr_bowtie -1  $hiseq_trim/R1_paired.gz -2 $hiseq_trim/R2_paired.gz --gene-gtf $dmel/gtf/Drosophila_melanogaster.BDGP6.84.gtf -o $hiseq_aln/mapsplice/


# bbmap

# load bbmap
module add UHTS/Aligner/BBMap/32.15
# build idx
cd $dmel/genome/ && { bbmap.sh ref=Drosophila_melanogaster.BDGP6.31.dna.genome.fa; cd - }
# prepare outdir
mkdir -p $hiseq_aln/bbmap/ 
# alignment
cd $hiseq_aln/bbmap/ && { bbmap.sh in=$hiseq_trim/R1_paired.gz in2=$hiseq_trim/R2_paired.gz out=bbmap.sam; cd - }


# 3) assembly
