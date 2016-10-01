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

# 3) assembly
