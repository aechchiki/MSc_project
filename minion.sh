## 0) original data

# original directory
minion_od=/home/jroux/archive/MinION/

# create wdir
minion_wd=/scratch/beegfs/monthly/aechchik/MSc/minion/
mkdir -p $minion_wd
# unzip files to wdir 
mkdir -p $minion_wd/r7
tar -xvf $minion_od/run_MinION_2015_09_29.tar -C $minion_wd/r7
mkdir -p $minion_wd/r9
tar -xvzf $minion_od/run_MinION_2016_06_23.tar.gz -C $minion_wd/r9

## 1) data extraction: fast5 to fasta/q
minion_raw_r7=$minion_wd/r7/DownloadsFromMetrichor/
minion_raw_r9=$minion_wd/r9/DownloadsFromMetrichor/
# load poretools
module add UHTS/Analysis/poretools/0.5.1
# prepare extraction directory
minion_reads=$minion_wd/reads
mkdir -p minion_reads
# extract
poretools fastq $minion_raw_r7/*.fast5 > $minion_reads/r7_reads.fastq
poretools fastq $minion_raw_r9/*.fast5 > $minion_reads/r9_reads.fastq

# select only 2d reads
# convert to fasta 
cat $minion_reads/r7_reads.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $minion_reads/r7_reads.fasta
cat $minion_reads/r9_reads.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $minion_reads/r9_reads.fasta
# print header and sequence for 2d reads
cat $minion_reads/r7_reads.fasta | awk '/_Basecall_2D_000_2d/{nr[NR]; nr[NR+1]}; NR in nr' > $minion_reads/r7_2d.fasta
cat $minion_reads/r9_reads.fasta | awk '/_Basecall_2D_000_2d/{nr[NR]; nr[NR+1]}; NR in nr' > $minion_reads/r9_2d.fasta


## 2) alignment

# note: for splice unaware aligners, map to transcriptome; for splice-aware aligners, map to genome
# note: alignment should be done using 2d reads only

# prepare directory
minion_aln=/scratch/beegfs/monthly/aechchik/MSc/minion/alignment
mkdir -p $minion_aln

# bbmap 
# note: on genome because splice-aware
# prepare directory
mkdir -p $minion_aln/bbmap
# load bbmap
module add UHTS/Aligner/BBMap/32.15
# build idx
cd $dmel/genome/ && { bbmap.sh ref=Drosophila_melanogaster.BDGP6.31.dna.genome.fa; cd - }
# alignment
cd $dmel/genome/ && { bbmap.sh in=$minion_reads/r7_2d.fasta out=$minion_aln/bbmap/bbmap_r7.sam; cd - }
cd $dmel/genome/ && { bbmap.sh in=$minion_reads/r9_2d.fasta out=$minion_aln/bbmap/bbmap_r9.sam; cd - }


# blasr
# note; on transcriptome because splice-unaware
# prepare directory
mkdir -p $minion_aln/blasr
# load blasr
# note: pacbio specific software
module add UHTS/PacBio/blasr/20140829
# alignment
blasr $minion_reads/r7_2d.fasta $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -nproc 24 -sam -out $minion_aln/blasr/blasr_r7.sam
blasr $minion_reads/r9_2d.fasta $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -nproc 24 -sam -out $minion_aln/blasr/blasr_r9.sam

# bwamem
# note: on genome because splice-aware
# prepare directory
mkdir -p $minion_aln/bwamem
# load bwa
module add UHTS/Aligner/bwa/0.7.13
# alignment
bwa mem -x ont2d $dmel/genome/Drosophila_melanogaster.BDGP6.31.dna.genome.fa $minion_reads/r7_2d.fasta > $minion_aln/bwamem/bwamem_r7.sam
bwa mem -x ont2d $dmel/genome/Drosophila_melanogaster.BDGP6.31.dna.genome.fa $minion_reads/r9_2d.fasta > $minion_aln/bwamem/bwamem_r9.sam

# gmap
# note: on genome because splice-aware
# note: not installed on cluster. ref to my home: need installation guidelines? on draft
# prepare directory
mkdir -p $minion_aln/gmap
# create database
/home/aechchik/bin/gmap-2016-06-30/bin/gmap_build -d GDB -D $dmel/chromosomes/ $dmel/chromosomes/*.fa 
# alignment 
/home/aechchik/bin/gmap-2016-06-30/bin/gmap -d GDB -D $dmel/chromosomes/ $minion_reads/r7_2d.fasta -f samse -n 0 -t 16 > $minion_aln/gmap/gmap_r7.sam
/home/aechchik/bin/gmap-2016-06-30/bin/gmap -d GDB -D $dmel/chromosomes/ $minion_reads/r9_2d.fasta -f samse -n 0 -t 16 > $minion_aln/gmap/gmap_r9.sam
# note: readme http://research-pub.gene.com/gmap/src/README

# graphmap
# note: on trasncriptome because non-splice aware
# note: not installed on cluster. ref to my /home: need installation guidelines?on draft
# prepare directory
mkdir -p $minion_aln/graphmap
# alignment
/home/aechchik/bin/graphmap/bin/Linux-x64/graphmap align -r $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -d $minion_reads/r7_2d.fasta -o $minion_aln/graphmap/graphmap_r7.sam
/home/aechchik/bin/graphmap/bin/Linux-x64/graphmap align -r $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -d $minion_reads/r9_2d.fasta -o $minion_aln/graphmap/graphmap_r9.sam

# error correction and remap? not sure. check venome snake publication!

# *) assembly??
