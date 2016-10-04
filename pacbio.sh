# 0) original data
# from portal: http://smrt-lgtf.vital-it.ch/user/listanalysis

# note: tmp done with subreads (fastq/a from bax.h5)
# note2: where did I get ROI from? not on smrt.vital-it.ch (?)

# original directory
pacbio_roi=/scratch/beegfs/monthly/aechchik/PacBioData/ROI
# note: 3 fractions! 
# reference
dmel=/scratch/beegfs/monthly/aechchik/MSc/dmel_files

# convert to fasta
cat $pacbio_roi/ROI_1-2k.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $pacbio_roi/ROI_1-2k.fasta
cat $pacbio_roi/ROI_2-3k.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $pacbio_roi/ROI_2-3k.fasta
cat $pacbio_roi/ROI_3-7k.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $pacbio_roi/ROI_3-7k.fasta

# 1) alignment

## roi alignment

# note: for splice unaware aligners, map to transcriptome; for splice-aware aligners, map to genome

# prepare directory
pacbio_aln=/scratch/beegfs/monthly/aechchik/MSc/pacbio/alignment/roi
mkdir -p $pacbio_aln

# bbmap 
# note: on genome because splice-aware
# prepare directory
mkdir -p $pacbio_aln/bbmap
# load bbmap
module add UHTS/Aligner/BBMap/32.15
# build idx
cd $dmel/genome/ && { bbmap.sh ref=Drosophila_melanogaster.BDGP6.31.dna.genome.fa; cd - ; }
# alignment
cd $dmel/genome/ && { bbmap.sh in=$pacbio_roi/ROI_1-2k.fasta out=$pacbio_aln/bbmap/bbmap_ROI_1-2k.sam; cd - ; }
cd $dmel/genome/ && { bbmap.sh in=$pacbio_roi/ROI_2-3k.fasta out=$pacbio_aln/bbmap/bbmap_ROI_2-3k.sam; cd - ; }
cd $dmel/genome/ && { bbmap.sh in=$pacbio_roi/ROI_3-7k.fasta out=$pacbio_aln/bbmap/bbmap_ROI_3-7k.sam; cd - ; }

# blasr
# note; on transcriptome because splice-unaware
# prepare directory
mkdir -p $pacbio_aln/blasr
# load blasr
# note: pacbio specific software
module add UHTS/PacBio/blasr/20140829
# alignment
blasr $pacbio_roi/ROI_1-2k.fasta $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -nproc 24 -sam -out $pacbio_aln/blasr/blasr_ROI_1-2k.sam
blasr $pacbio_roi/ROI_2-3k.fasta $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -nproc 24 -sam -out $pacbio_aln/blasr/blasr_ROI_2-3k.sam
blasr $pacbio_roi/ROI_3-7k.fasta $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -nproc 24 -sam -out $pacbio_aln/blasr/blasr_ROI_3-7k.sam

# bwamem
# note: on genome because splice-aware
# prepare directory
mkdir -p $pacbio_aln/bwamem
# load bwa
module add UHTS/Aligner/bwa/0.7.13
# alignment
bwa mem -x pacbio $dmel/genome/Drosophila_melanogaster.BDGP6.31.dna.genome.fa $pacbio_roi/ROI_1-2k.fasta > $pacbio_aln/bwamem/bwamem_ROI_1-2k.sam
bwa mem -x pacbio $dmel/genome/Drosophila_melanogaster.BDGP6.31.dna.genome.fa $pacbio_roi/ROI_2-3k.fasta > $pacbio_aln/bwamem/bwamem_ROI_2-3k.sam
bwa mem -x pacbio $dmel/genome/Drosophila_melanogaster.BDGP6.31.dna.genome.fa $pacbio_roi/ROI_3-7k.fasta > $pacbio_aln/bwamem/bwamem_ROI_3-7k.sam

# gmap
# note: on genome because splice-aware
# note: not installed on cluster. ref to my home: need installation guidelines? on draft
# prepare directory
mkdir -p $pacbio_aln/gmap
# load module 
module add UHTS/Analysis/gmap/2016.08.10
# create database
gmap_build -d GDB -D $dmel/chromosomes/ $dmel/chromosomes/*.fa 
# alignment 
gmap -d GDB -D $dmel/chromosomes/ $pacbio_roi/ROI_1-2k.fasta -f samse -n 0 -t 16 > $pacbio_aln/gmap/gmap_ROI_1-2k.sam
gmap -d GDB -D $dmel/chromosomes/ $pacbio_roi/ROI_2-3k.fasta -f samse -n 0 -t 16 > $pacbio_aln/gmap/gmap_ROI_2-3k.sam
gmap -d GDB -D $dmel/chromosomes/ $pacbio_roi/ROI_3-7k.fasta -f samse -n 0 -t 16 > $pacbio_aln/gmap/gmap_ROI_3-7k.sam

# graphmap
# note: on trasncriptome because non-splice aware
# note: not installed on cluster. ref to my /home: need installation guidelines?on draft
# prepare directory
mkdir -p $pacbio_aln/graphmap
# alignment
/home/aechchik/bin/graphmap align -r $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -d $pacbio_roi/ROI_1-2k.fasta -o $pacbio_aln/graphmap/graphmap_ROI_1-2k.sam
/home/aechchik/bin/graphmap align -r $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -d $pacbio_roi/ROI_2-3k.fasta -o $pacbio_aln/graphmap/graphmap_ROI_2-3k.sam
/home/aechchik/bin/graphmap align -r $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -d $pacbio_roi/ROI_3-7k.fasta -o $pacbio_aln/graphmap/graphmap_ROI_3-7k.sam
# note: readme http://research-pub.gene.com/gmap/src/README


### isoseq alignment

# note: for splice unaware aligners, map to transcriptome; for splice-aware aligners, map to genome

# prepare directory
pacbio_aln=/scratch/beegfs/monthly/aechchik/MSc/pacbio/alignment/isoseq
mkdir -p $pacbio_aln

# bbmap 
# note: on genome because splice-aware
# prepare directory
mkdir -p $pacbio_aln/bbmap
# load bbmap
module add UHTS/Aligner/BBMap/32.15
# build idx
cd $dmel/genome/ && { bbmap.sh ref=Drosophila_melanogaster.BDGP6.31.dna.genome.fa; cd - ; }
# alignment
cd $dmel/genome/ && { bbmap.sh in=$pacbio_isoseq/isoseq_1-2k.fasta out=$pacbio_aln/bbmap/bbmap_isoseq_1-2k.sam; cd - ; }
cd $dmel/genome/ && { bbmap.sh in=$pacbio_isoseq/isoseq_2-3k.fasta out=$pacbio_aln/bbmap/bbmap_isoseq_2-3k.sam; cd - ; }
cd $dmel/genome/ && { bbmap.sh in=$pacbio_isoseq/isoseq_3-7k.fasta out=$pacbio_aln/bbmap/bbmap_isoseq_3-7k.sam; cd - ; }

# blasr
# note; on transcriptome because splice-unaware
# prepare directory
mkdir -p $pacbio_aln/blasr
# load blasr
# note: pacbio specific software
module add UHTS/PacBio/blasr/20140829
# alignment
blasr $pacbio_isoseq/isoseq_1-2k.fasta $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -nproc 24 -sam -out $pacbio_aln/blasr/blasr_isoseq_1-2k.sam
blasr $pacbio_isoseq/isoseq_2-3k.fasta $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -nproc 24 -sam -out $pacbio_aln/blasr/blasr_isoseq_2-3k.sam
blasr $pacbio_isoseq/isoseq_3-7k.fasta $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -nproc 24 -sam -out $pacbio_aln/blasr/blasr_isoseq_3-7k.sam

# bwamem
# note: on genome because splice-aware
# prepare directory
mkdir -p $pacbio_aln/bwamem
# load bwa
module add UHTS/Aligner/bwa/0.7.13
# alignment
bwa mem -x pacbio $dmel/genome/Drosophila_melanogaster.BDGP6.31.dna.genome.fa $pacbio_isoseq/isoseq_1-2k.fasta > $pacbio_aln/bwamem/bwamem_isoseq_1-2k.sam
bwa mem -x pacbio $dmel/genome/Drosophila_melanogaster.BDGP6.31.dna.genome.fa $pacbio_isoseq/isoseq_2-3k.fasta > $pacbio_aln/bwamem/bwamem_isoseq_2-3k.sam
bwa mem -x pacbio $dmel/genome/Drosophila_melanogaster.BDGP6.31.dna.genome.fa $pacbio_isoseq/isoseq_3-7k.fasta > $pacbio_aln/bwamem/bwamem_isoseq_3-7k.sam

# gmap
# note: on genome because splice-aware
# note: not installed on cluster. ref to my home: need installation guidelines? on draft
# prepare directory
mkdir -p $pacbio_aln/gmap
# load module 
module add UHTS/Analysis/gmap/2016.08.10
# create database
gmap_build -d GDB -D $dmel/chromosomes/ $dmel/chromosomes/*.fa 
# alignment 
gmap -d GDB -D $dmel/chromosomes/ $pacbio_isoseq/isoseq_1-2k.fasta -f samse -n 0 -t 16 > $pacbio_aln/gmap/gmap_isoseq_1-2k.sam
gmap -d GDB -D $dmel/chromosomes/ $pacbio_isoseq/isoseq_2-3k.fasta -f samse -n 0 -t 16 > $pacbio_aln/gmap/gmap_isoseq_2-3k.sam
gmap -d GDB -D $dmel/chromosomes/ $pacbio_isoseq/isoseq_3-7k.fasta -f samse -n 0 -t 16 > $pacbio_aln/gmap/gmap_isoseq_3-7k.sam

# graphmap
# note: on trasncriptome because non-splice aware
# note: not installed on cluster. ref to my /home: need installation guidelines?on draft
# prepare directory
mkdir -p $pacbio_aln/graphmap
# alignment
/home/aechchik/bin/graphmap -r $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -d $pacbio_isoseq/isoseq_1-2k.fasta -o $pacbio_aln/graphmap/graphmap_isoseq_1-2k.sam
/home/aechchik/bin/graphmap align -r $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -d $pacbio_isoseq/isoseq_2-3k.fasta -o $pacbio_aln/graphmap/graphmap_isoseq_2-3k.sam
/home/aechchik/bin/graphmap align -r $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa -d $pacbio_isoseq/isoseq_3-7k.fasta -o $pacbio_aln/graphmap/graphmap_isoseq_3-7k.sam
# note: readme http://research-pub.gene.com/gmap/src/README


# 2) assembly?
