### MinION

## new 
# data in: /home/jroux/archive/MinION/run_MinION_2015_09_29.tar 
#97G archive!

# extract fastq from fast5
# module add UHTS/Analysis/poretools/0.5.1

cd /scratch/beegfs/monthly/aechchik/MSc/minion/data
bsub -J unzip_tar "tar -xvf /home/jroux/archive/MinION/run_MinION_2015_09_29.tar -C ."

# how may files are there? 
cd ./DownloadsFromMetrichor
# sample name: HiSeq_cDNAJRtrial2_4103_1_ch99_file436_strand.fast5
bsub -J fast5list "module add UHTS/Analysis/poretools/0.5.1; ls | grep fast5 > filelist"


bsub -J  "module add UHTS/Analysis/poretools/0.5.1; poretools fastq *.fast5 > fast5.fastq"

poretools fastq filelist > fast5.fastq


cat filelist | xargs cat >> fast5.fast5



echo !(*.fast5) | xargs ls > fast5.fastq

# extract fastq from fast5

bsub -J  "module add UHTS/Analysis/poretools/0.5.1; poretools fastq *.fast5 > fast5.fastq"



bsub -J poretoo "module add UHTS/Analysis/poretools/0.5.1; poretools fastq *.fast5 > fast5.fastq"



# read extraction: fast5 to fastq 

# check quality
 bsub -J porestat "module add UHTS/Analysis/poretools/0.5.1; poretools qualdist *.fast5 > quality.txt"

# build stats on reads
bsub -J porestat "module add UHTS/Analysis/poretools/0.5.1; poretools stats *.fast5 > stats.txt"

# convert fast5 to fasta 
bsub -J poretoo "module add UHTS/Analysis/poretools/0.5.1; poretools fastq *.fast5 > fast5.fastq"


# error correction
# from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4662598

# hybrid method: correct with previous data generated on hiseq (paired, 2x 100bp)
module add UHTS/Analysis/proovread/2.12.13;
# use proovread-flex, which is optimised for the uneven sequencing coverage seen in metagenomes and transcriptomes. 
# input: long reads (fastq), short reads (fastq)

https://github.com/BioInf-Wuerzburg/proovread/blob/master/README.pdf?raw=true
# not sure yet of how it works


# For de novo error correction: nanocorrect NO: deprecateed
# better use nanopolish To complete our de novo assembly paper we wrote a second, much more powerful, software package called nanopolish that uses the nanopore signal data to improve the accuracy of the assembly. 

# The reads that are input into the HMM must be output as a .fa file by poretools.
cd /scratch/beegfs/monthly/aechchik/MSc/minion/fast5_bc/fast5
module add UHTS/Analysis/poretools/0.5.1; poretools fasta *.fast5 > fast5_merged.fasta
# uoutput in /scratch/beegfs/monthly/aechchik/MSc/minion/fast5_bc/fasta
#??? how does it work 
# tried:
make -f /home/aechchik/bin/nanopolish/scripts/consensus.make READS=fast5_merged.fasta ASSEMBLY=fast5_merged_polished.fasta
# error:
make: *** No rule to make target `fast5_merged_polished.fasta', needed by `fast5_merged_polished.fasta.bwt'.  Stop.`


# We can either use BWA or LAST. In our experience BWA tends to be fast but less sensitive, and LAST tends to be slow but more sensitive.

# bwa index genome
module add UHTS/Aligner/bwa/0.7.13; 
bwa index /scratch/beegfs/monthly/aechchik/MSc/dmel_files/dmel_bdgp6_fa/Drosophila_melanogaster.BDGP6.31.dna.genome.fa
# create fasta index for samtoolsa
module add UHTS/Analysis/samtools/1.3;
samtools faidx /scratch/beegfs/monthly/aechchik/MSc/dmel_files/dmel_bdgp6_fa/Drosophila_melanogaster.BDGP6.31.dna.genome.fa
# run bwa and pipe straight to samtools to create BAM
bwa mem -x ont2d /scratch/beegfs/monthly/aechchik/MSc/dmel_files/dmel_bdgp6_fa/Drosophila_melanogaster.BDGP6.31.dna.genome.fa fast5_merged.fastq | samtools view -T /scratch/beegfs/monthly/aechchik/MSc/dmel_files/dmel_bdgp6_fa/Drosophila_melanogaster.BDGP6.31.dna.genome.fa -bS - | samtools sort -T 2D_vs_NC_Drosophila_melanogaster.BDGP6.31.bwa -o 2D_vs_NC_Drosophila_melanogaster.BDGP6.31.bam -


LAST straight to sorted BAM:

# create a LAST db of the reference genome
module add SequenceAnalysis/SequenceAlignment/last/531;
lastdb BDGP631 /scratch/beegfs/monthly/aechchik/MSc/dmel_files/dmel_bdgp6_fa/Drosophila_melanogaster.BDGP6.31.dna.genome.fa

# align high quaklity reads to reference genome with LAST
lastal -q 1 -a 1 -b 1 BDGP631 /scratch/beegfs/monthly/aechchik/MSc/minion/fast5_bc/fasta/fast5_merged.fasta > /scratch/beegfs/monthly/aechchik/MSc/minion/fast5_bc/fasta/fast5_merged.fasta.maf

# convert the MAF to BAM with complete CIGAR (matches and mismatches)

# get nanopore scripts 
git clone https://github.com/arq5x/nanopore-scripts.git
module add UHTS/Analysis/samtools/1.3;
python /home/aechchik/bin/nanopore-scripts/maf-convert.py sam /scratch/beegfs/monthly/aechchik/MSc/minion/fast5_bc/fasta/fast5_merged.fasta.maf | samtools view -T /scratch/beegfs/monthly/aechchik/MSc/dmel_files/dmel_bdgp6_fa/Drosophila_melanogaster.BDGP6.31.dna.genome.fa -bS - | samtools sort -T 2D_vs_NC_Drosophila_melanogaster.BDGP6.31.last -o 2D_vs_NC_Drosophila_melanogaster.BDGP6.31.last.bam -


