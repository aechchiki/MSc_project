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
minion_fq=$minion_wd/reads
mkdir -p minion_fq
# extract
poretools fastq $minion_raw_r7/*.fast5 > $minion_fq/r7_reads.fq
poretools fastq $minion_raw_r9/*.fast5 > $minion_fq/r9_reads.fq


# 1) alignment


# 2) assembly?
