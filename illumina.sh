# new 

# rename files: chekc old script

# data in /scratch/beegfs/monthly/aechchik/MSc/illumina/reads

# merge all reads in one R1 and one R2 
bsub -J R1all "cat *_R1.fastq.gz > R1_all.fastq.gz"; bsub -J R2all "cat *_R2.fastq.gz > R2_all.fastq.gz"

# fastqc on raw
bsub -J R1fastqc "module add UHTS/Quality_control/fastqc/0.11.2; fastqc R1_all.fastq.gz -o ."; bsub -J R2fastqc " module add UHTS/Quality_control/fastqc/0.11.2; fastqc R2_all.fastq.gz -o ."

# trimming 
bsub -J R1trimmomatic -q dee-hugemem  -M 20971520 "module add UHTS/Analysis/trimmomatic/0.33; trimmomatic PE -phred33 -trimlog TrimLog R1_all.fastq.gz R2_all.fastq.gz R1_paired.fastq.gz R1_unpaired.fastq.gz R2_paired.fastq.gz R2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:33"

# fastqc on trimmed 
bsub -J R1fastqc "module add UHTS/Quality_control/fastqc/0.11.2; fastqc R1_paired.fastq.gz -o ."; bsub -J R2fastqc " module add UHTS/Quality_control/fastqc/0.11.2; fastqc R2_paired.fastq.gz -o ."


