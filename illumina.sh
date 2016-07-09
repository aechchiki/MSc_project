# new 

# rename files: chekc old script

# data in /scratch/beegfs/monthly/aechchik/MSc/illumina/reads

# merge all reads in one R1 and one R2 
bsub -J R1all "cat *_R1.fastq.gz > R1_all.fastq.gz"; bsub -J R2all "cat *_R2.fastq.gz > R2_all.fastq.gz"

# fastqc on raw
bsub -J R1fastqc "module add UHTS/Quality_control/fastqc/0.11.2; fastqc R1_all.fastq.gz -o ."; bsub -J R2fastqc " module add UHTS/Quality_control/fastqc/0.11.2; fastqc R2_all.fastq.gz -o ."

# get trimmomatic TruSeq3-PE
curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
unzip Trimmomatic-0.33.zip
cd Trimmomatic-0.33/adapters


# trimming 
bsub -J R1trimmomatic -q dee-hugemem  -M 20971520 "module add UHTS/Analysis/trimmomatic/0.33; trimmomatic PE -phred33 -trimlog TrimLog R1_all.fastq.gz R2_all.fastq.gz R1_paired.fastq.gz R1_unpaired.fastq.gz R2_paired.fastq.gz R2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:33"
# man: http://www.usadellab.org/cms/?page=trimmomatic

# note
wc -l R1_paired.fastq.gz ; wc -l R2_paired.fastq.gz 
54561422 R1_paired.fastq.gz
54915439 R2_paired.fastq.gz

wc -l R1_paired.fastq.gz ; wc -l R2_paired.fastq.gz 
54614615 R1_paired.fastq.gz
55108556 R2_paired.fastq.gz



# print IDs list
bsub  "zcat R1_paired.fastq.gz | grep ^@ > R1_paired_IDs"; bsub "zcat R2_paired.fastq.gz | grep ^@ > R2_paired_IDs"
# print universal IDs 
bsub "sed 's/ .*//g' R1_paired_IDs > R1_IDs" ; bsub "sed 's/ .*//g' R2_paired_IDs > R2_IDs"
# print shared IDs
bsub "comm -12 R1_IDs_sorted R2_IDs_sorted > comm_IDs"

# extract files for procesing 
bsub "gunzip R1_paired.fastq.gz"; bsub "gunzip R2_paired.fastq.gz"

# get reads with shared IDs from R1
cat comm_IDs | while read line; do awk '/$line/ {for(i=1; i<=3; i++) {getline; print}}' R1_paired.fastq > R1_paired_kept.fastq ; done
# get reads with shared IDs from R2
cat comm_IDs | while read line; do awk '/$line/ {for(i=1; i<=3; i++) {getline; print}}' R2_paired.fastq > R2_paired_kept.fastq ; done


# fastqc on trimmed 
bsub -J R1fastqc "module add UHTS/Quality_control/fastqc/0.11.2; fastqc R1_paired.fastq.gz -o ."; bsub -J R2fastqc " module add UHTS/Quality_control/fastqc/0.11.2; fastqc R2_paired.fastq.gz -o ."


