# 0) original data
# from portal: http://smrt-lgtf.vital-it.ch/user/listanalysis

# note: tmp done with subreads (fastq/a from bax.h5)
# note2: where did I get ROI from? not on smrt.vital-it.ch (?)

# original directory
pacbio_roi=/scratch/beegfs/monthly/aechchik/PacBioData/ROI/
# note: 3 fractions! 

# convert to fasta
cat $pacbio_roi/ROI_1-2k.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $pacbio_roi/ROI_1-2k.fasta
cat $pacbio_roi/ROI_2-3k.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $pacbio_roi/ROI_2-3k.fasta
cat $pacbio_roi/ROI_3-7k.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $pacbio_roi/ROI_3-7k.fasta

# 1) alignment

# 2) assembly?
