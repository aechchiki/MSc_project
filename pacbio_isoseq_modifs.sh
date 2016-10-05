# to run the mapping scripts: need to have the same header. some mappers are wriitng additional mapping parameters to the read name, which makes the reads not recognized as being in the fasta file, thus not mapped. need to fix this

# pacbio isoseq reads are coming from two different pipelines. need to uniform format and cat them before launching the mapping script. 

# full-length reads are in fasta format, but the sequence is 60nt\n format
# consensus polished are in fastq format

# cat the file to be analyzed 
pacbio_reads_isoseq=/scratch/beegfs/monthly/aechchik/MSc/pacbio/reads/isoseq
mkdir -p $pacbio_reads_isoseq

pacbio_lgtf_isoseq=/scratch/beegfs/monthly/aechchik/PacBioData/isoseq
pacbio_lgtf_12=$pacbio_lgtf_isoseq/cDNA_1-2k_Isoseq
pacbio_lgtf_23=$pacbio_lgtf_isoseq/cDNA_2-3k_Isoseq
pacbio_lgtf_37=$pacbio_lgtf_isoseq/cDNA_3-7k_Isoseq

# convert fastq to fasta 

cat $pacbio_lgtf_12/Cluster/polished_high_qv_consensus_isoforms.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $pacbio_lgtf_12/Cluster/polished_high_qv_consensus_isoforms.fasta
cat $pacbio_lgtf_23/Cluster/polished_high_qv_consensus_isoforms.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $pacbio_lgtf_23/Cluster/polished_high_qv_consensus_isoforms.fasta
cat $pacbio_lgtf_37/Cluster/polished_high_qv_consensus_isoforms.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $pacbio_lgtf_37/Cluster/polished_high_qv_consensus_isoforms.fasta

# define files to cat

cluster=/Cluster/polished_high_qv_consensus_isoforms.fasta
classify=/Classify/isoseq_flnc.fasta

# define outdir 

pacbio_isoforms=/scratch/beegfs/monthly/aechchik/MSc/pacbio/reads/isoseq
mkdir -p $pacbio_isoforms

# cat files from Cluster and Classify to isoforms

cat $pacbio_lgtf_12$classify $pacbio_lgtf_12$cluster > $pacbio_isoforms/isoforms_1-2k.fasta
cat $pacbio_lgtf_23$classify $pacbio_lgtf_23$cluster > $pacbio_isoforms/isoforms_2-3k.fasta
cat $pacbio_lgtf_37$classify $pacbio_lgtf_37$cluster > $pacbio_isoforms/isoforms_3-7k.fasta

# conform fasta: remove 60nt\n format

cat $pacbio_isoforms/isoforms_1-2k.fasta | awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' > $pacbio_isoforms/isoforms_clean_1-2k.fasta
cat $pacbio_isoforms/isoforms_2-3k.fasta | awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' > $pacbio_isoforms/isoforms_clean_2-3k.fasta
cat $pacbio_isoforms/isoforms_3-7k.fasta | awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' > $pacbio_isoforms/isoforms_clean_3-7k.fasta

# cut non-meaningful read name in fasta header

cat $pacbio_isoforms/isoforms_clean_1-2k.fasta | sed 's/ isoform=c.*/ /g' | sed 's/CCS.*/CCS/g' > $pacbio_isoforms/isoseq_1-2k.fasta
cat $pacbio_isoforms/isoforms_clean_2-3k.fasta | sed 's/ isoform=c.*/ /g' | sed 's/CCS.*/CCS/g' > $pacbio_isoforms/isoseq_2-3k.fasta
cat $pacbio_isoforms/isoforms_clean_3-7k.fasta | sed 's/ isoform=c.*/ /g' | sed 's/CCS.*/CCS/g' > $pacbio_isoforms/isoseq_3-7k.fasta

# changes in the mapping file to conform to the fasta header of the corresponding fraction 

# define mapping directory

pacbio_aln=/scratch/beegfs/monthly/aechchik/MSc/pacbio/alignment/isoseq

# change header: cut after the first space and the first tab in mapping file

cat $pacbio_aln/bbmap/bbmap_isoseq_1-2k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/bbmap/bbmap_isoseqclean_1-2k.sam
cat $pacbio_aln/bbmap/bbmap_isoseq_2-3k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/bbmap/bbmap_isoseqclean_2-3k.sam
cat $pacbio_aln/bbmap/bbmap_isoseq_3-7k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/bbmap/bbmap_isoseqclean_3-7k.sam

cat $pacbio_aln/blasr/blasr_isoseq_1-2k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/blasr/blasr_isoseqclean_1-2k.sam
cat $pacbio_aln/blasr/blasr_isoseq_2-3k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/blasr/blasr_isoseqclean_2-3k.sam
cat $pacbio_aln/blasr/blasr_isoseq_3-7k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/blasr/blasr_isoseqclean_3-7k.sam

cat $pacbio_aln/bwamem/bwamem_isoseq_1-2k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/bwamem/bwamem_isoseqclean_1-2k.sam
cat $pacbio_aln/bwamem/bwamem_isoseq_2-3k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/bwamem/bwamem_isoseqclean_2-3k.sam
cat $pacbio_aln/bwamem/bwamem_isoseq_3-7k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/bwamem/bwamem_isoseqclean_3-7k.sam

cat $pacbio_aln/gmap/gmap_isoseq_1-2k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/gmap/gmap_isoseqclean_1-2k.sam
cat $pacbio_aln/gmap/gmap_isoseq_2-3k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/gmap/gmap_isoseqclean_2-3k.sam
cat $pacbio_aln/gmap/gmap_isoseq_3-7k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/gmap/gmap_isoseqclean_3-7k.sam

cat $pacbio_aln/graphmap/graphmap_isoseq_1-2k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/graphmap/graphmap_isoseqclean_1-2k.sam
cat $pacbio_aln/graphmap/graphmap_isoseq_2-3k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/graphmap/graphmap_isoseqclean_2-3k.sam
cat $pacbio_aln/graphmap/graphmap_isoseq_3-7k.sam | grep -v ^@ | sed 's/ [^\t]*\t/\t/' > $pacbio_aln/graphmap/graphmap_isoseqclean_3-7k.sam

