# how to generate the FinalTables for mapping

# 1) follow instructions on pacbio_isoseq_modifs.sh to generate a suitable sam file to be parsed by aechchiki/misc_genomics/MapTable.sh

# 2) run  aechchiki/misc_genomics/MapTable.sh

# 3) generate the table:

# define alignment directory
pacbio_aln=/scratch/beegfs/monthly/aechchik/MSc/pacbio/alignment/isoseq

# concatenate all "FinalTable" files (one per alignment stats launched)
find $pacbio_aln -name "FinalTable" -type f -exec cat {} + > $pacbio_aln/CatFinalTable

# generate header
cat $pacbio_aln/CatFinalTable | head -1 > $pacbio_aln/Head

# one row per mapping stats
cat $pacbio_aln/CatFinalTable | grep -v ^Mapp > $pacbio_aln/Body

# cat header and body
cat $pacbio_aln/Head $pacbio_aln/Body > $pacbio_aln/FinalTable.txt
