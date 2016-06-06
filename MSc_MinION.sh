### MinION

# check quality
 bsub -J porestat "module add UHTS/Analysis/poretools/0.5.1; poretools qualdist *.fast5 > quality.txt"

# build stats on reads
bsub -J porestat "module add UHTS/Analysis/poretools/0.5.1; poretools stats *.fast5 > stats.txt"

# convert fast5 to fasta 
bsub -J poretoo "module add UHTS/Analysis/poretools/0.5.1; poretools fasta *.fast5 > fast5.fasta"



 
 
