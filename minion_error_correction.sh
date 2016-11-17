# tried to use nanocorr, nanopolish, racon (depends on miniasm), but not successful. 
# report Canu step -correction

# r7 2D only

/home/aechchik/wgets/canu/Linux-amd64/bin/canu -correct -p r7 -d . genomeSize=61000000 -nanopore-raw r7_2d.fasta useGrid=false 

# r9 2D only

/home/aechchik/wgets/canu/Linux-amd64/bin/canu -correct -p r9 -d . genomeSize=61000000 -nanopore-raw r9_2d.fasta useGrid=false stopOnReadQuality=false
# canu was installed on mh home. required the latest version
# parameters: 
# -correct: to compute only read correction, no trimming or assembly 
# -p: prefix of canu output
# -d: directory to canu output
# genomeSize: estimated assembly size, in my case transcriptome size
# -nanopore-raw: error profile tuning for minion raw reads
# useGrid=False: disable default LSF
# stopOnReadQuality=false: not sure of what this means, but had to put it as parameter to avoid program crash (suggested in log, not necessary for r7 data)

