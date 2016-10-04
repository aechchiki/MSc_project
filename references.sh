## Get genome, annotation, extract transcriptome

# reference folder
dmel=/scratch/beegfs/monthly/aechchik/MSc/dmel_files

# get chromosomes 
mkdir -p $dmel/chromosomes/
cd $dmel/chromosomes/ && { curl -O ftp://ftp.ensemblgenomes.org/pub/metazoa/release-31/fasta/drosophila_melanogaster/dna/*.dna.chromosome*; gzip -d *.gz; cd -; }

# get full genome
mkdir -p $dmel/genome/
cd $dmel/genome/ && { curl -O ftp://ftp.ensemblgenomes.org/pub/metazoa/release-31/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.31.dna.genome.fa.gz; gzip -d *.gz; cd -; }

# get gtf annotation
mkdir -p $dmel/gtf/
cd $dmel/gtf/ && { curl -O ftp://ftp.ensembl.org/pub/release-84/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.84.gtf.gz; gzip -d *.gz; cd -}

# extract transcripts from full genome
mkdir -p $dmel/transcripts/
# load cufflinks
module add UHTS/Assembler/cufflinks/2.2.1
# extract annotated regions
gffread -g $dmel/genome/Drosophila_melanogaster.BDGP6.31.dna.genome.fa -x $dmel/transcripts/Drosophila_melanogaster.BDGP6.31.dna.transcripts.fa $dmel/gtf/Drosophila_melanogaster.BDGP6.84.gtf

