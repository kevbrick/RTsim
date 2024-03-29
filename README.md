## Simulating DNA replication in mammalian genomes (RT-Sim)

### Sample execution (to test script):
Rscript run_allCSmodel.R -i output.stats -t 200 -c 3 -n 0.15 -b mm10_ESC.sampleRT.bedgraph -o mm10.sample_origins.bedgraph

##### NOTE
Sample RT (mm10_ESC.sampleRT.bedgraph) and origins files (mm10.sample_origins.bedgraph) are included in the repo.

### Command line args:
arg | required | detail
--- | --- | --------------------------------
-b  |  *  | Replication timing BEDGRAPH file (for exact format, see below)
-c  |     | Total number of haploid genomes to simulate (default = 500)
-d  |     | DSB hotspots BEDGRAPH (NOT USED CURRENTLY)
-e  |     | Modelling run time (default: not used)
-f  |     | Percent replicating cells (default: not used)
-g  |     | Genome (default = mm10)
-i  |  *  | Output ID 
-m  |     | Save model (default = FALSE)
-n  |  *  | Origins per Mb (default = 0.2)
-o  |  *  | Origins BEDGRAPH file (column 4 in BEDGRAPH should be origin efficiency)
-r  |     | "Recycle" replication machinery (default = TRUE)
-s  |     | Use origin efficiency to weight random selection (default = TRUE)
-t  |  *  | Modelling cycles (default = 8000)
-u  |     | Output prefix
-v  |     | Exclude chromosome (NOT TESTED)
-w  |     | Random seed (all sims will give the same result)
-x  |     | Single chromosome (NOT TESTED)
-y  |     | Image output folder
-z  |     | Output RDS file (default = FALSE)

### Bash environment variables: 
$GENOMES  \
Path to a genomes folder. Must contain a sub-folder for the genome passed via -g to run_allCSmodel.R. \
Subfolder should contain two files : 
1) genome.fa - a FASTA file of the genome 
2) genome.fa.fai - a FASTA index file (use bedtools faidx) 

### Requirements: 
R version 3.5.2 

### R packages: 
animation \
data.table \
dplyr \
extrafont \
factoextra \
gganimate \
ggplot2 \
ggpmisc \
ggpubr \
grid \
gridExtra \
lsr \
Metrics \
numform  \
optparse \
plyr \
preprocessCore \
reshape2 \
tictoc \
zoo 



