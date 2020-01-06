Code to simulate replication timing in the mouse or human genome. 


## Simulating DNA replication in mammalian genomes (RT-Sim)

### Sample execution:

Rscript run_allCSmodel.R -i output.stats -t 1000 -c 3 -n 0.15 -b experimentalRT.bedgraph -o origins.bedgraph -g hg38 


### Bash environment variables: 
GENOMES : Path to a genomes folder. Must contain a sub-folder for the genome passed via -g to run_allCSmodel.R. Subfolder should contain two files 1) genome.fa - a FASTA file of the genome; 2) genome.fa.fai - a FASTA index file (use bedtools faidx)


### Requirements: 

R version 3.5.2


### R packages: 

animation

data.table

dplyr

extrafont

factoextra

gganimate

ggplot2

ggpmisc

ggpubr

grid

gridExtra

lsr

Metrics

numform 

optparse

plyr

preprocessCore

reshape2

tictoc

zoo




