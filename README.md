# MicrobDetect
MicrobDetect is a pipeline to identify and to profile the microbiome content from the whole genome shotgun sequencing data of the clinical samples. The clinical samples, e.g. the biopsy samples from the human stomach mucosal, might contain mostly human DNA and very low abundance of microbial DNA. The generic tools developed for the metagenomic whole shortgun sequencing data profiling tools may not be suitable for those clinical microbial analysis. We therefore developed "MicroDetect" pipeline.

# Applications
MicroDetect is developed for the analysis of:
  - detect possible existed pathongens in the clinical samples.
  - generate microbial profile from the whole genome sequencing data.
  - estimate the abudance of microbiome at different taxonomy levels, from strains to phylum.
  - multiple testing on the paired samples combining metadata revealing the potential outcome related microbes. 

# How to run 

a) Download and install MicrobDetect:
   
   git clone https://github.com/gk-zhang/MicrobDetect.git
   
b) Download clinical whole genome sequencing samples:

c) filter out the host reads

d) prepare the microbial/pathogen reference sequence from Patric database

e) map samples against reference database

d) evaluate and estimate the microbial abundance

f) profiling the microbial abundance for all the samples

g) downstream statistical analysis combining clinical data

