# FOCS
Codebase for FOCS, a new method for inferring enhancer-promoter (E-P) links and for integrative analysis of E-P links

Additional information about running FOCS can be found <a href="http://acgt.cs.tau.ac.il/focs/tutorial.html">here</a>

# Getting FOCS database
You can download the database <a href="http://acgt.cs.tau.ac.il/focs/download.html" target="_blank">here</a>.
The database includes preprocessed RData files of ENCODE, Roadmap Epigenomics, FANTOM5, and GRO-seq:
1) Predicted E-P links
2) Genomic positions of enhancers and promoters
3) Enhancer/promoter count/RPKM profiles
4) Sample annotations

In addition the database includes external data for integrative analysis of E-P links:
1) ChIA-PET interactions mediated by POL2 of 4 cell types (MCF7, HCT-116, K562 and Hela-S3) taken from <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39495" target="_blank">GSE39495</a>
2) Genomic positions of GTEx eQTL data (V6 and V6p builds) taken from <a href="https://www.gtexportal.org/home/" target="_blank">GTEx Portal</a>

# Dependencies
Please make sure you have the following pre-installed R packages:
* pscl
* MASS
* parallel
* AUC
* glmnet
* Matching
* RColorBrewer
* GenomicRanges
* ggplot2

# Running FOCS to infer E-P links
Please note that you can run FOCS only under linux platform
First create directories:
1) mkdir focs
2) cd focs
3) mkdir data
4) mkdir scripts
5) mkdir tmp

Copy the scripts in the R folder to ../focs/scripts/ folder
Additional instructions are given within tutorial.R script file

