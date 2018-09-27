# Nanopore Pipeline for Automated Analysis of cDNA Sequencing Data


This script is under construction, we are very grateful for suggestive contributions and ideas.

The pipeline was designed to integrate following tasks:

Combine fastq files from nanopore outputs, trim sequences and seperate barcoded samples, alignment and analysis



## How to install the package from GitHub

Download Script from Github, required packages are:
```
library(DESeq2)
library(rapport)
library(RColorBrewer)
library(AutoPipe)
```
Also required is installed Porechop, Minimap2, Sam tools. Be sure that all tools are avaiable within your system.

The Pipeline is working steps:

```
# Easy Workflow for Nanopore Alignment

Samples_discription=read.csv("samp_desc.csv", row.names=1)
#Example of the sample description is given


#Set up a object from "Poreseq" class
Set=Poreseq(Samples_discription)

# path of the folder where the "fastq" folder is located
Set@path=c("/media/daka/Data/RNA_Seq/Slices")

# Combine fastq files
GET_FASQ(Set)

# Multiplex Samples
TRIM_Barcodes(Set)

# Alignment der Sequences
Set=NANOPORE_Aligner(Set)

#Analyze Data
Set=Analyzer(Set,Write_exp=T,filter_genes=50,MA_PLOT=T)


```

Good luck...





### Authors

D. H. Heiland, Translational Research Group, Medcal-Center Freiburg, University of Freiburg

