#!/usr/bin/env Rscript
library(tidyverse)
library(spgs)

remotes::install_github("cran/seqRFLP")
remotes::install_github("benjjneb/dada2")
install.packages("Rcpp")
library(Rcpp)
library(seqRFLP)
library(optparse)
library(dada2)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

install.packages("Biostrings")
library(Biostrings)
# load libs funs
#keep this folder always like this
source(here::here("/home/genomics/icornelis/03_RawScripts/funs-libs.R"))
#source(here::here("/Volumes/genomics/icornelis/03_RawScripts/funs-libs.R"))

# make paths
proj.path <- here("/home/genomics/ndukan/NJ_2022_output/")


# upload barcodes and their reverse complement sequences
plates2 <- readxl::read_excel(paste0(proj.path,"20221122_Overview_DNAextractions_NJ2022_eDNA_OWF.xlsx"),sheet = "Data Demultiplexing")

#F_primer_tag <- DNAString(x="plates2$`F-primer tag`)
#create reverse complement of the sequences
plates2$FwdRevComp <- reverseComplement(plates2$`F-primer tag`, content = "dna", case = "upper")
plates2$RevRevComp <- reverseComplement(plates2$`R-primer tag`, content = "dna", case = "upper")

#save the files in fasta format into your folder
bcF <- dataframe2fas(data.frame(plates2$`Demultiplex name`, plates2$`F-primer tag`), (file=paste0(proj.path,"ZEROimpact_NJ2022_barcodes_F-tag.fas")))
bcR <- dataframe2fas(data.frame(plates2$`Demultiplex name`, plates2$`R-primer tag`), (file=paste0(proj.path,"ZEROimpact_NJ2022_barcodes_R-tag.fas")))

#for the reverse complements
bcRF <- dataframe2fas(data.frame(plates2$`Demultiplex name`, plates2$FwdRevComp), (file=paste0(proj.path,"ZEROimpact_NJ2022_barcodes_rcF-tag.fas")))
bcRR <- dataframe2fas(data.frame(plates2$`Demultiplex name`, plates2$RevRevComp), (file=paste0(proj.path,"ZEROimpact_NJ2022_Rbarcodes_rcR-tag.fas")))

