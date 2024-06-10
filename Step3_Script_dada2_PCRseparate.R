#set OS type for paths
if(.Platform$OS.type == "unix"){home<-"/home/"} else{
  home<-"//192.168.236.131/"}

#!/usr/bin/env Rscript
library(dada2)
library(vegan)
library(seqRFLP)
library(dplyr)
library(tibble)
library(stringr)

# report
writeLines("\nRunning dada2 denoising ...\n")

############## SET PARAMS ##############
############## SET PARAMS ##############

# load libs and funs
source("/home/genomics/icornelis/03_RawScripts/funs-libs.R")

# get args
option_list <- list( 
    make_option(c("-p","--primer"), type="character"),
    make_option(c("-l","--lib"), type="character")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))

# make paths
proj.path <- here("/home/genomics/icornelis/02_ZEROimpact/01_12S/VJ2022",paste0("MiFish_UE-S_concatenated"))
dir.sense <- paste0(proj.path,"/sense")
dir.antisense <- paste0(proj.path,"/antisense")

mkdir <- paste0(proj.path,"/results")

# confirm proj path 
writeLines(paste0("\n...\nOutput directory set to:\n",proj.path,"\n"))

# trucLens
# for Miya MiFish - truncLen 105 gives at least 29 bp overlap for the longest amplicons (e.g. Raja clavata @ 181 bp), and 40 bp for the regular 170 bp
trucVal <- c(105,105)


############## QUALITY TRIM TRUNCATE ##############
############## QUALITY TRIM TRUNCATE ##############

# report
writeLines("\n...\nQuality trimming and truncating\n")
Sys.sleep(3)

# quality trim Ns and truncate
#out.sense <- filterAndTrim(fwd=cpath("sense","trimmed","R1"), filt=cpath("sense","filtered","R1"), rev=cpath("sense","trimmed", "R2"), filt.rev=cpath("sense","filtered", "R2"), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=trucVal, multithread=TRUE, verbose=TRUE, matchIDs=TRUE)
#out.antisense <- filterAndTrim(fwd=cpath("antisense","trimmed","R1"), filt=cpath("antisense","filtered","R1"), rev=cpath("antisense","trimmed","R2"), filt.rev=cpath("antisense","filtered","R2"), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=trucVal, multithread=TRUE, verbose=TRUE, matchIDs=TRUE)
out.sense <- filterAndTrim(fwd=cpath("sense","trimmed","R1"), filt=cpath("sense","filtered","R1"), rev=cpath("sense","trimmed", "R2"), filt.rev=cpath("sense","filtered", "R2"), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, multithread=TRUE, verbose=TRUE, matchIDs=TRUE)
out.antisense <- filterAndTrim(fwd=cpath("antisense","trimmed","R1"), filt=cpath("antisense","filtered","R1"), rev=cpath("antisense","trimmed","R2"), filt.rev=cpath("antisense","filtered","R2"), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, multithread=TRUE, verbose=TRUE, matchIDs=TRUE)


############## LEARN ERRORS ##############
############## LEARN ERRORS ##############

# report
writeLines("\n...\nLearning errors\n")
Sys.sleep(3)

# learn errors
set.seed(42)
sense.filt.R1.errs <- learnErrors(cpath("sense","filtered","R1"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)
set.seed(42)
sense.filt.R2.errs <- learnErrors(cpath("sense","filtered","R2"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)
set.seed(42)
antisense.filt.R1.errs <- learnErrors(cpath("antisense","filtered","R1"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)
set.seed(42)
antisense.filt.R2.errs <- learnErrors(cpath("antisense","filtered","R2"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)


############## DENOISE ##############
############## DENOISE ##############

# report
writeLines("\n...\ndada2 denoising\n")
Sys.sleep(3)

# run dada denoising - takes time with pool=TRUE
sense.filt.R1.dada <- dada(cpath("sense","filtered","R1"), err=sense.filt.R1.errs, multithread=TRUE, pool=TRUE)
sense.filt.R2.dada <- dada(cpath("sense","filtered","R2"), err=sense.filt.R2.errs, multithread=TRUE, pool=TRUE)
antisense.filt.R1.dada <- dada(cpath("antisense","filtered","R1"), err=antisense.filt.R1.errs, multithread=TRUE, pool=TRUE)
antisense.filt.R2.dada <- dada(cpath("antisense","filtered","R2"), err=antisense.filt.R2.errs, multithread=TRUE, pool=TRUE)


############## DEREPLICATE ##############
############## DEREPLICATE ##############

# report
writeLines("\n...\nDereplication\n")
Sys.sleep(3)

# derep
sense.filt.R1.derep <- derepFastq(cpath("sense","filtered","R1"))
sense.filt.R2.derep <- derepFastq(cpath("sense","filtered","R2"))
antisense.filt.R1.derep <- derepFastq(cpath("antisense","filtered","R1"))
antisense.filt.R2.derep <- derepFastq(cpath("antisense","filtered","R2"))


############## MERGE ##############
############## MERGE ##############

# report
writeLines("\n...\nRead merging\n")
Sys.sleep(3)

# merge the R1 and R2
sense.merged <- mergePairs(dadaF=sense.filt.R1.dada, derepF=sense.filt.R1.derep, dadaR=sense.filt.R2.dada, derepR=sense.filt.R2.derep, verbose=TRUE, maxMismatch=0)
antisense.merged <- mergePairs(dadaF=antisense.filt.R1.dada, derepF=antisense.filt.R1.derep, dadaR=antisense.filt.R2.dada,  derepR=antisense.filt.R2.derep, verbose=TRUE, maxMismatch=0)

# make an OTU table
sense.seqtab <- makeSequenceTable(sense.merged)
antisense.seqtab <- makeSequenceTable(antisense.merged)

# reverse comp the antisense
colnames(antisense.seqtab) <- dada2::rc(colnames(antisense.seqtab))

# fix the names before merging
rownames(sense.seqtab) <- str_split_fixed(rownames(sense.seqtab),"\\.",4)[,1]
rownames(antisense.seqtab) <- str_split_fixed(rownames(antisense.seqtab),"\\.",4)[,1]

# merge the tables
merged.seqtab <- mergeSequenceTables(table1=sense.seqtab, table2=antisense.seqtab, repeats="sum")


############## REMOVE CHIMAERAS ##############
############## REMOVE CHIMAERAS ##############

# report
writeLines("\n...\nDetecting chimaeras\n")
Sys.sleep(3)

# remove chimaeras
merged.seqtab.nochim <- removeBimeraDenovo(merged.seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim <- as.data.frame(merged.seqtab.nochim)
seqtab.nochim.t <- as.data.frame(t(seqtab.nochim))

############## SAVE FILES ##############
############## SAVE FILES ##############

# report
writeLines("\n...\nSaving raw ASV files\n")
Sys.sleep(3)

# make df and fasta for IDs
otus.df <- tibble(names=paste0("MiFish_UE-","asv",str_pad(seq_along(colnames(merged.seqtab.nochim)),width=4,side="left",pad="0")), dnas=colnames(merged.seqtab.nochim)) %>% mutate(len=str_length(dnas))

# write out
write.FASTA(tab2fas(df=otus.df, seqcol="dnas", namecol="names"), file=paste0(proj.path,"/results/asvs.fa"))
write_tsv(as_tibble(otus.df, rownames="asv"), file=paste0(proj.path,"/results/asvs.tsv"))

# save the OTU table with sample names as df
colnames(seqtab.nochim) <- paste0("MiFish_UE-","asv",str_pad(seq_along(colnames(seqtab.nochim)),width=4,side="left",pad="0"))
write_tsv(as_tibble(t(seqtab.nochim), rownames="asv"), file=paste0(proj.path,"/results/asv-table.tsv"))

#count number of sequences per length category
library(Biostrings)
s <- readDNAStringSet(paste0(proj.path,"/results/asvs.fa"))
sum(width(s))/length(s)
sd(width(s))
sum(width(s)>186)
sum(width(s)<185)
sum(width(s)==185)
sum(width(s)<163)
min(width(s))
sort(width(s))

# save runinfo as df
sum(merged.seqtab.nochim)/sum(merged.seqtab)
getN <- function(x) sum(getUniques(x))
trackdenoise.sense <- cbind(sapply(sense.filt.R1.dada, getN), sapply(sense.filt.R2.dada, getN))
trackdenoise.antisense <- cbind(sapply(antisense.filt.R1.dada, getN),sapply(antisense.filt.R2.dada, getN))
trackmerge.sense <- cbind(sapply(sense.merged, getN))
trackmerge.antisense <- cbind(sapply(antisense.merged, getN))
trackmerge <- cbind(rowSums(merged.seqtab), rowSums(merged.seqtab.nochim))

dada2_track <- merge(out.sense, trackdenoise.sense, by=0, all.x=TRUE)
colnames(dada2_track) <- c("samples","input sense", "filter sense", "denoised F", "denoised R")
dada2_track <- merge(dada2_track, trackmerge.sense, by.x=1, by.y=0, all.x=TRUE)
colnames(dada2_track) <- c("samples","input sense", "filter sense","denoised F", "denoised R",
                           "merge F+R")
dada2_track <- merge(dada2_track, out.antisense, by.x=1, by=0, all.x=TRUE)
colnames(dada2_track) <- c("samples","input sense", "filter sense", "denoised F", "denoised R",
                           "merge F+R",
                           "input antisense","filter antisense")
dada2_track <- merge(dada2_track, trackdenoise.antisense, by.x=1, by.y=0, all.x=TRUE)
colnames(dada2_track) <- c("samples","input sense", "filter sense", "denoised F", "denoised R",
                           "merge F+R",
                           "input antisense","filter antisense","denoised F AS", "denoised R AS")
dada2_track <- merge(dada2_track, trackmerge.antisense, by.x=1, by.y=0, all.x=TRUE)
colnames(dada2_track) <- c("samples","input sense", "filter sense","denoised F", "denoised R",
                           "merge F+R",
                           "input antisense","filter antisense", "denoised F AS", "denoised R AS",
                           "merge F+R AS")
dada2_track$samples <- gsub(".R1.fastq.gz", "", dada2_track$samples)
dada2_track <- merge(dada2_track, trackmerge, by.x=1, by.y=0, all.x=TRUE)
colnames(dada2_track) <- c("samples","input sense", "filter sense", "denoised F", "denoised R",
                           "merge F+R",
                           "input antisense","filter antisense", "denoised F AS", "denoised R AS",
                           "merge F+R AS", 
                           "remove chimeras", "#ASV")
dada2_track[is.na(dada2_track)] <- 0

write_tsv(as_tibble(dada2_track, rownames="names"), file=paste0(proj.path,"/results/dada2_track_v2.tsv"))

# report
writeLines(paste0("\n...\nASVs written to:\n",paste0(proj.path,"/results/asvs.fna"),"\n"))
writeLines(paste0("\n...\nASV table written to:\n",paste0(proj.path,"/results/asv-table.tsv"),"\n"))

#save environment after merging in DADA2
save.image(paste0(proj.path,"/results/R_Environment/REnvironment_DADA2_merging.RData"))

#remove data that will not be used for further analysis
rm(antisense.filt.R1.dada, antisense.filt.R1.derep, antisense.filt.R1.errs, 
   antisense.filt.R2.dada, antisense.filt.R2.derep, antisense.filt.R2.errs, antisense.merged,
   sense.filt.R1.dada, sense.filt.R1.derep, sense.filt.R1.errs,
   sense.filt.R2.dada, sense.filt.R2.derep, sense.filt.R2.errs,
   out.antisense, out.sense, antisense.seqtab, sense.seqtab, 
   trackdenoise.sense, trackdenoise.antisense, trackmerge.sense, trackmerge.antisense, trackmerge)

#count number of sequences per sample an sort from lowest to highest
NumberOfSequences <-as.data.frame(sort.default(rowSums(seqtab.nochim[1:nrow(seqtab.nochim),]))) #the total number of reads per sample
NumberOfRawReads <- as.data.frame(colSums(seqtab.nochim)) # the total number of reads per ASV
NumberOfASVs <- as.data.frame(rowSums(seqtab.nochim !=0)) # the total number of ASVs per sample
write_tsv(as_tibble(NumberOfASVs, rownames="asv"), file=paste0(proj.path,"/results/NumberOfASVs.tsv"))

#Create rarecurve before concatenation of PCR replicates
rarecurve(seqtab.nochim, ylab = "ASVs", main = "Rarecurve of unconcatenated samples before taxonomic assignment", label = FALSE, step =1000) 
abline(10000)

# add taxa
set.seed(100)
taxa <- assignTaxonomy(merged.seqtab.nochim, "/home/genomics/icornelis/03_RawScripts/Step4_12S_ReferenceDB_TaxonomicAssignment/12S_references.fas", minBoot = 80, multithread = TRUE)
write_tsv(as_tibble(taxa, rownames="ASV"), file=paste0(proj.path,"/results/taxonomic-assignment_unrarefied.tsv"))
taxa <- data.frame(taxa)
#tax_rare <- assignTaxonomy(seqtab.concatenated.rarefied, "/home/genomics/icornelis/03_RawScripts/Step4_12S_ReferenceDB_TaxonomicAssignment/12S_references_20220906.fas", minBoot = 80, multithread = TRUE)
#tax_rare2 <- addSpecies(tax_rare, "/home/genomics/MiFish-UE_run2/Step3_DADA2/icornelis/12S references_species-assignment.fa")
#write_tsv(as_tibble(tax_rare, rownames="ASV"), file=paste0(proj.path,"/results/taxonomic-assignment_rarefied.tsv"))
#write_tsv(as_tibble(tax_rare2, rownames="ASV"), file=paste0(proj.path,"/results/taxonomic-assignment_rarefied_2.tsv"))
#taxa_rarefied <- data.frame(tax_rare)

#for unrarefied data
table_unrarefied<-as.data.frame(t(merged.seqtab.nochim))

#remove rownames
rownames(table_unrarefied)<-NULL

#add extra columns to table to add taxonomy, for this first transform the tax object to data frame
table_unrarefied$Kingdom<-taxa$Kingdom
table_unrarefied$Phylum<-taxa$Phylum
table_unrarefied$Class<-taxa$Class
table_unrarefied$Order<-taxa$Order
table_unrarefied$Family<-taxa$Family
table_unrarefied$Genus<-taxa$Genus
table_unrarefied$Species<-taxa$Species
table_unrarefied$ASV<-otus.df$names
write_tsv(as_tibble(table_unrarefied), file=paste0(proj.path,"/results/table_unrarefied.tsv"))
write.xlsx(as.data.frame(table_unrarefied),
           paste0(proj.path,"/results/table_unrarefied_raw.xlsx"), sheetName = "table_unrarefied_raw", colNames = TRUE, rowNames = FALSE, append = FALSE)
