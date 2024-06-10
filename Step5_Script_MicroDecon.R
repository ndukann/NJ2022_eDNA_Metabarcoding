if(.Platform$OS.type == "unix"){home<-"/home/"} else{
  home<-"//192.168.236.131/"}

library(ggplot2); packageVersion("ggplot2")
library(dplyr)
library(vegan)
library(openxlsx)
library(tibble)
library(tidyverse)
library(ade4)
library(here)
library(devtools)
devtools::install_github("donaldtmcknight/microDecon")
library(microDecon)

proj.path <- here("/home/genomics/ndukan/01_Zero-Impact/Temporal_Analysis")

table_raw <- readxl::read_excel(paste0(proj.path, "/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/table_unrarefied_raw_Full_TaxAss_WoRMS.xlsx"))

#see what species are on the list
table_raw_noneg <- table_raw[,!grepl("Neg", colnames(table_raw))]
table_raw_fish <- table_raw_noneg[table_raw_noneg$Class %in% fish_classes,]
table_raw_fish <- table_raw_fish[!rowSums(table_raw_fish[,1:(ncol(table_raw_fish)-11)])==0,]
table_raw_fish <- table_raw_fish[!(table_raw_fish$Species %in% freshwater_species), ]
table_raw_fish <- table_raw_fish[!table_raw_fish$Species == "NA", ]

a <- unique(table_raw_fish$Species) #71 taxa -4 higher level -7 freshwater = 60 species level, 64 taxa
b <- rownames(merged_table_unrarefied_fish)
setdiff(a,b)
  
#arrange the microDecon table. First column ASVs, starting from second column control samples, then biological samples and finally taxa (Species column)
table_raw_noT <- as.data.frame(table_raw[,1:(ncol(table_raw)-11)])

table_ASV <- as.data.frame(table_raw$ASV)
names(table_ASV) <- "OTU ID"
table_neg <- as.data.frame(table_raw[,grepl("Neg", colnames(table_raw))])
# table_neg_noField <- as.data.frame(table_neg[, !grepl("Neg_Field", colnames(table_neg))])
table_sample <- as.data.frame(table_raw_noT[,!grepl("Neg", colnames(table_raw_noT))])
table_taxa <- as.data.frame(table_raw$Species)
names(table_taxa) <- "Taxa"

table_decon <- cbind(table_ASV, table_neg,table_sample, table_taxa)
table_decon1 <- table_decon[, 2:(ncol(table_decon)-1)]
table_decon1 <- table_decon1[, !colSums(table_decon1)==0]
table_decon1 <- cbind(table_ASV, table_decon1, table_taxa)
table_decon2 <- table_decon1[,!grepl("Field", colnames(table_decon1))]
colnames(table_decon1)[1] <- "ASV"

# summ <- as.data.frame(colSums(table_decon1))
# summrow <- as.data.frame(rowSums(table_decon1))
# setdiff(colnames(table_decon),colnames(table_decon1))

# env <- readRDS(file = here("/home/genomics/ndukan/01_Zero-Impact/Temporal_Analysis/02_NJ2022_Miseq_5119441/R_Environment/env_unrarefied_ordered_BRZR.rds"))
samples_decon <- colnames(table_decon2[,2:(ncol(table_decon2)-1)])
samples_decon <- samples_decon[!grepl("Neg", samples_decon)]
Locations <- c(unique(gsub("\\_.*", "", samples_decon)))
neg_numbers <- length(c(colnames(table_decon2[, grepl("Neg", colnames(table_decon2))])))

sample_numbers <- NULL

for (i in 1:length(Locations)){
  sample_numbers[i] <- sum(str_count(gsub("\\_.*", "", samples_decon), paste0("\\b", Locations[i], "\\b"))) 
}
sample_numbers
# 30 30 12 15 15 30 27 25 30 15 27 30 30 15 30 30 27


table_raw_decontaminated <- decon(data = table_decon1, numb.blanks=neg_numbers, numb.ind=sample_numbers, taxa=T,
                             runs=2, thresh = 0.7, prop.thresh = 0.00005, regression=0, low.threshold=40, up.threshold=400)

table_decontaminated <- table_raw_decontaminated$decon.table #887 ASVs but 391 ASVs totally removed from the dataset
reads_removed <- table_raw_decontaminated$reads.removed 
difference_sum <- table_raw_decontaminated$sum.per.group
difference_mean <- table_raw_decontaminated$mean.per.group
ASVs_removed <- table_raw_decontaminated$OTUs.removed

#find which ASVs are totally removed and how many of them are fish ASVs
ASV_all <- c(table_raw$ASV)
ASV_decontaminated <- c(table_decontaminated$ASV)
ASV_names_removed <- setdiff(ASV_all, ASV_decontaminated)

table_removed <- table_raw[table_raw$ASV %in% ASV_names_removed,]
table_removed_fish <- table_removed[table_removed$Class %in% fish_classes,]
unique(table_removed_fish$Species)

#find which ASvs are identified as contaminants and how many of them are fish ASVs
ASV_contaminant <- c(ASVs_removed$ASV)
table_contaminant <- table_raw[table_raw$ASV %in% ASV_contaminant,]
table_contaminant_fish <- table_contaminant[table_contaminant$Class %in% fish_classes,]
unique(table_contaminant_fish$Species)

#check Zeus Faber
zeus <- table_raw[table_raw$Species=="Zeus faber",]
zeus <- zeus[,1:(ncol(zeus)-11)]
zeus <- zeus[, !colSums(zeus)==0]

#check Microstomus kitt
kitt <- table_raw[table_raw$Species=="Microstomus kitt",]
kitt <- kitt[,1:(ncol(kitt)-11)]
kitt <- kitt[, !colSums(kitt)==0]

spondy <- table_raw[table_raw$Species== "Spondyliosoma cantharus",]
spondy <- spondy[,1:(ncol(spondy)-11)]
spondy <- spondy[, !colSums(spondy)==0]

#without the field negatives
table_raw_decontaminated_noField <- decon(data = table_decon2, numb.blanks=neg_numbers, numb.ind=sample_numbers, taxa=T,
                                  runs=2, thresh = 0.7, prop.thresh = 0.00005, regression=0, low.threshold=40, up.threshold=400)

table_decontaminated_noField <- table_raw_decontaminated_noField$decon.table #247 ASVs contaminant but 121 ASVs are totally removed.
ASVs_removed_noField <- table_raw_decontaminated_noField$OTUs.removed


list_of_datasets <- list("Contaminant_All_Negs" = ASVs_removed, "Contaminant_noField_Negs" = ASVs_removed_noField)
write.xlsx(list_of_datasets, paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Contaminant_Tables_All.xlsx"), colNames = T)


#raw decontaminated table with field negatives
table_clean <- table_decontaminated[, !grepl("Taxa|blank", colnames(table_decontaminated))]
table_raw_clean <- merge(table_clean, table_raw[(ncol(table_raw)-10):ncol(table_raw)],
                         by.x = 1, by.y = "ASV", all.x = T)
colnames(table_raw_clean)[1] <- "ASV"
table_raw_clean <- table_raw_clean %>% relocate("ASV", .before = 'DADA2')

#raw decontaminated table without field negatives
table_clean_noField <- table_decontaminated_noField[, !grepl("Taxa|blank", colnames(table_decontaminated_noField))]
table_raw_clean_noField <- merge(table_clean_noField, table_raw[(ncol(table_raw)-10):ncol(table_raw)],
                         by.x = 1, by.y = "ASV", all.x = T)
colnames(table_raw_clean_noField)[1] <- "ASV"
table_raw_clean_noField <- table_raw_clean_noField %>% relocate("ASV", .before = 'DADA2')

table_raw<-as.data.frame(table_raw)
rownames(table_raw) <- table_raw$ASV

table_raw_clean <- as.data.frame(table_raw_clean)
rownames(table_raw_clean) <- table_raw_clean$ASV

totally_removed <- as.data.frame(setdiff(rownames(table_raw), rownames(table_raw_clean) ))
names(totally_removed) <-"Totally removed ASVs"

write.xlsx(totally_removed, 
           paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Totally_Removed_ASVs.xlsx"), 
           sheetName = "Sheet1", colNames = TRUE, rowNames = FALSE, append = FALSE)

#save the raw tables#save the raw tablestable_raw_clean
write.xlsx(table_raw_clean, 
           paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/table_unrarefied_raw_Full_TaxAss_WoRMS_MicroDecon.xlsx"), 
           sheetName = "Sheet1", colNames = TRUE, rowNames = FALSE, append = FALSE)

write.xlsx(table_raw_clean_noField, 
           paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/table_unrarefied_raw_Full_TaxAss_WoRMS_MicroDecon_noField.xlsx"), 
           sheetName = "Sheet1", colNames = TRUE, rowNames = FALSE, append = FALSE)

#concatenate the samples
seqtab_raw <- as.data.frame(t(table_raw_clean[,1:(ncol(table_raw_clean)-11)]))
colnames(seqtab_raw) <- table_raw_clean$ASV
seqtab_raw$names <- str_sub(rownames(seqtab_raw), end=-4) # removes PCR number (_S1, _S2 or _S3) from samplenames
seqtab_concatenated <- aggregate(seqtab_raw[,1:ncol(seqtab_raw)-1], by= list(seqtab_raw$names),FUN=sum)
rownames(seqtab_concatenated) <- seqtab_concatenated$Group.1
seqtab_concatenated <- seqtab_concatenated[,2:ncol(seqtab_raw)]  
NumberOfSequencesConcat <- as.data.frame(sort.default(rowSums(seqtab_concatenated[1:nrow(seqtab_concatenated),]))) #number of sequences per sample after concatenation, 8 samples <300 reads
NumberOfASVsConcat <- as.data.frame(sort.default(rowSums(seqtab_concatenated !=0))) #number of ASVs per sample after concatenation
write_tsv(as_tibble(NumberOfASVsConcat, rownames="asv"), file=paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/NumberOfASVsPerSample_Concat_unrarefied.tsv"))
write_tsv(as_tibble(NumberOfSequencesConcat, rownames="samples"), file=paste0(proj.path, "/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/NumberOfSequencesPerSamples_Concat_unrarefied.tsv"))

table_clean_concatenated_noT <- as.data.frame(t(seqtab_concatenated))
table_clean_concatenated_noT <- table_clean_concatenated_noT[, c(env_order$Niskin.sample[which(!env_order$Zone == "neg_control")])]
table_unrarefied_concatenated <- merge(table_clean_concatenated_noT, table_raw[(ncol(table_raw)-10):ncol(table_raw)],
                                       by.x = 0, by.y = "ASV", all.x = T)
colnames(table_unrarefied_concatenated)[1] <- "ASV"
table_unrarefied_concatenated <- table_unrarefied_concatenated %>% relocate("ASV", .before = 'DADA2')

write.xlsx(table_unrarefied_concatenated, 
           paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean_MicroDecon.xlsx"), 
           sheetName = "Sheet1", colNames = TRUE, rowNames = FALSE, append = FALSE)

#save the ASV contaminant table
write.xlsx(difference_sum, 
           paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Contaminant_reads_removed.xlsx"), 
           sheetName = "Sheet1", colNames = TRUE, rowNames = FALSE, append = FALSE)

#Create fish ASV table 
fish_classes <- readRDS(file = paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/Fish_classes.rds"))
table_unrarefied_fish <- as.data.frame(table_unrarefied_concatenated[table_unrarefied_concatenated$Class %in% fish_classes,])

#remove the freshwater species
freshwater_species <- c("Alburnus alburnus", "Barbatula barbatula", "Gobio gobio", "Rutilus rutilus", "Squalius cephalus", "Oreochromis niloticus", "Cyprinus carpio haematopterus x Megalobrama amblycephala")
table_unrarefied_fish <- table_unrarefied_fish[!(table_unrarefied_fish$Species %in% freshwater_species), ]


taxo <- "Species"
merged_table_unrarefied_fish <- aggregate(table_unrarefied_fish[,1:(ncol(table_unrarefied_fish)-11)], by= list(as.factor(table_unrarefied_fish[,taxo])),FUN=sum)
rownames(merged_table_unrarefied_fish) <- merged_table_unrarefied_fish$Group.1
merged_table_unrarefied_fish$Group.1 <- NULL

merged_table_unrarefied_fish <- merged_table_unrarefied_fish[!(rownames(merged_table_unrarefied_fish) == "NA"), ]
merged_table_unrarefied_fish <- merged_table_unrarefied_fish[, !colSums(merged_table_unrarefied_fish)==0]
merged_table_unrarefied_fish <- merged_table_unrarefied_fish[!rowSums(merged_table_unrarefied_fish)==0,]

#compare the fish species with the manuscript
merged_table_unrarefied_fish_old <- readxl::read_excel(paste0(proj.path, "/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Table_S10_Species_TotalReadCount.xlsx"))
merged_table_unrarefied_fish_old <- merged_table_unrarefied_fish_old[-(1:2), ]
merged_table_unrarefied_fish_old <- as.data.frame(merged_table_unrarefied_fish_old)
rownames(merged_table_unrarefied_fish_old) <- merged_table_unrarefied_fish_old[,1]
setdiff(rownames(merged_table_unrarefied_fish_old), rownames(merged_table_unrarefied_fish))

#check the read numbers of removed species and how they distributed across the samples and compare it with trawl data
raja_clavata <- table_raw[table_raw$Species=="Raja clavata",]
microstomus_kitt <- table_raw[table_raw$Species=="Microstomus kitt",]

#check the read numbers of fishes (new data)
Fish_reads <- as.data.frame(rowSums(merged_table_unrarefied_fish))
names(Fish_reads) <- "Reads"

#save the fish table and the read counts
table_unrarefied_fish <- table_unrarefied_fish[,!colSums(table_unrarefied_fish[,1:(ncol(table_unrarefied_fish)-11)])==0]

write.xlsx(table_unrarefied_fish, 
           paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Table_FishASVs_MicroDecon.xlsx"), 
           sheetName = "Sheet1", colNames = TRUE, rowNames = TRUE, append = FALSE)

write.xlsx(merged_table_unrarefied_fish, 
           paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Table_FishSpecies_MicroDecon.xlsx"), 
           sheetName = "Sheet1", colNames = TRUE, rowNames = TRUE, append = FALSE)

write.xlsx(Fish_reads, 
           paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/FishReadCounts_MicroDecon.xlsx"), 
           sheetName = "Sheet1", colNames = TRUE, rowNames = TRUE, append = FALSE)

#create environmental data
# seqtab_concatenated <- as.data.frame(table_unrarefied_concatenated[,1:(ncol(table_unrarefied_concatenated)-11)])
env_unrarefied <- as.data.frame(colnames(merged_table_unrarefied_fish))
colnames(env_unrarefied) <- 'Niskin.sample'
Coast <- readRDS(file = here("/home/genomics/ndukan/01_Zero-Impact/Temporal_Analysis/02_NJ2022_Miseq_5119441/R_Environment/Coast.rds"))
Coast <- Coast[Coast != "ftBRZRbis"]
Transition <- readRDS(file = here("/home/genomics/ndukan/01_Zero-Impact/Temporal_Analysis/02_NJ2022_Miseq_5119441/R_Environment/Transition.rds"))
Transition
Offshore <- readRDS(file = here("/home/genomics/ndukan/01_Zero-Impact/Temporal_Analysis/02_NJ2022_Miseq_5119441/R_Environment/Offshore.rds"))
Offshore

saveRDS(Coast, file = paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/Coast.rds"))
saveRDS(Transition, file = paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/Transition.rds"))
saveRDS(Offshore, file = paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/Offshore.rds"))

env_unrarefied$Location <- gsub("\\_.*", "", env_unrarefied$Niskin.sample)
env_unrarefied$Zone <- ifelse(env_unrarefied$Location %in% Coast, "Coast",
                              ifelse(env_unrarefied$Location %in% Transition, "Transition",
                                     ifelse(env_unrarefied$Location %in% Offshore, "Offshore", "neg_control")))
env_unrarefied$Zone_color <- ifelse(env_unrarefied$Zone=="Coast","seagreen3",
                                           ifelse(env_unrarefied$Zone=="Transition", "steelblue3",
                                                  ifelse(env_unrarefied$Zone=="Offshore","darkorange","red")))

zone_order <- c("Coast", "Transition", "Offshore")
env_unrarefied <- env_unrarefied[order(factor(env_unrarefied$Zone, levels = zone_order)), ]

saveRDS(env_unrarefied, file = paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/env_unrarefied_MicroDecon.rds"))


#species accumulation curves
seqtab_fish <- as.data.frame(table_unrarefied_fish[,1:(ncol(table_unrarefied_fish)-11)])

seqtab_115bis <- seqtab_fish[,grepl("115bis", colnames(seqtab_fish))]
seqtab_115bis <- seqtab_115bis[!rowSums(seqtab_115bis)==0,]
spec_accum_115bis <- iNEXT(seqtab_115bis, q=0, datatype="abundance")
plot(spec_accum_115bis, type=3)

seqtab_7105 <- seqtab_fish[,grepl("7105", colnames(seqtab_fish))]
seqtab_7105 <- seqtab_7105[!rowSums(seqtab_7105)==0,]
spec_accum_7105 <- iNEXT(seqtab_7105, q=0, datatype="abundance")
plot(spec_accum_7105, type=1)


coast_samples <- env$Niskin.sample[env$Zone == "Coast"]
table_coast <- seqtab_fish[, colnames(seqtab_fish) %in% coast_samples]
table_coast <- table_coast[!rowSums(table_coast)==0,]
spec_accum_coast <- iNEXT(table_coast, q=0, datatype="abundance")
plot(spec_accum_coast, type=1)
ggiNEXT(spec_accum_coast, type=1, se=TRUE, grey=FALSE, color.var = "Assemblage")

save.image((paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/R_Environment_MicroDecon/R_Environment_MicroDecon_Concatenation.RData")))
