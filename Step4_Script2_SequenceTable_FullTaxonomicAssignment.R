save.image((paste0(proj.path,"/Temporal_Analysis/02_NJ2022_Miseq_5119441/R_Environment/R_Environment_Step5TaxonomicAssignment_BLAST_Worms.RData")))

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
libraries <- c("here", "stringr", "taxize", "curl", "readr", "xlsx", "tibble", "openxlsx", "tidyverse", "data.table")

for (opties in libraries){
  
  if (opties %in% installed.packages()){
    
    library(opties,character.only = TRUE)
    
  } else {install.packages(opties,repos = "http://cran.us.r-project.org")
    
    library(opties,character.only = TRUE)
  }
}
#library(here, dplyr, stringr, taxize, curl, readr, xlsx, tibble, tidyverse, data.table) => didn't work like this. 

# load libs and funs
source("/home/genomics/icornelis/03_RawScripts/funs-libs.R")

proj.path <- here("/home/genomics/ndukan/01_Zero-Impact")
table_unrarefied <- readxl::read_excel(paste0(proj.path,"/NJ2022/table_unrarefied_DADA2assigned.xlsx"),sheet = 
                                            "table_unrarefied_DADA2aasigned")
blastn_Genbank <- readxl::read_excel(paste0(proj.path,"/NJ2022/blastn_GenBank_20230323.xlsx"),
                                     sheet = "blastn_GenBank_20230323")


#Define blast genbank => takes too much time, previous GenBank results will be used.
{
  blastn_gen_com<-paste(
    "/usr/local/bin/blastn" 
    ,"-num_threads 5" 
    ,"-db /home/genomics/bioinf_databases/genbank/nt2/nt" 
    ,"-query ", paste0(proj.path,"/NJ2022/asvs_run2.fa")
    ,"-out "  , paste0(proj.path,"/NJ2022/blastn_GenBank_",format(Sys.time(), "%Y%m%d"),".tsv")
    ,"-max_hsps 1 -max_target_seqs 1 -qcov_hsp_perc 75 -perc_identity 90" 
    ,"-outfmt '6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen evalue bitscore qcovs'"
  )
  # excute blastn
  system(blastn_gen_com)
}

#Define blast reference db
{
  blastn_ref_com<-paste(
    "/usr/local/bin/blastn" 
    ,"-num_threads 5" 
    ,"-db /home/genomics/icornelis/03_RawScripts/Step4_12S_ReferenceDB_TaxonomicAssignment/12S_references.fa" 
    ,"-query ", paste0(proj.path,"/NJ2022/asvs_run2.fa")
    ,"-out "  , paste0(proj.path,"/NJ2022/blastn_own_reference_tophits_",format(Sys.time(), "%Y%m%d"),".tsv")
    ,"-max_hsps 1 -max_target_seqs 1 -qcov_hsp_perc 75 -perc_identity 90" 
    ,"-outfmt '6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen evalue bitscore qcovs'"
  )
  # excute blastn
  system(blastn_ref_com)
}

{
  blastn_ref_com<-paste(
    "/usr/local/bin/blastn" 
    ,"-num_threads 5" 
    ,"-db /home/genomics/icornelis/03_RawScripts/Step4_12S_ReferenceDB_TaxonomicAssignment/12S_references.fa" 
    ,"-query ", paste0(proj.path,"/NJ2022/asvs_run2.fa")
    ,"-out "  , paste0(proj.path,"/NJ2022/blastn_own_reference_top10hits_",format(Sys.time(), "%Y%m%d"),".tsv")
    ,"-max_hsps 1 -max_target_seqs 10 -qcov_hsp_perc 75 -perc_identity 90" 
    ,"-outfmt '6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen evalue bitscore qcovs'"
  )
  # excute blastn
  system(blastn_ref_com)
}

# Read results
blastn_ref <- read_tsv(file.path(proj.path, "NJ2022", "blastn_own_reference_tophits_20230712.tsv"), 
                       col_names = c("ASV", "sseqid", "staxids", "Taxonomy", "pident", "qlen", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs"))

# Separate tax in each level
blastn_ref <- separate(blastn_ref, Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")

# in the blastn data select ASV and species column 
TA_DADA2 <- data.frame(table_unrarefied$ASV, table_unrarefied$Species) 
blastn_Genbank_select <- blastn_Genbank[!blastn_Genbank$pident < 97.000,]
blastn_Genbank <- data.frame(blastn_Genbank_select$ASV, blastn_Genbank_select$Species)
blastn_Genbank_select$sseqidlong <- blastn_Genbank_select$sseqid
blastn_Genbank_select$sseqid <- str_extract(blastn_Genbank_select$sseqidlong, "([A-Z])\\w+\\d")
blastn_Genbank_select$gi_id <- str_extract(blastn_Genbank_select$sseqidlong, "\\d{7,10}")

#blastn_ref <- data.frame(blastn_ref$ASV, blastn_ref$Species)
blastn_ref_select <- blastn_ref[!blastn_ref$pident < 97.000,]
blastn_ref <- data.frame(blastn_ref_select$ASV, blastn_ref_select$Species)

#scientific name from blastn_Genbank
unique_sseqid <- unique(blastn_Genbank_select$sseqid)
# taxize::use_entrez()
# usethis::edit_r_environ()
ENTREZ_KEY <- "b05d96cb7cce251e9b81cdf04abecc36d409"
Classification_ncbi <- classification(genbank2uid(id = unique(blastn_Genbank_select$gi_id), 
                                                  key = ENTREZ_KEY), db = 'ncbi')
Classification_NCBI_1 <-cbind(Classification_ncbi)
Classification_NCBI_1$gi_id <- unique(blastn_Genbank_select$gi_id)

# select just which lineage level to keep
Classification_NCBI <- select(Classification_NCBI_1, gi_id, query, superkingdom, 
                              phylum, class, order, family, genus, species)
colnames(Classification_NCBI) <- c("gi_id", "query", "Kingdom", "Phylum", "Class", 
                                   "Order", "Family", "Genus" , "Species")

idx <- which(blastn_Genbank_select$gi_id %in% Classification_NCBI$gi_id)
blastn_taxa_top_lineage <-
  blastn_Genbank_select |> left_join(Classification_NCBI,
                                     by ="gi_id")
blastn_Genbank$blastn_Genbank_select.Species <- blastn_taxa_top_lineage$Species.y

Classification_NCBI <- Classification_NCBI[!duplicated(Classification_NCBI$Species),]

#merge 3 datasets according to ASV name
table_Species <- merge(TA_DADA2, blastn_Genbank, by.x='table_unrarefied.ASV', by.y='blastn_Genbank_select.ASV', all.x=TRUE)
table_Species_2 <- merge(table_Species, blastn_ref, by.x='table_unrarefied.ASV', by.y='blastn_ref_select.ASV', all.x=TRUE)
colnames(table_Species_2) <- c("ASV", "DADA2", "blastn_GenBank", "blastn_ref")
#table_Species_3 <- table_Species_2[!duplicated(table_Species_3$table_Species_2.ASV),]

#write full taxonomic assignment in separate column Full
table_Species_2$Full <- str_extract(table_Species_2$DADA2, "^.*(?=(_))") #remove GenBank identifier from species name
table_Species_2$blastn_ref <- str_extract(table_Species_2$blastn_ref, "^.*(?=(_))") #remove GenBank identifier from species name
table_Species_2[is.na(table_Species_2)] <- "NA"
table_Species_2$Full  <- ifelse(table_Species_2$DADA2=="NA", table_Species_2$blastn_ref, table_Species_2$Full) #add blastn results against custom reference database
table_Species_2$Full  <- ifelse(table_Species_2$Full=="Skeletonema pseudocostatum", table_Species_2$blastn_GenBank, 
                                ifelse(table_Species_2$Full=="NA", table_Species_2$blastn_GenBank, table_Species_2$Full)) #add blastn results against GenBank
table_Species_2[is.na(table_Species_2)] <- "NA"

# fix names so it will be recognized by worms
#Species_FULL<-data.frame("Original_full"=unique(table_Species_2$Full))
openxlsx::write.xlsx(as.data.frame(unique(table_Species_2$Full)), paste0(proj.path,"/NJ2022/Species_FULL_WORMS_Input.xlsx"))

# use https://www.marinespecies.org/aphia.php?p=match to match full taxonomy
worms_names <- readxl::read_excel(paste0(proj.path,"/NJ2022/Species_FULL_WORMS_Output.xlsx"),sheet = 'WoRMS match')
worms_names_exact <- worms_names[worms_names$`Match type` %in% "exact",]
worms_names_exact <- worms_names_exact %>% filter(str_detect(str_trim(ScientificName_accepted), "\\s+"))
#Species <- unique(table_Species_2$Full)
#Species <- Species[!Species == 'NA']
#write_tsv(as.data.frame(Species),paste0(proj.path,"/MiFish_UE-S_concatenated/results/Species_FULL_WORMS_Input.tsv"),col_names = F)
#Worms_id <- get_wormsid(Species, marine_only = F)
Worms_id <- unique(worms_names$AphiaID_accepted)
Worms_id <- na.omit(Worms_id)
Classification_worms <- classification(Worms_id, db = 'worms')
Classification_worms <- cbind(Classification_worms)
Classification_WoRMS <- select(Classification_worms, species_id, kingdom, phylum, class, order, family, genus, species)
colnames(Classification_WoRMS) <- c("species_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus" , "Species")
Classification_WoRMS <- Classification_WoRMS[!duplicated(Classification_WoRMS$species),]

#set unaccepted names to accepted names
for (r in 1:nrow(table_Species_2)){
  if(table_Species_2$Full[r] %in% worms_names_exact$`unique(table_Species_2$Full)`){
    table_Species_2$Full[r] <- worms_names_exact$ScientificName_accepted[which(worms_names_exact$`unique(table_Species_2$Full)`==table_Species_2$Full[r])]
  }
}

#remove species with identical 12S sequences 
Identical_12S <- c("Eutrigla gurnardus", "Chelidonichthys lucerna","Chelidonichthys cuculus","Chelidonichthys spinosus",
                   "Hyperoplus lanceolatus","Hyperoplus immaculatus", "Ammodytes tobianus", "Ammodytes marinus",
                   "Alosa alosa", "Alosa fallax")
for (r in 1:nrow(table_Species_2)){
  if(table_Species_2$Full[r] %in% Identical_12S){
    table_Species_2$Full[r] <- "NA"
  }
}

select_columns<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus" , "Species")

#add blastn results to the unconcatenated table
table_unrarefied$DADA2 <- table_unrarefied$Species
table_unrarefied$blastn_GenBank <- table_Species_2$blastn_GenBank
table_unrarefied$blastn_ref <- table_Species_2$blastn_ref
table_unrarefied$Full <- table_Species_2$Full
table_unrarefied[which(table_unrarefied$DADA2=="Skeletonema pseudocostatum_(MK372941)"),select_columns]<-NA

# Find rows where all columns match the specified values to remove inaccurate classification during the dada2 assignment
matching_rows <- which(table_unrarefied$Kingdom == "Animalia" &
                         (is.na(table_unrarefied$Phylum) | table_unrarefied$Phylum == "Chordata") &
                         (is.na(table_unrarefied$Class) | table_unrarefied$Class == "Actinopteri") &
                         is.na(table_unrarefied$Order) &
                         is.na(table_unrarefied$Family) &
                         is.na(table_unrarefied$Genus) &
                         is.na(table_unrarefied$Species))

# Set all matching rows in all columns to "NA"
table_unrarefied[matching_rows, select_columns] <- NA

# Set new classification based on NCBI and WoRMS
for (r in 1:nrow(table_unrarefied)){
  if(table_unrarefied$Full[r] %in% Classification_NCBI$Species){
    table_unrarefied[r,select_columns] <- Classification_NCBI[which(Classification_NCBI$Species==table_unrarefied$Full[r]),select_columns]
  }
}

for (r in 1:nrow(table_unrarefied)){
  if(table_unrarefied$Full[r] %in% Classification_WoRMS$Species){
    table_unrarefied[r,select_columns] <- Classification_WoRMS[which(Classification_WoRMS$Species==table_unrarefied$Full[r]),select_columns]
  }
}

#Set Eukaryota to Animalia, Plantae, Chromista, ... 
Animalia_Phyla <- c("Arthropoda", "Mollusca", "Chordata", "Chordata",
                    "Nematoda", "Annelida", "Cnidaria", "Porifera", 
                    "Echinodermata", "Bryozoa", "Rotifera", "Nemertea", "Tardigrada")
for (r in 1:nrow(table_unrarefied)){
  if(table_unrarefied$Phylum[r] %in% Animalia_Phyla){
    table_unrarefied$Kingdom[r] <- "Animalia"
  }
}

Chromista_Phyla <- c("Bacillariophyta", "Ochrophyta", "Cryptophyta", "Haptophyta",
                     "Heliozoa", "Chromeridophyta", "Ciliophora", "Myzozoa",
                     "Bigyra", "Ochrophyta", "Oomycota", "Cercozoa", "Foraminifera", "Radiozoa")
for (r in 1:nrow(table_unrarefied)){
  if(table_unrarefied$Phylum[r] %in% Chromista_Phyla){
    table_unrarefied$Kingdom[r] <- "Chromista"
  }
}

Plantae_Phyla <- c("Chlorophyta", "Glaucophyta", "Rhodelphidia", "Rhodophyta",
                   "Prasinodermophyta", "Anthocerotophyta", "Bryophyta",
                   "Charophyta", "Marchantiophyta", "Tracheophyta",
                   "Streptophyta")
for (r in 1:nrow(table_unrarefied)){
  if(table_unrarefied$Phylum[r] %in% Plantae_Phyla){
    table_unrarefied$Kingdom[r] <- "Plantae"
  }
}

#For values where the full taxonomy is empty but there's still a higher classification add the higher classification with SP to the full column
find_first_non_na_column <- function(df) {
  non_na_columns <- c()
  for (i in 1:nrow(df)) {
    if (df$Full[i] %in% c("NA",NA)) {
      for (column in names(df)[(ncol(df)-10):ncol(df)]) {
        if(is.na(df$Family[i])){next}
        if (column != "Full" && df[i, column] %in% c("NA",NA)) {
          NOT_NA_column<-which(colnames(df)==column)-1
          df[i,"Full"]<-paste0(df[i,NOT_NA_column]," sp.")
          print(paste0("On row ",i," Full was changed to the ",colnames(df)[NOT_NA_column]," column + sp."))
          break
        }
      }
    }
  }
  return(df)
}

table_unrarefied <- find_first_non_na_column(table_unrarefied)
table_unrarefied[is.na(table_unrarefied)] <- "NA"
table_unrarefied$Species <- table_unrarefied$Full
table_unrarefied$Full <- NULL

openxlsx::write.xlsx(as.data.frame(table_unrarefied), 
                     paste0(proj.path,"/5119441_NJ2022_Miseq/table_unrarefied_raw_Full_TaxAss_WoRMS.xlsx"))

#source https://www.britannica.com/animal/fish/Annotated-classification
Non_Fish_Chordata <- c("Appendicularia", "Ascidiacea", "Ascidiacea", "Thaliacea",
                       "Amphibia", "Aves", "Aves", "Reptilia", "Reptilia", "Mammalia")

table_unrarefied_fish <- table_unrarefied[table_unrarefied$Kingdom %in% "Animalia",]
table_unrarefied_fish <- table_unrarefied_fish[table_unrarefied_fish$Phylum %in% "Chordata",]
table_unrarefied_fish <- table_unrarefied_fish[!table_unrarefied_fish$Class %in% Non_Fish_Chordata,]

fish_classes <- c(unique(table_unrarefied_fish$Class))
saveRDS(fish_classes, file = paste0(proj.path,"/NJ2022/R_Environment/Fish_classes.rds"))

fish_classes_Joran<-c( "Agnatha" #superclass
                       ,"Myxini" #class
                       ,"Petromyzonti" #class
                       ,"Elasmobranchii"#class
                       ,"Neoselachii" #subclass
                       ,"Batoidea" #infraclass
                       ,"Selachii" #infraclass
                       ,"Holocephali" #class
                       ,"Osteichthyes" #superclass
                       ,"Actinopterygii" #gigaclass
                       ,"Actinopteri" #superclass
                       ,"Chondrostei" #class
                       ,"Holostei"    #class   
                       ,"Teleostei"  #class
                       ,"Sarcopterygii" #gigaclass
                       ,"Coelacanthimorpha" #class
                       ,"Dipneusti" #class          
)

saveRDS(fish_classes_Joran, file = paste0(proj.path,"/NJ2022/R_Environment/fish_classes_Joran.rds"))
 
#Assign Alosa sp.
table_unrarefied$Species <- ifelse(table_unrarefied$Genus=="Alosa", "Alosa sp.", table_unrarefied$Species)

table_assignedASVs <- as.data.frame(table_unrarefied[!table_unrarefied$Species %in% c("NA"),])
table_fish <-  as.data.frame(table_unrarefied[table_unrarefied$Class %in% fish_classes,])
table_animalia <- as.data.frame(table_unrarefied[table_unrarefied$Kingdom == "Animalia",])
table_plantae <- as.data.frame(table_unrarefied[table_unrarefied$Kingdom == "Plantae",])
table_chromista <- as.data.frame(table_unrarefied[table_unrarefied$Kingdom == "Chromista",])
table_bacteria <- as.data.frame(table_unrarefied[table_unrarefied$Kingdom == "Bacteria",])
table_eukaryota <- as.data.frame(table_unrarefied[table_unrarefied$Kingdom == "Eukaryota",])
table_archaea <- as.data.frame(table_unrarefied[table_unrarefied$Kingdom == "Archaea",])

#check which species assigned with BLAST
DADA2_NA <- table_fish[table_fish$DADA2=="NA",]
BLAST_ref <- DADA2_NA[!DADA2_NA$blastn_ref=="NA",] #blast against ref db identifies 81 additional ASVS.
BLAST_ncbi <- DADA2_NA[DADA2_NA$blastn_ref=="NA",]

unique(table_unrarefied$Kingdom)
table_rest <- table_unrarefied %>%
  filter(!(Kingdom %in% c("Chromista", "Plantae", "Animalia", "Bacteria","Eukaryota","Archaea")))
table_rest_species <- table_rest %>%
  filter(!(Species %in% "NA"))


seqtab_chromista <- as.data.frame(table_chromista[,1:(ncol(table_chromista)-11)])
TotalRead_perASV_ch <- data.frame(rowSums(seqtab_chromista)) 
colnames(TotalRead_perASV_ch)<- "Read"
TotalRead_perASV_ch <- cbind(TotalRead_perASV_ch, table_chromista$Species)
colnames(TotalRead_perASV_ch)[ncol(TotalRead_perASV_ch)] <- "Species"
sum(TotalRead_perASV_ch$Read) #1629486

seqtab_animalia <- as.data.frame(table_animalia[,1:(ncol(table_animalia)-11)])
TotalRead_perASV_ani <- data.frame(rowSums(seqtab_animalia)) 
colnames(TotalRead_perASV_ani)<- "Read"
TotalRead_perASV_ani <- cbind(TotalRead_perASV_ani, table_animalia$Species)
colnames(TotalRead_perASV_ani)[ncol(TotalRead_perASV_ani)] <- "Species"
sum(TotalRead_perASV_ani$Read) #2300504

seqtab_plantae <- as.data.frame(table_plantae[,1:(ncol(table_plantae)-11)])
TotalRead_perASV_pla <- data.frame(rowSums(seqtab_plantae)) 
colnames(TotalRead_perASV_pla)<- "Read"
TotalRead_perASV_pla <- cbind(TotalRead_perASV_pla, table_plantae$Species)
colnames(TotalRead_perASV_pla)[ncol(TotalRead_perASV_pla)] <- "Species"
sum(TotalRead_perASV_pla$Read) #38131

seqtab_bacteria <- as.data.frame(table_bacteria[,1:(ncol(table_bacteria)-11)])
TotalRead_perASV_bac <- data.frame(rowSums(seqtab_bacteria)) 
colnames(TotalRead_perASV_bac)<- "Read"
TotalRead_perASV_bac <- cbind(TotalRead_perASV_bac, table_bacteria$Species)
colnames(TotalRead_perASV_bac)[ncol(TotalRead_perASV_bac)] <- "Species"
sum(TotalRead_perASV_bac$Read) #1414524

seqtab_eukaryota <- as.data.frame(table_eukaryota[,1:(ncol(table_eukaryota)-11)])
TotalRead_perASV_euk <- data.frame(rowSums(seqtab_eukaryota)) 
colnames(TotalRead_perASV_euk)<- "Read"
TotalRead_perASV_euk <- cbind(TotalRead_perASV_euk, table_eukaryota$Species)
colnames(TotalRead_perASV_euk)[ncol(TotalRead_perASV_euk)] <- "Species"
sum(TotalRead_perASV_euk$Read) #651803

unique(table_eukaryota$Family)
# [1] "NA"                 "Chattonellaceae"    "Katablepharidaceae" "Hemiselmidaceae"   
# [5] "Heterocapsaceae"    "Dinophysaceae"      "Geminigeraceae"     "Pyrenomonadaceae"  
# [9] "Aspergillaceae"(fungus)     "Cladosporiaceae"(fungus)    "Paramoebidae"(protozoa)  rest is chromista

table_eukaryota_NA <- as.data.frame(table_eukaryota[table_eukaryota$Family == "NA",])
table_eukaryota_notNA <- as.data.frame(table_eukaryota[!table_eukaryota$Family == "NA",])
seqtab_eukaryota_notNA <- as.data.frame(table_eukaryota_notNA[,1:(ncol(table_eukaryota_notNA)-11)])
TotalRead_perASV_euk_notNA <- data.frame(rowSums(seqtab_eukaryota_notNA)) 
colnames(TotalRead_perASV_euk_notNA)<- "Read"
TotalRead_perASV_euk_notNA <- cbind(TotalRead_perASV_euk_notNA, table_eukaryota_notNA$Species)
colnames(TotalRead_perASV_euk_notNA)[ncol(TotalRead_perASV_euk_notNA)] <- "Species"
sum(TotalRead_perASV_euk_notNA$Read) #90097

#how many total read counts of fish ASVs
seqtab_fish <- as.data.frame(table_fish[,1:(ncol(table_fish)-11)])
fish_ASV <- data.frame(Reads = rowSums(seqtab_fish))
sum(fish_ASV)
unique(table_fish$Class) #3
unique(table_fish$Order) #26
unique(table_fish$Family) #41
unique(table_fish$Genus) #61
unique(table_fish$Species) #72 with 4 elements assigned in family and genus level

#Assess actinopteri class
table_actinopteri <- as.data.frame(table_fish[table_fish$Class=="Actinopteri",])
seqtab_actinopteri <- as.data.frame(table_actinopteri[,1:(ncol(table_actinopteri)-11)])
actinopteri_ASV <- data.frame(Reads=rowSums(seqtab_actinopteri))
sum(actinopteri_ASV$Reads) #2150509
unique(table_actinopteri$Order) #21
unique(table_actinopteri$Family) #36
unique(table_actinopteri$Genus) #55
unique(table_actinopteri$Species) #66

#Assess cartilaginous
table_cartil <- as.data.frame(table_fish[!table_fish$Class=="Actinopteri",])
seqtab_cartil <- as.data.frame(table_cartil[,1:(ncol(table_cartil)-11)])
cartil_ASV <- data.frame(Reads=rowSums(seqtab_cartil))
sum(cartil_ASV$Reads) #2150509
unique(table_cartil$Order) #4
unique(table_cartil$Family) #5
unique(table_cartil$Genus) #5
unique(table_cartil$Species) #6

#analyze negative controls for ASVs and sequence counts
rownames(table_fish) <- table_fish$ASV
table_fish_neg <- as.data.frame(table_fish[,grepl("^Neg_", colnames(table_fish))])
neg_ASV <- data.frame(Reads = rowSums(table_fish_neg))
neg_ASV <- filter(neg_ASV, Reads!=0)
neg_ASV <- merge(neg_ASV, table_fish, by = "row.names", all.x = TRUE)
neg_ASV <- select(neg_ASV, "Row.names", "Reads", "Species")
sum(neg_ASV$Reads) #30424
unique(neg_ASV$Species)
openxlsx::write.xlsx(as.data.frame(neg_ASV), 
                     paste0(proj.path,"/NJ2022/table_negative_controls_ASV.xlsx"))

table_fish_neg_filter <- as.data.frame(table_fish_neg[,grepl("^Neg_Filter", colnames(table_fish_neg))])
neg_filter_ASV <- data.frame(Reads = rowSums(table_fish_neg_filter))
neg_filter_ASV <- filter(neg_filter_ASV, Reads!=0)
neg_filter_ASV <- merge(neg_filter_ASV, table_fish, by = "row.names", all.x = TRUE)
neg_filter_ASV <- select(neg_filter_ASV, "Row.names", "Reads", "Species")
sum(neg_filter_ASV$Reads) #3409
unique(neg_filter_ASV$Species)
# [1] "Sardina pilchardus"     "Limanda limanda"        "Solea solea"            "Pleuronectes platessa" 
# [5] "Engraulis encrasicolus" "Pomatoschistus minutus" "Merlangius merlangus"   "Platichthys flesus"    
# [9] "Raja brachyura"   

table_fish_neg_field <- as.data.frame(table_fish_neg[,grepl("^Neg_Field", colnames(table_fish_neg))])
neg_field_ASV <- data.frame(Reads = rowSums(table_fish_neg_field))
neg_field_ASV <- filter(neg_field_ASV, Reads!=0)
neg_field_ASV <- merge(neg_field_ASV, table_fish, by = "row.names", all.x = TRUE)
neg_field_ASV <- select(neg_field_ASV, "Row.names", "Reads", "Species")
sum(neg_field_ASV$Reads) #13845
unique(neg_field_ASV$Species)
# [1] "Sardina pilchardus"      "Limanda limanda"         "Solea solea"             "Pleuronectes platessa"  
# [5] "Ammodytidae sp."         "Scomber scombrus"        "Trachurus trachurus"     "Mullus surmuletus"      
# [9] "Merlangius merlangus"    "Triglidae sp."           "Dicentrarchus labrax"    "Echiichthys vipera"     
# [13] "Syngnathus rostellatus"  "Scyliorhinus canicula"   "Spondyliosoma cantharus" "Buglossidium luteum"    
# [17] "Microstomus kitt"        "Raja clavata"            "Callionymus lyra"        "Zeus faber" 

table_fish_neg_DNA <- as.data.frame(table_fish_neg[,grepl("^Neg_DNA", colnames(table_fish_neg))])
neg_DNA_ASV <- data.frame(Reads = rowSums(table_fish_neg_DNA))
neg_DNA_ASV <- filter(neg_DNA_ASV, Reads!=0)
neg_DNA_ASV <- merge(neg_DNA_ASV, table_fish, by = "row.names", all.x = TRUE)
neg_DNA_ASV <- select(neg_DNA_ASV, "Row.names", "Reads", "Species")
sum(neg_DNA_ASV$Reads) #13170
unique(neg_DNA_ASV$Species)
# [1] "Limanda limanda"           "Pleuronectes platessa"     "Merlangius merlangus"     
# [4] "Trisopterus luscus"        "Lepidorhombus wiffiagonis"

table_fish_neg_PCR <- as.data.frame(table_fish_neg[,grepl("^Neg_PCR", colnames(table_fish_neg))])
neg_PCR_ASV <- data.frame(Reads = rowSums(table_fish_neg_PCR))
neg_PCR_ASV <- filter(neg_PCR_ASV, Reads!=0)
sum(neg_PCR_ASV$Reads) #0
unique(neg_PCR_ASV$Species) #no fish ASV

