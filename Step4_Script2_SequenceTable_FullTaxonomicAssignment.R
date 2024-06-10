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
                   "Alosa alosa", "Alosa fallax", "Chelon ramada", "Chelon auratus", "Chelon labrosus")
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
