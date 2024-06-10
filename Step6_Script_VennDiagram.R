if(.Platform$OS.type == "unix"){home<-"/home/"} else{
  home<-"//192.168.236.131/"}

library(VennDiagram)
library(here)
library(openxlsx)

proj.path <- here("/home/genomics/ndukan/01_Zero-Impact/Temporal_Analysis")

#community matrix
community_NJ2022 <- readxl::read_excel(paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Table_FishSpecies_MicroDecon.xlsx"))
community_NJ2022 <- as.data.frame(community_NJ2022)
rownames(community_NJ2022) <- community_NJ2022[,1]
community_NJ2022[,1] <- NULL

data_morph <- readxl::read_excel(paste0(proj.path,"/02_NJ2022_Miseq_5119441/Morpholgical_abundance_matrix_STANDARDIZED.xlsx"),sheet = "Abundance matrix")
data_morph <- as.data.frame(data_morph)
rownames(data_morph) <- data_morph[,1]
data_morph[, 1] <- NULL

env <- readRDS(paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/env_unrarefied_MicroDecon.rds"))

#Venn diagram for all surface_bottom and morphology samples
env_keep_surface <- env[env$Depth=="Surface",]
env_keep_bottom <- env[env$Depth=="Bottom",]

samples_s <- c(env_keep_surface$Niskin.sample)
samples_b <- c(env_keep_bottom$Niskin.sample)

#subset eDNA data
samples_fish_s <- community_NJ2022[,colnames(community_NJ2022) %in% samples_s]
samples_fish_ss <- samples_fish_s[!rowSums(samples_fish_s)==0,]

samples_fish_b <- community_NJ2022[,colnames(community_NJ2022) %in% samples_b]
samples_fish_bb <- samples_fish_b[!rowSums(samples_fish_b)==0,]

#Create the venn diagram for surface_bottom_morphology data
eDNA_species_s <- rownames(samples_fish_ss)
eDNA_species_b <- rownames(samples_fish_bb)

# morph_species <- c(data_morph$SAM.SpeciesNameScientific)
colnames(data_morph)[colnames(data_morph) == "Pomatoschistus"] <- "Pomatoschistus sp."
morph_species <- c(colnames(data_morph))

counts1 <- c(length(eDNA_species_s), length(eDNA_species_b), length(unique(morph_species)))
counts1

pdf("Venn_numbers_NJ2022_MicroDecon.pdf",width=20,height=18,bg="white")

v <- venn.diagram(x=list(eDNA_species_s, eDNA_species_b, morph_species), category.names = c("eDNA Surface", "eDNA Bottom", "Morphology"),
                  fill = c("#41a48a", "#418393", "#d96648"), filename = NULL,
                  cex = 4, cat.cex = 3, fontfamily = "sans", cat.fontfamily = "sans", lwd=1, compression= "lzw",cat.fontface = "bold", radius=1,
                  main= "",main.fontfamily = "sans", main.cex=2, ext.text=TRUE, ellipse=FALSE)

grid.newpage()
grid.draw(v)

dev.off()

lapply(v, names)
lapply(v, function(i) i$label)

v[[8]]$label <- paste(setdiff(eDNA_species_s, union(eDNA_species_b, morph_species)), collapse="\n")
v[[11]]$label <- paste(setdiff(intersect(eDNA_species_s, eDNA_species_b), morph_species), collapse="\n")
v[[12]]$label <- paste(setdiff(eDNA_species_b, union(eDNA_species_s, morph_species)), collapse="\n")
v[[10]]$label <- paste(intersect(eDNA_species_s, intersect(eDNA_species_b, morph_species)), collapse="\n")
v[[9]]$label <- paste(setdiff(intersect(eDNA_species_b, morph_species), eDNA_species_s), collapse="\n")
v[[7]]$label <- paste(setdiff(morph_species, union(eDNA_species_s, eDNA_species_b)), collapse="\n")

# v[[10]]$label <- paste(setdiff(intersect(eDNA_species_s, morph_species), eDNA_species_b), collapse="\n")

grid.newpage(v)
grid.draw(v)


save.image((paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/R_Environment_MicroDecon/R_Environment_VennDiagram_MicroDecon.RData")))

species <- c("Cyclopterus lumpus", "Diplecogaster bimaculata", "Gasterosteus aculeatus", "Lophius piscatorius", "Pomatoschistus pictus", "Zeugopterus regius", "Anguilla anguilla", "Neogobius melanostomus", "Salmo salar", "Scophthalmus rhombus", "Trisopterus minutus", "Atherina presbyter", "Buglossidium luteum", "Trachinus draco")
species_sb <- as.data.frame(community_NJ2022[rownames(community_NJ2022) %in% species,])
a <- as.data.frame(community_NJ2022[grep("anguilla", rownames(community_NJ2022)),])
a <- as.data.frame(a[,!a==0])

ReadCounts <- as.data.frame(rowSums(species_sb[1:nrow(species_sb),])) #number of reads
names(ReadCounts)<- "Read_Counts"
Samples <- as.data.frame(rowSums(species_sb !=0)) #number of samples
names(Samples)<- "Sample_Counts"

table_species_sb <- cbind(ReadCounts, Samples)

write.xlsx(table_species_sb, 
           paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Species_Only_Surface_Bottom.xlsx"), 
           sheetName = "Sheet1", colNames = TRUE, rowNames = TRUE, append = FALSE)

