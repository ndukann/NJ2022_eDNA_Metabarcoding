if(.Platform$OS.type == "unix"){home<-"/home/"} else{
  home<-"//192.168.236.131/"}
source("/home/genomics/icornelis/03_RawScripts/funs-libs.R")

proj.path <- here("/home/genomics/ndukan/01_Zero-Impact/Temporal_Analysis")

library(vegan)
library(seqRFLP)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(openxlsx)
library(pairwiseAdonis)
library(MASS)
library(lmtest)

#community matrix
community_NJ2022 <- readxl::read_excel(paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Table_FishSpecies_MicroDecon.xlsx"))
community_NJ2022 <- as.data.frame(community_NJ2022)
rownames(community_NJ2022) <- community_NJ2022[,1]
community_NJ2022[,1] <- NULL
community_NJ2022 <- as.data.frame(t(community_NJ2022))

env <- readRDS(paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/env_unrarefied_MicroDecon.rds"))

colorder <- c(env$Niskin.sample)
community_NJ2022 <-community_NJ2022[colorder,]

#transform the data
community_NJ2022_transformed <- decostand(community_NJ2022, method="total")
community_NJ2022_transformed <- decostand(community_NJ2022_transformed, method="max")

##Permanova

#calculate bray curtis
community_NJ2022_transformed_bray <- vegdist(community_NJ2022_transformed, method = "bray",binary=TRUE) 
community_NJ2022_bray <- vegdist(community_NJ2022, method = "bray",binary=TRUE) 
bray_matrix <- as.matrix(community_NJ2022_transformed_bray)

Env_p <- env$Zone
Depth_p <- env$Depth
Description_p <- env$Description

adonis2(community_NJ2022_transformed_bray~Depth*Zone*Depth/Zone, env,  permutations = 9999, sqrt.dist=FALSE, by="terms")
pairwise.adonis2(community_NJ2022_transformed_bray~Zone, p.adjust.m="BH", reduce=NULL, env, perm=9999)
pairwise.adonis2(community_NJ2022_transformed_bray~Depth, p.adjust.m="BH", reduce=NULL, env, perm=9999)

par(mfrow=c(1,1))
Combined_betadisper_transformed <- betadisper(community_NJ2022_transformed_bray, Env_p, type="centroid")
Combined_betadisper_transformed <- betadisper(community_NJ2022_transformed_bray, Depth_p, type="centroid")
Combined_betadisper_transformed <- betadisper(community_NJ2022_transformed_bray, Description_p, type="centroid")

Combined_permdisp_transformed <- permutest(Combined_betadisper_transformed, permutations=9999)
Combined_permdisp_transformed 

par(mar=c(5, 6, 4, 4)) #to make the plot area smaller so y-lab is visible

png("Betadisper_plot_MicroDecon.png",width=3000,height=3000,units="px",res=400,bg="white")
plot(Combined_betadisper_transformed, hull=FALSE, ellipse=TRUE, cex.axis=1.5, cex.lab=1.5, main="")
dev.off()

save.image((paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/R_Environment_Anova_Permanova_MicroDecon.RData")))

