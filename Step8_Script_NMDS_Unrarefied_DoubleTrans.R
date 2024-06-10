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
library(phyloseq)
library(car)
library(lsmeans)
library(multcomp)
library(ggpubr)


libraries <- c("BiocManager"
               , "tidyr"
               , "taxonomizr"
               , "naturalsort"
               ,"stringr"
               , "phylotools"
               , "scales"
               , "ggpattern"
               , "ggh4x"
               , "reshape2"
)

for (opties in libraries){
  
  if (opties %in% installed.packages()){
    
    library(opties,character.only = TRUE)
    
  } else {install.packages(opties,repos = "http://cran.us.r-project.org")
    
    library(opties,character.only = TRUE)
  }
}

#community matrix
community_NJ2022 <- readxl::read_excel(paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Table_FishSpecies_MicroDecon.xlsx"))
community_NJ2022 <- as.data.frame(community_NJ2022)
rownames(community_NJ2022) <- community_NJ2022[,1]
community_NJ2022[,1] <- NULL
community_NJ2022 <- as.data.frame(t(community_NJ2022))
# community_NJ2022 <- community_NJ2022[, !colSums(community_NJ2022)==0]
# community_NJ2022 <- community_NJ2022[!rowSums(community_NJ2022)==0,]


env <- readRDS(paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/env_unrarefied_MicroDecon.rds"))
env$Depth <- ifelse(grepl("bottom", env$Niskin.sample), "Bottom", "Surface")
env$Description <- paste(env$Zone, env$Depth, sep = "_")
env$Description_color <- ifelse(env$Description=="Coast_Surface", "seagreen3",
                                ifelse(env$Description=="Coast_Bottom", "seagreen4",
                                       ifelse(env$Description=="Transition_Surface", "steelblue3",
                                              ifelse(env$Description=="Transition_Bottom", "steelblue4",
                                                     ifelse(env$Description=="Offshore_Surface", "darkorange",
                                                            ifelse(env$Description=="Offshore_Bottom", "darkorange4", env$Description))))))

saveRDS(env, file = paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/env_unrarefied_MicroDecon.rds"))

colorder <- c(env$Niskin.sample)

community_NJ2022 <- community_NJ2022[colorder,]

data.nmds <- decostand(community_NJ2022, method="total")
data.nmds <- decostand(data.nmds, method="max")



ord.NMDS=metaMDS(community_NJ2022, k=2, distace ="bray", trymax=100) #stress: 0.1654
ord.NMDS=metaMDS(data.nmds, k=2, distace ="bray", trymax=100) #stress: 0.1654

Scores_nmds <- scores(ord.NMDS)
Scores_nmds <- Scores_nmds[["sites"]]
Scores_nmds <- data.frame(Scores_nmds)
Scores_nmds$Niskin.sample <- rownames(Scores_nmds)
Scores_nmds_joined <- inner_join(Scores_nmds, env, by="Niskin.sample")

library(ggforce)
pdf("NMDS_MicroDecon_Unrarefied_NotTrans_Bray.pdf",width=20,height=15,bg="white")
# png("NMDS_Temporal_3Campaign_Miseq.png",width=4000,height=2000,units="px",res=300,bg="white")

NMDS_eDNA <- ggplot(Scores_nmds_joined, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(fill= Description), pch=21, size=6) +
  geom_mark_hull(aes(fill=factor(Description)), concavity=10, color="transparent", alpha= 0.2) +
  scale_fill_manual(values=c("seagreen3", "seagreen4","steelblue3","steelblue4" ,"darkorange","darkorange4"), drop=FALSE,
                    limits=c("Coast_Surface", "Coast_Bottom", "Transition_Surface", "Transition_Bottom","Offshore_Surface", "Offshore_Bottom"),
                    labels = c("Coast-Surface", "Coast-Bottom", "Transition-Surface", "Transition-Bottom","Offshore-Surface", "Offshore-Bottom")) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=30, colour="black"),
        axis.text.y = element_text(colour = "black", size = 26), 
        axis.text.x = element_text(colour = "black", size = 26), 
        legend.text = element_text(size = 30,  colour ="black"), 
        legend.position = "right", 
        legend.title = element_text(size = 30, colour="black", face="bold"),
        axis.title.y = element_text( size = 26, colour="black", angle=90), 
        axis.title.x = element_text( size = 26, colour = "black"), 
        panel.background = element_blank(), 
        # panel.background = element_rect(fill = "white", colour = "white", size = 1.2),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
        # legend.key=element_blank() +
        legend.key = element_rect("white")) +  
  labs(x = "NMDS1", fill="Zone-Depth", colour = "Description", y = "NMDS2")+
  xlim(-1.5, 1.5)+
  ylim(-0.75, 0.90)
  # geom_text(aes(label = Niskin.sample), size = 5, hjust = 0.5, vjust = -0.5, color = "black")
NMDS_eDNA

dev.off()

#indicator species analysis for eDNA analysis
library(indicspecies)
community_eDNA <- data.nmds
groups <- c(env$Zone)
indval = multipatt(community_eDNA, groups, duleg = TRUE,
                   control = how(nperm=9999)) 
indval_all <- indval[["sign"]]
summary(indval, indvalcomp = TRUE)

#indicator species analysis for trawl data
data_morph <- readxl::read_excel(paste0(proj.path,"/02_NJ2022_Miseq_5119441/Morpholgical_abundance_matrix_STANDARDIZED.xlsx"),sheet = "Abundance matrix")
data_morph <- as.data.frame(data_morph)
rownames(data_morph) <- data_morph[,1]
data_morph[, 1] <- NULL

groups2 <- c("Offshore", "Offshore", "Transition", "Transition", "Coast", "Transition", "Offshore", "Coast", "Transition", "Coast", "Coast", "Transition", "Transition", "Coast", "Transition", "Transition", "Coast")
indval2 = multipatt(data_morph, groups2, duleg = TRUE,
                   control = how(nperm=9999)) 
indval3 = multipatt(data_morph, groups2,
                    control = how(nperm=9999)) 
summary(indval3, indvalcomp = TRUE)
indval_all2 <- indval2[["sign"]]
indval_all2
summary(indval2, indvalcomp = TRUE)


save.image((paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/R_Environment_MicroDecon/R_Environment_NMDS_MicroDecon_Unrarefied_DoubleTrans.RData")))
 
