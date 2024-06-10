#set OS type for paths
if(.Platform$OS.type == "unix"){home<-"/home/"} else{
  home<-"//192.168.236.131/"}

#activate libraries
libraries <- c("BiocManager"
               , "ggplot2"
               , "tidyr"
               , "taxonomizr"
               , "naturalsort"
               , "stringr"
               , "phylotools"
               , "scales"
               , "ggpattern"
               , "ggh4x"
               , "reshape2"
               , "vegan"
               , "here"
               , "dplyr"
               , "phyloseq"
               , "ggpubr"
               , "metagMisc"
               , "iNEXT"
               , "ggClusterNet"
)

for (opties in libraries){
  
  if (opties %in% installed.packages()){
    
    library(opties,character.only = TRUE)
    
  } else {install.packages(opties,repos = "http://cran.us.r-project.org")
    
    library(opties,character.only = TRUE)
  }
}

proj.path <- here("/home/genomics/ndukan/01_Zero-Impact/Temporal_Analysis")

table_raw <- readxl::read_excel(paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean_MicroDecon.xlsx"))
env <- readRDS(paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/R_Environment_MicroDecon/env_unrarefied_MicroDecon.rds"))
colorder <- c(env$Niskin.sample)
env_fish <- readRDS(paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/R_Environment_MicroDecon/env_unrarefied_MicroDecon_Fish.rds"))
colorder_fish <- c(env_fish$Niskin.sample)
community_NJ2022 <- community_NJ2022[, colorder_fish]

# #check for the sequencing errors in samples
# seqtab_raw <- table_raw[,1:(ncol(table_raw)-11)]
# seqtab_raw <- seqtab_raw[,colorder]
# singleton_data <- as.data.frame(colnames(seqtab_raw))
# names(singleton_data) <- "Sample"
# singleton_data$Sample_size <- colSums(seqtab_raw)
# singleton_data$NumOfSingleton <- apply(seqtab_raw, 2, function(column) sum(column == 1))
# 
# ggplot(singleton_data, aes(x = Sample_size, y = NumOfSingleton)) +
#   geom_point() +  # This adds the scatter plot points
#   # geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Adds a linear regression line without the standard error
#   theme_minimal() +  # Uses a minimal theme for the plot
#   labs(
#     title = "Number of Singletons as a Function of Sample Size",
#     x = "Sample Size",
#     y = "Number of Singletons"
#   )


# ## Create a phyloseq object
# ps_eDNA <- table_raw[,colnames(table_raw) %in% env$Niskin.sample]
# smpl_eDNA <- env[env$Niskin.sample  %in%  colnames(ps_eDNA),]
# smpl_eDNA <- as.data.frame(smpl_eDNA)
# rownames(smpl_eDNA) <- colnames(ps_eDNA)
# ps_eDNA <- as.data.frame(ps_eDNA)
# rownames(ps_eDNA) <- table_raw$ASV
# Taxonomy_eDNA <- as.matrix(rownames(ps_eDNA))
# rownames(Taxonomy_eDNA) <- rownames(ps_eDNA)
# ps_unrarefied <- phyloseq(otu_table(ps_eDNA, taxa_are_rows = TRUE),
#                           sample_data(smpl_eDNA),
#                           tax_table(Taxonomy_eDNA))


###Fish ASV level Coverage

fish_classes <- readRDS(file = paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/Fish_classes.rds"))
table_unrarefied_fish <- as.data.frame(table_raw[table_raw$Class %in% fish_classes,])
freshwater_species <- c("Alburnus alburnus", "Barbatula barbatula", "Gobio gobio", "Rutilus rutilus", "Squalius cephalus", "Oreochromis niloticus", "Cyprinus carpio haematopterus x Megalobrama amblycephala")
table_unrarefied_fish <- table_unrarefied_fish[!(table_unrarefied_fish$Species %in% freshwater_species), ]
rownames(table_unrarefied_fish) <- table_unrarefied_fish$ASV
table_unrarefied_fish <- table_unrarefied_fish[, !colSums(table_unrarefied_fish[, 1:(ncol(table_unrarefied_fish)-11)])==0]
table_unrarefied_fish <- table_unrarefied_fish[!rowSums(table_unrarefied_fish[, 1:(ncol(table_unrarefied_fish)-11)])==0,]
seqtab_unrarefied_fish <- table_unrarefied_fish[,1:(ncol(table_unrarefied_fish)-11)]
singletons_ASV <- as.data.frame(rowSums(seqtab_unrarefied_fish)) #(no singleton fish ASV, min 3)
seqtab_unrarefied_fish <- seqtab_unrarefied_fish[,colorder_fish]

## Create a phyloseq object
ps_ASV <- seqtab_unrarefied_fish
smpl_ASV <- env[env$Niskin.sample  %in%  colnames(ps_ASV),]
smpl_ASV <- as.data.frame(smpl_ASV)
rownames(smpl_ASV) <- colnames(ps_ASV)
Taxonomy_ASV <- as.matrix(rownames(ps_ASV))
rownames(Taxonomy_ASV) <- rownames(ps_ASV)
ps_unrarefied_ASV <- phyloseq(otu_table(ps_ASV, taxa_are_rows = TRUE),
                          sample_data(smpl_ASV),
                          tax_table(Taxonomy_ASV))

inext_ASV_fish <- iNEXT(seqtab_unrarefied_fish, q=0, datatype="abundance", endpoint = NULL)
inext_coverage <- as.data.frame(inext_ASV_fish[["iNextEst"]][["coverage_based"]])

inext_ASV_fish_shannon <- iNEXT(seqtab_unrarefied_fish, q=1, datatype="abundance", endpoint = NULL)
inext_coverage_shannon <- as.data.frame(inext_ASV_fish_shannon[["iNextEst"]][["coverage_based"]])

##Calculate the Coverage
Coverage_ASV <- phyloseq_coverage(physeq = ps_unrarefied_ASV, correct_singletons = FALSE)

##Rarefy the data based on a coverage just below the minimum coverage using the function phyloseq_coverage_raref,
ps_rarefied_ASV <- phyloseq_coverage_raref(physeq = ps_unrarefied_ASV, 
                                            iter = 1, coverage = 0.983, drop_lowcoverage = T)
table_rarefied_ASV <- as.data.frame(ps_rarefied_ASV@otu_table)
table_rarefied_ASV$Species <- table_unrarefied_fish$Species

#prepare the community matrix of rarefied data
taxo <- "Species"
table_rarefied_merged_fish <- aggregate(table_rarefied_ASV[,1:(ncol(table_rarefied_ASV)-1)], 
                                   by= list(as.factor(table_rarefied_ASV[,taxo])),FUN=sum)
rownames(table_rarefied_merged_fish) <- as.character(table_rarefied_merged_fish$Group.1)
table_rarefied_merged_fish$Group.1 <- NULL
table_rarefied_merged_fish <- table_rarefied_merged_fish[!rowSums(table_rarefied_merged_fish) == 0,]
table_rarefied_merged_fish <- table_rarefied_merged_fish[,!colSums(table_rarefied_merged_fish) == 0]
table_rarefied_merged_fish <- table_rarefied_merged_fish[!rownames(table_rarefied_merged_fish)=="NA",]
table_rarefied_merged_fish <- as.data.frame(t(table_rarefied_merged_fish))

#remove the samples from env data that are not present in the rarefied data
keep_samples <- c(rownames(table_rarefied_merged_fish))
env_rarefied <- env_fish[env_fish$Niskin.sample %in% keep_samples,]
colorder_rarefied <- c(env_rarefied$Niskin.sample)
table_rarefied_merged_fish <- table_rarefied_merged_fish[colorder_rarefied,]

#Plot richness and Shannon with rarefied data
#transform the data
table_rarefied_merged_fish_transformed <- decostand(table_rarefied_merged_fish, method="total")
table_rarefied_merged_fish_transformed <- decostand(table_rarefied_merged_fish_transformed, method="max")

#calculate Shannon
shannon <- data.frame(diversity(t(table_rarefied_merged_fish), index = "shannon"))
colnames(shannon) <- "Shannon"

shannon_transformed <- data.frame(diversity(t(table_rarefied_merged_fish_transformed), index = "shannon")) 
colnames(shannon_transformed) <- "Shannon_transformed"

env_rarefied$Shannon <- shannon$Shannon
env_rarefied$Shannon_transformed <- shannon_transformed$Shannon_transformed

qqPlot(env_rarefied$Shannon)
qqPlot(env_rarefied$Shannon_transformed)

env_rarefied$Zone <- factor(env_rarefied$Zone, levels = c("Coast", "Transition", "Offshore"))

Plot_shannon <- ggplot(env_rarefied, aes(x = Zone, y = Shannon, fill= Depth)) +
  geom_boxplot(width=0.5, alpha=0.8, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.07), aes(color = Depth), size = 2) + # Adds points, aligned with boxes
  scale_fill_manual(values = c("steelblue3", "lightskyblue1")) +
  scale_color_manual(values = c("steelblue4", "lightskyblue")) + 
  labs(x = "", y = "Shannon Diversity Index") +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5, size=30),
        axis.title.y = element_text(color = "black", size=26, hjust = 0.5, vjust = 1.5),
        axis.text.y = element_text(size = 26), 
        axis.text.x = element_text(size = 26),
        legend.text = element_text(color = "black",size = 24),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA)) 
Plot_shannon

#calculate richness
env_rarefied$NumOfSp <- rowSums(table_rarefied_merged_fish>0)
qqPlot(env_rarefied$NumOfSp)

env_rarefied$Zone <- factor(env_rarefied$Zone, levels = c("Coast", "Transition", "Offshore"))

Plot_richness <- ggplot(env_rarefied, aes(x = Zone, y = NumOfSp, fill= Depth)) +
  geom_boxplot(width=0.5, alpha=0.8, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.07), aes(color = Depth), size = 2) + # Adds points, aligned with boxes
  scale_fill_manual(values = c("steelblue3", "lightskyblue1")) +
  scale_color_manual(values = c("steelblue4", "lightskyblue")) + 
  labs(x = "", y = "Species Richness") +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5, size=30),
        axis.title.y = element_text(color = "black", size=26, hjust = 0.5, vjust = 1.5),
        axis.text.y = element_text(size = 26), 
        axis.text.x = element_text(size = 26),
        legend.text = element_text(color = "black",size = 24),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA)) 
Plot_richness

png("Shannon_Richness_MicroDecon_CoverageBased_rarefied.png",width=4000,height=2000,units="px",res=300,bg="white")
ggarrange(Plot_shannon, Plot_richness,
          ncol = 2, nrow = 1,  align = "hv", 
          labels = c("A", "B"),
          widths = c(1, 1), heights = c(1, 1),
          common.legend = TRUE,
          # legend = "right",
          font.label = list(size=22))
dev.off()

#do the statistical analysis
#Shannon
par(mfrow = c(1,1))
qqPlot(env_rarefied$Shannon) #not transformed is better
qqPlot(env_rarefied$Shannon_transformed)

library(lmerTest)

# model_shannon <- lm(Shannon_transformed ~ Depth*Zone, data=env_rarefied)
model_shannon <- lm(Shannon ~ Depth*Zone, data=env_rarefied) 
model_shannon <- rlm(Shannon ~ Depth*Zone, data=env_rarefied)  #less sensitive to outliers
anova_shannon <- Anova(model_shannon)
anova_shannon

aov <- aov(Shannon ~ Zone, data=env_rarefied)
aov
summary(aov)
TukeyHSD(aov, "Zone", ordered=FALSE, conf.level = 0.95) #sample size should be equal

model_shannon_mixed <- lmer(Shannon ~ Depth * Zone + (1 | Location), data=env_rarefied)
anova_shannon <- Anova(model_shannon_mixed)
anova_shannon

par(mfrow = c(2,2))
plot(model_shannon)

kruskal.test(Shannon ~ Zone, data=env_rarefied)
kruskal.test(Shannon ~ Depth, data=env_rarefied)

# library(rstatix)
# dunn_test(Shannon ~ Zone, data=env_rarefied, p.adjust.method = "BH")

# wilcox.test(Shannon ~ Depth, data=env_rarefied) #factor needs exactly 2 levels

# model_shannon_glm <- glm(Shannon ~ Depth * Zone, family = poisson, data = env_rarefied)
# anova_shannon <- Anova(model_shannon_glm)
# anova_shannon

# aov <- aov(Shannon_transformed ~ Zone, data=env_rarefied)
# aov
# summary(aov)
# TukeyHSD(aov, "Zone", ordered=FALSE, conf.level = 0.95) #sample size should be equal

# # Run the Kruskal-Wallis test
# kruskal_result <- kruskal.test(Shannon ~ Depth*Zone, data = env)
# print(kruskal_result)

#Species richness
par(mfrow = c(1,1))
qqPlot(env_rarefied$NumOfSp)

# model_richness <- rlm(NumOfSp ~ Depth*Zone, data=env)  
model_richness <- lm(NumOfSp ~ Depth*Zone, data=env_rarefied) 
anova_richness <- Anova(model_richness)
anova_richness 
summary(model_richness)

par(mfrow = c(2,2))
plot(model_richness)

aov <- aov(NumOfSp ~ Zone, data=env_rarefied)
aov
summary(aov)
TukeyHSD(aov, "Zone", ordered=FALSE, conf.level = 0.95) #sample size should be equal


#NMDS for rarefied data
data.nmds <- decostand(table_rarefied_merged_fish_transformed, method="total")
data.nmds <- decostand(data.nmds, method="max")

ord.NMDS=metaMDS(data.nmds, k=2, distace ="jaccard", trymax=100) #stress: 0.1802 transformed
stressplot(ord.NMDS)

Scores_nmds <- scores(ord.NMDS)
Scores_nmds <- Scores_nmds[["sites"]]
Scores_nmds <- data.frame(Scores_nmds)
Scores_nmds$Niskin.sample <- rownames(Scores_nmds)
Scores_nmds_joined <- inner_join(Scores_nmds, env_rarefied, by="Niskin.sample")

library(ggforce)
pdf("NMDS_MicroDecon_CoverageBasedRarefied_DoubleTrans_Jaccard.pdf",width=20,height=15,bg="white")
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

# #TO DO
# 1) Merge the species table
# 2) Calculate Shannon and richness and plot
# 3) Make a plot for hill numbers
# 4) Do the NMDS with the rarefied data (bray curtis rarefied unrarefied, jaccard rarefied unrarefied)


###Species level coverage
community_NJ2022 <- readxl::read_excel(paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Table_FishSpecies_MicroDecon.xlsx"))
community_NJ2022 <- as.data.frame(community_NJ2022)
rownames(community_NJ2022) <- community_NJ2022[,1]
community_NJ2022[,1] <- NULL
singletons <- as.data.frame(rowSums(community_NJ2022))

## Create a phyloseq object
ps_fish <- community_NJ2022[,colnames(community_NJ2022) %in% env_fish$Niskin.sample]
smpl_fish <- env_fish[env_fish$Niskin.sample  %in%  colnames(ps_fish),]
smpl_fish <- as.data.frame(smpl_fish)
rownames(smpl_fish) <- colnames(ps_fish)
ps_fish <- as.data.frame(ps_fish)
Taxonomy_fish <- as.matrix(rownames(ps_fish))
rownames(Taxonomy_fish) <- rownames(ps_fish)
ps_unrarefied_fish <- phyloseq(otu_table(ps_fish, taxa_are_rows = TRUE),
                          sample_data(smpl_fish),
                          tax_table(Taxonomy_fish))

#q=0
inext_fish <- iNEXT(community_NJ2022, q=0, datatype="abundance", endpoint = NULL)
inext_coverage_fish <- as.data.frame(inext_fish[["iNextEst"]][["coverage_based"]])

#q=1 
inext_fish_shannon <- iNEXT(community_NJ2022, q=1, datatype="abundance", endpoint = NULL)
inext_coverage_fish_shannon <- as.data.frame(inext_fish_shannon[["iNextEst"]][["coverage_based"]])

##Calculate the Coverage
Coverage_fish <- phyloseq_coverage(physeq = ps_unrarefied_fish, correct_singletons = FALSE)

##Estimate the diversity
Fish_Diversity <- estimateD(ps_fish, datatype="abundance", base="coverage",
                                 level=0.97, q=1,
                                 conf=0.95) #select the minimum coverage

##Rarefy the data based on a coverage just below the minimum coverage using the function phyloseq_coverage_raref,
##(using the minimum coverage will remove the sample from the rarefied dataset)
ps_rarefied_fish <- phyloseq_coverage_raref(physeq = ps_unrarefied_fish, 
                                            iter = 1, coverage = 0.98, drop_lowcoverage = T)
table_rarefied_fish <- as.data.frame(ps_rarefied_fish@otu_table)

save.image((paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/R_Environment_MicroDecon/R_Environment_CoverageBasedRarefaction_AlphaDiversity.RData")))

