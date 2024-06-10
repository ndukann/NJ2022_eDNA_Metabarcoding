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

table_cleaned <- readxl::read_excel(paste0(proj.path, "/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean_MicroDecon.xlsx"))
seqtab_cleaned <- as.data.frame(table_cleaned[,1:(ncol(table_cleaned)-11)])
seqtab_cleaned <- as.data.frame(t(seqtab_cleaned))

setdiff(colnames(table_cleaned), colnames(seqtab_cleaned))
setdiff(rownames(NumberOfSequencesConcat), rownames(community_NJ2022))
NumberOfSequencesConcat <- as.data.frame(sort.default(rowSums(seqtab_cleaned[1:nrow(seqtab_cleaned),]))) #number of sequences per sample after concatenation, 8 samples <300 reads
NumberOfASVsConcat <- as.data.frame(sort.default(rowSums(seqtab_cleaned !=0))) #number of ASVs per sample after concatenation

#community matrix
community_NJ2022 <- readxl::read_excel(paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Table_FishSpecies_MicroDecon.xlsx"))
community_NJ2022 <- as.data.frame(community_NJ2022)
rownames(community_NJ2022) <- community_NJ2022[,1]
community_NJ2022[,1] <- NULL
community_NJ2022 <- as.data.frame(t(community_NJ2022))

env <- readRDS(paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/env_unrarefied_MicroDecon.rds"))

colorder <- c(env$Niskin.sample)
community_NJ2022 <-community_NJ2022[colorder,]

#transfor the data
community_NJ2022_transformed <- decostand(community_NJ2022, method="total")
community_NJ2022_transformed <- decostand(community_NJ2022_transformed, method="max")

#calculate Shannon
shannon <- data.frame(diversity(t(community_NJ2022), index = "shannon"))
colnames(shannon) <- "Shannon"

shannon_transformed <- data.frame(diversity(t(community_NJ2022_transformed), index = "shannon")) 
colnames(shannon_transformed) <- "Shannon_transformed"

env$Shannon <- shannon$Shannon
env$Shannon_transformed <- shannon_transformed$Shannon_transformed

qqPlot(env$Shannon)
qqPlot(env$Shannon_transformed)

env$Zone <- factor(env$Zone, levels = c("Coast", "Transition", "Offshore"))

Plot_shannon <- ggplot(env, aes(x = Zone, y = Shannon, fill= Depth)) +
  geom_boxplot(width=0.5, alpha=0.8, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.07), aes(color = Depth), size = 2) + # Adds points, aligned with boxes
  scale_fill_manual(values = c("turquoise4", "turquoise2")) +
  scale_color_manual(values = c("turquoise4", "turquoise3")) + 
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

#Calculate species richness

env$NumOfSp <- rowSums(community_NJ2022>0)
qqPlot(env$NumOfSp)

env$Zone <- factor(env$Zone, levels = c("Coast", "Transition", "Offshore"))

Plot_richness <- ggplot(env, aes(x = Zone, y = NumOfSp, fill= Depth)) +
  geom_boxplot(width=0.5, alpha=0.8, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.07), aes(color = Depth), size = 2) + # Adds points, aligned with boxes
  scale_fill_manual(values = c("turquoise4", "turquoise2")) +
  scale_color_manual(values = c("turquoise4", "turquoise3")) + 
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

png("Shannon_Richness_MicroDecon_turquoise.png",width=4000,height=2000,units="px",res=300,bg="white")
pdf("Shannon_Richness_MicroDecon_turquoise.pdf",width=24,height=11,bg="white")

ggarrange(Plot_shannon, Plot_richness,
          ncol = 2, nrow = 1,  align = "hv", 
          labels = c("A", "B"),
          widths = c(1, 1), heights = c(1, 1),
          common.legend = TRUE,
          # legend = "right",
          font.label = list(size=22))
dev.off()

#######Statistical analysis for alpha diversity

#Shannon
# model_shannon <- lm(Shannon_transformed ~ Depth*Zone, data=env)

model_shannon <- rlm(Shannon ~ Depth*Zone, data=env)  #less sensitive to outliers
anova_shannon <- Anova(model_shannon)
anova_shannon 
summary(model_shannon)
# 
# wald_test_results <- coeftest(model_shannon)
# print(wald_test_results)
# 
# library(robust)
# model_shannon <- lmRob(Shannon ~ Depth*Zone, data=env)  #less sensitive to outliers
# anova(model_shannon)

par(mfrow = c(2,2))
plot(model_shannon)

# # Run the Kruskal-Wallis test
# kruskal_result <- kruskal.test(Shannon ~ Depth*Zone, data = env)
# print(kruskal_result)

#Species richness
par(mfrow = c(1,1))
qqPlot(env$NumOfSp)

model_richness <- rlm(NumOfSp ~ Depth*Zone, data=env)  
model_richness <- lm(NumOfSp ~ Depth*Zone, data=env) 
anova_richness <- Anova(model_richness)
anova_richness 
summary(model_richness)

aov <- aov(NumOfSp ~ Zone, data=env)
aov
summary(aov)
TukeyHSD(aov, "Zone", ordered=FALSE, conf.level = 0.95) #sample size should be equal

#######Permanova

#calculate bray curtis
community_NJ2022_transformed_bray <- vegdist(community_NJ2022_transformed, method = "bray",binary=TRUE) 
community_NJ2022_bray <- vegdist(community_NJ2022, method = "bray",binary=TRUE) 
bray_matrix <- as.matrix(community_NJ2022_transformed_bray)

Env_p <- env$Zone
Depth_p <- env$Depth
Description_p <- env$Description

adonis2(community_NJ2022_transformed_bray~Depth*Zone*Depth/Zone, env,  permutations = 9999, sqrt.dist=FALSE, by="terms")
# adonis2(community_NJ2022_bray~Depth*Zone*Depth/Zone, env,  permutations = 9999, sqrt.dist=FALSE, by="terms")
pairwise.adonis2(community_NJ2022_transformed_bray~Zone, p.adjust.m="BH", reduce=NULL, env, perm=9999)
pairwise.adonis2(community_NJ2022_transformed_bray~Depth, p.adjust.m="BH", reduce=NULL, env, perm=9999)


par(mfrow=c(1,1))
#factor Combined ; doesn't accept models. Groups Can consist of a factor with a single level (i.e., one group).
Combined_betadisper_transformed <- betadisper(community_NJ2022_transformed_bray, Env_p, type="centroid")
Combined_betadisper_transformed <- betadisper(community_NJ2022_transformed_bray, Depth_p, type="centroid")

Combined_permdisp_transformed <- permutest(Combined_betadisper_transformed, permutations=9999)
Combined_permdisp_transformed #SIGN

#mar => A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. 
#The default is c(5, 4, 4, 2) + 0.1.
par(mar=c(5, 6, 4, 4)) #to make the plot area smaller so y-lab is visible

png("Betadisper_plot_MicroDecon.png",width=3000,height=3000,units="px",res=400,bg="white")
pdf("Betadisper_plot_MicroDecon.pdf",width=8,height=6,bg="white")
plot(Combined_betadisper_transformed, hull=FALSE, ellipse=TRUE, cex.axis=1.5, cex.lab=1.5, main="")
dev.off()

save.image((paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/R_Environment_Anova_Permanova_MicroDecon.RData")))

