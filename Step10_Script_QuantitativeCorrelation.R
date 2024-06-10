if(.Platform$OS.type == "unix"){home<-"/home/"} else{
  home<-"//192.168.236.131/"}

library(vegan)

source("/home/genomics/icornelis/03_RawScripts/funs-libs.R")
proj.path <- here("/home/genomics/ndukan/01_Zero-Impact/Temporal_Analysis")

community_abundance_NJ2022 <- readxl::read_excel(paste0(proj.path, "/02_NJ2022_Miseq_5119441/Morpholgical_abundance_matrix_STANDARDIZED.xlsx"))
community_abundance_NJ2022 <- as.data.frame(community_abundance_NJ2022)
rownames(community_abundance_NJ2022) <- community_abundance_NJ2022[,1]
community_abundance_NJ2022[,1] <- NULL

community_abundance_NJ2022_transformed <- decostand(community_abundance_NJ2022, method= "total")
community_abundance_NJ2022_transformed <- decostand(community_abundance_NJ2022_transformed, method= "total")

community_biomass_NJ2022 <- readxl::read_excel(paste0(proj.path, "/02_NJ2022_Miseq_5119441/Morpholgical_biomass_matrix_STANDARDIZED.xlsx"))
community_biomass_NJ2022 <- as.data.frame(community_biomass_NJ2022)
rownames(community_biomass_NJ2022) <- community_biomass_NJ2022[,1]
community_biomass_NJ2022[,1] <- NULL

community_biomass_NJ2022_transformed <- decostand(community_biomass_NJ2022, method= "total")
community_biomass_NJ2022_transformed <- decostand(community_biomass_NJ2022_transformed, method= "max")


community_eDNA_NJ2022 <- readxl::read_excel(paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Table_FishSpecies_MicroDecon.xlsx"))
community_eDNA_NJ2022 <- as.data.frame(community_eDNA_NJ2022)
rownames(community_eDNA_NJ2022) <- community_eDNA_NJ2022[,1]
community_eDNA_NJ2022[,1] <- NULL
community_eDNA_NJ2022 <- as.data.frame(t(community_eDNA_NJ2022))

community_eDNA_NJ2022_transformed <- decostand(community_eDNA_NJ2022, method= "total")
community_eDNA_NJ2022_transformed <- decostand(community_eDNA_NJ2022, method= "max")

community_eDNA_NJ2022_transformed <- as.data.frame(t(community_eDNA_NJ2022_transformed))
community_eDNA_NJ2022 <- as.data.frame(t(community_eDNA_NJ2022))

env <- readRDS(paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/R_Environment_MicroDecon/env_unrarefied_fish_MicroDecon.rds"))
colorder <- c(env$Niskin.sample)
community_eDNA_NJ2022_transformed <- community_eDNA_NJ2022_transformed[,colorder]

#number of stations that the species were detected
DetectionTrawl <- as.data.frame(colSums(community_abundance_NJ2022 !=0)) 
names(DetectionTrawl)<- "Location_Counts"

## Limanda limanda
#choose read count data for the species
read_asv <- filter(community_eDNA_NJ2022_transformed, rownames(community_eDNA_NJ2022_transformed) == "Limanda limanda")
read_asv <- data.frame(t(read_asv))

read_asv$Location <- env$Location

#avarage the read counts of replicates
mean_read_asv <- aggregate(Limanda.limanda ~ Location, data = read_asv, FUN = mean)
rownames(mean_read_asv) <- mean_read_asv[,1]
names(mean_read_asv)[2] <- "MeanRead_count"
mean_read_asv[1]<- NULL

#sum the read counts of replicates
sum_read_asv <- aggregate(Limanda.limanda ~ Location, data = read_asv, FUN = sum)
rownames(sum_read_asv) <- sum_read_asv[,1]
names(sum_read_asv)[2] <- "SumRead_count"
sum_read_asv[1] <- NULL

#choose abundance data for the species
abundance <- select(community_abundance_NJ2022, "Limanda limanda")
names(abundance)[1]<- "abundance"

#choose biomass data for the species
biomass <- select(community_biomass_NJ2022, "Limanda limanda")
names(biomass)[1]<- "biomass"

#merge the tables
merged_data <- merge(mean_read_asv, sum_read_asv, by = "row.names", all = TRUE)
rownames(merged_data) <-merged_data$Row.names
merged_data$Row.names <- NULL
merged_data <- merge(merged_data, abundance, by = "row.names", all = TRUE)
rownames(merged_data) <-merged_data$Row.names
merged_data$Row.names <- NULL
merged_data_l <- merge(merged_data, biomass, by = "row.names", all = TRUE)
rownames(merged_data_l) <-merged_data_l$Row.names
merged_data_l$Row.names <- NULL
merged_data_l$Species <- "Limanda limanda"

## Pleuronectes platessa
#choose read count data for the species
read_asv <- filter(community_eDNA_NJ2022_transformed, rownames(community_eDNA_NJ2022_transformed) == "Pleuronectes platessa")
read_asv <- data.frame(t(read_asv))

read_asv$Location <- env$Location

#avarage the read counts of replicates
mean_read_asv <- aggregate(Pleuronectes.platessa ~ Location, data = read_asv, FUN = mean)
rownames(mean_read_asv) <- mean_read_asv[,1]
names(mean_read_asv)[2] <- "MeanRead_count"
mean_read_asv[1]<- NULL

#sum the read counts of replicates
sum_read_asv <- aggregate(Pleuronectes.platessa ~ Location, data = read_asv, FUN = sum)
rownames(sum_read_asv) <- sum_read_asv[,1]
names(sum_read_asv)[2] <- "SumRead_count"
sum_read_asv[1] <- NULL

#choose abundance data for the species
abundance <- select(community_abundance_NJ2022, "Pleuronectes platessa")
names(abundance)[1]<- "abundance"

#choose biomass data for the species
biomass <- select(community_biomass_NJ2022, "Pleuronectes platessa")
names(biomass)[1]<- "biomass"

#merge the tables
merged_data <- merge(mean_read_asv, sum_read_asv, by = "row.names", all = TRUE)
rownames(merged_data) <-merged_data$Row.names
merged_data$Row.names <- NULL
merged_data <- merge(merged_data, abundance, by = "row.names", all = TRUE)
rownames(merged_data) <-merged_data$Row.names
merged_data$Row.names <- NULL
merged_data_p <- merge(merged_data, biomass, by = "row.names", all = TRUE)
rownames(merged_data_p) <-merged_data_p$Row.names
merged_data_p$Row.names <- NULL

merged_data_p$Species <- "Pleuronectes platessa"

## Solea solea
#choose read count data for the species
read_asv <- filter(community_eDNA_NJ2022_transformed, rownames(community_eDNA_NJ2022_transformed) == "Solea solea")
read_asv <- data.frame(t(read_asv))

read_asv$Location <- env$Location

#avarage the read counts of replicates
mean_read_asv <- aggregate(Solea.solea ~ Location, data = read_asv, FUN = mean)
rownames(mean_read_asv) <- mean_read_asv[,1]
names(mean_read_asv)[2] <- "MeanRead_count"
mean_read_asv[1]<- NULL

#sum the read counts of replicates
sum_read_asv <- aggregate(Solea.solea ~ Location, data = read_asv, FUN = sum)
rownames(sum_read_asv) <- sum_read_asv[,1]
names(sum_read_asv)[2] <- "SumRead_count"
sum_read_asv[1] <- NULL

#choose abundance data for the species
abundance <- select(community_abundance_NJ2022, "Solea solea")
names(abundance)[1]<- "abundance"

#choose biomass data for the species
biomass <- select(community_biomass_NJ2022, "Solea solea")
names(biomass)[1]<- "biomass"

#merge the tables
merged_data <- merge(mean_read_asv, sum_read_asv, by = "row.names", all = TRUE)
rownames(merged_data) <-merged_data$Row.names
merged_data$Row.names <- NULL
merged_data <- merge(merged_data, abundance, by = "row.names", all = TRUE)
rownames(merged_data) <-merged_data$Row.names
merged_data$Row.names <- NULL
merged_data_s <- merge(merged_data, biomass, by = "row.names", all = TRUE)
rownames(merged_data_s) <-merged_data_s$Row.names
merged_data_s$Row.names <- NULL

merged_data_s$Species <- "Solea solea"

## Echiichthys vipera
#choose read count data for the species
read_asv <- filter(community_eDNA_NJ2022_transformed, rownames(community_eDNA_NJ2022_transformed) == "Echiichthys vipera")
read_asv <- data.frame(t(read_asv))

read_asv$Location <- env$Location

#avarage the read counts of replicates
mean_read_asv <- aggregate(Echiichthys.vipera ~ Location, data = read_asv, FUN = mean)
rownames(mean_read_asv) <- mean_read_asv[,1]
names(mean_read_asv)[2] <- "MeanRead_count"
mean_read_asv[1]<- NULL

#sum the read counts of replicates
sum_read_asv <- aggregate(Echiichthys.vipera ~ Location, data = read_asv, FUN = sum)
rownames(sum_read_asv) <- sum_read_asv[,1]
names(sum_read_asv)[2] <- "SumRead_count"
sum_read_asv[1] <- NULL

#choose abundance data for the species
abundance <- select(community_abundance_NJ2022, "Echiichthys vipera")
names(abundance)[1]<- "abundance"

#choose biomass data for the species
biomass <- select(community_biomass_NJ2022, "Echiichthys vipera")
names(biomass)[1]<- "biomass"

#merge the tables
merged_data <- merge(mean_read_asv, sum_read_asv, by = "row.names", all = TRUE)
rownames(merged_data) <-merged_data$Row.names
merged_data$Row.names <- NULL
merged_data <- merge(merged_data, abundance, by = "row.names", all = TRUE)
rownames(merged_data) <-merged_data$Row.names
merged_data$Row.names <- NULL
merged_data_e <- merge(merged_data, biomass, by = "row.names", all = TRUE)
rownames(merged_data_e) <-merged_data_e$Row.names
merged_data_e$Row.names <- NULL

merged_data_e$Species <- "Echiichthys vipera"

## Trachurus trachurus
#choose read count data for the species
read_asv <- filter(community_eDNA_NJ2022_transformed, rownames(community_eDNA_NJ2022_transformed) == "Trachurus trachurus")
read_asv <- data.frame(t(read_asv))

read_asv$Location <- env$Location

#avarage the read counts of replicates
mean_read_asv <- aggregate(Trachurus.trachurus ~ Location, data = read_asv, FUN = mean)
rownames(mean_read_asv) <- mean_read_asv[,1]
names(mean_read_asv)[2] <- "MeanRead_count"
mean_read_asv[1]<- NULL

#sum the read counts of replicates
sum_read_asv <- aggregate(Trachurus.trachurus ~ Location, data = read_asv, FUN = sum)
rownames(sum_read_asv) <- sum_read_asv[,1]
names(sum_read_asv)[2] <- "SumRead_count"
sum_read_asv[1] <- NULL

#choose abundance data for the species
abundance <- select(community_abundance_NJ2022, "Trachurus trachurus")
names(abundance)[1]<- "abundance"

#choose biomass data for the species
biomass <- select(community_biomass_NJ2022, "Trachurus trachurus")
names(biomass)[1]<- "biomass"

#merge the tables
merged_data <- merge(mean_read_asv, sum_read_asv, by = "row.names", all = TRUE)
rownames(merged_data) <-merged_data$Row.names
merged_data$Row.names <- NULL
merged_data <- merge(merged_data, abundance, by = "row.names", all = TRUE)
rownames(merged_data) <-merged_data$Row.names
merged_data$Row.names <- NULL
merged_data_t <- merge(merged_data, biomass, by = "row.names", all = TRUE)
rownames(merged_data_t) <-merged_data_e$Row.names
merged_data_t$Row.names <- NULL

merged_data_t$Species <- "Trachurus trachurus"

## Merlangius merlangus
#choose read count data for the species
read_asv <- filter(community_eDNA_NJ2022_transformed, rownames(community_eDNA_NJ2022_transformed) == "Merlangius merlangus")
read_asv <- data.frame(t(read_asv))

read_asv$Location <- env$Location

#avarage the read counts of replicates
mean_read_asv <- aggregate(Merlangius.merlangus ~ Location, data = read_asv, FUN = mean)
rownames(mean_read_asv) <- mean_read_asv[,1]
names(mean_read_asv)[2] <- "MeanRead_count"
mean_read_asv[1]<- NULL

#sum the read counts of replicates
sum_read_asv <- aggregate(Merlangius.merlangus ~ Location, data = read_asv, FUN = sum)
rownames(sum_read_asv) <- sum_read_asv[,1]
names(sum_read_asv)[2] <- "SumRead_count"
sum_read_asv[1] <- NULL

#choose abundance data for the species
abundance <- select(community_abundance_NJ2022, "Merlangius merlangus")
names(abundance)[1] <- "abundance"

#choose biomass data for the species
biomass <- select(community_biomass_NJ2022, "Merlangius merlangus")
names(biomass)[1]<- "biomass"

#merge the tables
merged_data <- merge(mean_read_asv, sum_read_asv, by = "row.names", all = TRUE)
rownames(merged_data) <-merged_data$Row.names
merged_data$Row.names <- NULL
merged_data <- merge(merged_data, abundance, by = "row.names", all = TRUE)
rownames(merged_data) <-merged_data$Row.names
merged_data$Row.names <- NULL
merged_data_m <- merge(merged_data, biomass, by = "row.names", all = TRUE)
rownames(merged_data_m) <-merged_data_m$Row.names
merged_data_m$Row.names <- NULL

merged_data_m$Species <- "Merlangius merlangus"

merged_data_all <- rbind(merged_data_e, merged_data_l, merged_data_m, merged_data_p, merged_data_s, merged_data_t)

library(ggthemes)
library(ggsci)

# Create the correlation plot for biomass with linear trend + confidence interval
bio <- ggplot(merged_data_all, aes(x=biomass, y=MeanRead_count)) +
  geom_point(aes(color=Species),size=4, show.legend = FALSE) +
  # geom_smooth(method=lm , color="orangered2", fill="slategray2", se=TRUE) +
  geom_smooth(aes(color= Species, fill= Species),  method=lm , se=TRUE, show.legend = FALSE) +
  theme_minimal()+
  # ggtitle("") +
  theme_bw()+
  # scale_color_manual(values = ggthemes("excel_new", 6, type = "continuous")) +  # Set color palette
  # scale_fill_manual(values = ggthemes("excel_new", 6, type = "continuous")) +  # Set fill palette
  scale_color_manual(values = pal_nejm(palette = "default", alpha = 1)(6)) +  # Set color palette with alpha
  scale_fill_manual(values = pal_nejm(palette = "default", alpha = 1)(6)) +  # Set fill palette with alpha
  theme(plot.title = element_text(hjust = 0.5, size = 32))+
  xlab("Biomass") +
  ylab("eDNA read count") +
  theme(axis.title.x = element_text(size = 28, margin = margin(t = 15)),
        axis.title.y = element_text(size = 28, margin = margin(r = 15)),
        axis.text = element_text(size= 20),
        strip.text = element_text(size = 24))+
  facet_wrap(~Species, ncol=3, nrow=2 )+
  stat_cor( label.sep="\n", method="spearman", cor.coef.name= "rho", size=10)
# geom_text(aes(label = paste("rho =", round(cor(Relative_biomass, Relative_MeanRead, method = "spearman"), 2))),
#           x = Inf, y = -Inf, hjust = 1, vjust = 0, size = 12)  # doesn't make sense
bio

# png("Correlation_RelativeBiomass_NJ2022_new.png",width=3000,height=3000,units="px",res=400,bg="white")
pdf("Correlation_Biomass_NJ2022_mean_MicroDecon_relabunBiomass.pdf",width=22,height=14,bg="white")

ggarrange(bio,
          ncol = 1, nrow = 1,  align = "hv", 
          # labels = c("A", "B"),
          widths = c(1, 1), heights = c(1, 1))

dev.off()

# Create the correlation plot for abundance with linear trend + confidence interval
abun <- ggplot(merged_data_all, aes(x=abundance, y=MeanRead_count)) +
  geom_point(aes(color=Species),size=4, show.legend = FALSE) +
  # geom_smooth(method=lm , color="orangered2", fill="slategray2", se=TRUE) +
  geom_smooth(aes(color= Species, fill= Species),  method=lm , se=TRUE, show.legend = FALSE) +
  theme_minimal()+
  # ggtitle("") +
  theme_bw()+
  # scale_color_manual(values = ggthemes("excel_new", 6, type = "continuous")) +  # Set color palette
  # scale_fill_manual(values = ggthemes("excel_new", 6, type = "continuous")) +  # Set fill palette
  scale_color_manual(values = pal_nejm(palette = "default", alpha = 1)(6)) +  # Set color palette with alpha
  scale_fill_manual(values = pal_nejm(palette = "default", alpha = 1)(6)) +  # Set fill palette with alpha
  theme(plot.title = element_text(hjust = 0.5, size = 32))+
  xlab("Abundance") +
  ylab("eDNA read count") +
  theme(axis.title.x = element_text(size = 28, margin = margin(t = 15)),
        axis.title.y = element_text(size = 28, margin = margin(r = 15)),
        axis.text = element_text(size= 20),
        strip.text = element_text(size = 24))+
  facet_wrap(~Species, ncol=3, nrow=2 )+
  stat_cor( label.sep="\n", method="spearman", cor.coef.name= "rho", size=10)
# geom_text(aes(label = paste("rho =", round(cor(Relative_biomass, Relative_MeanRead, method = "spearman"), 2))),
#           x = Inf, y = -Inf, hjust = 1, vjust = 0, size = 12)  # doesn't make sense
abun

# png("Correlation_RelativeAbundance_NJ2022_MicroDecon.png",width=3000,height=3000,units="px",res=400,bg="white")
pdf("Correlation_Abundance_NJ2022_mean_MicroDecon_relabunBiomass.pdf",width=22,height=14,bg="white")

ggarrange(abun,
          ncol = 1, nrow = 1,  align = "hv", 
          # labels = c("A", "B"),
          widths = c(1, 1), heights = c(1, 1))

dev.off()

save.image((paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/R_Environment_MicroDecon/R_Environment_QuantitativeCorrelation_MicroDecon.RData")))

remove(merged_data_all, merged_data, merged_data_e, merged_data_l, merged_data_m, merged_data_p, merged_data_s, merged_data_t)