if(.Platform$OS.type == "unix"){home<-"/home/"} else{
  home<-"//192.168.236.131/"}

libraries <- c("BiocManager"
               , "ggplot2"
               , "tidyr"
               , "taxonomizr"
               , "naturalsort"
               , "stringr"
               , "phylotools"
               , "vegan"
               , "scales"
               , "ggpattern"
               , "ggh4x"
               , "reshape2"
               , "grid"
               , "ggdendro"
               , "ggpubr"
)

for (opties in libraries){
  
  if (opties %in% installed.packages()){
    
    library(opties,character.only = TRUE)
    
  } else {install.packages(opties,repos = "http://cran.us.r-project.org")
    
    library(opties,character.only = TRUE)
  }
}

proj.path <- here("/home/genomics/ndukan/01_Zero-Impact/Temporal_Analysis")

env <- readRDS(paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/env_unrarefied_MicroDecon.rds"))
community_NJ2022 <- readxl::read_excel(paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Table_FishSpecies_MicroDecon.xlsx"))
community_NJ2022 <- as.data.frame(community_NJ2022)
rownames(community_NJ2022) <- community_NJ2022[,1]
community_NJ2022[,1] <- NULL
community_NJ2022 <- as.data.frame(t(community_NJ2022))

colorder <- c(env$Niskin.sample)
community_NJ2022 <-community_NJ2022[colorder,]

#transform the data
community_NJ2022_transformed <- decostand(community_NJ2022, method="total")
community_NJ2022_transformed <- decostand(community_NJ2022_transformed, method="max")

community_NJ2022_transformed <- as.data.frame(t(community_NJ2022_transformed))
community_NJ2022 <- as.data.frame(t(community_NJ2022))

Fish_dendro_nozero <- as.dendrogram(hclust(d=dist(x = community_NJ2022_transformed), method = "ward.D"))
# Fish_dendro_nozero <- as.dendrogram(hclust(d=dist(x = relative_abundance), method = "ward.D"))
dendro_plot_nozero <- ggdendrogram(data = Fish_dendro_nozero, rotate=T)
print(dendro_plot_nozero)

#Order Fish Species according to dendrogram
Fish_order_nozero <- order.dendrogram(Fish_dendro_nozero)
OrderFish_nozero <- community_NJ2022[Fish_order_nozero, ]

relative_abundance <- data.frame(matrix(nrow=nrow(OrderFish_nozero), ncol=ncol(OrderFish_nozero)
                                        , dimnames = list(c(rownames(OrderFish_nozero)), c(colnames(OrderFish_nozero)))))
for(i in 1:ncol(OrderFish_nozero)){relative_abundance[,i] <- (OrderFish_nozero[,i]/colSums(OrderFish_nozero[1:ncol(OrderFish_nozero)])[i])*100}

total_rel <- as.data.frame(rowSums(relative_abundance))
write.xlsx(total_rel, 
           paste0(proj.path,"/02_NJ2022_Miseq_5119441/NJ2022_Taxonomic_Assignment_Tables/Results_MicroDecon/Total_RelativeAbundances.xlsx"), 
           sheetName = "Sheet1", colNames = TRUE, rowNames = TRUE, append = FALSE)

#adding this step orders the species better in this dataset
Fish_dendro_nozero2 <- as.dendrogram(hclust(d=dist(x = relative_abundance), method = "ward.D"))
dendro_plot_nozero2 <- ggdendrogram(data = Fish_dendro_nozero2, rotate=T)
print(dendro_plot_nozero2)

Fish_order_nozero2 <- order.dendrogram(Fish_dendro_nozero2)
OrderFish_nozero2 <- OrderFish_nozero[Fish_order_nozero2, ]

relative_abundance <- tibble::rownames_to_column(relative_abundance, "Species") #put the row names in to column
colnames(relative_abundance) <- gsub("X","",colnames(relative_abundance)) #replace X with empty on the colnames
relative_abundance_location <- melt(relative_abundance, na.rm = FALSE, value.name = "value", id='Species') #makes the data set long format, 3 colulns one for species, one for samples, one for values. the length is 61(3species)x157(#samples)=9577
colnames(relative_abundance_location) <- c("Species", "Location", "Relative_Abundance")
relative_abundance_location[relative_abundance_location == 0] <- NA

Plot_Heat <-ggplot(data=relative_abundance_location
                   , aes(y=factor(Species, levels = rev(rownames(OrderFish_nozero2))),x=Location, fill= Relative_Abundance)
)+ 
  geom_tile(
  )+
  scale_fill_gradientn(colours=c( # For all the bad values <5%
    "darkgreen","green","yellow") # For the passed values
    , limits=c(0,100) # sets abundance range from 0 to 100 to get consitency when comparing samples
    # # , values=rescale(c(0,(5-(10^-10)) # For all the bad values <5%
    # # ,5,50,100)) # For the passed values
    , na.value="lightgrey"
    , breaks=c(25,50,75,100)
    , xlab("% \nRead\nAbundance\n")
  )+
  theme(legend.position = "right"
        , legend.title=element_text(size=18) #title text size
        , axis.text.y = element_text( face="italic",size=13, colour = "black") #species size
        # , plot.title = element_text(size=16) #size of the title of the plot
        , legend.text = element_text(size = 16)
        , strip.text.x = element_text(size = 11)
        , axis.text.x = element_text(size=12, colour = c(env$Description_color), angle=90, vjust=0.5)
        , plot.title = element_blank()
        , axis.title = element_blank()
  )+
  guides(pattern = guide_legend(override.aes = list(fill = "green"),
                                title="")
  )
Plot_Heat

pdf("HeatMap_MicroDecon_RelAbun_DoubleTrans_Final.pdf",width=20,height=10,bg="white")


ggarrange(Plot_Heat,
          ncol = 1, nrow = 1,  align = "hv", 
          # labels = c("A", "B"),
          widths = c(1, 1), heights = c(1, 1))

dev.off()

pdf("HeatMap_dendrogram_MicroDecon_final.pdf",width=12,height=8,bg="white")

dendro_plot_nozero2

dev.off()

save.image((paste0(proj.path,"/02_NJ2022_Miseq_5119441/R_Environment/R_Environment_MicroDecon/REnvironment_HeatMap_MicroDecon.RData")))
