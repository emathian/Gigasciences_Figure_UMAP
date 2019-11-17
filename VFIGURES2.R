# Librairies --------------------------------------------------------------
library(umap)
library(ggplot2)
library(caTools)
library(ade4)
library(umap)
library(dplyr)
library(data.table)
library(rspatial)
library(gridExtra)
library(ggpubr)
library(viridis)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(cluster)
library(eulerr)
library(ade4)
source('DR_Method/Dimensionality_reduction_comparison_metrics.R')
library(gganimate)
library(latex2exp)
library(gifski) # for gganimate
library(png) # for gganimate
set.seed(1564404882)
theme_set(theme_bw())

# Dimensionality reduction projection -------------------------------------

## Import data -------------------------------------------------------------

Attributes_UMAP_TCACLCNECSCL <- read.table("data/Attributes_UMAP_TCACLCNECSCLC.tsv",sep='\t' , header = T)
vst50_TCACLCNECSCLC<- read.table("data/VST_nosex_TCACLCNECSCLC.txt", header = T)
vst50_TCACLCNECSCLC <- data.frame(t(vst50_TCACLCNECSCLC ))
vst50_TCACLCNECSCLC_designRD <- vst50_TCACLCNECSCLC # Structure match with the one require for PCA and UMAP
vst50_TCACLCNECSCLC_designRD <- vst50_TCACLCNECSCLC_designRD[- which(rownames(vst50_TCACLCNECSCLC_designRD) == "S00716_A"|rownames(vst50_TCACLCNECSCLC_designRD) == "S02322_B"),]
rownames(vst50_TCACLCNECSCLC_designRD)[which(rownames(vst50_TCACLCNECSCLC_designRD) == "S02322_A")] <- "S02322.R1"

## UMAP Projections --------------------------------------------------------

### NN = 121 

#### Legends
umap121 = umap(vst50_TCACLCNECSCLC_designRD, n_neighbors = 121)
umap121_res_df <- as.data.frame(umap121$layout)
umap121_res_df  = setDT(umap121_res_df , keep.rownames = TRUE)[]
colnames(umap121_res_df)[1] <- "Sample_ID"
umap121_res_df <- umap121_res_df[order(umap121_res_df$Sample_ID),]
p2_standard_nn121 <- ggplot(umap121_res_df, aes(x=V1, y=V2,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha =.8) + scale_color_brewer(palette="Spectral") +labs(title = "121")+ theme(legend.title = element_blank(), legend.text = element_text(size = 16))
p2_standard_nn121_legend <- get_legend(p2_standard_nn121)
p2_standard_legend <-get_legend(p2_standard_nn121)

p2Kmeans_standard_nn121 <- ggplot(umap121_res_df, aes(x=V1, y=V2,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha =.8) + scale_color_brewer(palette="Spectral") +
  theme(legend.text = element_text(size = 11) ,
        legend.title = element_blank()) +  guides(col = guide_legend( ncol = 2)) # Y axis text
p2Kmeans_standard_nn121_legend <- get_legend(p2Kmeans_standard_nn121)

p2_standard_nn121 <- ggplot(umap121_res_df, aes(x=V1, y=V2,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha =.8) + scale_color_brewer(palette="Spectral")
p2_standard_nn121 <- p2_standard_nn121 +  labs(title="", y="dim2", x="dim1")  +
  theme(legend.position = "none", legend.text = element_text(size = 11) ,legend.title = element_blank())# +  guides(col = guide_legend( ncol = 4)) # Y axis text

p2_standard_nn121

### NN 15

umap15 = umap(vst50_TCACLCNECSCLC_designRD)
umap15_res_df <- as.data.frame(umap15$layout)
umap15_res_df  = setDT(umap15_res_df , keep.rownames = TRUE)[]
colnames(umap15_res_df)[1] <- "Sample_ID"


### NN 112

umap112 = umap(vst50_TCACLCNECSCLC_designRD, n_neighbors = 112)
umap112_res_df <- as.data.frame(umap112$layout)
umap112_res_df  = setDT(umap112_res_df , keep.rownames = TRUE)[]
colnames(umap112_res_df)[1] <- "Sample_ID"


### NN 135

umap135 = umap(vst50_TCACLCNECSCLC_designRD, n_neighbors = 135  )
umap135_res_df <- as.data.frame(umap135$layout)
umap135_res_df  = setDT(umap135_res_df , keep.rownames = TRUE)[]
colnames(umap135_res_df)[1] <- "Sample_ID"

### NN 150

umap150 = umap(vst50_TCACLCNECSCLC_designRD)
umap150_res_df <- as.data.frame(umap150$layout)
umap150_res_df  = setDT(umap150_res_df , keep.rownames = TRUE)[]
colnames(umap150_res_df)[1] <- "Sample_ID"

### NN 208

#umap208 = umap(vst50_TCACLCNECSCLC_designRD, n_neighbors = 208)
#umap208_res_df <- as.data.frame(umap208$layout)
#umap208_res_df  = setDT(umap208_res_df , keep.rownames = TRUE)[]
#colnames(umap208_res_df)[1] <- "Sample_ID"
umap208_res_df <- read.table("data/Coords_umap_nn208.tsv", header = F)
colnames(umap208_res_df) <- c("Sample_ID", "V1", "V2")
umap208_res_df <- umap208_res_df[order(umap208_res_df$Sample_ID),]

p2_standard_nn208 <- ggplot(umap208_res_df, aes(x=V1, y=V2,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha =.8) + scale_color_brewer(palette="Spectral") +labs(title= "208")+ theme(legend.title = element_blank()) + theme(legend.position = "none")
p2_standard_nn208 <- p2_standard_nn208 +  labs(title="", y="dim2", x="dim1")  +
  theme(legend.position = "none",
        plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
        plot.subtitle =element_text(size=14, hjust=0.5),
        plot.caption =element_text(size=12,  hjust=0.5),
        axis.title.x=element_text(size=16),  # X axis title
        axis.title.y=element_text(size=16),  # Y axis title
        axis.text.x=element_text(size=14),  # X axis text
        axis.text.y=element_text(size=14),
        legend.text = element_text(size = 11) ,
        legend.title = element_blank())# +  guides(col = guide_legend( ncol = 4)) # Y axis text

## PCA ---------------------------------------------------------------------

acp_2D <- dudi.pca(vst50_TCACLCNECSCLC_designRD, center = T , scale = F , scannf = F , nf=2) #
acp_fig1_li_df =  as.data.frame(acp_2D$li)
acp_fig1_li_df = setDT(acp_fig1_li_df, keep.rownames = TRUE)[]
colnames(acp_fig1_li_df)[1]<-'Sample_ID'
acp_fig1_li_df <- acp_fig1_li_df[order(acp_fig1_li_df$Sample_ID),]

acp_5D <- dudi.pca(vst50_TCACLCNECSCLC_designRD, center = T , scale = F , scannf = F , nf=5) #
acp_5D_li_df =  as.data.frame(acp_5D$li)
acp_5D_li_df  = setDT(acp_5D_li_df , keep.rownames = TRUE)[]
colnames(acp_5D_li_df )[1]<-'Sample_ID'
acp_5D_li_df <- acp_5D_li_df[order(acp_5D_li_df$Sample_ID),]
#rownames(vst50_TCACLCNECSCLC_designRD) == acp_fig1_li_df$Sample_ID

p1 <- ggplot(acp_5D_li_df, aes(x=Axis1 , y=Axis2,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha=0.8)+
  scale_color_brewer(palette="Spectral")  + scale_shape_manual(values = c(16, 7, 15,17,8 ))#wes_palette(n=5, name="FantasticFox1")
p1 <- p1 +  labs(title="",   y="PC2 (9.7%)", x="PC1 (27%)") +
  theme( legend.position="none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=16),  # X axis title
         axis.title.y=element_text(size=16),  # Y axis title
         axis.text.x=element_text(size=14),  # X axis text
         axis.text.y=element_text(size=14),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text # Y axis text
p1


## Projections -------------------------------------------------------------


p5D1d2d <- ggplot(acp_5D_li_df, aes(x=Axis1 , y=Axis2,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha=0.8)+
  scale_color_brewer(palette="Spectral")  + scale_shape_manual(values = c(16, 7, 15,17,8 ))#wes_palette(n=5, name="FantasticFox1")
p5D1d2d <- p5D1d2d +  labs(title="",   y="PC2 (9.7%)", x="PC1 (26.6%)") +
  theme( legend.position="none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=12),  # X axis title
         axis.title.y=element_text(size=12),  # Y axis title
         axis.text.x=element_text(size=12),  # X axis text
         axis.text.y=element_text(size=12),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text # Y axis text

p5D1d3d <- ggplot(acp_5D_li_df, aes(x=Axis1 , y=Axis3,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha=0.8)+
  scale_color_brewer(palette="Spectral")  + scale_shape_manual(values = c(16, 7, 15,17,8 ))#wes_palette(n=5, name="FantasticFox1")
p5D1d3d <- p5D1d3d +  labs(title="",   y="PC3 (4.4%)", x="PC1 (26.6%)") +
  theme( legend.position="none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=12),  # X axis title
         axis.title.y=element_text(size=12),  # Y axis title
         axis.text.x=element_text(size=12),  # X axis text
         axis.text.y=element_text(size=12),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text # Y axis text

p5D1d4d <- ggplot(acp_5D_li_df, aes(x=Axis1 , y=Axis4,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha=0.8)+
  scale_color_brewer(palette="Spectral")  + scale_shape_manual(values = c(16, 7, 15,17,8 ))#wes_palette(n=5, name="FantasticFox1")
p5D1d4d  <- p5D1d4d  +  labs(title="",   y="PC4 (3.6%)", x="PC1 (26.6%)") +
  theme( legend.position="none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=12),  # X axis title
         axis.title.y=element_text(size=12),  # Y axis title
         axis.text.x=element_text(size=12),  # X axis text
         axis.text.y=element_text(size=12),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text # Y axis text

p5D1d5d <- ggplot(acp_5D_li_df, aes(x=Axis1 , y=Axis5,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha=0.8)+
  scale_color_brewer(palette="Spectral")  + scale_shape_manual(values = c(16, 7, 15,17,8 ))#wes_palette(n=5, name="FantasticFox1")
p5D1d5d  <- p5D1d5d  +  labs(title="",   y="PC5 (3.3%)", x="PC1 (26.6%)") +
  theme( legend.position="none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=12),  # X axis title
         axis.title.y=element_text(size=12),  # Y axis title
         axis.text.x=element_text(size=12),  # X axis text
         axis.text.y=element_text(size=12),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text # Y axis text

p5D2d3d <- ggplot(acp_5D_li_df, aes(x=Axis2 , y=Axis3,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha=0.8)+
  scale_color_brewer(palette="Spectral")  + scale_shape_manual(values = c(16, 7, 15,17,8 ))#wes_palette(n=5, name="FantasticFox1")
p5D2d3d  <- p5D2d3d  +  labs(title="",   y="PC3 (4.4%)", x="PC2 (9.7%)") +
  theme( legend.position="none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=12),  # X axis title
         axis.title.y=element_text(size=12),  # Y axis title
         axis.text.x=element_text(size=12),  # X axis text
         axis.text.y=element_text(size=12),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text # Y axis text


p5D2d4d <- ggplot(acp_5D_li_df, aes(x=Axis2 , y=Axis4,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha=0.8)+
  scale_color_brewer(palette="Spectral")  + scale_shape_manual(values = c(16, 7, 15,17,8 ))#wes_palette(n=5, name="FantasticFox1")
p5D2d4d  <- p5D2d4d  +  labs(title="",   y="PC4 (3.6%)", x="PC2 (9.7%)") +
  theme( legend.position="none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=12),  # X axis title
         axis.title.y=element_text(size=12),  # Y axis title
         axis.text.x=element_text(size=12),  # X axis text
         axis.text.y=element_text(size=12),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text # Y axis text


p5D2d5d <- ggplot(acp_5D_li_df, aes(x=Axis2 , y=Axis5,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha=0.8)+
  scale_color_brewer(palette="Spectral")  + scale_shape_manual(values = c(16, 7, 15,17,8 ))#wes_palette(n=5, name="FantasticFox1")
p5D2d5d  <- p5D2d5d  +  labs(title="",   y="PC5 (3.3%)", x="PC2 (9.7%)") +
  theme( legend.position="none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=12),  # X axis title
         axis.title.y=element_text(size=12),  # Y axis title
         axis.text.x=element_text(size=12),  # X axis text
         axis.text.y=element_text(size=12),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text # Y axis text


p5D3d4d <- ggplot(acp_5D_li_df, aes(x=Axis3 , y=Axis4,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha=0.8)+
  scale_color_brewer(palette="Spectral")  + scale_shape_manual(values = c(16, 7, 15,17,8 ))#wes_palette(n=5, name="FantasticFox1")
p5D3d4d  <- p5D3d4d  +  labs(title="",   y="PC4 (3.6%)", x="PC3 (4.4%)") +
  theme( legend.position="none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=12),  # X axis title
         axis.title.y=element_text(size=12),  # Y axis title
         axis.text.x=element_text(size=12),  # X axis text
         axis.text.y=element_text(size=12),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text # Y axis text


p5D3d5d <- ggplot(acp_5D_li_df, aes(x=Axis3 , y=Axis5,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha=0.8)+
  scale_color_brewer(palette="Spectral")  + scale_shape_manual(values = c(16, 7, 15,17,8 ))#wes_palette(n=5, name="FantasticFox1")
p5D3d5d  <- p5D3d5d  +  labs(title="",   y="PC5 (3.3%)", x="PC3 (4.4%)") +
  theme( legend.position="none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=12),  # X axis title
         axis.title.y=element_text(size=12),  # Y axis title
         axis.text.x=element_text(size=12),  # X axis text
         axis.text.y=element_text(size=12),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text # Y axis text


p5D4d5d <- ggplot(acp_5D_li_df, aes(x=Axis4 , y=Axis5,  color=Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)) +  geom_point(size=4, alpha=0.8)+
  scale_color_brewer(palette="Spectral")  + scale_shape_manual(values = c(16, 7, 15,17,8 ))#wes_palette(n=5, name="FantasticFox1")
p5D4d5d  <- p5D4d5d  +  labs(title="",   y="PC5 (3.3%)", x="PC4 (3.6%)") +
  theme( legend.position="none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=12),  # X axis title
         axis.title.y=element_text(size=12),  # Y axis title
         axis.text.x=element_text(size=12),  # X axis text
         axis.text.y=element_text(size=12),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text # Y axis text


ggarrange(ggarrange(p2_standard_nn208,p2_standard_nn121_legend , ncol = 2, widths = c(0.8,.2)),  ggarrange(p5D1d2d, p5D1d3d, p5D1d4d, p5D1d5d, p5D2d3d, p5D2d4d, p5D2d5d, p5D3d4d, p5D3d5d, p5D4d5d, ncol = 2, nrow = 5) , nrow = 2, heights = c(0.3,0.7) , labels = c("A", "B"))


# Kmeans clustering-------------------------------------------------------
## Kmeans Clustering on UMAP -----------------------------------------------


# We remove the small group, that is to say the cluster the supra-carcinoids and the LCNEC/NA
umap208_res_df_V2 <- umap208_res_df
umap208_res_df_V2 <- cbind(umap208_res_df_V2, "Molecular_clusters" = Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)
umap208_res_df_V2 <- umap208_res_df_V2[-which(umap208_res_df_V2$Molecular_clusters == "LCNEC/NA" |  umap208_res_df_V2$Molecular_clusters == "Supra_carcinoid"),]

# Kmeans
umap_nn208_C6 <- kmeans(umap208_res_df_V2 [,2:3],centers=6,nstart=30)
umap_nn208_C6$cluster


### Kmeans Projections ------------------------------------------------------

pKmeans <- ggplot(umap208_res_df_V2, aes(x=V1, y=V2,  color= as.factor(umap_nn208_C6$cluster ))) +  geom_point(size=4, alpha =.8) + scale_color_brewer(palette="Spectral")
pKmeans <- pKmeans +  labs(title="", y="", x="")  +
  theme( legend.position = "bottom",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=16),  # X axis title
         axis.title.y=element_text(size=16),  # Y axis title
         axis.text.x=element_text(size=14),  # X axis text
         axis.text.y=element_text(size=14),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank()) +  guides(col = guide_legend( ncol = 3)) # Y axis text

pKmeans_legend <- get_legend(pKmeans)

pKmeans <- ggplot(umap208_res_df_V2, aes(x=V1, y=V2,  color= as.factor(umap_nn208_C6$cluster ))) +  geom_point(size=4, alpha =.8) + scale_color_brewer(palette="Spectral")
pKmeans <- pKmeans +  labs(title="", y="", x="")  +
  theme( legend.position = "none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=16),  # X axis title
         axis.title.y=element_text(size=16),  # Y axis title
         axis.text.x=element_text(size=14),  # X axis text
         axis.text.y=element_text(size=14),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank()) +  guides(col = guide_legend( ncol = 3)) # Y axis text


### Kmeans confusion matrix -------------------------------------------------

umap208_res_df_V2$Molecular_clusters <- as.character(umap208_res_df_V2$Molecular_clusters)
Molecular_clusters_6g <- unlist(lapply( 1:length( umap208_res_df_V2$Molecular_clusters), function(i){
  
  if (umap208_res_df_V2$Molecular_clusters[i] == "SCLC/LCNEC-like"){
    "LCNEC/TypeI"
  }
  else if (umap208_res_df_V2$Molecular_clusters[i] == "LCNEC/SCLC-like"){
    "SCLC/SCLC-like"
  }
  else{
    umap208_res_df_V2$Molecular_clusters[i]
  }
  
}))

x = table(umap_nn208_C6$cluster,Molecular_clusters_6g ); x
x = prop.table(x, 2)
x2 = as.matrix(x)

melted_prop <- melt(x2)
head(melted_prop)

kmean_HM <- ggplot(data = melted_prop, aes(x=Molecular_clusters_6g , y=as.factor(Var1), fill=value))  +geom_tile()+    scale_fill_gradient(low = "#FFFFCC", high = "#41B6C4")+
  geom_text(aes(label=round(value,2)), color="black",  face = "bold",size=4)
kmean_HM <- kmean_HM +  labs(title="", x="Molecular Clusters", y="Kmeans Cluster")  +
  theme( legend.position = "none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=16),  # X axis title
         axis.title.y=element_text(size=16),  # Y axis title
         axis.text.x=element_text(size=14, angle=45,  hjust=1),  # X axis text
         axis.text.y=element_text(size=14),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text


kmean_HM


### Silhouette conefficient --------------------------------------------------

silhouette(umap_nn208_C6$cluster, dist(umap208_res_df_V2[,2:3]))
summary(silhouette(umap_nn208_C6$cluster, dist(umap208_res_df_V2[,2:3])), FUN = mean)


### Fisher test -------------------------------------------------------------

# Carcinoids A1 cluster 4
fisher.test( matrix(c(33,0,3,166), ncol = 2 , nrow = 2))

# Carcinoids A2 cluster 6
fisher.test( matrix(c(32,0,0,170), ncol = 2 , nrow = 2))

# Carcinoids B cluster 1
fisher.test( matrix(c(20,0,0,182), ncol = 2 , nrow = 2))


# LCNEC T1 Cluster 2
fisher.test( matrix(c(28,9,11,154), ncol = 2 , nrow = 2))


# LCNEC T2 Cluster 3
fisher.test( matrix(c(20,12,3,167), ncol = 2 , nrow = 2))


# SCLC Cluster 5
fisher.test( matrix(c(45,1,8,148), ncol = 2 , nrow = 2), simulate.p.value = TRUE)


## Kmeans on PCA ------------------------------------------------------------
### Kmeans Projections ------------------------------------------------------

sample_umap <- data.frame("Sample_ID" =  umap208_res_df_V2$Sample_ID, "Molecular_clusters" = umap208_res_df_V2$Molecular_clusters)
acp_5D_li_df_vkmeans <- merge(sample_umap, acp_5D_li_df, by = "Sample_ID")


acp_5D_li_df_vkmeans_C6<- kmeans(acp_5D_li_df_vkmeans[,3:dim(acp_5D_li_df_vkmeans)[2]],centers=6,nstart=30)
acp_5D_li_df_vkmeans_C6$cluster

pKmeans_ACPL <- ggplot(acp_5D_li_df_vkmeans, aes(x=Axis1 , y=Axis2,  color= as.factor(acp_5D_li_df_vkmeans_C6$cluster ))) +  geom_point(size=4, alpha =.8) + scale_color_brewer(palette="Spectral")
#pKmeans_legend <- get_legend(pKmeans_ACPL)

pKmeans_ACP <- ggplot(acp_5D_li_df_vkmeans, aes(x=Axis1 , y=Axis2,  color= as.factor(acp_5D_li_df_vkmeans_C6$cluster ))) +  geom_point(size=4, alpha =.8) + scale_color_brewer(palette="Spectral")
pKmeans_ACP <- pKmeans_ACP +  labs(title="", y="", x="")  +
  theme( legend.position = "none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=16),  # X axis title
         axis.title.y=element_text(size=16),  # Y axis title
         axis.text.x=element_text(size=14),  # X axis text
         axis.text.y=element_text(size=14),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text


### Kmeans confusion matrix -------------------------------------------------
acp_5D_li_df_vkmeans$Molecular_clusters <- as.character(acp_5D_li_df_vkmeans$Molecular_clusters)
Molecular_clusters_6g_ACP <- unlist(lapply( 1:length(acp_5D_li_df_vkmeans$Molecular_clusters), function(i){
  
  if (acp_5D_li_df_vkmeans$Molecular_clusters[i] == "SCLC/LCNEC-like"){
    "LCNEC/TypeI"
  }
  else if (acp_5D_li_df_vkmeans$Molecular_clusters[i] == "LCNEC/SCLC-like"){
    "SCLC/SCLC-like"
  }
  else{
    acp_5D_li_df_vkmeans$Molecular_clusters[i]
  }
  
}))


x_acp = table(acp_5D_li_df_vkmeans_C6$cluster,Molecular_clusters_6g_ACP); x_acp
x_acp = prop.table(x_acp, 2)
x2_acp = as.matrix(x_acp)

melted_prop_acp <- melt(x2_acp)
head(melted_prop_acp)

### Silhouette conefficient --------------------------------------------------
silhouette(acp_5D_li_df_vkmeans_C6$cluster, dist(acp_5D_li_df_vkmeans[,3:dim(acp_5D_li_df_vkmeans)[2]]))
summary(silhouette(acp_5D_li_df_vkmeans_C6$cluster, dist(acp_5D_li_df_vkmeans[,3:dim(acp_5D_li_df_vkmeans)[2]])), FUN = mean)

kmean_HM_ACP <- ggplot(data = melted_prop_acp, aes(x=Molecular_clusters_6g_ACP , y=as.factor(Var1), fill=value))  +geom_tile()+    scale_fill_gradient(low = "#FFFFCC", high = "#41B6C4")+
  geom_text(aes(label=round(value,2)), color="black",  face = "bold",size=4)
kmean_HM_ACP <- kmean_HM_ACP+  labs(title="", x="Molecular Clusters", y="Kmeans Cluster")  +
  theme( legend.position = "none",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=16),  # X axis title
         axis.title.y=element_text(size=16),  # Y axis title
         axis.text.x=element_text(size=14, angle=45,  hjust=1),  # X axis text
         axis.text.y=element_text(size=14),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text


kmean_HM_ACP

### Fisher test -------------------------------------------------------------

# Carcinoids A1 cluster 4
fisher.test( matrix(c(33,0,2,167), ncol = 2 , nrow = 2))

# Carcinoids A2 cluster 6
fisher.test( matrix(c(26,6,16,154), ncol = 2 , nrow = 2))
chisq.test(matrix(c(26,6,16,154), ncol = 2 , nrow = 2))

# Carcinoids B cluster 1
fisher.test( matrix(c(10,10,26,156), ncol = 2 , nrow = 2))


# LCNEC T1 Cluster 2
fisher.test( matrix(c(23,16,20,153), ncol = 2 , nrow = 2))


# LCNEC T2 Cluster 3
fisher.test( matrix(c(18,7,25,152), ncol = 2 , nrow = 2))


# SCLC Cluster 5
fisher.test( matrix(c(40,0,13,149), ncol = 2 , nrow = 2), simulate.p.value = TRUE)

fisher.test( matrix(c(20,15,23,149), ncol = 2 , nrow = 2), simulate.p.value = TRUE)


ggarrange(          ggarrange(p2_standard_nn208, pKmeans,  p1, pKmeans_ACP,  p2Kmeans_standard_nn121_legend, pKmeans_legend, ncol = 2, nrow = 3,  heights = c(2,2,1), labels = c("A","B","D", "E")),
                    ggarrange( kmean_HM ,  kmean_HM_ACP,  nrow =2,  labels = c("C","F")),
                    ncol = 2, widths = c(2,1)
)


######################################################################################################@
# Sequence difference View
# Sequence difference view ----------------------------------------------
## Function for SD layout --------------------------------------------------
Seq_graph_by_k_vfig  <- function (data_Seq, Names=NULL, list_col=NULL, data_diff_mean_K = NULL, log = FALSE, my_legend_position, norm = T ){
  if (is.null(data_diff_mean_K) == TRUE) {
    data_diff_mean_k <- data.frame("k" =  unique(data_Seq$K))
    for (j in seq(from = 3, to = dim(data_Seq)[2], by = 1)) {
      mean_by_k <- tapply(data_Seq[, j], data_Seq$K , mean)
      data_diff_mean_k <- cbind(data_diff_mean_k, mean_by_k)
    }
    colnames(data_diff_mean_k)[2:length(colnames(data_diff_mean_k))] <- colnames(data_Seq)[3:dim(data_Seq)[2]]
    
    if (is.null(Names) == FALSE){
      if (length(Names) != (dim(data_Seq)[2] - 3)){
        warning("The list of names gave as input doesn't match with the number of curve.")
      }
      else{
        colnames(data_diff_mean_k)[2:length(colnames(data_diff_mean_k))] <- Names
      }
    }
    
  }
  else{
    data_diff_mean_k <- data_diff_mean_K
    
  }
  # Normalisation
  if (norm == T){
    data_diff_mean_kRef  <- data_diff_mean_k
    data_diff_mean_kRef[,2] <-  (data_diff_mean_k[,3] -  data_diff_mean_k[,2])/( data_diff_mean_k[,3]- data_diff_mean_k[,2])
    data_diff_mean_kRef[,3] <-  ( data_diff_mean_k[,3]- data_diff_mean_k[,3] )/( data_diff_mean_k[,3]- data_diff_mean_k[,2])
    data_diff_mean_kRef[,4] <-  ( data_diff_mean_k[,3] - data_diff_mean_k[,4] )/( data_diff_mean_k[,3]- data_diff_mean_k[,2])
    data_diff_mean_kRef[,5] <-  (data_diff_mean_k[,3] -  data_diff_mean_k[,5])/( data_diff_mean_k[,3]- data_diff_mean_k[,2])
    data_diff_mean_kRef[,6] <-  (data_diff_mean_k[,3] -  data_diff_mean_k[,6])/( data_diff_mean_k[,3]- data_diff_mean_k[,2])
    print( data_diff_mean_kRef)
    myylab <- TeX("Normalized $\\bar{SD}_k$")
  }
  else{
    data_diff_mean_kRef <- data_diff_mean_k
    myylab <- TeX("$\\bar{SD}_k$")
    
  }
  if (log == FALSE){
    data_diff_mean_kRef_graph <- data.frame('k' = data_diff_mean_kRef$k , 'diff_seq' =( data_diff_mean_kRef[, 2]), 'Method' = rep(as.character(colnames(data_diff_mean_kRef)[2]), length(data_diff_mean_kRef$k)))
    if (dim(data_diff_mean_kRef)[2]>3){
      for (i in 3:(dim(data_diff_mean_kRef)[2])){
        c_df <- data.frame('k' = data_diff_mean_kRef$k , 'diff_seq' = (data_diff_mean_kRef[, i]), 'Method' = rep(as.character(colnames(data_diff_mean_kRef)[i]), length(data_diff_mean_kRef$k)))
        data_diff_mean_kRef_graph <- rbind(data_diff_mean_kRef_graph, c_df)
      }
    }
    
    theme_set(theme_bw())
    p <- ggplot(data_diff_mean_kRef_graph, aes(x=k, y=diff_seq,  color=Method)) + geom_line() + geom_point()+
      scale_color_viridis(discrete=TRUE) 
    p <- p +  labs(title="", 
                   y=myylab, x= "k") +
      theme( legend.position= my_legend_position ,#
             plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
             plot.subtitle =element_text(size=14, hjust=0.5),
             plot.caption =element_text(size=12,  hjust=0.5),
             axis.title.x=element_text(size=16, face = "italic"),  # X axis title
             axis.title.y=element_text(size=16),  # Y axis title
             axis.text.x=element_text(size=14),  # X axis text
             axis.text.y=element_text(size=14),
             legend.text = element_text(size = 14) ,
             legend.title = element_blank())#+  guides(col = guide_legend( ncol = 2))  # Y axis text
    print(p)
    return(p)
  }
  else {
    data_diff_mean_kRef_graph <- data.frame('k' = data_diff_mean_kRef$k , 'diff_seq' = log( data_diff_mean_kRef[, 2],10), 'Method' = rep(as.character(colnames(data_diff_mean_kRef)[2]), length(data_diff_mean_kRef$k)))
    if (dim(data_diff_mean_kRef)[2]>3){
      for (i in 3:(dim(data_diff_mean_kRef)[2])){
        c_df <- data.frame('k' = data_diff_mean_kRef$k , 'diff_seq' = log(data_diff_mean_kRef[, i],10), 'Method' = rep(as.character(colnames(data_diff_mean_kRef)[i]), length(data_diff_mean_kRef$k)))
        data_diff_mean_kRef_graph <- rbind(data_diff_mean_kRef_graph, c_df)
      }
    }
    theme_set(theme_bw())
    p <- ggplot(data_diff_mean_kRef_graph, aes(x=k, y=diff_seq,  color=Method)) + geom_line() + geom_point()+
      scale_color_viridis(discrete=TRUE) 
    p <- p +  labs(title="", 
                   y=TeX("$log(\\bar{SD}_k)$"), x="k") +
      theme( legend.position= my_legend_position ,#
             plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
             plot.subtitle =element_text(size=14, hjust=0.5),
             plot.caption =element_text(size=12,  hjust=0.5),
             axis.title.x=element_text(size=16, face = "italic"),  # X axis title
             axis.title.y=element_text(size=16),  # Y axis title
             axis.text.x=element_text(size=14),  # X axis text
             axis.text.y=element_text(size=14),
             legend.text = element_text(size = 14) ,
             legend.title = element_blank())#+  guides(col = guide_legend( ncol = 2))  # Y axis text
    print(p)
    return(p)
  } 
}



## Pre processing ----------------------------------------------------------

acp_fig1_li_df$Sample_ID <- as.character(acp_fig1_li_df$Sample_ID)
acp_fig1_li_df <- acp_fig1_li_df[order(acp_fig1_li_df$Sample_ID),]
head(acp_fig1_li_df)

acp_5D_li_df$Sample_ID <- as.character(acp_5D_li_df$Sample_ID)
acp_5D_li_df <- acp_5D_li_df[order(acp_5D_li_df$Sample_ID),]
head(acp_5D_li_df)

umap15_res_df$Sample_ID <- as.character(umap15_res_df$Sample_ID)
umap15_res_df <- umap15_res_df[order(umap15_res_df$Sample_ID),]
head(umap15_res_df)


umap112_res_df$Sample_ID <- as.character(umap112_res_df$Sample_ID)
umap112_res_df <- umap112_res_df[order(umap112_res_df$Sample_ID),]
head(umap112_res_df)

umap135_res_df$Sample_ID <- as.character(umap135_res_df$Sample_ID)
umap135_res_df <- umap135_res_df[order(umap135_res_df$Sample_ID),]
head(umap135_res_df)

umap121_res_df$Sample_ID <- as.character(umap121_res_df$Sample_ID)
umap121_res_df <- umap121_res_df[order(umap121_res_df$Sample_ID),]
head(umap121_res_df)

umap208_res_df$Sample_ID <- as.character(umap208_res_df$Sample_ID)
umap208_res_df <- umap208_res_df[order(umap208_res_df$Sample_ID),]
head(umap208_res_df)

vst50_TCACLCNECSCLC_designRDSAMPLE <- rownames(vst50_TCACLCNECSCLC_designRD)
vst50_TCACLCNECSCLC_VSamp <- cbind("Sample_ID" = vst50_TCACLCNECSCLC_designRDSAMPLE, vst50_TCACLCNECSCLC_designRD)
vst50_TCACLCNECSCLC_VSamp <- vst50_TCACLCNECSCLC_VSamp[order(vst50_TCACLCNECSCLC_VSamp$Sample_ID),]
vst50_TCACLCNECSCLC_VSamp[1:10,1:10]
List_projection <- list(data.frame(acp_fig1_li_df), data.frame(acp_5D_li_df),data.frame(umap15_res_df[,1:3]),data.frame(umap112_res_df[,1:3]) ,data.frame(umap121_res_df[,1:3]), data.frame(umap135_res_df[,1:3]),data.frame(umap208_res_df[,1:3]))
str(List_projection)

List_projection2 <- list(data.frame(acp_fig1_li_df), data.frame(umap15_res_df[,1:3]),data.frame(umap208_res_df[,1:3]))

## SD calculation ---------------------------------------------------------
gloabal_seq_list = Seq_calcul( List_projection2,  data.frame(vst50_TCACLCNECSCLC_VSamp), listK =seq(from= 10, to = 208, by = 15) )
#Main_SQ_res_NN <- Seq_main(l_data = List_projection , dataRef = data.frame(vst50_TCACLCNECSCLC_VSamp) , listK = seq(from= 1, to = 208, by = 5), colnames_res_df = c("pca_2D", "pca_5D","umap_md=0.1_nn=15", "umap_md=0.1_nn=112","umap_md=0.1_nn=121", "umap_md=0.1_nn=135","umap_md=0.1_nn=208")  , filename = NULL , graphics = TRUE, stats = TRUE) #
#saveRDS(Main_SQ_res_NN, "Main_SQ_res_NN.rds")
Main_SQ_res_NN <- readRDS("Main_SQ_res_NN.rds")
head(Main_SQ_res_NN[[1]])
colnames(Main_SQ_res_NN[[1]]) <- c("Sample_ID", "K", "PCA 2D", "PCA 4D", "UMAP nn = 15", "UMAP nn = 121", "UMAP nn = 208")
