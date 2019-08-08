# V Git 
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
library(gifski) # for gganimate
library(png) # for gganimate
set.seed(1564404882)
theme_set(theme_bw())

# Dimensionality reduction projection -------------------------------------

## Import data -------------------------------------------------------------

Attributes_UMAP_TCACLCNECSCL <- read.table("/data/Attributes_UMAP_TCACLCNECSCLC.tsv",sep='\t' , header = T)
vst50_TCACLCNECSCLC<- read.table("/data/VST_nosex_50pc_TCACLCNECSCLC.txt", header = T)
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

### NN 208
umap208_res_df <- read.table("data/Coords_umap_nn208.tsv")
colnames(umap208_res_df) <- c("Sample_ID", "V1", "V2") 
#umap208 = umap(vst50_TCACLCNECSCLC_designRD, n_neighbors = 208)
#umap208_res_df <- as.data.frame(umap208$layout)
#umap208_res_df  = setDT(umap208_res_df , keep.rownames = TRUE)[]
#colnames(umap208_res_df)[1] <- "Sample_ID"
#umap208_res_df <- umap208_res_df[order(umap208_res_df$Sample_ID),]

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
List_projection <- list(data.frame(acp_fig1_li_df), data.frame(acp_5D_li_df),data.frame(umap15_res_df[,1:3]), data.frame(umap121_res_df[,1:3]), data.frame(umap208_res_df[,1:3]))
str(List_projection)


## SD calculation ---------------------------------------------------------

#Main_SQ_res_NN <- Seq_main(l_data = List_projection , dataRef = data.frame(vst50_TCACLCNECSCLC_VSamp) , listK = seq(from= 1, to = 208, by = 5), colnames_res_df = c("pca_2D", "pca_5D","umap_md=0.1_nn=15", "umap_md=0.1_nn=121", "umap_md=0.1_nn=208")  , filename = NULL , graphics = TRUE, stats = TRUE) #
#saveRDS(Main_SQ_res_NN, "Main_SQ_res_NN.rds")
Main_SQ_res_NN <- readRDS("/RData/Main_SQ_res_NN.rds")
head(Main_SQ_res_NN[[1]])
colnames(Main_SQ_res_NN[[1]]) <- c("Sample_ID", "K", "PCA 2D", "PCA 4D", "UMAP nn = 15", "UMAP nn = 121", "UMAP nn = 208")

## SD main graphics --------------------------------------------------------

sd_main <- Seq_graph_by_k_vfig(data_Seq = Main_SQ_res_NN[[1]], Names=NULL, list_col=NULL, data_diff_mean_K = NULL, log =FALSE, "none", norm = FALSE)
sd_main
sd_main_norm <- Seq_graph_by_k_vfig(data_Seq = Main_SQ_res_NN[[1]], Names=NULL, list_col=NULL, data_diff_mean_K = NULL, log =FALSE, "none", norm = TRUE)

sd_main_normL <- Seq_graph_by_k_vfig(data_Seq = Main_SQ_res_NN[[1]], Names=NULL, list_col=NULL, data_diff_mean_K = NULL, log =FALSE, "left", norm = TRUE)
sd_main_norm_legend <- get_legend(sd_main_normL )

head(Main_SQ_res_NN[[1]])
colnames(Main_SQ_res_NN[[1]]) <- c("Sample_ID", "K", "PCA 2D", "PCA 4D", "UMAP nn = 15", "UMAP nn = 121", "UMAP nn = 208")
sd_main <- Seq_graph_by_k_vfig(data_Seq = Main_SQ_res_NN[[1]], Names=NULL, list_col=NULL, data_diff_mean_K = NULL, log =FALSE, "none", norm = FALSE)
sd_main
sd_main_norm <- Seq_graph_by_k_vfig(data_Seq = Main_SQ_res_NN[[1]], Names=NULL, list_col=NULL, data_diff_mean_K = NULL, log =FALSE, "none", norm = TRUE)
length(intersect(SdMap_Umap$Sample_ID[1:30] , SdMap_Umap_refACP$Sample_ID[1:30]))

## SD Main statistics ----------------------------------------------------

pwt_df <- data.frame('mean_seq' = Main_SQ_res_NN[[2]][, 2], 'method'= rep(paste(colnames(Main_SQ_res_NN[[2]])[2], 2, sep = ""), dim(Main_SQ_res_NN[[2]])[1]))
for (i in 3:dim(Main_SQ_res_NN[[2]])[2]){
  c_df <- data.frame('mean_seq' = Main_SQ_res_NN[[2]][, i], 'method'=rep(paste(colnames(Main_SQ_res_NN[[2]])[i], i, sep = ""), dim(Main_SQ_res_NN[[2]])[1]))
  pwt_df <- rbind(pwt_df, c_df )
}
paired_test_m_BILAT  <- pairwise.t.test(pwt_df$mean_seq, p.adj = "holm", pwt_df$method,    paired = TRUE, alternative = 'two.sided')$p.value ;       paired_test_m
paired_test_m_BILAT  <- as.matrix(paired_test_m_BILAT)
colnames(paired_test_m_BILAT ) <- c("pca 2D","pca 5D", "umap nn = 15","umap nn = 121")
rownames(paired_test_m_BILAT ) <- c("pca 5D", "umap nn = 15","umap nn = 121", "umap nn = 208" )
paired_test_m_BILAT  <-  format(paired_test_m_BILAT  , scientific=TRUE,digits=3)
tablemain_BILAT <- ggtexttable(paired_test_m_BILAT , rows =rownames(paired_test_m_BILAT),theme = ttheme("classic"))
tablemain_BILAT


pwt_df <- data.frame('mean_seq' = Main_SQ_res_NN[[2]][, 2], 'method'= rep(paste(colnames(Main_SQ_res_NN[[2]])[2], 2, sep = ""), dim(Main_SQ_res_NN[[2]])[1]))
for (i in 3:dim(Main_SQ_res_NN[[2]])[2]){
  c_df <- data.frame('mean_seq' = Main_SQ_res_NN[[2]][, i], 'method'=rep(paste(colnames(Main_SQ_res_NN[[2]])[i], i, sep = ""), dim(Main_SQ_res_NN[[2]])[1]))
  pwt_df <- rbind(pwt_df, c_df )
}
paired_test_m_GREATER  <- pairwise.t.test(pwt_df$mean_seq, p.adj = "holm", pwt_df$method,    paired = TRUE, alternative = 'greater')$p.value ;     
paired_test_m_GREATER  <- as.matrix(paired_test_m_GREATER)
colnames(paired_test_m_GREATER ) <- c("pca 2D","pca 5D", "umap nn = 15","umap nn = 121")
rownames(paired_test_m_GREATER ) <- c("pca 5D", "umap nn = 15","umap nn = 121", "umap nn = 208" )
paired_test_m_GREATER  <-  format(paired_test_m_GREATER  , scientific=TRUE,digits=3)
tablemain_GREATER <- ggtexttable(paired_test_m_GREATER , rows =rownames(paired_test_m_GREATER),theme = ttheme("classic"))
tablemain_GREATER



pwt_df <- data.frame('mean_seq' = Main_SQ_res_NN[[2]][, 2], 'method'= rep(paste(colnames(Main_SQ_res_NN[[2]])[2], 2, sep = ""), dim(Main_SQ_res_NN[[2]])[1]))
for (i in 3:dim(Main_SQ_res_NN[[2]])[2]){
  c_df <- data.frame('mean_seq' = Main_SQ_res_NN[[2]][, i], 'method'=rep(paste(colnames(Main_SQ_res_NN[[2]])[i], i, sep = ""), dim(Main_SQ_res_NN[[2]])[1]))
  pwt_df <- rbind(pwt_df, c_df )
}
paired_test_m_LESS  <- pairwise.t.test(pwt_df$mean_seq, p.adj = "holm", pwt_df$method,    paired = TRUE, alternative = 'less')$p.value ;     
paired_test_m_LESS  <- as.matrix(paired_test_m_LESS)
colnames(paired_test_m_LESS ) <- c("pca 2D","pca 5D", "umap nn = 15","umap nn = 121")
rownames(paired_test_m_LESS ) <- c("pca 5D", "umap nn = 15","umap nn = 121", "umap nn = 208" )
paired_test_m_LESS  <-  format(paired_test_m_LESS  , scientific=TRUE,digits=3)
tablemain_LESS <- ggtexttable(paired_test_m_LESS , rows =rownames(paired_test_m_LESS),theme = ttheme("classic"))
tablemain_LESS

Main_stat_table <- matrix(ncol =5, nrow = 5)
colnames(Main_stat_table) = c("PCA 2D", "PCA 5D", "UMAP nn = 15", "UMAP nn = 121", "UMAP nn = 208")
rownames(Main_stat_table) = c("PCA 2D", "PCA 5D", "UMAP nn = 15", "UMAP nn = 121", "UMAP nn = 208")
Main_stat_table[ 1,] = c("-", "1.00e+00", "1.76e-04", "1.06e-02", "2.84e-01")
Main_stat_table[ 2,] = c("5.99e-08", "-", "1.42e-05", "1.16e-05", "4.60e-06")
Main_stat_table[ 3,] = c("1.00e+00", "1.00e+00", "-", "1.00e+00", "1.00e+00")
Main_stat_table[ 4,] = c("1.00e+00", "1.00e+00", "3.18e-05", "-", "1.00e+00")
Main_stat_table[ 5,] = c("1.00e+00", "1.00e+00", "4.09e-05", "2.88e-0.4", "-")
Main_stat_table
ggarrange(tablemain_BILAT, tablemain_GREATER, tablemain_LESS, nrow = 3)

tab <- ggtexttable(Main_stat_table, 
                   theme = ttheme(
                     colnames.style = colnames_style(face = "bold"),
                     rownames.style = rownames_style(face = "bold"))
                   )
tab
tab <- table_cell_bg(tab, row = 2, column = 2,
                       fill="darkgray")
tab <- table_cell_bg(tab, row = 3, column = 3,
                     fill="darkgray")
tab <- table_cell_bg(tab, row = 4, column = 4,
                     fill="darkgray")
tab <- table_cell_bg(tab, row = 5, column = 5,
                     fill="darkgray")
tab <- table_cell_bg(tab, row = 6, column = 6,
                     fill="darkgray")


tab

# Main Graphic   #

ggarrange(ggarrange(sd_main, sd_main_norm, sd_main_norm_legend, ncol = 3, widths = c(0.45,0.4,0.2), labels = c("A","B","")) ,ggarrange(tab, labels = c("C")) , nrow= 2)


# SD local analysis -------------------------------------------------------
## SD clculation (Local) -----------------------------------------------------------

#Main_SQ_res_NN_20 <- Seq_main(l_data = List_projection , dataRef = data.frame(vst50_TCACLCNECSCLC_VSamp) , listK = seq(from= 1, to = 20, by = 1), colnames_res_df = c("pca_2D", "pca_5D","umap_md=0.1_nn=15", "umap_md=0.1_nn=121", "umap_md=0.1_nn=208")  , filename = NULL , graphics = TRUE, stats = TRUE) #
#saveRDS(Main_SQ_res_NN_20, "Main_SQ_res_NN_20.rds")
Main_SQ_res_NN_20 <- readRDS("/RData/Main_SQ_res_NN_20.rds")

Main_SQ_res_NN_20v2 <- Main_SQ_res_NN_20[[1]][,1:2]
Main_SQ_res_NN_20v2 <- cbind(Main_SQ_res_NN_20v2,  Main_SQ_res_NN_20[[1]][,4:dim( Main_SQ_res_NN_20[[1]])[2]] )
seq_local <- Seq_graph_by_k_vfig(data_Seq = Main_SQ_res_NN_20[[1]], Names=NULL, list_col=NULL, data_diff_mean_K = NULL, log = FALSE, c(0.2,0.8), norm = F)
seq_local



## Statistic tests (Local)---------------------------------------------------------


pwt_df <- data.frame('mean_seq' = Main_SQ_res_NN_20[[2]][, 2], 'method'= rep(paste(colnames(Main_SQ_res_NN_20[[2]])[2], 2, sep = ""), dim(Main_SQ_res_NN_20[[2]])[1]))
for (i in 3:dim(Main_SQ_res_NN_20[[2]])[2]){
   c_df <- data.frame('mean_seq' = Main_SQ_res_NN_20[[2]][, i], 'method'=rep(paste(colnames(Main_SQ_res_NN_20[[2]])[i], i, sep = ""), dim(Main_SQ_res_NN_20[[2]])[1]))
   pwt_df <- rbind(pwt_df, c_df )
}
paired_test_m_BILAT  <- pairwise.wilcox.test(pwt_df$mean_seq, p.adj = "holm", pwt_df$method,    paired = TRUE, alternative = 'two.sided')$p.value ;    
paired_test_m_BILAT  <- as.matrix(paired_test_m_BILAT)
colnames(paired_test_m_BILAT ) <- c("pca 2D","pca 5D", "umap nn = 15","umap nn = 121")
rownames(paired_test_m_BILAT ) <- c("pca 5D", "umap nn = 15","umap nn = 121", "umap nn = 208" )
paired_test_m_BILAT  <-  format(paired_test_m_BILAT  , scientific=TRUE,digits=3)
tableloc_BILAT <- ggtexttable(paired_test_m_BILAT , rows =rownames(paired_test_m_BILAT),theme = ttheme("classic"))
tableloc_BILAT


paired_test_m_GREATER  <- pairwise.wilcox.test(pwt_df$mean_seq, p.adj = "holm", pwt_df$method,    paired = TRUE, alternative = 'greater')$p.value ;       
paired_test_m_GREATER  <- as.matrix(paired_test_m_GREATER)
colnames(paired_test_m_GREATER ) <- c("pca 2D","pca 5D", "umap nn = 15","umap nn = 121")
rownames(paired_test_m_GREATER ) <- c("pca 5D", "umap nn = 15","umap nn = 121", "umap nn = 208" )
paired_test_m_GREATER  <-  format(paired_test_m_GREATER  , scientific=TRUE,digits=3)
tableloc_GREATER <- ggtexttable(paired_test_m_GREATER , rows =rownames(paired_test_m_GREATER),theme = ttheme("classic"))
tableloc_GREATER



paired_test_m_LESS  <- pairwise.wilcox.test(pwt_df$mean_seq, p.adj = "holm", pwt_df$method,    paired = TRUE, alternative = 'less')$p.value ;       
paired_test_m_LESS  <- as.matrix(paired_test_m_LESS)
colnames(paired_test_m_LESS ) <- c("pca 2D","pca 5D", "umap nn = 15","umap nn = 121")
rownames(paired_test_m_LESS ) <- c("pca 5D", "umap nn = 15","umap nn = 121", "umap nn = 208" )
paired_test_m_LESS  <-  format(paired_test_m_LESS  , scientific=TRUE,digits=3)
tableloc_LESS <- ggtexttable(paired_test_m_LESS , rows =rownames(paired_test_m_LESS),theme = ttheme("classic"))
tableloc_LESS


# Overview
ggarrange(tableloc_BILAT, tableloc_GREATER, tableloc_LESS, nrow = 3)



## SD Maps --------------------------------------------------------------

Sd_map_df <- data.frame("Sample_ID" = Main_SQ_res_NN$Seq_df$Sample_ID, "k" = Main_SQ_res_NN$Seq_df$K , "UMAP208" = Main_SQ_res_NN$Seq_df$`UMAP nn = 208` )
sd_map_lnen <- SD_map_f(Sd_map_df , umap208_res_df, "bottom")
sd_map_lnen[[1]]
sd_map_lnen[[2]]


## Individual analysis of the neighborhood ---------------------------------

### Projection of S02297 nearest neighbors according UMAP Coordinates -----------------------

dist1 <- as.matrix(dist( umap208_res_df[, 2:dim(umap208_res_df)[2]], method = "euclidian", diag = TRUE, upper = TRUE))
umap208_Sample <- as.character(umap208_res_df$Sample_ID)
rownames(dist1) <- umap208_Sample
colnames(dist1) <- umap208_Sample
dist1 <- as.data.frame(dist1)
dist1S02297 <-  dist1[order(dist1$S02297),]
rownames(dist1S02297)
dist1S02297 <- data.frame("Sample_ID" = rownames(dist1S02297), "rank" = seq(1:length(rownames(dist1S02297))))

UmapS02297 <- merge(umap208_res_df, dist1S02297,  by= "Sample_ID")

pS02297L <- ggplot(UmapS02297, aes(x=V1, y=V2,  color=rank, label =UmapS02297$Sample_ID)) +  geom_point(size=4, alpha =.8) + scale_color_distiller(palette = "Spectral")
pS02297L <- pS02297L + geom_point(aes(x = UmapS02297$V1[which(UmapS02297$Sample_ID == "S02297")] , y =  UmapS02297$V2[which(UmapS02297$Sample_ID == "S02297")] ), size=4, col = "black")   +
  geom_text(aes( UmapS02297$V1[ UmapS02297$Sample_ID== "S02297"],UmapS02297$V2[ UmapS02297$Sample_ID== "S02297"], label = UmapS02297$Sample_ID[ UmapS02297$Sample_ID== "S02297"]) , col = 'black' , hjust=0, vjust=0 )+ 
    labs(title="", 
          y=TeX("dim2"), x="dim1") +
  theme( legend.position = "bottom",
    plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
    plot.subtitle =element_text(size=14, hjust=0.5),
    plot.caption =element_text(size=12,  hjust=0.5),
    axis.title.x=element_text(size=16),  # X axis title
    axis.title.y=element_text(size=16),  # Y axis title
    axis.text.x=element_text(size=14),  # X axis text
    axis.text.y=element_text(size=14),
    legend.text = element_text(size = 10) ,
    legend.title = element_blank())#+  guides(col = guide_legend( ncol = 2))  # Y axis text
pS02297L


### Projection of S02297nearest neighbors  according the original space -----------------------


### Calculation of SD with the PCA as reference -----------------------------

List_projection <- list(data.frame(umap208_res_df[,1:3]))
str(List_projection)

#Main_SQ_res_NN_refACP <- Seq_main(l_data = List_projection , dataRef = data.frame(acp_5D_li_df) , listK = seq(from= 1, to = 208, by = 5), colnames_res_df = c( "umap_md=0.1_nn=208")  , filename = NULL , graphics = TRUE, stats = TRUE) #
#saveRDS(Main_SQ_res_NN_refACP, "Main_SQ_res_NN_refACP.rds")

Main_SQ_res_NN_refACP <- readRDS("/RData/Main_SQ_res_NN_refACP.rds")

Sd_map_df_refACP <- data.frame("Sample_ID" = Main_SQ_res_NN_refACP$Seq_df$Sample_ID, "k" = Main_SQ_res_NN_refACP$Seq_df$K , "UMAP208" = Main_SQ_res_NN_refACP$Seq_df$`umap_md=0.1_nn=208` )
sd_map_lnen_pca <- SD_map_f(Sd_map_df_refACP , umap208_res_df, "bottom")
sd_map_lnen_pca[[1]]
sd_map_lnen_pca[[2]]

### Projection of S02297nearest neighbors  according PCA 5D --------

dist2 <- as.matrix(dist(acp_5D_li_df[, 2:dim(acp_5D_li_df)[2]], method = "euclidian", diag = TRUE, upper = TRUE))
acp_5D_li_df_Sample <- as.character(acp_5D_li_df$Sample_ID)
rownames(dist2) <- acp_5D_li_df_Sample 
colnames(dist2) <- acp_5D_li_df_Sample 
dist2 <- as.data.frame(dist2)
dist2S02297ACP <-  dist2[order(dist2$S02297),]
rownames(dist2S02297ACP)
dist2S02297ACP <- data.frame("Sample_ID" = rownames(dist2S02297ACP), "rank" = seq(1:length(rownames(dist2S02297ACP))))

UmapS02297ACPH <- merge(umap208_res_df, dist2S02297ACP,  by= "Sample_ID")

pS02297ACPH <- ggplot(UmapS02297ACPH, aes(x=V1, y=V2,  color=rank, label =UmapS02297ACPH$Sample_ID)) +  geom_point(size=4, alpha =.8) + scale_color_distiller(palette = "Spectral")
pS02297ACPH <- pS02297ACPH + geom_point(aes(x = UmapS02297ACPH$V1[which(UmapS02297ACPH$Sample_ID == "S02297")] , y =  UmapS02297ACPH$V2[which(UmapS02297ACPH$Sample_ID == "S02297")] ), size=4, col = "black")   +
  geom_text(aes( UmapS02297ACPH$V1[ UmapS02297ACPH$Sample_ID== "S02297"],UmapS02297ACPH$V2[ UmapS02297ACPH$Sample_ID== "S02297"], label = UmapS02297ACPH$Sample_ID[ UmapS02297ACPH$Sample_ID== "S02297"]) , col = 'black' , hjust=0, vjust=0 )+
  labs(title="", 
       y=TeX("dim2"), x="dim1") +
  theme( legend.position = "bottom",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=16),  # X axis title
         axis.title.y=element_text(size=16),  # Y axis title
         axis.text.x=element_text(size=14),  # X axis text
         axis.text.y=element_text(size=14),
         legend.text = element_text(size = 10) ,
         legend.title = element_blank())#+  guides(col = guide_legend( ncol = 2))  # Y axis text



# MAp main graphic #
ggarrange(pSD_Main_refACP,pS02297L, pS02297ACPH,  labels = c("A", "B", "C"), nrow = 1)



# Centraliy preservation --------------------------------------------------

List_projection <- list(data.frame(acp_fig1_li_df), data.frame(acp_5D_li_df),data.frame(umap15_res_df[,1:3]), data.frame(umap121_res_df[,1:3]), data.frame(umap208_res_df[,1:3]))
str(List_projection)
length(List_projection)
#Main_CP_res_NN <- CP_main(l_data = List_projection , list_K = c(seq(from= 1, to = 208, by = 5),208) , dataRef =  data.frame(vst50_TCACLCNECSCLC_VSamp) , colnames_res_df = c("pca_2D", "pca_5D","umap_md=0.1_nn=15", "umap_md=0.1_nn=121", "umap_md=0.1_nn=208") , filename = NULL , graphics = TRUE, stats = TRUE)
#saveRDS(Main_CP_res_NN , "Main_CP_res_NN.rds")
Main_CP_res_NN <- readRDS('/RData/Main_CP_res_NN.rds')
Main_CP_res_NN$Graphic


CP_graph_by_k_VFIG  <-function (data_CP,  ref_CP_data, Names=NULL, list_col=NULL, log=FALSE, my_legend_position = "bottom"){
  if (dim(data_CP)[1] != dim(ref_CP_data)[1]){
    warning("Input data frames don't have the same number of line.")
    break
  }
  else{
    colnames(data_CP)[1] <- "Sample_ID" ; colnames(ref_CP_data)[1] <- "Sample_ID"
    colnames(data_CP)[2] <- "K" ; colnames(ref_CP_data)[2] <- "K"
    data_CP <- cbind(data_CP, ref_CP_data[, 3])
    data_CP[,3:dim(data_CP)[2]] <- data_CP[,3:dim(data_CP)[2]]
    data_diff_mean_k <- data.frame("k" =  unique(data_CP$K))
    if (is.null(Names)==TRUE){
      Names <- colnames(data_CP)[3:(length(colnames(data_CP))-1)]
    }
    else if (length(Names) != dim(data_CP)[2]-3){
      warning("The length of the list of Names is inconsistant in regard of the length of the input data.")
      Names <- colnames(data_CP)[3:(length(colnames(data_CP))-1)]
    }
    for (I in seq(from = 3, to = (dim(data_CP)[2]-1), by = 1)) {
      abs_diff <- abs(scale(as.numeric(data_CP[, I])) - scale(as.numeric(data_CP[, dim(data_CP)[2]])))
      c_abs_diff_df <- data.frame("abs_diff" = abs_diff, 'K' = data_CP$K)
      abs_diff_k <- tapply(c_abs_diff_df$abs_diff, c_abs_diff_df$K, mean)
      data_diff_mean_k <- cbind(data_diff_mean_k, abs_diff_k)
    }
    colnames(data_diff_mean_k)[2:length(colnames(data_diff_mean_k))] <- Names
  }
  data_diff_mean_k_graph <- data.frame('k' = data_diff_mean_k$k , 'diff_cp' = data_diff_mean_k[, 2], 'Method' = rep(as.character(colnames(data_diff_mean_k)[2]), length(data_diff_mean_k$k)))
  if (dim(data_diff_mean_k)[2] >= 3){
    for (i in 3:(dim(data_diff_mean_k)[2])){
      c_df <- data.frame('k' = data_diff_mean_k$k , 'diff_cp' = data_diff_mean_k[, i], 'Method' = rep(as.character(colnames(data_diff_mean_k)[i]), length(data_diff_mean_k$k)))
      data_diff_mean_k_graph <- rbind(data_diff_mean_k_graph, c_df)
    }
  }
  print("head(data_diff_mean_k_graph )")
  print(head(data_diff_mean_k_graph ))
  print(summary(as.factor(data_diff_mean_k_graph$Method)))
  theme_set(theme_bw())
  if (log==FALSE){
    p <- ggplot(data_diff_mean_k_graph, aes(x=k, y=diff_cp,  color=Method)) + geom_line() + geom_point()+
      scale_color_viridis(discrete=TRUE) 
    p <- p +  labs(title="", caption = "",
                   y=TeX("$CP_k$"), x="K") +theme( legend.position= my_legend_position ,#
                                                   plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
                                                   plot.subtitle =element_text(size=14, hjust=0.5),
                                                   plot.caption =element_text(size=12,  hjust=0.5),
                                                   axis.title.x=element_text(size=16),  # X axis title
                                                   axis.title.y=element_text(size=16),  # Y axis title
                                                   axis.text.x=element_text(size=14),  # X axis text
                                                   axis.text.y=element_text(size=14),
                                                   legend.text = element_text(size = 14) ,
                                                   legend.title = element_blank())
  }
  else {
    p <- ggplot(data_diff_mean_k_graph, aes(x=k, y=log(diff_cp),  color=Method)) + geom_line() + geom_point()+
      scale_color_viridis(discrete=TRUE) 
    p <- p +  labs(title="", caption = "",
                   y=TeX("$log(CP_k)$"), x="K") +theme( legend.position= my_legend_position ,#
                                                        plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
                                                        plot.subtitle =element_text(size=14, hjust=0.5),
                                                        plot.caption =element_text(size=12,  hjust=0.5),
                                                        axis.title.x=element_text(size=16),  # X axis title
                                                        axis.title.y=element_text(size=16),  # Y axis title
                                                        axis.text.x=element_text(size=14),  # X axis text
                                                        axis.text.y=element_text(size=14),
                                                        legend.text = element_text(size = 14) ,
                                                        legend.title = element_blank())
      }
  return(p)
}

CP_df_res = Main_CP_res_NN[[1]][, 1:(dim(Main_CP_res_NN[[1]])[2]-1 )]
CP_ref_df_res  = Main_CP_res_NN[[1]][,1:2]
CP_ref_df_res = cbind(CP_ref_df_res, Main_CP_res_NN[[1]]$REF)
CP_Main_graph <- CP_graph_by_k_VFIG(CP_df_res,  CP_ref_df_res, Names=NULL, list_col=NULL, log=FALSE, my_legend_position = c(.2,.9))
CP_Main_graph


# CP statistics -----------------------------------------------------------

pwt_df_208 <- data.frame('mean_cp' = Main_CP_res_NN[[2]][, 2], 'method'= rep(paste(colnames(Main_CP_res_NN[[2]])[2], 2, sep = ""), dim(Main_CP_res_NN[[2]])[1]))
for (i in 3:dim(Main_CP_res_NN[[2]])[2]){
  c_df <- data.frame('mean_cp' = Main_CP_res_NN[[2]][, i], 'method'=rep(paste(colnames(Main_CP_res_NN[[2]])[i], i, sep = ""), dim(Main_CP_res_NN[[2]])[1]))
  pwt_df_208 <- rbind(pwt_df_208, c_df )
}

paired_Ttest_mg_BILAT <-pairwise.t.test(pwt_df_208$mean_cp, pwt_df_208$method,  p.adj = "holm",  paired = TRUE , alternative = "two.sided")$p.value 
paired_Ttest_m_matg_BILAT <- as.matrix(paired_Ttest_mg_BILAT)
colnames(paired_Ttest_m_matg_BILAT) <- c("pca 2D","pca 5D", "umap nn = 15","umap nn = 121")
rownames(paired_Ttest_m_matg_BILAT) <- c("pca 5D", "umap nn = 15","umap nn = 121", "umap nn = 208" )
paired_Ttest_m_matg_BILAT <-  format(paired_Ttest_m_matg_BILAT , scientific=TRUE,digits=3)
table_globale_BILAT <- ggtexttable(paired_Ttest_m_matg_BILAT , rows =rownames(paired_Ttest_m_matg_BILAT),
                                   theme = ttheme("classic"))+ labs(title = "Alternative less")
table_globale_BILAT 

paired_Ttest_mg_GREATER <-pairwise.t.test(pwt_df_208$mean_cp, pwt_df_208$method,  p.adj = "holm",  paired = TRUE , alternative = "greater")$p.value
paired_Ttest_m_matg_GREATER <- as.matrix(paired_Ttest_mg_GREATER)
colnames(paired_Ttest_m_matg_GREATER) <- c("pca 2D","pca 5D", "umap nn = 15","umap nn = 121")
rownames(paired_Ttest_m_matg_GREATER) <- c("pca 5D", "umap nn = 15","umap nn = 121", "umap nn = 208" )
paired_Ttest_m_matg_GREATER <-  format(paired_Ttest_m_matg_GREATER , scientific=TRUE,digits=3)
table_globale_GREATER <- ggtexttable(paired_Ttest_m_matg_GREATER , rows =rownames(paired_Ttest_m_matg_GREATER),
                                     theme = ttheme("classic"))+ labs(title = "Alternative less")
table_globale_GREATER


paired_Ttest_mg_LESS <-pairwise.t.test(pwt_df_208$mean_cp, pwt_df_208$method,  p.adj = "holm",  paired = TRUE , alternative = "less")$p.value 
paired_Ttest_m_matg_LESS <- as.matrix(paired_Ttest_mg_LESS)
colnames(paired_Ttest_m_matg_LESS) <- c("pca 2D","pca 5D", "umap nn = 15","umap nn = 121")
rownames(paired_Ttest_m_matg_LESS) <- c("pca 5D", "umap nn = 15","umap nn = 121", "umap nn = 208" )
paired_Ttest_m_matg_LESS <-  format(paired_Ttest_m_matg_LESS , scientific=TRUE,digits=3)
table_globale_LESS <- ggtexttable(paired_Ttest_m_matg_LESS , rows =rownames(paired_Ttest_m_matg_LESS),
                                     theme = ttheme("classic"))+ labs(title = "Alternative less")
table_globale_LESS

ggarrange(table_globale_BILAT , table_globale_GREATER, table_globale_LESS, nrow = 3)


# Cp  Maps ----------------------------------------------------------------

CP_df_res_K208 <- Main_CP_res_NN[[1]][Main_CP_res_NN[[1]]$K == 208,]
CP_df_res_K208 <- CP_df_res_K208[order(CP_df_res_K208$Sample_ID),]
pUMAP_CP208 <- ggplot(umap208_res_df, aes(x=V1, y=V2,  color=(CP_df_res_K208$`umap_md=0.1_nn=208`/max(CP_df_res_K208$`umap_md=0.1_nn=208`)))) +  geom_point(size=4, alpha =.8) + scale_color_distiller(palette = "Spectral")
pUMAP_CP208 <- pUMAP_CP208 +  labs(title="", 
                             y=TeX("dim2"), x="dim1") +
  theme( legend.position = "right",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=16),  # X axis title
         axis.title.y=element_text(size=16),  # Y axis title
         axis.text.x=element_text(size=14),  # X axis text
         axis.text.y=element_text(size=14),
         legend.text = element_text(size = 8) ,
         legend.title = element_blank())#+  guides(col = guide_legend( ncol = 2))  # Y axis text

pUMAP_CP208



CP_df_res_K208 <- Main_CP_res_NN[[1]][Main_CP_res_NN[[1]]$K == 208,]
CP_df_res_K208 <- CP_df_res_K208[order(CP_df_res_K208$Sample_ID),]

pgenes_CP208 <- ggplot(umap208_res_df, aes(x=V1, y=V2,  color=CP_df_res_K208$REF/max(CP_df_res_K208$REF))) +  geom_point(size=4, alpha =.8) + scale_color_distiller(palette = "Spectral")
pgenes_CP208 <- pgenes_CP208  +  labs(title="", 
                                   y=TeX("dim2"), x="dim1") +
  theme( legend.position ="left",
         plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=16),  # X axis title
         axis.title.y=element_text(size=16),  # Y axis title
         axis.text.x=element_text(size=14),  # X axis text
         axis.text.y=element_text(size=14),
         legend.text = element_text(size = 8) ,
         legend.title = element_blank())#+  guides(col = guide_legend( ncol = 2))  # Y axis text

pgenes_CP208


Main_stat_table_cp <- matrix(ncol =5, nrow = 5)
colnames(Main_stat_table_cp) = c("PCA 2D", "PCA 5D", "UMAP nn = 15", "UMAP nn = 121", "UMAP nn = 208")
rownames(Main_stat_table_cp) = c("PCA 2D", "PCA 5D", "UMAP nn = 15", "UMAP nn = 121", "UMAP nn = 208")
Main_stat_table_cp[ 1,] = c("-", "1.00e+00", "5.41e-08", "1.00e+00", "1.00e+00")
Main_stat_table_cp[ 2,] = c("3.29e-19", "-", "7.87e-09", "8.96e-07", "4.46e-04")
Main_stat_table_cp[ 3,] = c("1.00e+00", "1.00e+00", "-", "1.00e+00", "1.00e+00")
Main_stat_table_cp[ 4,] = c("2.57e+01", "1.00e+00", "1.28e-07", "-", "1.00e+00")
Main_stat_table_cp[ 5,] = c("7.55e+04", "1.00e+00", "1.30e-07", "6.68e-0.4", "-")
Main_stat_table_cp

tab <- ggtexttable(Main_stat_table_cp, 
                   theme = ttheme(
                     colnames.style = colnames_style(face = "bold"),
                     rownames.style = rownames_style(face = "bold"))
)
tab
tab <- table_cell_bg(tab, row = 2, column = 2,
                     fill="darkgray")
tab <- table_cell_bg(tab, row = 3, column = 3,
                     fill="darkgray")
tab <- table_cell_bg(tab, row = 4, column = 4,
                     fill="darkgray")
tab <- table_cell_bg(tab, row = 5, column = 5,
                     fill="darkgray")
tab <- table_cell_bg(tab, row = 6, column = 6,
                     fill="darkgray")


tab

ggarrange(ggarrange(CP_Main_graph,tab, nrow = 2, heights = c(0.6,0.4), labels = c("A","B")), ggarrange(pUMAP_CP208, pgenes_CP208, nrow = 2, common.legend = T, labels = "C"), ncol = 2, widths = c(0.6,.4))




# Moran Index ------------------------------------------------------------

# Mooran Index calculation ------------------------------------------------

List_coords <- list('PCA 5D' = data.frame(acp_5D_li_df), 'umap_nn=208' =  data.frame(umap208_res_df[,1:3]), "gene_expr_mat" = data.frame(vst50_TCACLCNECSCLC_VSamp) )
vst50_TCACLCNECSCLC_VSamp <- vst50_TCACLCNECSCLC_VSamp[order(vst50_TCACLCNECSCLC_VSamp$Sample_ID),]

# start_time <- Sys.time()
# no_cores <- detectCores() ; no_cores
# cl <- makeCluster(no_cores)
# registerDoParallel(cl)
# MI_foreach <- foreach(i=1:6398,.combine=c, .packages=c("spdep",'FNN', 'base') ) %dopar%{
#   c_vst50_TCACLCNECSCLC_VSamp <- data.frame("Sample_ID"  = vst50_TCACLCNECSCLC_VSamp[,1],  vst50_TCACLCNECSCLC_VSamp[,i])
#   colnames(c_vst50_TCACLCNECSCLC_VSamp)[2] <- colnames(vst50_TCACLCNECSCLC_VSamp)[i]
#   moran_I_main(l_coords_data = List_coords  , spatial_data =data.frame(c_vst50_TCACLCNECSCLC_VSamp), listK= seq(2,208,15), nsim = 10, Stat=FALSE, Graph = FALSE, methods_name = c('pca5D','umap121','gene_expr_mat'))
# }
# stopCluster(cl)
# end_time <- Sys.time()
# start_time - end_time
#saveRDS(MI_foreach, 'MI_foreach108.rds')
# 3h30 !

MI_foreachR <- readRDS("/RData/MI_foreach0730.rds")


# Moran index overlapping between lists  -----------------------------------


## Main graphics -----------------------------------------------------------

MI_mean_gene_PCA5D <- unlist(lapply(1:6398, function(i){
  c_gene_PCA5D <- mean(unlist(lapply(1:14, function(j){
    MI_foreachR[[i]][1,1,j]
  })) )
}))


MI_mean_gene_UMAP208 <- unlist(lapply(1:6398, function(i){
  c_gene_PCA5D <- mean(unlist(lapply(1:14, function(j){
    MI_foreachR[[i]][2,1,j]
  })) )
}))

MI_mean_gene_ExprMat <- unlist(lapply(1:6398, function(i){
  c_gene_PCA5D <- mean(unlist(lapply(1:14, function(j){
    MI_foreachR[[i]][3,1,j]
  })) )
}))

MI_gene_name <- unlist(lapply(1:6398, function(i){
  colnames(MI_foreachR[[i]])
}))

MI_mean_genes_df <- data.frame('Gene_names' = MI_gene_name, 'PCA5D' = MI_mean_gene_PCA5D , "UMAP208" = MI_mean_gene_UMAP208,  "GeneExprMat" = MI_mean_gene_ExprMat)


PCA5D_order <- MI_mean_genes_df[order(MI_mean_genes_df$PCA5D, decreasing = T),]
PCA5D_order_genes <- as.character(PCA5D_order$Gene_names)
PCA5D_order_genes[1:10]

UMAP208_order <- MI_mean_genes_df[order(MI_mean_genes_df$UMAP208, decreasing = T),]
UMAP208_order_genes <- as.character(UMAP208_order$Gene_names)
UMAP208_order_genes[1:10]

GeneExpr_order <- MI_mean_genes_df[order(MI_mean_genes_df$GeneExprMat, decreasing = T),]
GeneExpr_order_genes <-  as.character(GeneExpr_order$Gene_names)
GeneExpr_order_genes[1:10]


PCA5D_GeneExprMat_overlap <-unlist(lapply(1:6398, function(i){
  length(intersect(PCA5D_order_genes[1:i], GeneExpr_order_genes[1:i])) /i 
}))

UMAP208_GeneExprMat_overlap <-unlist(lapply(1:6398, function(i){
  length(intersect(UMAP208_order_genes[1:i], GeneExpr_order_genes[1:i])) /i 
}))


PCA5D_UMAP208_overlap <- unlist(lapply(1:6398, function(i){
  length(intersect(UMAP208_order_genes[1:i], PCA5D_order_genes[1:i])) /i 
}))

PCA5D_UMAP208_GeneExprMat_overlap <- unlist(lapply(1:6398, function(i){
  length(intersect( intersect(UMAP208_order_genes[1:i], PCA5D_order_genes[1:i]), GeneExpr_order_genes[1:i])) /i 
}))



Expected <- unlist(lapply(1:6398, function(i){
  (i *(i/6398)^2) / i
}))

cols <- c("Expected" =  "#FDE725FF", "PCA 5D - GeneExprMat"="#440154FF","PCA 5D - UMAP 208 - GeneExprMat" =  "#5DC863FF",  "PCA5D - UMAP208"= "#21908CFF", "UMAP 208 - GeneExprMat"="#3B528BFF")

p_overlap <- ggplot() 
p_overlap <- p_overlap + geom_point(aes(x=seq(1:6398), y=PCA5D_GeneExprMat_overlap , colour="PCA 5D - GeneExprMat"), size = 0.1)
p_overlap <- p_overlap + geom_point(aes(seq(1:6398),  y=UMAP208_GeneExprMat_overlap, colour="UMAP 208 - GeneExprMat"), size = 0.1)
p_overlap <- p_overlap + geom_point(aes(seq(1:6398),  y=PCA5D_UMAP208_overlap , colour="PCA5D - UMAP208"), size = 0.1)
p_overlap <- p_overlap + geom_point(aes(seq(1:6398),  y=PCA5D_UMAP208_GeneExprMat_overlap, colour="PCA 5D - UMAP 208 - GeneExprMat"), size = 0.1)
p_overlap <- p_overlap + geom_line(aes(seq(1:6398),  y=Expected, colour="Expected"), size = 1)
p_overlap <- p_overlap +  labs(title="", 
                               y="Overlapping %", x="N", caption = "") + 
   scale_colour_manual(values=cols, labels =c("Expected", "PCA 5D - GeneExprMat","PCA 5D - UMAP 208 - GeneExprMat",  "PCA5D - UMAP208", "UMAP 208 - GeneExprMat")) +
  theme( legend.position =  "bottom",
    plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
    plot.subtitle =element_text(size=14, hjust=0.5),
    plot.caption =element_text(size=12,  hjust=0.5),
    axis.title.x=element_text(size=16),  # X axis title
    axis.title.y=element_text(size=16),  # Y axis title
    axis.text.x=element_text(size=14),  # X axis text
    axis.text.y=element_text(size=14),
    legend.text = element_text(size = 14) ,
    legend.title = element_blank())  + geom_vline(xintercept = 1000, linetype="dotted") +  guides(col = guide_legend( ncol = 2))  # Y axis text # Y axis text


p_ovelap_zoom <- ggplot() 
p_ovelap_zoom <- p_ovelap_zoom + geom_point(aes(x=seq(1:1050), y=PCA5D_GeneExprMat_overlap[1:1050] , colour="PCA 5D - GeneExprMat"), size = 0.1)
p_ovelap_zoom <- p_ovelap_zoom + geom_point(aes(seq(1:1050),  y=UMAP208_GeneExprMat_overlap[1:1050], colour="UMAP 208 - GeneExprMat"), size = 0.1)
p_ovelap_zoom <- p_ovelap_zoom + geom_point(aes(seq(1:1050),  y=PCA5D_UMAP208_overlap[1:1050] , colour="PCA5D - UMAP208"), size = 0.1)
p_ovelap_zoom <- p_ovelap_zoom + geom_point(aes(seq(1:1050),  y=PCA5D_UMAP208_GeneExprMat_overlap[1:1050], colour="PCA 5D - UMAP 208 - GeneExprMat"), size = 0.1)
p_ovelap_zoom <- p_ovelap_zoom + geom_line(aes(seq(1:1050),  y=Expected[1:1050], colour="Expected"), size = 1)
p_ovelap_zoom <- p_ovelap_zoom +  labs(title="", 
                               y="", x="", caption = "") + 
  scale_colour_manual(values=cols) +
  theme( legend.position = "none",
    plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
    plot.subtitle =element_text(size=14, hjust=0.5),
    plot.caption =element_text(size=12,  hjust=0.5),
    axis.title.x=element_text(size=16),  # X axis title
    axis.title.y=element_text(size=16),  # Y axis title
    axis.text.x=element_text(size=14),  # X axis text
    axis.text.y=element_text(size=14),
    legend.text = element_text(size = 14) ,
    legend.title = element_blank())  + 
  geom_vline(xintercept = 1000, linetype="dashed")# Y axis text


#print(p_ovelap_zoom) 


t.test(PCA5D_GeneExprMat_overlap, UMAP208_GeneExprMat_overlap, paired = T, alternative = 'greater')

# Euler diagram ------------------------------------------------------------

first_cor_list = list("PCA5D" = PCA5D_order_genes[1:1000], "UMAP 208" = UMAP208_order_genes[1:1000], "GenesExprMat" =GeneExpr_order_genes[1:1000])

  print(venn(first_cor_list, show.plot=FALSE))


  library(VennDiagram)
  FirstList <- list(PCA5D_order_genes1000, UMAP208_order_genes1000, GeneExpr_order_genes1000 )
  venn(FirstList, show.plot=TRUE)
  
  library(eulerr)
  
  
  MyVenn<- c(c(A = 11, B = 83, C = 37, 
               "A&B" = 35, "A&C" = 9, "B&C"= 81,"A&B&C" = 873))
  PEul2 <- plot(euler(MyVenn), labels = c("PCA 5D ", "UMAP \n nn = 208 ", "GenesExprMat"), fills = c("#5DC863FF", "#FDE725FF", "#21908CFF"),  quantities = TRUE, alpha= .6)
  PEul2
  
ggarrange(p_overlap, NULL, PEul2, labels = c("A", "B",""), ncol = 3, widths = c(.7, .15 , .3))


# MI_OTP Example ----------------------------------------------------------

Embl_OTP ="ENSG00000171540.6"
otp_df <- data.frame("Sample_ID" = vst50_TCACLCNECSCLC_VSamp$Sample_ID, "OTP" =  vst50_TCACLCNECSCLC_VSamp[,which(colnames(vst50_TCACLCNECSCLC_VSamp)== Embl_OTP)] )
otp_df <- otp_df[order(otp_df$Sample_ID),]
MI_OTP <- moran_I_main(l_coords_data = List_coords, spatial_data = data.frame(otp_df), listK= seq(1,207,1), nsim = 10, Stat=FALSE, Graph = TRUE, methods_name = c('PCA 5D','UMAP nn = 208', "GeneExprMat"))
MI_OTP_by_k = moran_I_scatter_plot_by_k(data = as.array(MI_OTP)[[1]], Xlab = NULL, Ylab=NULL, Title= NULL)
MI_otp_plot <- MI_OTP_by_k[[2]]
#c(viridisLite::viridis(3), viridisLite::viridis(4)[2]))
p_otp1 <- ggplot(MI_otp_plot , aes(x=K_level, y=moranI, color= Methods )) + geom_line(aes( color=Methods )) + geom_point()+
  scale_color_manual(values =c("PCA 5D" = "#440154FF",  "UMAP nn = 208"  ="#FDE725FF" , "GeneExprMat" = "#21908CFF") , breaks= c("PCA 5D", "UMAP nn = 208", 'GeneExprMat' ) ,labels = c("PCA 5D", "UMAP nn = 208", "GenesExprMat"))# + scale_shape_manual(values = c(16, 7, 15,17,8 ))#wes_palette(n=5, name="FantasticFox1")
p_otp1 <- p_otp1 +  labs(title="",   y="Moran Index", x="K") +theme( legend.position =  c(0.8,0.9),
                                                                     plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
                                                                     plot.subtitle =element_text(size=14, hjust=0.5),
                                                                     plot.caption =element_text(size=12,  hjust=0.5),
                                                                     axis.title.x=element_text(size=16),  # X axis title
                                                                     axis.title.y=element_text(size=16),  # Y axis title
                                                                     axis.text.x=element_text(size=14),  # X axis text
                                                                     axis.text.y=element_text(size=14),
                                                                     legend.text = element_text(size = 14) ,
                                                                     legend.title = element_blank())  # Y axis text


OTPMap_df <- merge(umap208_res_df,otp_df, by = "Sample_ID")

pOTPMap <- ggplot(OTPMap_df, aes(x=V1, y=V2,  color=OTP)) +  geom_point(size=4, alpha =.8) + scale_color_distiller(palette = "Spectral")
pOTPMap <- pOTPMap+  labs(title="", 
                              y=TeX("dim2"), x="dim1") +
  theme( 
    plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
    plot.subtitle =element_text(size=14, hjust=0.5),
    plot.caption =element_text(size=12,  hjust=0.5),
    axis.title.x=element_text(size=16),  # X axis title
    axis.title.y=element_text(size=16),  # Y axis title
    axis.text.x=element_text(size=14),  # X axis text
    axis.text.y=element_text(size=14),
    legend.text = element_text(size = 14) ,
    legend.title = element_blank())#+  guides(col = guide_legend( ncol = 2))  # Y axis text





# MI MKI67 ----------------------------------------------------------------


Embl_MKI67 = as.character(Ref_gene$V1[Ref_gene$V7 == "MKI67"])
MKI67_df <- data.frame("Sample_ID" = vst50_TCACLCNECSCLC_VSamp$Sample_ID, "MKI67" =  vst50_TCACLCNECSCLC_VSamp[,which(colnames(vst50_TCACLCNECSCLC_VSamp)== Embl_MKI67)] )
MKI67_df <- MKI67_df[order(MKI67_df$Sample_ID),]
MI_MKI67 <- moran_I_main(l_coords_data = List_coords, spatial_data = data.frame(MKI67_df), listK= seq(1,207,1), nsim = 10, Stat=FALSE, Graph = TRUE, methods_name = c('PCA 5D','UMAP nn = 208', "GeneExprMat"))
MI_MKI67_by_k = moran_I_scatter_plot_by_k(data = as.array(MI_MKI67)[[1]], Xlab = NULL, Ylab=NULL, Title= NULL)
MI_MKI67_plot <- MI_MKI67_by_k[[2]]
MI_MKI67_plot$Methods <- as.character(MI_MKI67_plot$Methods)

p_MKI671 <- ggplot(MI_MKI67_plot , aes(x=K_level, y=moranI, color= Methods )) + geom_line(aes( color=Methods )) + geom_point()+
  scale_color_manual(values =c("PCA 5D" = "#440154FF",  "UMAP nn = 208"  ="#FDE725FF" , "GeneExprMat" = "#21908CFF") , breaks= c("PCA 5D", "UMAP nn = 208", 'GeneExprMat' ) ,labels = c("PCA 5D", "UMAP nn = 208", "GenesExprMat"))# + scale_shape_manual(values = c(16, 7, 15,17,8 ))#wes_palette(n=5, name="FantasticFox1")
p_MKI671 <- p_MKI671 +  labs(title="",   y="Moran Index", x="u") +theme( legend.position =  c(0.8,0.9),
                                                                     plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
                                                                     plot.subtitle =element_text(size=14, hjust=0.5),
                                                                     plot.caption =element_text(size=12,  hjust=0.5),
                                                                     axis.title.x=element_text(size=16),  # X axis title
                                                                     axis.title.y=element_text(size=16),  # Y axis title
                                                                     axis.text.x=element_text(size=14),  # X axis text
                                                                     axis.text.y=element_text(size=14),
                                                                     legend.text = element_text(size = 14) ,
                                                                     legend.title = element_blank())  # Y axis text


MKI67Map_df <- merge(umap208_res_df, MKI67_df, by = "Sample_ID")

pMKI67Map <- ggplot(MKI67Map_df, aes(x=V1, y=V2,  color=MKI67)) +  geom_point(size=4, alpha =.8) + scale_color_distiller(palette = "Spectral")
pMKI67Map <- pMKI67Map+  labs(title="", 
                             y=TeX("dim2"), x="dim1") +
  theme( legend.position = c(0.2,0.2),
    plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
    plot.subtitle =element_text(size=14, hjust=0.5),
    plot.caption =element_text(size=12,  hjust=0.5),
    axis.title.x=element_text(size=16),  # X axis title
    axis.title.y=element_text(size=16),  # Y axis title
    axis.text.x=element_text(size=14),  # X axis text
    axis.text.y=element_text(size=14),
    legend.text = element_text(size = 14) ,
    legend.title = element_blank())#+  guides(col = guide_legend( ncol = 2))  # Y axis text

pMKI67Map


ggarrange(p_MKI671, pMKI67Map, ncol = 2, labels = c("A", "B"))



# MoranIndex theoric analysis ---------------------------------------------

# umap208_Moran_df <- data.frame()
# umap208_Moran_df <- cbind("umap208_res_df"= umap208_res_df, "Molecular_cluster" =Attributes_UMAP_TCACLCNECSCL$Molecular_clusters )
# MidX <- ifelse(umap208_Moran_df$umap208_res_df.V1 < mean(umap208_Moran_df$umap208_res_df.V1), 0, 1)
# MidY <- ifelse(umap208_Moran_df$umap208_res_df.V2 < mean(umap208_Moran_df$umap208_res_df.V2), 0, 1)
# Uniform <- rep(1, length(umap208_Moran_df$umap208_res_df.V2) -1)
# Uniform <- c(Uniform , 2)
# alea <- rnorm(208,mean = 3,sd=207)
# spa_MI_test <- data_frame("Sample_ID" = umap208_Moran_df$umap208_res_df.Sample_ID,  "MidX" = MidX, "MidY"= MidY, "Uniform" = Uniform, "Alea" = alea )
# 
# List_coords <- list('umap_nn=208' =  data.frame(umap208_res_df[,1:3]) ) 
# 
# MI_MidX <- moran_I_main(l_coords_data = List_coords, spatial_data = data.frame(spa_MI_test), listK= seq(1,207,1), nsim = 10, Stat=FALSE, Graph = TRUE, methods_name = c('UMAP nn = 208'))
# 
# MI_MidY <- moran_I_main(l_coords_data = List_coords, spatial_data = data.frame("Sample_ID" = umap208_Moran_df$umap208_res_df.Sample_ID, "MidY" = MidY  ), listK= seq(1,207,1), nsim = 10, Stat=FALSE, Graph = TRUE, methods_name = c('UMAP nn = 208'))
# MI_MidY = moran_I_scatter_plot_by_k(data = as.array(MI_MidY)[[1]], Xlab = NULL, Ylab=NULL, Title= NULL)
# 
# attX <- rep(colnames(as.array(MI_MidX)[[1]])[1],dim(as.array(MI_MidX)[[1]])[3])
# attY <- rep(colnames(as.array(MI_MidX)[[1]])[2],dim(as.array(MI_MidX)[[1]])[3])
# attU <- rep(colnames(as.array(MI_MidX)[[1]])[3],dim(as.array(MI_MidX)[[1]])[3])
# attA <- rep(colnames(as.array(MI_MidX)[[1]])[4],dim(as.array(MI_MidX)[[1]])[3])
# 
# 
# att <- data.frame(c(attX, attY, attU, attA))
# att
# mix <- c()
# miy <- c()
# miu <- c()
# mia <- c()
# 
# for(i in 1:207){
#   mix <- c(mix ,as.array(MI_MidX)[[1]][1,1,i] )
#   miy <- c(miy ,as.array(MI_MidX)[[1]][1,2,i] )
#   miu <- c(miu ,as.array(MI_MidX)[[1]][1,3,i] )
#   miu <- c(miu ,as.array(MI_MidX)[[1]][1,4,i] )
#   
#   
# }
# 
# att <- cbind(att, c(mix,miy,miu, mia))
# klevel <- rep(seq(1:207),4)
# att  <- cbind(att,klevel )
# colnames(att) = c("Att" ,"Mi", "klevel")
# 
# p_Mi_exam <- ggplot(att , aes(x=klevel, y=Mi, color= Att )) + geom_line()+ geom_point()
# p_Mi_exam
# 
# 
# p2_standard_nn208X<- ggplot(umap208_res_df, aes(x=V1, y=V2,  color=as.factor(MidX))) +  geom_point(size=4, alpha =.8) + scale_color_brewer(palette="Spectral") +labs(title= "208")+ theme(legend.title = element_blank()) + theme(legend.position = "none")
# p2_standard_nn208X <- p2_standard_nn208X +  labs(title="", y="dim2", x="dim1")
# p2_standard_nn208X
# 
# p2_standard_nn208Y<- ggplot(umap208_res_df, aes(x=V1, y=V2,  color=as.factor(MidY))) +  geom_point(size=4, alpha =.8) + scale_color_brewer(palette="Spectral") +labs(title= "208")+ theme(legend.title = element_blank()) + theme(legend.position = "none")
# p2_standard_nn208Y <- p2_standard_nn208Y +  labs(title="", y="dim2", Y="dim1")
# p2_standard_nn208Y
# 
# p2_standard_nn208U<- ggplot(umap208_res_df, aes(x=V1, y=V2,  color=as.factor(Uniform))) +  geom_point(size=4, alpha =.8) + scale_color_brewer(palette="Spectral") +labs(title= "208")+ theme(legend.title = element_blank()) + theme(legend.position = "none")
# p2_standard_nn208U <- p2_standard_nn208U +  labs(title="", y="dim2", U="dim1")
# p2_standard_nn208U
# 
# p2_standard_nn208A<- ggplot(umap208_res_df, aes(x=V1, y=V2,  color=alea)) +  geom_point(size=4, alpha =.8) + scale_color_distiller(palette="Spectral") +labs(title= "208")+ theme(legend.title = element_blank()) + theme(legend.position = "none")
# p2_standard_nn208A <- p2_standard_nn208A +  labs(title="", Y="dim2", A="dim1")
# p2_standard_nn208A
# 
# 
# ggarrange(p_Mi_exam, ggarrange(p2_standard_nn208X, p2_standard_nn208Y, p2_standard_nn208U, p2_standard_nn208A, nrow = 4), ncol = 2)
# 
# 



# UMAP gganimate ----------------------------------------------------------

ggplot(gapminder, aes(gdpPercap, lifeExp, size = pop, colour = country)) +
  geom_point(alpha = 0.7, show.legend = FALSE) +
  scale_colour_manual(values = country_colors) +
  scale_size(range = c(2, 12)) +
  scale_x_log10() +
  facet_wrap(~continent) +
  # Here comes the gganimate specific bits
  labs(title = 'Year: {frame_time}', x = 'GDP per capita', y = 'life expectancy') +
  transition_time(year) +
  ease_aes('linear')


umap_res_df <- data.frame()
for (i in seq(10,208,5)){
  print(i)
  umap_c_nn = umap(vst50_TCACLCNECSCLC_designRD, n_neighbors = i)
  umap_c_nn_res_df <- as.data.frame(umap_c_nn$layout)
  umap_c_nn_res_df  = setDT(umap_c_nn_res_df , keep.rownames = TRUE)[]
  colnames(umap_c_nn_res_df)[1] <- "Sample_ID"
  umap_c_nn_res_df <- umap_c_nn_res_df[order(umap_c_nn_res_df$Sample_ID),]
  umap_c_nn_res_df <- cbind(umap_c_nn_res_df , "nn" = rep(i, dim(umap_c_nn_res_df)[1]))
  umap_c_nn_res_df <- cbind(umap_c_nn_res_df , "Molecular_clusters" = Attributes_UMAP_TCACLCNECSCL$Molecular_clusters)
  umap_res_df <- rbind(umap_res_df, umap_c_nn_res_df)
}

umap_res_dfV2 <- umap_res_df[umap_res_df$nn > 20]
umap_res_dfV2 <- umap_res_dfV2[complete.cases(umap_res_dfV2),]
p <- ggplot(umap_res_dfV2 , aes(x=V1, y=V2,  color=Molecular_clusters)) +  geom_point()+
  scale_color_brewer(palette="Spectral") 
p <- p + transition_time(nn)+labs(title = 'n_neighbors:{frame_time}')  + ease_aes('linear')


