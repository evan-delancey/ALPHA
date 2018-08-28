####################################################################
#------------------------------------------------------------------#
#------------------------------------------------------------------#
#Bean plot chart preparation for data selection evaluation         #
#Filename: "InputDataExploration_Beanplot_V7.R"                    #
#Written and developed by Alex O. Onojeghuo (PhD)                  #
#Alberta Biodiversity Monitoring Institute, August 18 2018         #
#------------------------------------------------------------------#
#------------------------------------------------------------------#
####################################################################

##Load libraries
library(ggplot2)
library(dplyr)
library(raster)
library(rgdal)
library(gtable)
library(grid)
library(zoo)
library(reshape2)
library(gridExtra)
library(beanplot)
##install.packages(c("raster", "sp"))

rasterOptions(maxmemory = 1e+09,tmpdir = "D:/RtmpRasterDump")

train.pts <- readOGR("D:/Wetland_Mapping_Grassland_AB/7_VariableSelection/ViolinPlot", "Training_data_v2")

setwd("D:/Wetland_Mapping_Grassland_AB/3_Input_layers")

### List tif files in directory
list.files(pattern = "*.tif$")

###List required files
fls <- c("VHmean.tif", "VHsd.tif", "VVmean.tif", "VVsd.tif", "Dissimilarity_VHmean.tif", "Entropy_VHmean.tif", 
         "DPOL_S1.tif", "ARI.tif", "EVI.tif", "IRECI.tif", "MNDWI.tif", "NDVI.tif", "NDVI705.tif", "NDWI.tif", 
         "PC1.tif", "PC1_10.tif", "PC2.tif", "PC2_10.tif", "PC3.tif", "PC3_10.tif", "REIP.tif", "C10_LIDAR.tif",
         "HAND1_SRTM.tif", "HAND2_SRTM.tif", "SWI_LIDAR.tif", "TPI50_LIDAR.tif", "TPI500_LIDAR.tif", "TRI_LIDAR.tif", 
         "VBF_LIDAR.tif")

##Read in training data	
train <- raster("D:/Wetland_Mapping_Grassland_AB/6_GEE_inputs/training_v2.tif")
#Set output directory
OUTPUTS <- "D:/Wetland_Mapping_Grassland_AB/7_VariableSelection/BeanPlot_v2/Result"

##Extract training points from raster
df <- c(1:length(train.pts))
for (i in 1:length(fls)){
	r <- raster(fls[i])
	ext <- extract(r, train.pts)
	df <- cbind(df, ext)
	print(i)
}
##Extract training points from training data
ext <- extract(train, train.pts)
##Create dataframe between extracted points and training points
df <- cbind(df, ext)
##exclude first column from dataframe
df <- df[,-1]
##Remove all NA values
df <- na.omit(df)
##Define column names for dataframe
colnames(df) <- c("VH", "VHsd", "VV", "VVsd", "Dis_VH", "Ent_VH", "DPOL", "ARI", "EVI", "IRECI", "MNDWI", "NDVI", "NDVI705", "NDWI", 
                  "PC1", "PC1_10", "PC2", "PC2_10", "PC3", "PC3_10", "REIP", "C10", "HAND1", "HAND2", "SWI", "TPI50", "TPI500", "TRI", 
                  "VBF", "LC_class")
##If df is not a dataframe yet, convert to dataframe
df <- as.data.frame(df)
##check correlation of columns in dataframe
cor(df)
##Reclassify data with new class numbers
df$class <- ifelse(df$LC_class == 1, "Upland", ifelse(df$LC_class ==2, "Wetland", "Open water"))
##Set working directory to store outputs
setwd(OUTPUTS)

##Split screen for plotting & output results as pdf
Funbean<-function(x){
X1<-cbind(df$class,x)
X1a<-X1[ which (X1$class == "Wetland"), ] 
X1b<-X1[ which (X1$class == "Upland"), ] 
X2<-as.data.frame(X1b[,2],X1a[2,])}
X5<-list(df[,c(1,31)],df[c(2,31)],df[c(3,31)],df[c(4,31)],df[c(5,31)],df[c(6,31)],df[c(7,31)],df[c(8,31)],df[c(9,31)],df[c(10,31)],df[c(11,31)],
         df[c(12,31)],df[c(13,31)],df[c(14,31)],df[c(15,31)],df[c(16,31)],df[c(17,31)], df[c(18,31)], df[c(19,31)], df[c(20,31)], df[c(21,31)], 
         df[c(22,31)], df[c(23,31)], df[c(24,31)], df[c(25,31)], df[c(26,31)], df[c(27,31)], df[c(28,31)], df[c(29,31)])

names(X5)<-c("VH", "VHsd", "VV", "VVsd", "Dis_VH", "Ent_VH", "DPOL", "ARI", "EVI", "IRECI", "MNDWI", "NDVI", "NDVI705", "NDWI", 
             "PC1", "PC1_10", "PC2", "PC2_10", "PC3", "PC3_10", "REIP", "C10", "HAND1", "HAND2", "SWI", "TPI50", "TPI500", "TRI", 
             "VBF")

Funbean<-function(X1){
  X1a<-X1[ which (X1$class == "Wetland"), ] 
  X1b<-X1[ which (X1$class == "Upland"), ] 
beanplot(X1b[,1], X1a[,1],side = "both", main = names(X1[1]),ll = 0.03, beanlinewd = 1.3, alpha = 0.8, xlab=NA, xaxt="n" ,
         col =list(c("darksalmon","gray50"), c("deepskyblue1", "gray50")), border = NA) }


pdf(file="Variables9.pdf",height=15 , width=10)
par(mfrow=c(6,5))
lapply(X5[1:29],Funbean)
legend("topright",fill = c("darksalmon", "deepskyblue1"), legend = c("Upland", "Wetland"),bty = "n", cex=0.85)
dev.off()


