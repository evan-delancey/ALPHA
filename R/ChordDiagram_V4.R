####################################################################
#------------------------------------------------------------------#
#------------------------------------------------------------------#
#Chord diagram                                                     #
#Filename: "ChordDiagram_V4.R"                                     #
#Written and developed by Alex O. Onojeghuo (PhD)                  #
#Alberta Biodiversity Monitoring Institute, August 20 2018         #
#------------------------------------------------------------------#
#------------------------------------------------------------------#
####################################################################

##Load library
library(xlsx)
library(circlize)
library(igraph)

##Specify output folder
OUTPUT <- "F:/8_ChordDiagram"

##Set working directory
setwd("F:/8_ChordDiagram/data")
##Read in .csv file
##The csv file contained r values of the input varibles as matrix 
dat = read.csv("correlation_excel_V5_29Variables.csv", sep=",", row.names = 1)
matrix.dat <- as.matrix(dat)

## Set grid colour based on type of input variables (i.e. Sentinel-1, Sentinel-2, and Topographic)
grid.col = c(VH = "red", VHsd = "red", VV = "red", VVsd = "red", Dis_VH = "red", Ent_VH = "red", DPOL = "red", ARI = "green", EVI = "green", IRECI = "green", 
             MNDWI = "green", NDVI = "green", NDVI705 = "green", NDWI = "green", PC1 = "green", PC1_10 = "green", PC2 = "green", PC2_10 = "green", PC3 = "green", 
             PC3_10 = "green", REIP = "green", C10 = "blue", HAND1 = "blue", HAND2 = "blue", SWI = "blue", TPI50 = "blue", TPI500 = "blue", TRI = "blue", VBF = "blue")

## Make the circular plot
chordDiagram(matrix.dat, grid.col = grid.col, transparency = 0.5)
