#########################################################################
# WETLAND INDICES version 3.3 #

# Created by Liam Beaudette, GIS Technician
# Email: lbeaudet@ualberta.ca
# Alberta Biodiversity Monitoring Institute
# May 2018

#########################################################################

# This script is designed to calculate several topographic indices from 
# multiple digital elevation models. The indices were chosen for use in wetland
# mapping throughout the province of Alberta and include the SAGA Wetness Index,
# TOpographic Position Index (250m, 500m, 750m window sizes), 
# Multiresolution Index of Valley Bottom Flatness, Topographic Ruggedness Index, 
# and Concavity.

#########################################################################

#Import libraries
library(RSAGA)
library(raster)
library(sp)
library(rgdal)
library(tools)

rasterOptions(maxmemory = 1e+09,tmpdir = "J:/RtmpRasterDump")
myenv <- rsaga.env(path = "C:/Program Files/saga-6.3.0_x64", modules="C:/PROGRA~1/SAGA-6~1.0_X/tools", workspace = ".", parallel = TRUE)
myenv

rsaga.get.libraries(path=myenv$modules)

#Input Parameters
#Clip DEM
DEM_clip = FALSE

#SAGA Wetness Index (SWI)
SWI <- FALSE

#Topographic Position Index 250 (TPI250)
TPI250 <- FALSE

#Topographic Position Index 500 (TPI500)
TPI500 <- TRUE

#Topographic Position Index 750 (TPI750)
TPI750 <- FALSE

#Multiresolution Index of Valley Bottom Flatness (VBF)
VBF <- FALSE

#Terrain Ruggedness Index (TRI)
TRI <- FALSE

#Concavity 10 (C10)
C10 <- FALSE

#Concavity 15 (C15)
C15 <- FALSE

#Concavity 20 (C20)
C20 <- FALSE

#Set the workspace to the working folder
workdir <- "V:/Wetlands/SAGA_wetlands"
setwd(workdir)
getwd()

#Results location
results_folder <- "V:/Wetlands/SAGA_wetlands/Results"

#Folder containing index layers (Should be in the same location as your working directory and shapefiles)
index_folder <- "Index_layers"

#Filepath to the DEMS in .tif format
DEMS_dir <- "J:/LandCover/DEMs/Clip"


########################   DEM CLIPPPING   ##############################


#Loops through DEMs folder

dem_list <- list.files(DEMS_dir, pattern = ".tif$")
for (d in dem_list){
  dem_name<- file_path_sans_ext(d)
  print(paste0("Working with DEM: ", dem_name))
  DEM <- paste0(DEMS_dir, "/", d)
  dem_type <- paste0("_" , dem_name)
  dem_crops <- paste0(results_folder, "/Indexed_DEMs")

  #Clip the DEM to the index areas if TRUE
  if (DEM_clip == TRUE){
  
    #Creates a folder 
    dir.create(dem_crops)
    print("Folder added to Results")
  
    #Loops through index folder
    fclist <- ogrListLayers(index_folder)
    for (i in fclist){
      fclass <- readOGR(index_folder, i)
      r <- raster(DEM)
      ext_r <- extent(r)
      ext_fc <- extent(fclass)
    
      #Clips each index layer to the DEM
      print("Cropping...")
      setwd(dem_crops)
      filename <- paste0(dem_crops, "/", i, dem_type)
      #crop(r, fclass, format="GTiff", filename=filename, overwrite=TRUE)
      print("Raster written")
      setwd(workdir)
    }
  } else {
    print("Clip DEM not requested")
  }
print("Clipping complete")
}  


########################   SAGA WETNESS INDEX  ##########################


#If SWI = TRUE, performs SAGA Wetness Index
if (SWI == TRUE){
  print("Performing SAGA wetness index...")
  
  #Creates a folder to store results
  folder <- print(paste0(results_folder,"/SWI"))
  dir.create(folder)
  print("SWI folder added to Results")
  
  #Obtains specific names
  for (d in dem_list){
    dem_name<- file_path_sans_ext(d)
    print(paste0("Working with DEM: ", dem_name))
    DEM <- paste0(DEMS_dir, "/", dem_name)
    dem_type <- paste0("_", dem_name)
    
    #Creates a temporary folder to store intermediate data
    swi_folder <- print(paste0(folder,"/SWI_tmp", dem_type))
    dir.create(swi_folder)
    
    #Loops through clipped DEM folder
    swi_list <- list.files(dem_crops, pattern = dem_type)
    for (i in swi_list){
      DEM_tile <- print(paste0(dem_crops, "/", i))
      outname <- file_path_sans_ext(i)
      setwd(swi_folder)
      swi_output <- print(paste0(outname, "_tmp"))
      
      #Generates SWI as SAGA Grid
      print("Generating SWI...")
      rsaga.wetness.index(DEM_tile, swi_output, env=myenv)
      setwd(workdir)
      gc()
    }

    #Converts SAGA grids to TIFF format
    tif_list <- list.files(swi_folder, full.names = FALSE, pattern = ".sdat$")
    for (i in tif_list){ 
      base <- basename(i)
      swi_tif <- file_path_sans_ext(base)
      s <- raster(i)
      filename <- paste0(swi_folder, "/", swi_tif, ".tif")
      writeRaster(s, filename=filename, format="GTiff", overwrite=TRUE)
      print(paste0(swi_tif, " saved as tif"))
    
      #Deletes intermediate files
      SAGA_files <- paste0(swi_folder, "/", swi_tif)
      mgrd <- paste0(SAGA_files, ".mgrd")
      prj <- paste0(SAGA_files, ".prj")
      xml <- paste0(SAGA_files, ".sdat.aux.xml")
      sgrd <- paste0(SAGA_files, ".sgrd")
      file.remove(i)
      file.remove(mgrd)
      file.remove(prj)
      file.remove(xml)
      file.remove(sgrd)
      print("Intermediate files deleted")
      gc()
    }
  
    #Creates a vector of all the .tif files required for the mosaic
    slist <- list.files(swi_folder, pattern = ".tif$")
    slist
    a <- list()
    for (i in slist){
      s_file <- paste0(swi_folder, "/", i)
      s <- raster(s_file)
      outname <- file_path_sans_ext(i)
      assign(paste0("r", outname), s)
      a[[which(slist==i)]] <- s
    }
    x <- a[!sapply(a, is.null)]
    
    #Mosaics the rasters within the folder
    filename <- paste0(folder, "/SWI", dem_type, ".tif")
    names(x)[1:2] <- c('x', 'y')
    x$fun <- mean
    x$na.rm <- TRUE
    x$filename <- filename
    x$overwrite <- TRUE
    y <- do.call(mosaic, x)
    plot(y)
    print("Mosaic successful")
    
    #Deletes intermediate files
    unlink(swi_folder, recursive = TRUE)
  }
print ("SWI generation Complete")

} else {
  print("SWI not requested")
}


###############   TOPOGRAPHIC POSITION INDEX - 250   ####################


#If TPI250 = TRUE, performs Topographic Position Index 250 (TPI250)
if (TPI250 == TRUE){
  print("Performing topographic position index 250...")
  
  #Creates a folder to store results
  folder <- paste0(results_folder,"/TPI250")
  dir.create(folder)
  print("TPI250 folder added to Results")
  
  #Obtains specific names
  for (d in dem_list){
    dem_name<- file_path_sans_ext(d)
    print(paste0("Working with DEM: ", dem_name))
    DEM <- paste0(DEMS_dir, "/", dem_name)
    dem_type <- paste0("_", dem_name)
  
    #Creates a temporary folder to store intermediate data
    tpi250_folder <- print(paste0(folder,"/TPI250_tmp", dem_type))
    dir.create(tpi250_folder)

    #Loops through DEM folder
    tpi250_list <- list.files(dem_crops, pattern = dem_type)
    for (i in tpi250_list){
      DEM_tile <- print(paste0(dem_crops, "/", i))
      outname <- file_path_sans_ext(i)
      setwd(tpi250_folder)
      GRID_files <- print(paste0(outname, "_tmp"))
      tpi250_output <- print(paste0(GRID_files, ".sgrd"))
      
      #Generates TPI as SAGA Grid
      print("Generating TPI250...")
      DEM <- paste0(dem_crops, "/", i) 
      TPI <- paste0(tpi250_folder, "/", tpi250_output)
      RAD_MIN <- "0"
      RAD_MAX <- "250"
      rsaga.geoprocessor("ta_morphometry", 18, list(DEM=DEM, TPI=TPI, RADIUS_MIN=RAD_MIN, RADIUS_MAX=RAD_MAX), env=myenv)
      print("TPI250 created...")
      setwd(workdir)
      gc()
    }
  
    #Converts SAGA grids to TIFF format
    tif_list <- list.files(tpi250_folder, full.names = FALSE, pattern = ".sdat$")
    for (i in tif_list){ 
      sdat <- paste0(tpi250_folder, "/", i)
      tpi250_tif <- file_path_sans_ext(i)
      setwd(tpi250_folder)
      sdat_raster <- raster(sdat)
      TPIname <- (paste0(tpi250_folder, "/", tpi250_tif, ".tif"))
      writeRaster(sdat_raster, filename=TPIname, format="GTiff", overwrite=TRUE)
      print(paste0(tpi250_tif, "file saved as .tif"))
      
      #Deletes intermediate files
      SAGA_files <- paste0(tpi250_folder, "/", tpi250_tif)
      mgrd <- paste0(SAGA_files, ".mgrd")
      prj <- paste0(SAGA_files, ".prj")
      xml <- paste0(SAGA_files, ".sdat.aux.xml")
      sgrd <- paste0(SAGA_files, ".sgrd")
      file.remove(sdat)
      file.remove(mgrd)
      file.remove(prj)
      file.remove(xml)
      file.remove(sgrd)
      print("Intermediate files deleted")
      setwd(workdir)
      gc()
    }
  
    #Creates a vector of all the .tif files required for the mosaic
    t250list <- list.files(tpi250_folder, pattern = ".tif$")
    t250list
    w <- list()
    for (i in t250list){
      t_file <- paste0(tpi250_folder, "/", i)
      t <- raster(t_file)
      outname <- file_path_sans_ext(i)
      assign(paste0("r", outname), t)
      w[[which(t250list==i)]] <- t
    }
    z <- w[!sapply(w, is.null)]
    
    #Mosaics the rasters within the folder
    filename <- paste0(folder, "/TPI250", dem_type, ".tif")
    names(z)[1:2] <- c('x', 'y')
    z$fun <- mean
    z$na.rm <- TRUE
    z$filename <- filename
    z$overwrite <- TRUE
    y <- do.call(mosaic, z)
    plot(y)
    print("Mosaic successful")
    setwd(workdir)
  
    #Deletes intermediate files
    #unlink(tpi250_folder, recursive = TRUE)
  }
  print("TPI250 generation complete")
  
} else {
  print("TPI250 not requested")
}


###############   TOPOGRAPHIC POSITION INDEX - 500   ####################


#If TPI500 = TRUE, performs Topographic Position Index 500 (TPI500)
if (TPI500 == TRUE){
  print("Performing topographic position index 500...")
  
  #Creates a folder to store results
  folder <- paste0(results_folder,"/TPI500")
  dir.create(folder)
  print("TPI500 folder added to Results")
  
  #Obtains specific names
  for (d in dem_list){
    dem_name<- file_path_sans_ext(d)
    print(paste0("Working with DEM: ", dem_name))
    DEM <- paste0(DEMS_dir, "/", dem_name)
    dem_type <- paste0("_", dem_name)
    
    #Creates a temporary folder to store intermediate data
    tpi500_folder <- print(paste0(folder,"/TPI500_tmp", dem_type))
    dir.create(tpi500_folder)
    
    #Loops through DEM folder
    tpi500_list <- list.files(dem_crops, pattern = dem_type)
    for (i in tpi500_list){
      DEM_tile <- print(paste0(dem_crops, "/", i))
      outname <- file_path_sans_ext(i)
      setwd(tpi500_folder)
      GRID_files <- print(paste0(outname, "_tmp"))
      tpi500_output <- print(paste0(GRID_files, ".sgrd"))
      
      #Generates TPI as SAGA Grid
      print("Generating TPI500...")
      DEM <- paste0(dem_crops, "/", i) 
      TPI <- paste0(tpi500_folder, "/", tpi500_output)
      RAD_MIN <- "0"
      RAD_MAX <- "500"
      rsaga.geoprocessor("ta_morphometry", 18, list(DEM=DEM, TPI=TPI, RADIUS_MIN=RAD_MIN, RADIUS_MAX=RAD_MAX), env=myenv)
      print("TPI500 created...")
      setwd(workdir)
      gc()
    }
    
    #Converts SAGA grids to TIFF format
    tif_list <- list.files(tpi500_folder, full.names = FALSE, pattern = ".sdat$")
    for (i in tif_list){ 
      sdat <- paste0(tpi500_folder, "/", i)
      tpi500_tif <- file_path_sans_ext(i)
      setwd(tpi500_folder)
      sdat_raster <- raster(sdat)
      TPIname <- (paste0(tpi500_folder, "/", tpi500_tif, ".tif"))
      writeRaster(sdat_raster, filename=TPIname, format="GTiff", overwrite=TRUE)
      print(paste0(tpi500_tif, "file saved as .tif"))
      
      #Deletes intermediate files
      SAGA_files <- paste0(tpi500_folder, "/", tpi500_tif)
      mgrd <- paste0(SAGA_files, ".mgrd")
      prj <- paste0(SAGA_files, ".prj")
      xml <- paste0(SAGA_files, ".sdat.aux.xml")
      sgrd <- paste0(SAGA_files, ".sgrd")
      file.remove(sdat)
      file.remove(mgrd)
      file.remove(prj)
      file.remove(xml)
      file.remove(sgrd)
      print("Intermediate files deleted")
      setwd(workdir)
      gc()
    }
    
    #Creates a vector of all the .tif files required for the mosaic
    t500list <- list.files(tpi500_folder, pattern = ".tif$")
    t500list
    w <- list()
    for (i in t500list){
      t_file <- paste0(tpi500_folder, "/", i)
      t <- raster(t_file)
      outname <- file_path_sans_ext(i)
      assign(paste0("r", outname), t)
      w[[which(t500list==i)]] <- t
    }
    z <- w[!sapply(w, is.null)]
    
    #Mosaics the rasters within the folder
    filename <- paste0(folder, "/TPI500", dem_type, ".tif")
    names(z)[1:2] <- c('x', 'y')
    z$fun <- mean
    z$na.rm <- TRUE
    z$filename <- filename
    z$overwrite <- TRUE
    y <- do.call(mosaic, z)
    plot(y)
    print("Mosaic successful")
    
    #Deletes intermediate files
    #unlink(tpi500_folder, recursive = TRUE)
  }
  print("TPI500 generation complete")
  
} else {
  print("TPI500 not requested")
}


###############   TOPOGRAPHIC POSITION INDEX - 750   ####################


#If TPI750 = TRUE, performs Topographic Position Index 750 (TPI750)
if (TPI750 == TRUE){
  print("Performing topographic position index 750...")
  
  #Creates a folder to store results
  folder <- paste0(results_folder,"/TPI750")
  dir.create(folder)
  print("TPI750 folder added to Results")
  
  #Obtains specific names
  for (d in dem_list){
    dem_name<- file_path_sans_ext(d)
    print(paste0("Working with DEM: ", dem_name))
    DEM <- paste0(DEMS_dir, "/", dem_name)
    dem_type <- paste0("_", dem_name)
    
    #Creates a temporary folder to store intermediate data
    tpi750_folder <- print(paste0(folder,"/TPI750_tmp", dem_type))
    dir.create(tpi750_folder)
    
    #Loops through DEM folder
    tpi750_list <- list.files(dem_crops, pattern = dem_type)
    for (i in tpi750_list){
      DEM_tile <- print(paste0(dem_crops, "/", i))
      outname <- file_path_sans_ext(i)
      setwd(tpi750_folder)
      GRID_files <- print(paste0(outname, "_tmp"))
      tpi750_output <- print(paste0(GRID_files, ".sgrd"))
      
      #Generates TPI as SAGA Grid
      print("Generating TPI750...")
      DEM <- paste0(dem_crops, "/", i) 
      TPI <- paste0(tpi750_folder, "/", tpi750_output)
      RAD_MIN <- "0"
      RAD_MAX <- "750"
      rsaga.geoprocessor("ta_morphometry", 18, list(DEM=DEM, TPI=TPI, RADIUS_MIN=RAD_MIN, RADIUS_MAX=RAD_MAX), env=myenv)
      print("TPI750 created...")
      setwd(workdir)
      gc()
    }
    
    #Converts SAGA grids to TIFF format
    tif_list <- list.files(tpi750_folder, full.names = FALSE, pattern = ".sdat$")
    for (i in tif_list){ 
      sdat <- paste0(tpi750_folder, "/", i)
      tpi750_tif <- file_path_sans_ext(i)
      setwd(tpi750_folder)
      sdat_raster <- raster(sdat)
      TPIname <- (paste0(tpi750_folder, "/", tpi750_tif, ".tif"))
      writeRaster(sdat_raster, filename=TPIname, format="GTiff", overwrite=TRUE)
      print(paste0(tpi750_tif, "file saved as .tif"))
      
      #Deletes intermediate files
      SAGA_files <- paste0(tpi750_folder, "/", tpi750_tif)
      mgrd <- paste0(SAGA_files, ".mgrd")
      prj <- paste0(SAGA_files, ".prj")
      xml <- paste0(SAGA_files, ".sdat.aux.xml")
      sgrd <- paste0(SAGA_files, ".sgrd")
      file.remove(sdat)
      file.remove(mgrd)
      file.remove(prj)
      file.remove(xml)
      file.remove(sgrd)
      print("Intermediate files deleted")
      setwd(workdir)
      gc()
    }
    
    #Creates a vector of all the .tif files required for the mosaic
    t750list <- list.files(tpi750_folder, pattern = ".tif$")
    t750list
    w <- list()
    for (i in t750list){
      t_file <- paste0(tpi750_folder, "/", i)
      t <- raster(t_file)
      outname <- file_path_sans_ext(i)
      assign(paste0("r", outname), t)
      w[[which(t750list==i)]] <- t
    }
    z <- w[!sapply(w, is.null)]
    
    #Mosaics the rasters within the folder
    filename <- paste0(folder, "/TPI750", dem_type, ".tif")
    names(z)[1:2] <- c('x', 'y')
    z$fun <- mean
    z$na.rm <- TRUE
    z$filename <- filename
    z$overwrite <- TRUE
    y <- do.call(mosaic, z)
    plot(y)
    print("Mosaic successful")
    
    #Deletes intermediate files
    #unlink(tpi750_folder, recursive = TRUE)
  }
  print("TPI750 generation complete")
  
} else {
  print("TPI750 not requested")
}


###################   VALLEY BOTTOM FLATNESS   ##########################
  
  
#If VBF = TRUE, performs Multiresolution Index of Valley Bottom Flatness
if (VBF == TRUE){
  print("Performing Multiresolution Index of Valley Bottom Flatness...")
  
  #Creates a folder to store results
  folder <- print(paste0(results_folder,"/VBF"))
  dir.create(folder)
  print("VBF folder added to results")
  
  #Obtains specific names
  for (d in dem_list){
    dem_name<- file_path_sans_ext(d)
    print(paste0("Working with DEM: ", dem_name))
    DEM <- paste0(DEMS_dir, "/", dem_name)
    dem_type <- paste0("_", dem_name)
  
    #Creates a folder to store results
    VBF_folder <- print(paste0(folder,"/VBF_tmp", dem_type))
    dir.create(VBF_folder)

    #Loops through DEM folder
    VBF_list <- list.files(dem_crops, pattern = dem_type)
    for (i in VBF_list){
      DEM_tile <- print(paste0(dem_crops, "/", i))
      outname <- file_path_sans_ext(i)
      setwd(VBF_folder)
      GRID_files <- print(paste0(outname, "_tmp"))
      VBF_output <- print(paste0(GRID_files, ".sgrd"))
      
      #Generates VBF as SAGA Grid
      print("Generating Valley Bottom Flatness...")
      DEM <- paste0(dem_crops, "/", i) 
      MRVBF <- paste0(VBF_folder, "/", VBF_output)
      rsaga.geoprocessor("ta_morphometry", 8, list(DEM=DEM, MRVBF=MRVBF), env=myenv)
      print("MRIVBF created...")
    }
      
    #Converts SAGA grids to TIFF format
    tif_list <- list.files(VBF_folder, full.names = FALSE, pattern = ".sdat$")
    for (i in tif_list){ 
      sdat <- paste0(VBF_folder, "/", i)
      VBF_tif <- file_path_sans_ext(i)
      setwd(VBF_folder)
      s <- raster(sdat)
      filename <- paste0(VBF_folder, "/", VBF_tif, ".tif")
      writeRaster(s, filename=filename, format="GTiff", overwrite=TRUE)
      print(paste0(VBF_tif, "saved as tif"))
      
      #Deletes intermediate files
      SAGA_files <- paste0(VBF_folder, "/", VBF_tif)
      mgrd <- paste0(SAGA_files, ".mgrd")
      prj <- paste0(SAGA_files, ".prj")
      xml <- paste0(SAGA_files, ".sdat.aux.xml")
      sgrd <- paste0(SAGA_files, ".sgrd")
      file.remove(sdat)
      file.remove(mgrd)
      file.remove(prj)
      file.remove(xml)
      file.remove(sgrd)
      print("Intermediate files deleted")
      setwd(workdir)
      gc()
    }
      
    #Creates a vector of all the .tif files required for the mosaic
    VBFlist <- list.files(VBF_folder, pattern = ".tif$")
    w <- list()
    for (i in VBFlist){
      t_file <- paste0(VBF_folder, "/", i)
      t <- raster(t_file)
      outname <- file_path_sans_ext(i)
      assign(paste0("r", outname), t)
      w[[which(VBFlist==i)]] <- t
    }
    z <- w[!sapply(w, is.null)]
    
    #Mosaics the rasters within the folder
    filename <- paste0(folder, "/VBF", dem_type, ".tif")
    names(z)[1:2] <- c('x', 'y')
    z$fun <- mean
    z$na.rm <- TRUE
    z$filename <- filename
    z$overwrite <- TRUE
    y <- do.call(mosaic, z)
    plot(y)
    print("Mosaic successful")
    
    #Deletes intermediate files
    #unlink(VBF_folder, recursive = TRUE)
  }
  print ("Valley Bottom Flatness generation complete")
  
} else {
  print("Valley Bottom Flatness not requested")
}


###################   TERRAIN RUGGEDNESS INDEX   ##########################


#If TRI = TRUE, performs Terrain Ruggedness Index
if (TRI == TRUE){
  print("Performing Terrain Ruggedness Index...")
  
  #Creates a folder to store results
  folder <- print(paste0(results_folder,"/TRI"))
  dir.create(folder)
  print("TRI folder added to results")
  
  #Obtains specific names
  for (d in dem_list){
    dem_name<- file_path_sans_ext(d)
    print(paste0("Working with DEM: ", dem_name))
    DEM <- paste0(DEMS_dir, "/", dem_name)
    dem_type <- paste0("_", dem_name)
    
    #Creates a folder to store results
    TRI_folder <- print(paste0(folder,"/TRI_tmp", dem_type))
    dir.create(TRI_folder)
    
    #Loops through DEM folder
    TRI_list <- list.files(dem_crops, pattern = dem_type)
    for (i in TRI_list){
      DEM_tile <- print(paste0(dem_crops, "/", i))
      outname <- file_path_sans_ext(i)
      setwd(TRI_folder)
      GRID_files <- print(paste0(outname, "_tmp"))
      TRI_output <- print(paste0(GRID_files, ".sgrd"))
      
      #Generates TRI as SAGA Grid
      print("Generating Terrain Ruggedness Index...")
      DEM <- paste0(dem_crops, "/", i) 
      TRI <- paste0(TRI_folder, "/", TRI_output)
      rsaga.geoprocessor("ta_morphometry", 16, list(DEM=DEM, TRI=TRI, RADIUS="4"), env=myenv)
      print("TRI created...")
    }
    
    #Converts SAGA grids to TIFF format
    tif_list <- list.files(TRI_folder, full.names = FALSE, pattern = ".sdat$")
    for (i in tif_list){ 
      sdat <- paste0(TRI_folder, "/", i)
      TRI_tif <- file_path_sans_ext(i)
      setwd(TRI_folder)
      s <- raster(sdat)
      filename <- paste0(TRI_folder, "/", TRI_tif, ".tif")
      writeRaster(s, filename=filename, format="GTiff", overwrite=TRUE)
      print(paste0(TRI_tif, "saved as tif"))
      
      #Deletes intermediate files
      SAGA_files <- paste0(TRI_folder, "/", TRI_tif)
      mgrd <- paste0(SAGA_files, ".mgrd")
      prj <- paste0(SAGA_files, ".prj")
      xml <- paste0(SAGA_files, ".sdat.aux.xml")
      sgrd <- paste0(SAGA_files, ".sgrd")
      file.remove(sdat)
      file.remove(mgrd)
      file.remove(prj)
      file.remove(xml)
      file.remove(sgrd)
      print("Intermediate files deleted")
      setwd(workdir)
      gc()
    }
    
    #Creates a vector of all the .tif files required for the mosaic
    TRIlist <- list.files(TRI_folder, pattern = ".tif$")
    w <- list()
    for (i in TRIlist){
      t_file <- paste0(TRI_folder, "/", i)
      t <- raster(t_file)
      outname <- file_path_sans_ext(i)
      assign(paste0("r", outname), t)
      w[[which(TRIlist==i)]] <- t
    }
    z <- w[!sapply(w, is.null)]
    
    #Mosaics the rasters within the folder
    filename <- paste0(folder, "/TRI", dem_type, ".tif")
    names(z)[1:2] <- c('x', 'y')
    z$fun <- mean
    z$na.rm <- TRUE
    z$filename <- filename
    z$overwrite <- TRUE
    y <- do.call(mosaic, z)
    plot(y)
    print("Mosaic successful")
    
    #Deletes intermediate files
    #unlink(TRI_folder, recursive = TRUE)
  }
  print ("Terrain Ruggedness Index generation complete")
  
} else {
  print("Terrain Ruggedness Index not requested")
}




########################   CONCAVITY - 10   #############################


#If C10 = TRUE, performs concavity 10
if (C10 == TRUE){
  print("Performing Curvature 10...")
  
  #Creates a folder to store results
  folder <- print(paste0(results_folder,"/C10"))
  dir.create(folder)
  print("C10 folder added to Results")
  
  #Obtains specific names
  for (d in dem_list){
    dem_name<- file_path_sans_ext(d)
    print(paste0("Working with DEM: ", dem_name))
    DEM <- paste0(DEMS_dir, "/", dem_name)
    dem_type <- paste0("_", dem_name)
  
    #Creates a folder to store intermediate data
    C10_folder <- print(paste0(folder, "/C10_tmp", dem_type))
    dir.create(C10_folder)
    print("C10 folder added to Results")

    #Loops through DEM folder
    C10_list <- list.files(dem_crops, pattern = dem_type)
    for (i in C10_list){
      DEM_tile <- print(paste0(dem_crops, "/", i))
      outname <- file_path_sans_ext(i)
      setwd(C10_folder)
      GRID_files <- print(paste0(outname, "_tmp"))
      C10_output <- print(paste0(GRID_files, ".sgrd"))
      
      #Generates C10 as SAGA Grid
      print("Generating Curvature 10...")
      DEM <- paste0(dem_crops, "/", i) 
      CURVE <- paste0(C10_folder, "/", C10_output)
      rsaga.geoprocessor("ta_morphometry", 21, list(DEM=DEM, CONVEXITY = CURVE, TYPE="1", METHOD="0", SCALE= "10"), env=myenv)
      print("Creating Curvature 10...")
    }
    
    #Converts SAGA grids to TIFF format
    tif_list <- list.files(C10_folder, full.names = FALSE, pattern = ".sdat$")
    for (i in tif_list){ 
      sdat <- paste0(C10_folder, "/", i)
      C10_tif <- file_path_sans_ext(i)
      setwd(C10_folder)
      s <- raster(sdat)
      filename <- paste0(C10_folder, "/", C10_tif, ".tif")
      writeRaster(s, filename=filename, format="GTiff", overwrite=TRUE)
      print(paste0(C10_tif, "saved as tif"))
      
      #Deletes intermediate files
      SAGA_files <- paste0(C10_folder, "/", C10_tif)
      mgrd <- paste0(SAGA_files, ".mgrd")
      prj <- paste0(SAGA_files, ".prj")
      xml <- paste0(SAGA_files, ".sdat.aux.xml")
      sgrd <- paste0(SAGA_files, ".sgrd")
      file.remove(sdat)
      file.remove(mgrd)
      file.remove(prj)
      file.remove(xml)
      file.remove(sgrd)
      print("Intermediate files deleted")
      setwd(workdir)
      gc()
    }

    #Creates a vector of all the .tif files required for the mosaic
    C10list <- list.files(C10_folder, pattern = ".tif$")
    w <- list()
    for (i in C10list){
      t_file <- paste0(C10_folder, "/", i)
      t <- raster(t_file)
      outname <- file_path_sans_ext(i)
      assign(paste0("r", outname), t)
      w[[which(C10list==i)]] <- t
    }
    z <- w[!sapply(w, is.null)]
    
    #Mosaics the rasters within the folder
    filename <- paste0(folder, "/C10", dem_type, ".tif")
    names(z)[1:2] <- c('x', 'y')
    z$fun <- mean
    z$na.rm <- TRUE
    z$filename <- filename
    z$overwrite <- TRUE
    y <- do.call(mosaic, z)
    plot(y)
    print("Mosaic successful")
    
    #Deletes intermediate files
    #unlink(C10_folder, recursive = TRUE)
  }
  print ("Curvature 10 generation complete")
} else {
  print("Curvature 10 not requested")
}


########################   CONCAVITY - 15   #############################


#If C15 = TRUE, performs concavity 15
if (C15 == TRUE){
  print("Performing Curvature 15...")
  
  #Creates a folder to store results
  folder <- print(paste0(results_folder,"/C15"))
  dir.create(folder)
  print("C15 folder added to Results")
  
  #Obtains specific names
  for (d in dem_list){
    dem_name<- file_path_sans_ext(d)
    print(paste0("Working with DEM: ", dem_name))
    DEM <- paste0(DEMS_dir, "/", dem_name)
    dem_type <- paste0("_", dem_name)
    
    #Creates a folder to store intermediate data
    C15_folder <- print(paste0(folder, "/C15_tmp", dem_type))
    dir.create(C15_folder)
    print("C15 folder added to Results")
    
    #Loops through DEM folder
    C15_list <- list.files(dem_crops, pattern = dem_type)
    for (i in C15_list){
      DEM_tile <- print(paste0(dem_crops, "/", i))
      outname <- file_path_sans_ext(i)
      setwd(C15_folder)
      GRID_files <- print(paste0(outname, "_tmp"))
      C15_output <- print(paste0(GRID_files, ".sgrd"))
      
      #Generates C15 as SAGA Grid
      print("Generating Curvature 15...")
      DEM <- paste0(dem_crops, "/", i) 
      CURVE <- paste0(C15_folder, "/", C15_output)
      rsaga.geoprocessor("ta_morphometry", 21, list(DEM=DEM, CONVEXITY = CURVE, TYPE="1", METHOD="0", SCALE= "15"), env=myenv)
      print("Creating Curvature 15...")
    }
    
    #Converts SAGA grids to TIFF format
    tif_list <- list.files(C15_folder, full.names = FALSE, pattern = ".sdat$")
    for (i in tif_list){ 
      sdat <- paste0(C15_folder, "/", i)
      C15_tif <- file_path_sans_ext(i)
      setwd(C15_folder)
      s <- raster(sdat)
      filename <- paste0(C15_folder, "/", C15_tif, ".tif")
      writeRaster(s, filename=filename, format="GTiff", overwrite=TRUE)
      print(paste0(C15_tif, "saved as tif"))
      
      #Deletes intermediate files
      SAGA_files <- paste0(C15_folder, "/", C15_tif)
      mgrd <- paste0(SAGA_files, ".mgrd")
      prj <- paste0(SAGA_files, ".prj")
      xml <- paste0(SAGA_files, ".sdat.aux.xml")
      sgrd <- paste0(SAGA_files, ".sgrd")
      file.remove(sdat)
      file.remove(mgrd)
      file.remove(prj)
      file.remove(xml)
      file.remove(sgrd)
      print("Intermediate files deleted")
      setwd(workdir)
      gc()
    }
    
    #Creates a vector of all the .tif files required for the mosaic
    C15list <- list.files(C15_folder, pattern = ".tif$")
    w <- list()
    for (i in C15list){
      t_file <- paste0(C15_folder, "/", i)
      t <- raster(t_file)
      outname <- file_path_sans_ext(i)
      assign(paste0("r", outname), t)
      w[[which(C15list==i)]] <- t
    }
    z <- w[!sapply(w, is.null)]
    
    #Mosaics the rasters within the folder
    filename <- paste0(folder, "/C15", dem_type, ".tif")
    names(z)[1:2] <- c('x', 'y')
    z$fun <- mean
    z$na.rm <- TRUE
    z$filename <- filename
    z$overwrite <- TRUE
    y <- do.call(mosaic, z)
    plot(y)
    print("Mosaic successful")
    
    #Deletes intermediate files
    #unlink(C15_folder, recursive = TRUE)
  }
  print ("Curvature 15 generation complete")
} else {
  print("Curvature 15 not requested")
}

########################   CONCAVITY - 20   #############################


#If C20 = TRUE, performs concavity 20
if (C20 == TRUE){
  print("Performing Curvature 20...")
  
  #Creates a folder to store results
  folder <- print(paste0(results_folder,"/C20"))
  dir.create(folder)
  print("C20 folder added to Results")
  
  #Obtains specific names
  for (d in dem_list){
    dem_name<- file_path_sans_ext(d)
    print(paste0("Working with DEM: ", dem_name))
    DEM <- paste0(DEMS_dir, "/", dem_name)
    dem_type <- paste0("_", dem_name)
    
    #Creates a folder to store intermediate data
    C20_folder <- print(paste0(folder, "/C20_tmp", dem_type))
    dir.create(C20_folder)
    print("C20 folder added to Results")
    
    #Loops through DEM folder
    C20_list <- list.files(dem_crops, pattern = dem_type)
    for (i in C20_list){
      DEM_tile <- print(paste0(dem_crops, "/", i))
      outname <- file_path_sans_ext(i)
      setwd(C20_folder)
      GRID_files <- print(paste0(outname, "_tmp"))
      C20_output <- print(paste0(GRID_files, ".sgrd"))
      
      #Generates C20 as SAGA Grid
      print("Generating Curvature 20...")
      DEM <- paste0(dem_crops, "/", i) 
      CURVE <- paste0(C20_folder, "/", C20_output)
      rsaga.geoprocessor("ta_morphometry", 21, list(DEM=DEM, CONVEXITY = CURVE, TYPE="1", METHOD="0", SCALE= "20"), env=myenv)
      print("Creating Curvature 20...")
    }
    
    #Converts SAGA grids to TIFF format
    tif_list <- list.files(C20_folder, full.names = FALSE, pattern = ".sdat$")
    for (i in tif_list){ 
      sdat <- paste0(C20_folder, "/", i)
      C20_tif <- file_path_sans_ext(i)
      setwd(C20_folder)
      s <- raster(sdat)
      filename <- paste0(C20_folder, "/", C20_tif, ".tif")
      writeRaster(s, filename=filename, format="GTiff", overwrite=TRUE)
      print(paste0(C20_tif, "saved as tif"))
      
      #Deletes intermediate files
      SAGA_files <- paste0(C20_folder, "/", C20_tif)
      mgrd <- paste0(SAGA_files, ".mgrd")
      prj <- paste0(SAGA_files, ".prj")
      xml <- paste0(SAGA_files, ".sdat.aux.xml")
      sgrd <- paste0(SAGA_files, ".sgrd")
      file.remove(sdat)
      file.remove(mgrd)
      file.remove(prj)
      file.remove(xml)
      file.remove(sgrd)
      print("Intermediate files deleted")
      setwd(workdir)
      gc()
    }
    
    #Creates a vector of all the .tif files required for the mosaic
    C20list <- list.files(C20_folder, pattern = ".tif$")
    w <- list()
    for (i in C20list){
      t_file <- paste0(C20_folder, "/", i)
      t <- raster(t_file)
      outname <- file_path_sans_ext(i)
      assign(paste0("r", outname), t)
      w[[which(C20list==i)]] <- t
    }
    z <- w[!sapply(w, is.null)]
    
    #Mosaics the rasters within the folder
    filename <- paste0(folder, "/C20", dem_type, ".tif")
    names(z)[1:2] <- c('x', 'y')
    z$fun <- mean
    z$na.rm <- TRUE
    z$filename <- filename
    z$overwrite <- TRUE
    y <- do.call(mosaic, z)
    plot(y)
    print("Mosaic successful")
    
    #Deletes intermediate files
    #unlink(C20_folder, recursive = TRUE)
  }
  print ("Curvature 20 generation complete")
} else {
  print("Curvature 20 not requested")
}

print("Started July 26 2018 at 9:30am")
print("Script ran successfully")