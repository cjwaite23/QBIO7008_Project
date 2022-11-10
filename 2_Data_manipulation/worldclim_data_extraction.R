##### Description ######
# Extracting Climate Data for Trinidad rivers
# Data files are too large so are stored on my personal hard drive

##### Best guess of coordinates:
# Taylor: 10.7087479, -61.2746940
#         10°42'31.4928''N 61°16'28.8984''W
# Caigual: 10.7154825, -61.2745807
#         10°42'55.7388''N 61°16'28.4916''W
# Upper La Laja: 10.7193228, -61.2667921
#         10°43'9.5628''N 61°16'0.4512''W
# Lower La Laja: 10.7125326, -61.2635067
#         10°42'45.1188''N 61°15'48.6252''W

##### Libraries #####
library(raster)
library(tidyverse)

# Create RasterStack objects
prec_files <- list.files("C:/Users/callu/Documents/QBIO7008 Guppy Project/7008_large_files/1_Data/worldclim/precipitation00-18/", ".tif", full.names=TRUE)
prec <- stack(prec_files)
tmin_files <- list.files("C:/Users/callu/Documents/QBIO7008 Guppy Project/7008_large_files/1_Data/worldclim/tempmin00-18/", ".tif", full.names=TRUE)
tmin <- stack(tmin_files)
tmax_files <- list.files("C:/Users/callu/Documents/QBIO7008 Guppy Project/7008_large_files/1_Data/worldclim/tempmax00-18/", ".tif", full.names=TRUE)
tmax <- stack(tmax_files)

months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
date_names <- paste(months, rep(as.character(2000:2018), each = 12), sep = "_")
names(prec) <- date_names
names(tmin) <- date_names
names(tmax) <- date_names

# Create a data.frame with sample site coordinates
site <- c("Trinidad")
lon <- c(-61.27)
lat <- c(10.71)
samples <- data.frame(lon, lat)
# Extract data from RasterLayer
prec_data <- extract(prec, samples)
tmin_data <- extract(tmin, samples)
tmax_data <- extract(tmax, samples)

clim_data <- rbind(prec_data, tmin_data, tmax_data)
row.names(clim_data) <- c("prec", "tmin", "tmax")
save(clim_data, file = "2_Data_manipulation/data_files/climate_data.RData")
