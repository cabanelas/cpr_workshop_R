################################################################################
#############        Continuous Plankton Recorder     ##########################
#############    CPR-BEAMS Workshop -- November 2025  ##########################
## by: Alexandra Cabanelas 
## created NOV-2025
################################################################################

#https://website.whoi.edu/cpr-beams/
#https://doi.mba.ac.uk/data/3567
#https://drive.google.com/drive/u/0/folders/1FopWXQNcl4reXZtT5QPUmLUlFO911373

# translating Pierre Helaouet MATLAB code to R 
# script #3 = STEP4 Select areas to extract CPR data
# translating Prog_CPRBeam_WS_Step4_TrimToArea.m to R 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(dplyr)
library(ggplot2)
library(sf)
library(sp)
library(RColorBrewer)
library(viridis)
library(purrr)

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
# --- Set WD
# --- Load Rdata file
load("CPR_Data_CPRBeam.RData")

## ------------------------------------------ ##
#         Example 1: simple rectangle -----
## ------------------------------------------ ##
##        Extract data in specific area  ----

# --- Extract spatio-temporal coordinates ---
T_Coord <- Data_LargeZooplankton %>% 
  select(1:8) #metadata columns

# --- Define spatial and temporal filter ---
Logi_tot <- with(T_Coord,
                 Year      >= 1958 & Year      <= 2022 &
                 Longitude >= -76  & Longitude <= -61  &
                 Latitude  >= 36.9 & Latitude  <= 46.8)

# --- Create subset to check if area is correct ---
Test_Area <- T_Coord[Logi_tot, ]
Test_Area$YearCat <- factor(Test_Area$Year)

# --- Crude map ---
ggplot(Test_Area, aes(x = Longitude, y = Latitude, color = YearCat)) +
  borders("world", color = "gray80", fill = "gray95") +
  geom_point(alpha = 0.4, size = 2) +
  coord_cartesian(xlim = c(-76, -61), ylim = c(36, 46)) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude", color = "Year")

##        Efficiently re-organize dataset  ----
#Area_1_Ex1 <- list(
  # --- Spatio-temporal metadata ---
#  SpatioTemp = Data_LargeZooplankton[Logi_tot, 1:8],
  # --- Reorganize abundance data ---
#  Data = list(
#    LargeZoo = Data_LargeZooplankton[Logi_tot, 9:ncol(Data_LargeZooplankton)],
#    SmallZoo = Data_SmallZooplankton[Logi_tot, 9:ncol(Data_SmallZooplankton)],
#    Phyto    = Data_Phytoplankton[Logi_tot, 9:ncol(Data_Phytoplankton)]
#  ),
  # --- Reorganize reference lists ---
#  List = list(
#    LargeZoo = List_LargeZooplankton,
#    SmallZoo = List_SmallZooplankton,
#    Phyto    = List_Phytoplankton
#  ),
  # --- Reorganize taxonomic metadata ---
#  Taxo = list(
#    LargeZoo = Taxo_LargeZooplankton,
#    SmallZoo = Taxo_SmallZooplankton,
#    Phyto    = Taxo_Phytoplankton
#  )
#)

# --- Clearvars ---
rm(T_Coord, Logi_tot, Test_Area) 

## ------------------------------------------ ##
#    Example 2: Google Earth Polygon -----
## ------------------------------------------ ##
##        Extract data in specific area  ----

# Go to folder where *.kml file is stored
# --- Find *.kml file ---
S_Kml <- list.files(pattern = "\\.kml$", full.names = TRUE)[1]#first .kml file
# --- Extract *.kml file as structure ---
S_Area <- sf::st_read(S_Kml, quiet = TRUE) #sf polygon
# --- Extract coordinates from any data file ---
T_Coord <- Data_LargeZooplankton[, 1:8]

# --- Convert sampling points to spatial object ---
Points_sf <- T_Coord %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# --- Find all samples in polygon extracted from *.kml ---
sf::st_crs(Points_sf) == sf::st_crs(S_Area)
#Logi_Space <- st_within(Points_sf, S_Area, sparse = FALSE)[, 1]#this works for kml with only 1 polygon
Logi_Space <- lengths(sf::st_within(Points_sf, S_Area)) > 0

# --- Create logical for selected time period ---
Logi_Time <- T_Coord$Year >= 1958 & T_Coord$Year <= 2022

# --- Combine logicals ---
Logi_tot <- Logi_Time & Logi_Space

##        Efficiently re-organize dataset  ----
Area_1_Ex2 <- list(
  # Reorganize Spatio-temporal data ---
  SpatioTemp = T_Coord[Logi_tot, 1:8],
  
  # --- Reorganize abundance data ---
  Data = list(
    LargeZoo = Data_LargeZooplankton[Logi_tot, 9:ncol(Data_LargeZooplankton)],
    SmallZoo = Data_SmallZooplankton[Logi_tot, 9:ncol(Data_SmallZooplankton)],
    Phyto    = Data_Phytoplankton[Logi_tot, 9:ncol(Data_Phytoplankton)]
  ),
  
  # --- Reorganize reference lists ---
  List = list(
    LargeZoo = List_LargeZooplankton,
    SmallZoo = List_SmallZooplankton,
    Phyto    = List_Phytoplankton
  ),
  
  # --- Reorganize taxonomic metadata ---
  Taxo = list(
    LargeZoo = Taxo_LargeZooplankton,
    SmallZoo = Taxo_SmallZooplankton,
    Phyto    = Taxo_Phytoplankton
  )
)

# --- Crude map ---
Area_1_Ex2$SpatioTemp$YearCat <- as.factor(Area_1_Ex2$SpatioTemp$Year)

ggplot(Area_1_Ex2$SpatioTemp, 
       aes(x = Longitude, y = Latitude, color = YearCat)) +
  borders("world", color = "gray80", fill = "gray95") +
  geom_point(alpha = 0.4, size = 2) +
  coord_cartesian(xlim = c(-76, -61), ylim = c(32, 46)) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude", color = "Year")

# --- Clearvars ---
rm(S_Kml, S_Area, T_Coord, Points_sf, Logi_Space, Logi_Time, Logi_tot)

## ------------------------------------------ ##
#    Example 3: Several Polygons -----
## ------------------------------------------ ##
##        Extract data in specific area  ----

# --- Load shapefiles containing areas ---
S_Folder <- list.files(path = "raw/CPRBeam_ICES_areas", 
                       pattern = "\\.shp$", full.names = TRUE)[1]
S_ICES_Area_All <- st_read(S_Folder, quiet = TRUE)

##        Only selected areas  ----
# Example: Major_FA = 27, SubArea = 7, Division = e;f;g;h;j
# 1) Extract all info contained in 'Major_FA'
S_ICES_Area <- S_ICES_Area_All %>%
  dplyr::filter(
    # 2) "Major_FA" = 27 
    Major_FA == "27",
    # 3) "SubArea" = 7 
    SubArea  == "7",
    # 4) "Division" in {'e','f','g','h','j'}
    Division %in% c("e", "f", "g", "h", "j")
  )

##        Un-project data  ----
# info about projection - embedded in the sf object
st_crs(S_ICES_Area) #st_crs(S_ICES_Area)$proj4string  

# --- Transform to WGS84 (EPSG:4326) ---
S_ICES_Area_proj <- st_transform(S_ICES_Area, crs = 4326)

##        Find which sample in which area ----
# --- Extract coordinates from any data file
T_Coord <- Data_LargeZooplankton[, 1:8]

# --- Convert sampling points to sf object  ---
Points_sf <- T_Coord %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# Logical matrix: rows = samples, cols = selected ICES areas
# MATLAB: Logi_Space(:,ii) = inpolygon(...)
Logi_Space <- sf::st_within(Points_sf, S_ICES_Area_proj, sparse = FALSE)

##        Efficiently re-organize dataset ----
# --- Create logical for selected time period ---
Logi_Time <- T_Coord$Year >= 1958 & T_Coord$Year <= 2022

# Initialize list (MATLAB: Area_Ex3(ii))
Area_Ex3 <- vector("list", length = nrow(S_ICES_Area_proj))

# --- Loop through all selected areas --- 
for (ii in seq_len(nrow(S_ICES_Area_proj))) {
  
  # --- Combined logical ---
  # MATLAB: and(Logi_Time, Logi_Space(:,ii))
  Logi_tot <- Logi_Time & Logi_Space[, ii]
  
  # --- Reorganize data for current area ---
  Area_Ex3[[ii]] <- list(
    # --- Add area name for each dimension (to check) ---
    Name = S_ICES_Area_proj$Area_Full[ii],
    
    # --- Reorganize Spatio-temporal data ---
    SpatioTemp = T_Coord[Logi_tot, ],
    
    # --- Abundance data ---
    Data = list(
      LargeZoo = Data_LargeZooplankton[Logi_tot, 9:ncol(Data_LargeZooplankton)],
      SmallZoo = Data_SmallZooplankton[Logi_tot, 9:ncol(Data_SmallZooplankton)],
      Phyto    = Data_Phytoplankton[Logi_tot, 9:ncol(Data_Phytoplankton)]
    ),
    
    # --- List (same for all areas) ---
    List = list(
      LargeZoo = List_LargeZooplankton,
      SmallZoo = List_SmallZooplankton,
      Phyto    = List_Phytoplankton
    ),
    
    # --- Taxo (same for all areas) ---
    Taxo = list(
      LargeZoo = Taxo_LargeZooplankton,
      SmallZoo = Taxo_SmallZooplankton,
      Phyto    = Taxo_Phytoplankton
    )
  )
}

# --- Assign area names to list elements ---
names(Area_Ex3) <- S_ICES_Area_proj$Area_Full

##        Check map ----

# --- Define a unique color for each area ---
colormap_Spe <- viridis::viridis(length(Area_Ex3), option = "D")
#colormap_Spe <- RColorBrewer::brewer.pal(n = max(3, min(length(Area_Ex3), 9)), name = "Set1")
#colormap_Spe <- rep(colormap_Spe, length.out = length(Area_Ex3))
#colormap_Spe <- scales::hue_pal()(length(Area_Ex3))

# --- Combine all sample points into one data frame ---
# Add area name and color to each sample
Sample_All <- purrr::map2_dfr(Area_Ex3, seq_along(Area_Ex3), function(area, ii) {
  area$SpatioTemp %>%
    dplyr::mutate(
      AreaName = area$Name,
      Color = colormap_Spe[ii]
    )
})

# --- Convert to sf object for plotting ---
Sample_sf <- sf::st_as_sf(Sample_All, coords = c("Longitude", "Latitude"), crs = 4326)

# --- PLOT 1 ---
ggplot() +
  # ICES polygons
  geom_sf(data = S_ICES_Area_proj, aes(fill = Area_Full), 
          alpha = 0.25, color = "black") +
  # sample points
  geom_sf(data = Sample_sf, aes(color = AreaName), size = 2, 
          alpha = 0.5, shape = 21, stroke = 0.3) +
  borders("world", color = "gray80", fill = "gray70") +
  coord_sf(xlim = c(-12, -1), ylim = c(48, 53), expand = FALSE) +
  scale_color_manual(values = colormap_Spe) +
  scale_fill_manual(values = colormap_Spe) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude",
       color = "ICES Area", fill = "ICES Area")

# --- PLOT 2 ---
ggplot() +
  borders("world", fill = "grey80", colour = "grey70") +
  geom_sf(data  = S_ICES_Area_proj, aes(fill = Area_Full),
          alpha = 0.25, color = NA) +
  geom_sf(data  = Sample_sf, aes(color = AreaName),
          size  = 2, alpha = 0.5) +
  coord_sf(xlim = c(-12, -1), ylim = c(48, 53), expand = FALSE) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude", 
       fill = "ICES areas", color = "ICES areas")