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
# script #7 = STEP8 regularise in space
# translating Prog_CPRBeam_WS_Step8_RegulInSpace.m to R 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)             
library(rnaturalearth)  
library(scales)

# helper ported from f_CPRBeam_pos2dist.m (used by the IDW section)
source("f_CPRBeam_pos2dist.R")

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
# --- Set WD
# --- Load Rdata file
load("CPR_Data_CPRBeam.RData") # you could also do load("CPR_Data_CPRBeam.mat")

## ------------------------------------------ ##
#      Extract data in specific area  -----
## ------------------------------------------ ##
##     Example 1: Simple Rectangle     ----

# --- Extract coordinates from any data file ---
T_Coord <- Data_LargeZooplankton %>% 
  select(1:8) #metadata columns

# --- Define spatial and temporal filter ---
Logi_tot <- with(T_Coord,
                 Year      >= 1958 & Year      <= 2022 &
                   Longitude >= -4.5 & Longitude <= 10   &
                   Latitude  >= 50   & Latitude  <= 61)

## ------------------------------------------ ##
#     Efficiently reorganize dataset  -----
## ------------------------------------------ ##
Area_2 <- list(
  # --- Spatio-temporal metadata ---
  SpatioTemp = T_Coord[Logi_tot, ],
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

# --- Clearvars ---
rm(T_Coord, Logi_tot)

## ------------------------------------------ ##
#        Extract CFIN data  -----
## ------------------------------------------ ##
# --- Identify CFIN by ID = 40 ---
Logi_TaxaID <- Area_2$List$LargeZoo$Accepted_ID == 40

# --- inspect matching taxa ---
Area_2$List$LargeZoo[Logi_TaxaID, 1:2]

# --- Extract CFIN data column ---
T_CFIN <- Area_2$Data$LargeZoo[, Logi_TaxaID, drop = FALSE]
colnames(T_CFIN) <- "CFIN"

# --- Create a table with spatio-temporal coordinates and CFIN ---
T_Extract <- cbind(Area_2$SpatioTemp, T_CFIN)

## ------------------------------------------ ##
#        Calculate CFIN timeseries  -----
## ------------------------------------------ ##
T_Stats <- T_Extract %>%
  group_by(Year, Month) %>%
  summarise(
    mean_CFIN = mean(CFIN, na.rm = TRUE),
    .groups = "drop") %>%
  complete(
    Year = 1958:2022,
    Month = 1:12) %>%
  arrange(Year, Month)

## ------------------------------------------ ##
#      Visualize timeseries  -----
## ------------------------------------------ ##
ggplot(T_Stats, aes(x = Year, y = factor(Month), fill = mean_CFIN)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "D", na.value = "white", name = "Mean CFIN") +
  scale_x_continuous(
    breaks = seq(1958, 2022, by = 2),
    expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(
    limits = as.character(12:1),
    labels = rev(month.abb),
    expand = expansion(mult = c(0, 0))) +
  labs(title = "Monthly Mean Abundance of Calanus finmarchicus",
       x = "Year", y = "Month") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   color = "black", size = 11))

## ------------------------------------------ ##
#   Regularize in space with simple binning  -----
## ------------------------------------------ ##

# --- Set-up regularisation ---
V_YearWindow <- 1958:2022
V_Step <- 0.5
V_Lat <- seq(50, 61, by = V_Step)
V_Lon <- seq(-4, 10, by = V_Step)

# --- Create input variable ---
# T_Extract cols after dropping Sample: 1=Lat 2=Lon 3=Year 4=Month ... 8=CFIN
V_Interp <- as.matrix(T_Extract[, -1])  

# --- Pre-dimension output arrays ---
dims <- c(length(V_Lat), length(V_Lon), 12, length(V_YearWindow))
M_Abun <- array(NA_real_, dim = dims)
M_SampEff <- array(NA_real_, dim = dims)

# --- Create coordinate matrices ---
# M_Lat <- matrix(rep(V_Lat, each = length(V_Lon)), nrow = length(V_Lat), ncol = length(V_Lon))
# M_Lon <- matrix(rep(V_Lon, times = length(V_Lat)), nrow = length(V_Lat), ncol = length(V_Lon), byrow = TRUE)

M_Lat <- matrix(rep(V_Lat, times = length(V_Lon)),
                nrow = length(V_Lat), ncol = length(V_Lon))   # cols = V_Lat
M_Lon <- matrix(rep(V_Lon, each  = length(V_Lat)),
                nrow = length(V_Lat), ncol = length(V_Lon))   # rows = V_Lon

# --- Loop over years, months, and spatial bins ---
for (ii in seq_along(V_YearWindow)) {
  # Subset data for current year
  temp <- V_Interp[V_Interp[, 3] == V_YearWindow[ii], , drop = FALSE]
  
  # Loop over each month (1-12)
  for (jj in 1:12) {
    # subset for the current month
    temp2 <- temp[temp[, 4] == jj, , drop = FALSE]
    
    # Loop over each longitude bin
    for (kk in seq_along(V_Lon)) {
      # Loop over each latitude bin
      for (ll in seq_along(V_Lat)) {
        # Define spatial bin boundaries
        lon_min <- V_Lon[kk] - V_Step / 2
        lon_max <- V_Lon[kk] + V_Step / 2
        lat_min <- V_Lat[ll] - V_Step / 2
        lat_max <- V_Lat[ll] + V_Step / 2
        
        f <- which(
          temp2[, 2] >= lon_min & temp2[, 2] < lon_max &
            temp2[, 1] >= lat_min & temp2[, 1] < lat_max
        )
        
        if (length(f) > 0) {
          # Mean abundance (column 8) for this bin, month, and year
          M_Abun[ll, kk, jj, ii] <- mean(temp2[f, 8], na.rm = TRUE)
          # Sampling effort: number of observations in the bin
          M_SampEff[ll, kk, jj, ii] <- length(f)
        }
      }
    }
  }
}

# --- Store in a list ---
#similar to MATLAB struct
S_Grid <- list(
  M_Abun = M_Abun,
  M_SampEff = M_SampEff,
  M_Lat = M_Lat,
  M_Lon = M_Lon
)

## ------------------------------------------ ##
#   Graph (simple binning)  -----
## ------------------------------------------ ##

# Basemap (land polygons) for the North Sea region
world <- ne_countries(scale = "medium", returnclass = "sf")

# --- Collapse 4D array [Lat, Lon, Month, Year] down to 2D [Lat, Lon] ---
# MATLAB: mean over months (dim 3), then mean over years.
# apply() keeps the listed margins and averages across the rest.
M_Map3D_Abun <- apply(S_Grid$M_Abun, c(1, 2, 4), mean, na.rm = TRUE) # -> [Lat, Lon, Year]
M_Map2D_Abun <- apply(M_Map3D_Abun,  c(1, 2),    mean, na.rm = TRUE) # -> [Lat, Lon]

M_Map3D_SampEff <- apply(S_Grid$M_SampEff, c(1, 2, 4), mean, na.rm = TRUE)
M_Map2D_SampEff <- apply(M_Map3D_SampEff,  c(1, 2),    mean, na.rm = TRUE)
# could do sum instead^

# apply() returns NaN for all-NA cells; make them NA so they drop out of plots
M_Map2D_Abun[is.nan(M_Map2D_Abun)]       <- NA
M_Map2D_SampEff[is.nan(M_Map2D_SampEff)] <- NA

# --- Long data frame of grid nodes for plotting ---
# expand.grid(Lat, Lon) varies Lat fastest, matching column-major as.vector()
grid_df <- expand.grid(Lat = V_Lat, Lon = V_Lon)
grid_df$Abun    <- as.vector(M_Map2D_Abun)
grid_df$SampEff <- as.vector(M_Map2D_SampEff)

# Common map scaffolding
base_map <- function(title) {
  list(
    geom_sf(data = world, fill = "grey50", color = NA),
    coord_sf(xlim = c(-4.5, 10), ylim = c(50, 61), expand = FALSE),
    labs(title = title, x = NULL, y = NULL),
    theme_minimal(base_size = 14)
  )
}

# 1) Map the original (raw) data
ggplot() +
  geom_point(data = T_Extract,
             aes(Longitude, Latitude, color = log10(CFIN + 1)),
             size = 2, alpha = 0.5) +
  scale_color_viridis_c(limits = c(0, 2), oob = squish,
                        name = expression(log[10]*"(abundance + 1)")) +
  base_map("Original data")

# 2) Map the regularised (binned) data
ggplot() +
  geom_tile(data = grid_df,
            aes(Lon, Lat, fill = log10(Abun + 1)),
            width = V_Step, height = V_Step) +
  scale_fill_viridis_c(limits = c(0, 2), oob = squish, na.value = NA,
                       name = expression(log[10]*"(abundance + 1)")) +
  base_map("Regularised (0.5 deg binning)")

# 3) Map the number of samples per node
ggplot() +
  geom_tile(data = grid_df,
            aes(Lon, Lat, fill = SampEff),
            width = V_Step, height = V_Step) +
  scale_fill_viridis_c(limits = c(0, 50), oob = squish, na.value = NA,
                       name = "Samples per node") +
  base_map("Sampling effort per node")

## ------------------------------------------ ##
#   Regularise in space with IDW  -----
## ------------------------------------------ ##

# --- Set-up regularisation ---
V_YearWindow <- 2000:2022      # IDW uses the more recent, better-sampled period
V_Step <- 0.5
V_Lat <- seq(50, 61, by = V_Step)
V_Lon <- seq(-4, 10, by = V_Step)

# --- Input variable (same column layout as above) ---
V_Interp <- as.matrix(T_Extract[, -1])

# --- Grid-node coordinate matrices (same convention as the binning fix) ---
M_Lat <- matrix(rep(V_Lat, times = length(V_Lon)),
                nrow = length(V_Lat), ncol = length(V_Lon))
M_Lon <- matrix(rep(V_Lon, each  = length(V_Lat)),
                nrow = length(V_Lat), ncol = length(V_Lon))
node_lat <- as.vector(M_Lat)   # column-major, so node kk == matrix cell kk
node_lon <- as.vector(M_Lon)
n_nodes  <- length(node_lat)

# --- IDW settings ---
V_r     <- 50   # search radius in km
V_nSamp <- 0    # 0 = use ALL samples within the radius; >0 = nearest N only
power   <- 2    # inverse-distance power (d^-2)

# --- Pre-dimension output ---
M_Abun_IDW <- array(NA_real_,
                    dim = c(length(V_Lat), length(V_Lon), 12, length(V_YearWindow)))

# --- Loops ---
for (ii in seq_along(V_YearWindow)) {          # each year
  temp <- V_Interp[V_Interp[, 3] == V_YearWindow[ii], , drop = FALSE]
  
  for (jj in 1:12) {                           # each month
    temp2 <- temp[temp[, 4] == jj, , drop = FALSE]
    temp_M <- rep(NA_real_, n_nodes)
    
    if (nrow(temp2) > 0) {                      # skip empty year/month subsets
      samp_lat <- temp2[, 1]
      samp_lon <- temp2[, 2]
      samp_ab  <- temp2[, 8]
      
      for (kk in seq_len(n_nodes)) {           # each grid node
        # distance from this node to every sample (vectorised, method 2 = spherical)
        d <- f_CPRBeam_pos2dist(node_lat[kk], node_lon[kk],
                                samp_lat, samp_lon, method = 2)
        
        keep <- which(d < V_r & !is.na(samp_ab))   # inside radius & has a value
        if (length(keep) > 0) {
          dk <- d[keep]
          ak <- samp_ab[keep]
          dk[dk == 0] <- 1e-6                  # guard against a sample on the node
          
          if (V_nSamp > 0 && length(ak) >= V_nSamp) {
            o  <- order(dk)[1:V_nSamp]         # keep the V_nSamp nearest
            dk <- dk[o]
            ak <- ak[o]
          }
          
          w <- dk^(-power)                     # inverse-distance weights
          temp_M[kk] <- sum(ak * w) / sum(w)
        }
      }
    }
    # reshape node vector back to [Lat, Lon] (column-major, matches node ordering)
    M_Abun_IDW[, , jj, ii] <- matrix(temp_M, nrow = length(V_Lat), ncol = length(V_Lon))
  }
}

S_Interp <- list(M_Abun = M_Abun_IDW, M_Lat = M_Lat, M_Lon = M_Lon)

## ------------------------------------------ ##
#   Graph (IDW)  -----
## ------------------------------------------ ##

# Collapse 4D -> 2D exactly as before
M_Map3D_Abun <- apply(S_Interp$M_Abun, c(1, 2, 4), mean, na.rm = TRUE)
M_Map2D_Abun <- apply(M_Map3D_Abun,    c(1, 2),    mean, na.rm = TRUE)
M_Map2D_Abun[is.nan(M_Map2D_Abun)] <- NA

grid_idw <- expand.grid(Lat = V_Lat, Lon = V_Lon)
grid_idw$Abun <- as.vector(M_Map2D_Abun)

# 1) Original data (>= 2000 to match the IDW window)
ggplot() +
  geom_point(data = subset(T_Extract, Year >= 2000),
             aes(Longitude, Latitude, color = log10(CFIN + 1)),
             size = 2, alpha = 0.5) +
  scale_color_viridis_c(limits = c(0, 2), oob = squish,
                        name = expression(log[10]*"(abundance + 1)")) +
  base_map("Original data (2000-2022)")

# 2) IDW-regularised data
ggplot() +
  geom_tile(data = grid_idw,
            aes(Lon, Lat, fill = log10(Abun + 1)),
            width = V_Step, height = V_Step) +
  scale_fill_viridis_c(limits = c(0, 2), oob = squish, na.value = NA,
                       name = expression(log[10]*"(abundance + 1)")) +
  base_map("IDW interpolation (r = 50 km)")

## -----------------------------------------------------------------------------
## Performance note:
## The IDW block is a faithful translation of the MATLAB nested loops and can be slow. 
## could replace f_CPRBeam_pos2dist() with a fully vectorised call such as
## geosphere::distGeo(cbind(node_lon[kk], node_lat[kk]), cbind(samp_lon, samp_lat))
## which returns meters; divide by 1000 for km.
## -----------------------------------------------------------------------------