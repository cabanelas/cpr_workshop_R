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
V_Interp <- as.matrix(T_Extract[, -1])  

# --- Pre-dimension output arrays ---
dims <- c(length(V_Lat), length(V_Lon), 12, length(V_YearWindow))
M_Abun <- array(NA_real_, dim = dims)
M_SampEff <- array(NA_real_, dim = dims)

# --- Create coordinate matrices ---
M_Lat <- matrix(rep(V_Lat, each = length(V_Lon)), nrow = length(V_Lat), ncol = length(V_Lon))
M_Lon <- matrix(rep(V_Lon, times = length(V_Lat)), nrow = length(V_Lat), ncol = length(V_Lon), byrow = TRUE)

# --- Loop over years, months, and spatial bins ---
for (ii in seq_along(V_YearWindow)) {
  # Subset data for current year
  temp <- V_Interp[V_Interp[, 3] == V_YearWindow[ii], ]
  
  # Loop over each month (1-12)
  for (jj in 1:12) {
    # subset for the current month
    temp2 <- temp[temp[, 4] == jj, ]
    
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
#   Graph  -----
## ------------------------------------------ ##

# NEED TO ADD THIS PART
# MATLAB LINE 119 
