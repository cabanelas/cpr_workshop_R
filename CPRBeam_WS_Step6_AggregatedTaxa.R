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
# script #5 = STEP6 create timeseries
# translating Prog_CPRBeam_WS_Step6_AggregatedTaxa.m to R 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(dplyr)
library(ggplot2)
library(tidyr)

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
#          Quick map -----
## ------------------------------------------ ##
T_Map <- Area_2$SpatioTemp
T_Map$YearCat <- as.factor(T_Map$Year)

ggplot(T_Map, aes(x = Longitude, y = Latitude, color = YearCat)) +
  borders("world", colour = "gray80", fill = "gray70") +
  geom_point(alpha = 0.5, size = 2) +
  coord_quickmap(xlim = c(-10, 15), ylim = c(50, 61)) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude", color = "Year")

## ------------------------------------------ ##
# Create new group: Example 1: large copepods -----
## ------------------------------------------ ##

# --- Create logical vector for large copepods ---
Logi_LCop <- as.logical(Area_2$List$LargeZoo$Copepods)

# --- Subset data to only large copepods ---
V_LCop <- Area_2$Data$LargeZoo[, Logi_LCop]

# --- Compute row-wise mean abundance across large copepods ---
V_Mean_LCop <- rowMeans(V_LCop, na.rm = TRUE)

# --- Concatenate with spatio-temporal coordinates ---
Area_2$LCop$T_Data <- cbind(Area_2$SpatioTemp, Mean_LCop = V_Mean_LCop)

# --- Store taxa list for large copepods ---
Area_2$LCop$T_List <- Area_2$List$LargeZoo[Logi_LCop, ]

## ------------------------------------------ ##
#    Generate timeseries from new dataset -----
## ------------------------------------------ ##
# --- Summarize by Year and Month ---
Area_2$T_Stats_LCop <- Area_2$LCop$T_Data %>%
  group_by(Year, Month) %>%
  summarise(
    mean_mean = mean(Mean_LCop, na.rm = TRUE),
    .groups = "drop") %>%
  complete(Year = 1958:2022,
           Month = 1:12) %>%
  arrange(Year, Month)

## ------------------------------------------ ##
#        Large copepods timeseries -----
## ------------------------------------------ ##
# --- Create time index for plotting ---
Area_2$T_Stats_LCop$TimeIndex <- seq_len(nrow(Area_2$T_Stats_LCop))

# --- Plot time series ---
ggplot(Area_2$T_Stats_LCop, aes(x = TimeIndex, y = mean_mean)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  scale_x_continuous(
    breaks = seq(1, nrow(Area_2$T_Stats_LCop), by = 60),  
    labels = seq(1958, 2022, by = 5),
    expand = expansion(mult = c(0, 0))) +
  labs(title = "Monthly Abundances of Mean Large Copepods",
       x = "Year", y = "Abundance") +
  theme_minimal(base_size = 14)

## ------------------------------------------ ##
# Create new group: Example 2: all Calanidae -----
## ------------------------------------------ ##
# --- Create logical vector for large Calanidae ---
Logi_LCalano <- Area_2$Taxo$LargeZoo$Family == "Calanidae"
Logi_LCalano <- replace_na(Logi_LCalano, FALSE)

# --- Subset data to only large Calanidae ---
V_LCalano <- Area_2$Data$LargeZoo[, Logi_LCalano]

# --- Compute row-wise mean abundance across Calanidae taxa ---
V_Mean_LCalano <- rowMeans(V_LCalano, na.rm = TRUE)

# --- Concatenate with spatio-temporal coordinates ---
Area_2$LCalano$T_Data <- cbind(Area_2$SpatioTemp, Mean_LCalano = V_Mean_LCalano)

# --- Store taxa list for Calanidae ---
Area_2$LCalano$T_List <- Area_2$List$LargeZoo[Logi_LCalano, ]

## ------------------------------------------ ##
#    Generate timeseries from new dataset -----
## ------------------------------------------ ##
# --- Summarize by Year and Month ---
Area_2$T_Stats_LCalano <- Area_2$LCalano$T_Data %>%
  group_by(Year, Month) %>%
  summarise(
    mean_mean = mean(Mean_LCalano, na.rm = TRUE),
    .groups = "drop") %>%
  complete(Year = 1958:2022,
           Month = 1:12) %>%
  arrange(Year, Month)

## ------------------------------------------ ##
#        Large Calanidae timeseries -----
## ------------------------------------------ ##
# --- Create time index for plotting ---
Area_2$T_Stats_LCalano$TimeIndex <- seq_len(nrow(Area_2$T_Stats_LCalano))

# --- Plot time series ---
ggplot(Area_2$T_Stats_LCalano, aes(x = TimeIndex, y = mean_mean)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  scale_x_continuous(
    breaks = seq(1, nrow(Area_2$T_Stats_LCalano), by = 60),  
    labels = seq(1958, 2022, by = 5),
    expand = expansion(mult = c(0, 0))) +
  labs(title = "Monthly Abundances of Mean Large Calanidae",
       x = "Year", y = "Abundance") +
  theme_minimal(base_size = 14)