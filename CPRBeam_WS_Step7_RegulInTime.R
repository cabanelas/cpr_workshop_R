################################################################################
#############        Continuous Plankton Recorder     ##########################
#############    CPR-BEAMS Workshop -- November 2025  ##########################
## by: Alexandra Cabanelas 
## created NOV-2025
################################################################################

#https://website.whoi.edu/cpr-beams/
#https://drive.google.com/drive/u/0/folders/1FopWXQNcl4reXZtT5QPUmLUlFO911373

# translating Pierre Helaouet MATLAB code to R 
# script #6 = STEP7 regularise in time
# translating Prog_CPRBeam_WS_Step7_RegulInTime.m to R 

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
  theme_minimal(base_size = 12)

## ------------------------------------------ ##
#   Create timeseries with missing values  -----
## ------------------------------------------ ##
# --- Extract CFIN TS as vector ---
V_CFIN_Raw <- T_Stats$mean_CFIN

# --- Generate random indices for missing values (10% of total) ---
set.seed(123)
V_IdxMiss <- sample(seq_along(V_CFIN_Raw))  # same as MATLAB randperm

# --- Introduce missing values ---
V_CFIN_Missing <- V_CFIN_Raw
V_CFIN_Missing[V_IdxMiss[1:78]] <- NA  # 10% of missing values (= 78)

# --- Create time index vectors ---
V_X <- seq_along(V_CFIN_Raw)
V_X_Missing <- V_X
V_X_Missing[V_IdxMiss[1:78]] <- NA

rm(V_IdxMiss)

## ------------------------------------------ ##
#           Interpolation 1  -----
## ------------------------------------------ ##

V_CFIN_Interp <- V_CFIN_Missing
# Create interpolant using only non-missing points
missing_idx <- which(is.na(V_CFIN_Missing)) 

# --- Linear interpolation using approx() ---
V_CFIN_Interp[missing_idx] <- approx(
  x = V_X[!is.na(V_CFIN_Missing)],
  y = V_CFIN_Missing[!is.na(V_CFIN_Missing)],
  xout = V_X[missing_idx],
  rule = 2
)$y

# --- Calculate correlation ---
R_Corr <- cor(V_CFIN_Raw, V_CFIN_Interp, use = "complete.obs")

# --- Plot ---
# Create data frame for plotting
df_plot <- data.frame(
  Time = V_X,
  Original = V_CFIN_Raw,
  Interpolated = V_CFIN_Interp,
  Missing = V_CFIN_Missing
) 

df_long <- df_plot %>%
  pivot_longer(
    cols = c(Original, Interpolated, Missing),
    names_to = "Series",
    values_to = "Abundance")

series_colors <- c(
  Original     = "#0072CE",   
  Interpolated = "#A1132F",   
  Missing      = "#A8A8A8")

ggplot(df_long, aes(x = Time, y = Abundance, color = Series)) +
  geom_line(linewidth = 1.2, alpha = 0.6, na.rm = TRUE) +
  scale_color_manual(values = series_colors) +
  scale_x_continuous(
    breaks = seq(1, length(V_X), by = 60),
    labels = seq(1958, 2022, by = 5),
    expand = expansion(mult = c(0, 0))) +
  labs(title = paste0("TS interpolation / r² = ", round(R_Corr^2, 3)),
    x = "Year", y = "Abundance", color = NULL) +
  theme_minimal(base_size = 14)




