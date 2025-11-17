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
# script #6 = STEP7 regularise in time
# translating Prog_CPRBeam_WS_Step7_RegulInTime.m to R 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(dplyr)
library(ggplot2)
library(tidyr)
library(signal) # for interp with method = "pchip"
library(akima) 

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
# --- Linear Interpolation ---

# --- Duplicate time series with missing values ---
V_CFIN_Interp <- V_CFIN_Missing
# --- Identify missing indices ---
missing_idx <- which(is.na(V_CFIN_Missing)) 

# --- Linear interpolation ---
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
  Original     = "#A8A8A8",   
  Interpolated = "#A1132F",   
  Missing      = "#0072CE")

ggplot(df_long, aes(x = Time, y = Abundance, color = Series)) +
  geom_line(linewidth = 1.2, alpha = 0.6, na.rm = TRUE) +
  scale_color_manual(values = series_colors) +
  scale_x_continuous(
    breaks = seq(1, length(V_X), by = 60),
    labels = seq(1958, 2022, by = 5),
    expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(title = paste0("TS interpolation / r² = ", round(R_Corr^2, 3)),
    x = "Year", y = "Abundance", color = NULL) +
  theme_minimal(base_size = 14)

rm(V_CFIN_Interp, R_Corr, df_plot, df_long)

## ------------------------------------------ ##
#           Interpolation 2  -----
## ------------------------------------------ ##
# --- pchip ---
# shape-preserving piecewise cubic interpolation 

# --- Duplicate time series with missing values ---
V_CFIN_Interp <- V_CFIN_Missing
# --- Identify missing indices ---
missing_idx <- which(is.na(V_CFIN_Missing))

# --- Shape-preserving cubic interpolation ---
V_CFIN_Interp[missing_idx] <- interp1(
  x = V_X[!is.na(V_CFIN_Missing)],
  y = V_CFIN_Missing[!is.na(V_CFIN_Missing)],
  xi = V_X[missing_idx],
  method = "pchip"
)

# --- Calculate correlation ---
R_Corr <- cor(V_CFIN_Raw, V_CFIN_Interp, use = "complete.obs")

# --- Prepare data for plotting ---
df_plot <- data.frame(
  Time = V_X,
  Original = V_CFIN_Raw,
  Interpolated = V_CFIN_Interp,
  Missing = V_CFIN_Missing
)

df_long <- df_plot %>%
  pivot_longer(cols = c(Original, Interpolated, Missing),
               names_to = "Series",
               values_to = "Abundance")

# --- Plot ---
ggplot(df_long, aes(x = Time, y = Abundance, color = Series)) +
  geom_line(linewidth = 1.2, alpha = 0.6, na.rm = TRUE) +
  scale_color_manual(values = series_colors) +
  scale_x_continuous(
    breaks = seq(1, length(V_X), by = 60),
    labels = seq(1958, 2022, by = 5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(
    title = paste0("TS interpolation / r² = ", round(R_Corr^2, 3)),
    x = "Year", y = "Abundance", color = NULL
  ) +
  theme_minimal(base_size = 14)

rm(V_CFIN_Interp, R_Corr, df_plot, df_long)

## ------------------------------------------ ##
#           Interpolation 3  -----
## ------------------------------------------ ##
# --- makima ---
# piecewise function of polynomials with degree at most three 

# --- Duplicate time series with missing values ---
V_CFIN_Interp <- V_CFIN_Missing
# --- Identify missing indices ---
missing_idx <- which(is.na(V_CFIN_Missing))

# --- Akima interpolation (similar to MATLAB makima) ---
V_CFIN_Interp[missing_idx] <- aspline(
  x = V_X[!is.na(V_CFIN_Missing)],
  y = V_CFIN_Missing[!is.na(V_CFIN_Missing)],
  xout = V_X[missing_idx]
)$y

# --- Calculate correlation ---
R_Corr <- cor(V_CFIN_Raw, V_CFIN_Interp, use = "complete.obs")

# --- Prepare data for plotting ---
df_plot <- data.frame(
  Time = V_X,
  Original = V_CFIN_Raw,
  Interpolated = V_CFIN_Interp,
  Missing = V_CFIN_Missing
)

df_long <- df_plot %>%
  pivot_longer(cols = c(Original, Interpolated, Missing),
               names_to = "Series",
               values_to = "Abundance")

# --- Plot ---
ggplot(df_long, aes(x = Time, y = Abundance, color = Series)) +
  geom_line(linewidth = 1.2, alpha = 0.6, na.rm = TRUE) +
  scale_color_manual(values = series_colors) +
  scale_x_continuous(
    breaks = seq(1, length(V_X), by = 60),
    labels = seq(1958, 2022, by = 5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(
    title = paste0("TS interpolation / r² = ", round(R_Corr^2, 3)),
    x = "Year", y = "Abundance", color = NULL
  ) +
  theme_minimal(base_size = 14)

rm(V_CFIN_Interp, R_Corr, df_plot, df_long)

## ------------------------------------------ ##
#           Interpolation 4  -----
## ------------------------------------------ ##
# --- "Colebrook" interpolation method ---

# --- Reshape time series into 12 months × 65 years ---
M_CFIN_Missing <- matrix(V_CFIN_Missing, nrow = 12, ncol = 65)
M_CFIN_Interp <- M_CFIN_Missing

# --- Monthly and annual means ---
V_Mean_Year <- colMeans(M_CFIN_Missing, na.rm = TRUE)   # annual means
V_GlobMean_Year <- mean(V_Mean_Year, na.rm = TRUE)      # global annual mean
V_Mean_Month <- rowMeans(M_CFIN_Missing, na.rm = TRUE)  # monthly means

# --- Apply Colebrook interpolation ---
for (ii in seq_len(ncol(M_CFIN_Missing))) {
  temp <- M_CFIN_Missing[, ii]
  n_missing <- sum(is.na(temp))
  
  if (n_missing == 0) {
    next  # no missing values
  } else if (n_missing > 12) {
    temp2 <- temp
    temp2[is.na(temp2)] <- -9999
    M_CFIN_Interp[, ii] <- temp2
  } else {
    Idx_MissVal <- which(is.na(temp))
    if (!is.na(V_Mean_Year[ii]) && V_Mean_Year[ii] != 0) {
      M_CFIN_Interp[Idx_MissVal, ii] <- V_Mean_Month[Idx_MissVal] * (V_Mean_Year[ii] / V_GlobMean_Year)
    } else {
      M_CFIN_Interp[Idx_MissVal, ii] <- 0
    }
  }
}

# --- Replace placeholder with NA ---
M_CFIN_Interp[M_CFIN_Interp == -9999] <- NA
# --- Flatten back to vector ---
V_CFIN_Interp <- as.vector(M_CFIN_Interp)

# --- Calculate correlation ---
R_Corr <- cor(V_CFIN_Raw, V_CFIN_Interp, use = "complete.obs")

# --- Prepare data for plotting ---
df_plot <- data.frame(
  Time = V_X,
  Original = V_CFIN_Raw,
  Interpolated = V_CFIN_Interp,
  Missing = V_CFIN_Missing
)

df_long <- df_plot %>%
  pivot_longer(cols = c(Original, Interpolated, Missing),
               names_to = "Series",
               values_to = "Abundance")

# --- Plot ---
ggplot(df_long, aes(x = Time, y = Abundance, color = Series)) +
  geom_line(linewidth = 1.2, alpha = 0.6, na.rm = TRUE) +
  scale_color_manual(values = series_colors) +
  scale_x_continuous(
    breaks = seq(1, length(V_X), by = 60),
    labels = seq(1958, 2022, by = 5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(
    title = paste0("TS interpolation / r² = ", round(R_Corr^2, 3)),
    x = "Year", y = "Abundance", color = NULL
  ) +
  theme_minimal(base_size = 14)
