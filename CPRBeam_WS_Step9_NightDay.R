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
# script #8 = STEP9 Night and Day
# translating Prog_CPRBeam_WS_Step9_NightDay.m to R 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)      
library(patchwork)   

# solar-elevation helper, eqv to f_CPRBeam_SolarAzEl.m
source("f_CPRBeam_SolarAzEl.R")

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
# --- Set WD (via the .Rproj) ---
# --- Load Rdata file ---
load("CPR_Data_CPRBeam.RData") # you could also do load("CPR_Data_CPRBeam.mat")

## ------------------------------------------ ##
#      Extract data & calculate solar elevation  -----
## ------------------------------------------ ##
# Build a UTC timestamp per sample from the metadata columns
# cols: Sample, Latitude, Longitude, Year, Month, Day, Hour, Minute
V_DateTime <- ISOdatetime(
  Data_LargeZooplankton$Year,
  Data_LargeZooplankton$Month,
  Data_LargeZooplankton$Day,
  Data_LargeZooplankton$Hour,
  Data_LargeZooplankton$Minute,
  0,
  tz = "UTC"
)

# Solar azimuth & elevation for every sample
# Alt = 0; lat/lon come from the sample coordinates
Solar <- f_CPRBeam_SolarAzEl(
  V_DateTime,
  Data_LargeZooplankton$Latitude,
  Data_LargeZooplankton$Longitude,
  alt = 0
)
V_El <- Solar$El   # solar elevation (degrees)
# alternative: oce::sunAngle(V_DateTime,
#   longitude = Data_LargeZooplankton$Longitude,
#   latitude  = Data_LargeZooplankton$Latitude)$altitude

## ------------------------------------------ ##
#      Represent elevation values  -----
## ------------------------------------------ ##
# MATLAB drew a per-sample bar of elevation (gold above horizon, grey below).
# In R a histogram of the elevations is the clearer equivalent, with the
# +/- 6 deg day/night cut-offs marked.
elev_df <- data.frame(El = V_El)
elev_df$period <- with(elev_df, ifelse(
  is.na(El), NA,
  ifelse(El < -6, "Night (< -6 deg)",
         ifelse(El > 6, "Day (> 6 deg)", "Twilight"))
))

ggplot(subset(elev_df, !is.na(El)), aes(El, fill = period)) +
  geom_histogram(binwidth = 2, colour = NA) +
  geom_vline(xintercept = c(-6, 6), linetype = 2, colour = "grey30") +
  scale_fill_manual(values = c("Day (> 6 deg)"   = "#F2C300",
                               "Twilight"         = "grey70",
                               "Night (< -6 deg)" = "grey40"),
                    name = NULL) +
  labs(x = "Solar elevation (degrees)", y = "Number of samples",
       title = "Distribution of solar elevation at sampling") +
  theme_minimal(base_size = 13)

## ------------------------------------------ ##
#      Table of solar elevation  -----
## ------------------------------------------ ##
T_Elevation <- cbind(Data_LargeZooplankton[, 1:8],
                     Elevation = V_El)

## ------------------------------------------ ##
#      Year x Month figure helper  -----
## ------------------------------------------ ##
# Reusable panel: a Year x Month heatmap with annual (right) and monthly (bottom)
# marginal-mean bars, mirroring the MATLAB 4x4 subplot layout. Reused for the
# elevation panel and the night/day copepod panels.
# `df` must have columns Year, Month, value on the full 1958:2022 x 1:12 grid.
plot_year_month <- function(df, title, fill_name,
                            option = "viridis", probs = c(0.20, 0.99)) {

  df$value[is.nan(df$value)] <- NA
  lim <- quantile(df$value, probs, na.rm = TRUE)

  ann <- df %>% group_by(Year)  %>% summarise(m = mean(value, na.rm = TRUE), .groups = "drop")
  mon <- df %>% group_by(Month) %>% summarise(m = mean(value, na.rm = TRUE), .groups = "drop")

  fill_scale <- scale_fill_viridis_c(option = option, limits = lim,
                                     oob = squish, na.value = "grey90",
                                     name = fill_name)

  p_main <- ggplot(df, aes(Month, Year, fill = value)) +
    geom_tile() +
    fill_scale +
    scale_x_continuous(breaks = 1:12, labels = month.abb, expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(1960, 2020, 10), expand = c(0, 0)) +
    labs(x = NULL, y = "Year", title = title) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  p_ann <- ggplot(ann, aes(m, Year, fill = m)) +
    geom_col() +
    scale_fill_viridis_c(option = option, limits = lim, oob = squish, guide = "none") +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Annual", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_blank())

  p_mon <- ggplot(mon, aes(Month, m, fill = m)) +
    geom_col() +
    scale_fill_viridis_c(option = option, limits = lim, oob = squish, guide = "none") +
    scale_x_continuous(breaks = 1:12, labels = month.abb, expand = c(0, 0)) +
    labs(x = NULL, y = "Monthly") +
    theme_minimal(base_size = 11)

  # main (top-left), annual bar (right), monthly bar (bottom)
  p_main + p_ann + p_mon + plot_spacer() +
    plot_layout(design = "AB\nCD", widths = c(4, 1), heights = c(4, 1),
                guides = "collect") &
    theme(legend.position = "bottom")
}

## ------------------------------------------ ##
#      Solar elevation: stats + figure  -----
## ------------------------------------------ ##
Stats_Elev <- T_Elevation %>%
  group_by(Year, Month) %>%
  summarise(value = mean(Elevation, na.rm = TRUE), .groups = "drop") %>%
  complete(Year = 1958:2022, Month = 1:12) %>%
  arrange(Year, Month)

plot_year_month(Stats_Elev,
                title = "Mean solar elevation",
                fill_name = "Elevation (deg)",
                option = "viridis")

## ------------------------------------------ ##
#      Mean large copepods by Night and Day  -----
## ------------------------------------------ ##
# Abundance block and copepod columns
M_Abund     <- as.matrix(Data_LargeZooplankton[, 9:ncol(Data_LargeZooplankton)])
cop_logical <- as.logical(List_LargeZooplankton$Copepods)
stopifnot(length(cop_logical) == ncol(M_Abund))  # columns must align with the list

# Per-sample mean abundance across all large-copepod taxa
V_CopMean <- rowMeans(M_Abund[, cop_logical, drop = FALSE], na.rm = TRUE)

# Day / night masks from solar elevation (NA-safe)
is_night <- !is.na(V_El) & V_El < -6
is_day   <- !is.na(V_El) & V_El >  6

# Metadata + per-sample copepod mean, split by period
meta <- Data_LargeZooplankton[, 1:8]

# Helper: build the Year x Month copepod stats for a subset
ym_copepods <- function(mask) {
  data.frame(meta[mask, ], value = V_CopMean[mask]) %>%
    group_by(Year, Month) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    complete(Year = 1958:2022, Month = 1:12) %>%
    arrange(Year, Month)
}

S_DayNight <- list(
  Night = list(Name = "Night", Stats = ym_copepods(is_night)),
  Day   = list(Name = "Day",   Stats = ym_copepods(is_day))
)

## ------------------------------------------ ##
#      Night / Day figures  -----
## ------------------------------------------ ##
# Night (dark "mako" palette, echoing the MATLAB gray colormap for night)
plot_year_month(S_DayNight$Night$Stats,
                title = "Night - mean large copepods",
                fill_name = "Mean abundance",
                option = "mako")

# Day
plot_year_month(S_DayNight$Day$Stats,
                title = "Day - mean large copepods",
                fill_name = "Mean abundance",
                option = "viridis")
