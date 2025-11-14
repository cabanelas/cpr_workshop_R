################################################################################
#############        Continuous Plankton Recorder     ##########################
#############    CPR-BEAMS Workshop -- November 2025  ##########################
## by: Alexandra Cabanelas 
## created NOV-2025
################################################################################

#https://website.whoi.edu/cpr-beams/
#https://drive.google.com/drive/u/0/folders/1FopWXQNcl4reXZtT5QPUmLUlFO911373

# translating Pierre Helaouet MATLAB code to R 
# script #2 = STEP3 Check data loaded in R
# translating Prog_CPRBeam_WS_Step3_Check.m to R 
# heatmap of sampling effort

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(dplyr)
library(ggplot2)

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
# --- 1) Set WD

# --- 2) Load Rdata file
load("CPR_Data_CPRBeam.RData")

# --- 3) Check
## Confirm number of loaded objects matches number of csv files in the 
#CPRBeam_DataExtract folder
length(ls()) == length(list.files("raw/CPRBeam_DataExtract", pattern = "\\.csv$"))
## Check that all tables have the correct dimensions
dim(Data_LargeZooplankton)

## ------------------------------------------ ##
#           --- 4) Crude map -----
## ------------------------------------------ ##
# --- sampling locations colored by year ---
T_Map <- Data_LargeZooplankton %>%
  select(1:8) %>%  #metadata columns
  mutate(YearCat = as.factor(Year))

#rm(list = setdiff(ls(), "T_Map"))

ggplot(T_Map, aes(x = Longitude, y = Latitude, color = YearCat)) +
  borders("world", colour = "gray80", fill = "gray95") +
  #or   geom_sf(data = world, fill = "gray95", color = "gray80") +
  geom_point(alpha = 0.5, size = 2) +
  coord_cartesian(xlim = c(-110, 50), ylim = c(30, 90)) +
  #oe   coord_sf(xlim = c(-110, 50), ylim = c(30, 90), expand = FALSE) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude", color = "Year")

## ------------------------------------------ ##
#           --- 5) Heatmap -----
## ------------------------------------------ ##
# Heatmap of sampling effort by Year and Month
# Generate table with sampling effort 
T_SampEff <- Data_LargeZooplankton %>%
  select(Year, Month) %>%
  count(Year, Month, name = "GroupCount")

ggplot(T_SampEff, aes(x = factor(Year), y = factor(Month), 
                      fill = GroupCount)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(title = "Sampling Effort by Year and Month",
       x = "Year", y = "Month", fill = "Samples") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   color = "black", size = 10),
        legend.position = "right",
        legend.direction = "vertical")

ggplot(T_SampEff, aes(x = factor(Year), y = factor(Month), 
                      fill = GroupCount)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Samples") +
  scale_x_discrete(breaks = scales::pretty_breaks()) +
  labs(x = "Year", y = "Month",
       title = "Sampling effort by Year × Month") +
  theme_minimal(base_size = 11)




## ------------------------------------------ ##
#            Other map options -----
## ------------------------------------------ ##
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
T_sf  <- st_as_sf(T_Map, coords = c("Longitude","Latitude"), crs = 4326, remove = FALSE)
world_hi <- ne_countries(scale = "large", returnclass = "sf")
ggplot() +
  geom_sf(data = world_hi, fill = "grey95", color = "grey75", linewidth = 0.25) +
  geom_sf(data = T_sf, aes(color = YearCat), size = 1.6, alpha = 0.5) +
  coord_sf(xlim = c(-110, 50), ylim = c(30, 90), expand = FALSE) +
  labs(x = NULL, y = NULL, color = "Year") +
  theme_minimal(base_size = 11)

library(ggrastr)
bbox <- list(xmin = -110, xmax = 50, ymin = 30, ymax = 90)
world_lite <- map_data("world") |>
  dplyr::filter(long >= bbox$xmin, long <= bbox$xmax,
                lat  >= bbox$ymin, lat  <= bbox$ymax)
ggplot() +
  geom_polygon(data = world_lite,
               aes(long, lat, group = group),
               fill = "grey95", color = "grey80", linewidth = 0.2) +
  geom_point_rast(data = T_Map,
                  aes(x = Longitude, y = Latitude, color = YearCat),
                  raster.dpi = 150, size = 1.2, alpha = 0.6) +
  coord_quickmap(xlim = c(-110, 50), ylim = c(30, 90), expand = FALSE) +
  theme_minimal(base_size = 11)
