################################################################################
#############        Continuous Plankton Recorder     ##########################
#############    CPR-BEAMS Workshop -- November 2025  ##########################
## by: Alexandra Cabanelas 
## created NOV-2025
################################################################################

#https://website.whoi.edu/cpr-beams/
#https://drive.google.com/drive/u/0/folders/1FopWXQNcl4reXZtT5QPUmLUlFO911373

# translating Pierre Helaouet MATLAB code to R 
# script #4 = STEP5 Create timeseries
# translating Prog_CPRBeam_WS_Step5_Timeseries.m to R 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(dplyr)
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
                 Longitude >= -76  & Longitude <= -61  &
                 Latitude  >= 36.9 & Latitude  <= 46.8)

## ------------------------------------------ ##
#     Efficiently reorganize dataset  -----
## ------------------------------------------ ##
Area_1 <- list(
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
gc()

## ------------------------------------------ ##
#  Extracting data required to create ts in area_1 -----
## ------------------------------------------ ##
# Example with Calanus finmarchicus
# Find C. finmarchicus (CFIN) in the large zooplankton list

# --- By taxon name ---
Logi_TaxaName <- Area_1$List$LargeZoo$Taxon_Name == "Calanus finmarchicus"

# --- By CPR ID ---
Logi_TaxaID <- Area_1$List$LargeZoo$Accepted_ID == 40

# --- Check that the logical is identifying CFIN ---
Area_1$List$LargeZoo[Logi_TaxaName, 1:2]
Area_1$List$LargeZoo[Logi_TaxaID, 1:2]

# --- Extract CFIN abundance data (first column) ---
T_CFIN <- Area_1$Data$LargeZoo[, 1, drop = FALSE] #which(Logi_TaxaName)
colnames(T_CFIN) <- "CFIN"  # Rename column for clarity

# --- Create a table with spatio-temporal coordinates and CFIN ---
T_Extract <- cbind(Area_1$SpatioTemp, T_CFIN)

# --- Export as CSV ---
#write.csv(T_Extract, "output/CPRBeam_CFINinArea1.csv", row.names = FALSE)

# --- clean workspace ---
rm(Logi_TaxaName, Logi_TaxaID, T_CFIN)
gc()

## ------------------------------------------ ##
#     Create CFIN timeseries for area_1 -----
## ------------------------------------------ ##

# --- Summarize CFIN by Year and Month ---
T_Stats <- T_Extract %>%
  group_by(Year, Month) %>%
  summarise(
    GroupCount = dplyr::n(),              # like MATLAB GroupCount
    mean_CFIN  = mean(CFIN, na.rm = TRUE),# mean on CFIN
    .groups    = "drop"
  )

## ------------------------------------------ ##
# Problem: No sampling effort from 1978 to 1990 (included) -----
## ------------------------------------------ ##
T_Extract %>%
  count(Year) %>%
  filter(Year >= 1978 & Year <= 1990) # no samples in those years

## ------------------------------------------ ##
#       Force full Year-Month grid -----
## ------------------------------------------ ##
T_Stats_Full <- T_Stats %>%
  tidyr::complete( # MATLABs IncludeEmptyGroups= true; repmat + reshape
    Year  = 1958:2022,
    Month = 1:12,
    fill  = list(GroupCount = 0, mean_CFIN = NA_real_)
  ) %>%
  arrange(Year, Month)

## ------------------------------------------ ##
#           Visualize timeseries -----
## ------------------------------------------ ##
ggplot(T_Stats_Full, aes(x = Year, y = factor(Month), fill = mean_CFIN)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "D", na.value = "white", name = "Mean CFIN") + 
  #option = "plasma", na.value = "grey90"
  scale_x_continuous(
    breaks = seq(min(T_Stats_Full$Year), max(T_Stats_Full$Year), by = 2), 
    expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_discrete(
    limits = as.character(12:1), 
    labels = rev(month.abb), 
    expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = "Mean Abundance of Calanus finmarchicus",
       x = "Year", y = "Month", fill = "Mean CFIN") +
  theme_minimal() +
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   color = "black", size = 11))

## ------------------------------------------ ##
#       What happened in May 1977? -----
## ------------------------------------------ ##
# What is a "normal" average in may?

# --- Filter May data --- 
MayStats <- T_Stats_Full %>% filter(Month == 5)

# --- Calculate long-term May average --- 
V_MeanMay <- mean(MayStats$mean_CFIN, na.rm = TRUE)

# --- Plot mean CFIN abundance in May --- 
ggplot(MayStats, aes(x = Year, y = mean_CFIN)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = V_MeanMay, color = "darkorange", linewidth = 1.2) +
  labs(title = "Mean CFIN Abundance in May (1958–2022)",
       subtitle = paste("Long-term May average:", round(V_MeanMay, 2)),
       x = "Year", y = "Mean CFIN") +
  scale_x_continuous(expand = expansion(mult = c(0, 0)),
                     breaks = seq(1958, 2022, by = 2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_minimal() +
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   color = "black", size = 11))

## ------------------------------------------ ##
## What is a "normal" number of samples in may? ----
## ------------------------------------------ ##
# Inspect May 1977 sample data
T_Extract %>%
  filter(Year == 1977, Month == 5) %>%
  select(Year, Month, Day, CFIN)

## ------------------------------------------ ##
## Sample count in May across years ----
## ------------------------------------------ ##

# --- Calculate long-term May sample count average ---
V_MeanMay_Samp <- mean(MayStats$GroupCount, na.rm = TRUE)

# --- Plot sample count in May ---
ggplot(MayStats, aes(x = Year, y = GroupCount)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = V_MeanMay_Samp, color = "darkorange", linewidth = 1.2) +
  labs(
    title = "Sampling Effort in May (1958–2022)",
    subtitle = paste("Long-term May sample count:", round(V_MeanMay_Samp, 1)),
    x = "Year",
    y = "Number of Samples") +
  scale_x_continuous(expand = expansion(mult = c(0, 0)),
                     breaks = seq(1958, 2022, by = 2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_minimal() +
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   color = "black", size = 11))

## ------------------------------------------ ##
##        Visualise timeseries ----
## ------------------------------------------ ##
## Heatmap 1: mean_CFIN (0–200) ##
ggplot(T_Stats_Full, aes(x = Year, y = factor(Month), fill = mean_CFIN)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = viridis::viridis(100), 
                       limits = c(0, 200), na.value = "white", name = "Mean CFIN") +
  scale_x_continuous(breaks = seq(1958, 2022, by = 2),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(limits = as.character(12:1),
                   labels = rev(month.abb),
                   expand = expansion(mult = c(0, 0))) +
  labs(title = "Monthly Mean CFIN Abundance", x = "Year", y = "Month") +
  theme_minimal() +
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   color = "black", size = 11))

## Heatmap 2: GroupCount (0–50)  ##
ggplot(T_Stats_Full, aes(x = Year, y = factor(Month), fill = GroupCount)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = viridis::viridis(100),
                       limits = c(0, 50), 
                       na.value = "white", 
                       name = "Sample Count") +
  scale_x_continuous(breaks = seq(1958, 2022, by = 2),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(limits = as.character(12:1),
                   labels = rev(month.abb),
                   expand = expansion(mult = c(0, 0))) +
  labs(title = "Monthly Sampling Effort", x = "Year", y = "Month") +
  theme_minimal() +
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   color = "black", size = 11))

## ------------------------------------------ ##
##       Write the table to CSV ----
## ------------------------------------------ ##
#write.csv(T_Stats_Full, "output/CPRBeam_TS_CFINinArea1.csv", row.names = FALSE)

