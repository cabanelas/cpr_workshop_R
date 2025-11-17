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
# script #1 = STEP2 Load data in R
# translating Prog_CPRBeam_WS_Step2_Load.m to R 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(tidyverse) 

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
# --- 1) Set WD

# --- 2) List all csv files
csv_files <- list.files("raw/CPRBeam_DataExtract", 
                        pattern = "\\.csv$", 
                        full.names = TRUE)

# --- 3) Function to read and name each table 
read_and_name_table <- function(file_path) {
  file_name <- basename(file_path)
  name_parts <- str_split(file_name, "_", simplify = TRUE)
  
  # Create object name from parts 3 and 4
  if (ncol(name_parts) >= 4) {
    object_name <- paste0(name_parts[3], "_", tools::file_path_sans_ext(name_parts[4]))
  } else {
    object_name <- tools::file_path_sans_ext(file_name)
  }
  
  # Read CSV
  table_data <- read_csv(file_path, show_col_types = FALSE)
  
  # Assign to global environment
  assign(object_name, table_data, envir = .GlobalEnv)
  
  return(object_name)
}

# --- 4) Read and assign all tables
loaded_tables <- map_chr(csv_files, read_and_name_table)

# --- 5) Save all loaded tables to an RData file
save(list = loaded_tables, file = "CPR_Data_CPRBeam.RData")