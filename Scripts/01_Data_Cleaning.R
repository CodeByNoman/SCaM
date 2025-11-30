###############################################################################
#------------- Soil Texture Data Cleaning & Validation Workflow ---------------#
###############################################################################

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(readxl)

# 1. Read raw data -------------------------------------------------------------
soil <- read_xlsx("Data/Raw/Soil_Water_Data.xlsx")
print(soil)

# 2. Columns rename for clarity ------------------------------------------------
soil <- soil %>%
  rename(SAND = `Sand_%`, SILT = `Silt_%`, CLAY = `Clay_%`)

# 3. Missing values removal ----------------------------------------------------
soil <- soil %>%
  filter(!is.na(SAND) & !is.na(SILT) & !is.na(CLAY))

# 4. Validation of sum of soil texture components ------------------------------
soil <- soil %>%
  mutate(total = SAND + SILT + CLAY) %>%
  filter(total == 100 )

# 5. Save cleaned & seperated soil texture data --------------------------------
ST.data <- soil %>% select(Address_code, X_Coordinates, Y_Coordinates, 
                           SAND, SILT, CLAY)
write_csv(ST.data, "Data/Processed/Soil_Data_Cleaned.csv")
print(ST.data)

