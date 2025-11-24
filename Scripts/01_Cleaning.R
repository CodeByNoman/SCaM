library(tidyverse)
library(readxl)

# 1. Read raw data
soil <- read_xlsx("Data/Raw/Soil_Water_Data.xlsx")
print(soil)

# 2. Columns rename for clarity
soil <- soil %>%
  rename(SAND = `Sand_%`, SILT = `Silt_%`, CLAY = `Clay_%`, 
         OM = `O.M`, SS = `Soil_Saturation_%`)

# 3. Missing values removal
soil <- soil %>%
  filter(!is.na(SAND) & !is.na(SILT) & !is.na(CLAY) 
         & !is.na(OM) & !is.na(SS))

# 4. Validation of sum of soil texture components / 
soil <- soil %>%
  mutate(total = SAND + SILT + CLAY) %>%
  filter(total == 100 )

# 5. Save cleaned & seperated soil texture data
ST.data <- soil %>% select(Address_code, SAND, SILT, CLAY, OM, SS)

write_csv(ST.data, "Data/Processed/Soil_Data_Cleaned.csv")
print(ST.data)

