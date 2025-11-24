# Soil Texture Classification and Spatial Mapping Using Field Samples and World Soil Grid Data
This project focuses on classifying soil texture using field sample data and generating spatially continuous maps using interpolation techniques in R. The results are compared with the World Soil Grid predictions to understand spatial variability and model accuracy.
## Research Question
How do soil texture classes derived from field-measured soil samples differ from those estimated by the World Soil Grid, and what spatial patterns emerge when using interpolation techniques in R?
## Hypothesis
Field-derived soil texture classes will show finer spatial variation than World Soil Grid predictions, with higher disagreement in areas of micro-scale heterogeneity.
## Objective
To classify soil texture using the USDA texture triangle, map spatial variability using interpolation, and compare results with World Soil Grid data.
## Data Sources
#### 1. Field Soil Data
- Sand (%), Silt (%), Clay (%)
- Coordinates (WGS84)
- Format: CSV / XLSX
#### 2. World Soil Grid (ISRIC)
- Sand, Silt, Clay (0-5 cm depth)
- Resolution: 250 m
- Format: GeoTIFF
#### 3. Study Area
-  Shapefile defining the project boundary
## Folder Structure
```r
# Soil-Texture-Mapping/
│
├── README.md                    # Project overview
├── LICENSE                      # MIT license
├── .gitignore                   # Ignored files
│
├── Data/
│   ├── Raw/                     # Original datasets (CSV, shapefiles)
│   └── Processed/               # Cleaned/normalized data
│
├── Requirements.txt             # R package list
│
├── Scripts/
│   ├── 01_Data_Cleaning.R       # Data preprocessing
│   ├── 02_Classification.R      # USDA soil texture classification
│   └── 03_Mapping.R             # Spatial interpolation & mapping
│
├── Outputs/
│   ├── Figures/                 # Generated maps
│   └── Tables/                  # Classification results
│
└── Documents/
    ├── Final_Report.qmd         # Quarto analytical report
    └── Presentation.qmd         # Slides
```
## Methods
- Clean and standardize soil sample data
- Classify soil texture using soiltexture package
- Interpolate sand/silt/clay using IDW & Kriging
- Classify interpolated rasters
- Compare with World Soil Grid data
- Generate maps and summary tables
## How to Run This Project
**1.** Install required packages:
```r
# source("requirements.txt")
```
**2.** Run scripts in the following order:
- `01_Data_Cleaning.R`
- `02_Classification.R`
- `03_Mapping.R`

**3.** Render the Quarto report:
```r
# Quarto render Documents/Final_Report.qmd
```
## Author
**Noman Ahmad**  
Doctorate in Agricultural and Natural Sciences  
AGP3141 – Environmental Data Visualization in R  
Date: Spring 2025  
[Linkedin](https://www.linkedin.com/in/noman-ahmad-6960bb176/); [ResearchGate](https://www.researchgate.net/profile/Noman-Ahmad-11?ev=hdr_xprf); [ORCiD](https://orcid.org/0000-0002-2663-9190)


