library(soiltexture)

# 1. Read processed data
soil_data <- read.csv("Data/Processed/Soil_Data_Cleaned.csv")

# 2. Get ISSS/USDA/FAO texture matrix
tex_matrix <- TT.points.in.classes(tri.data = soil_data[
  c("SAND", "SILT", "CLAY")], class.sys = "USDA.TT")
print(tex_matrix)

# 3. Convert matrix (0/1) to a single texture class column
soil_data$Texture_Class <- apply(tex_matrix, 1, function(row) 
  {classes <- names(row)[row == 1]; if (length(classes) == 0) 
    return(NA); classes[1]})
SoilTC <- soil_data %>% select(Address_code, SAND, SILT, CLAY, Texture_Class)
print(SoilTC)
write.csv(SoilTC, "A:/Class_Assignments/SCaM/Outputs/Tables/Soiltexc.csv")

# 4. Soil Texture Triangle

tiff("A:/Class_Assignments/SCaM/Outputs/Figures/Soil_Texture_Triangle.tiff",
     width = 2000, height = 1800, res = 300, compression = "lzw")

par(family = "serif")

classes <- unique(soil_data$Texture_Class)
classes <- classes[!is.na(classes)]

palette_colors <- rainbow(length(classes))
names(palette_colors) <- classes

TT.plot(class.sys = "USDA.TT",  main = "USDA Soil Texture Classification",
  cex.lab = 1.1, cex.axis = 1.0, cex.main = 1.25, frame.bg.col = "white")

TT.points(tri.data = soil_data, geo = TT.geo.get("USDA.TT"),
  col = palette_colors[soil_data$Texture_Class], pch = 19, cex = 1.2)

legend("topleft", inset = c(1.04, 0.06), legend = classes, 
       col = palette_colors[classes], title = "Classes", title.cex = 1.0, 
       title.font = 2, pch = 19, pt.cex = 1.0, cex = 0.9, bty = "o")

dev.off()

# 5. 
