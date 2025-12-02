###############################################################################
#--------------- Soil Texture Interpolation (Kriging) Workflow ---------------#
###############################################################################

# 0. Libraries -----------------------------------------------------------------
library(sf)
library(sp)
library(terra)
library(gstat)
library(raster)
library(geodata)
library(ggplot2)
library(patchwork)
library(ggspatial)
library(soiltexture)

# 1. Paths & CRS ---------------------------------------------------------------
input_csv   <- "Data/Processed/Soil_Data_Cleaned.csv"
input_shp   <- "Data/Raw/study_area.shp"
out_fig_dir <- "Outputs/Figures"
out_tbl_dir <- "Outputs/Tables"

dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_tbl_dir, showWarnings = FALSE, recursive = TRUE)

# CRS definitions
epsg_input   <- 4326    # raw data (WGS84)
epsg_kriging <- 32642   # UTM (for variogram/kriging)
epsg_map     <- 4326    # final map (WGS84)

# 2. Load Data & Prepare -------------------------------------------------------
soil_df <- read.csv(input_csv, stringsAsFactors = FALSE)
req_cols <- c("X_Coordinates","Y_Coordinates","SAND","SILT","CLAY")
stopifnot(all(req_cols %in% names(soil_df)))

soil_sf <- st_as_sf(soil_df, coords = c("X_Coordinates","Y_Coordinates"),
                    crs = epsg_input, remove = FALSE)

aoi_sf <- st_read(input_shp, quiet = TRUE)
if (is.na(st_crs(aoi_sf))) st_crs(aoi_sf) <- epsg_input

# project to UTM for kriging
soil_utm <- st_transform(soil_sf, epsg_kriging)
aoi_utm  <- st_transform(aoi_sf,  epsg_kriging)

# convert to sp for gstat (sf -> sp conversion preserves CRS)
soil_sp <- as(soil_utm, "Spatial")   # gstat expects sp objects

# 3. Build Prediction Grid (m) -------------------------------------------------
res_m <- 250
bb <- st_bbox(aoi_utm)

grid_df <- expand.grid(
  x = seq(bb["xmin"], bb["xmax"], by = res_m),
  y = seq(bb["ymin"], bb["ymax"], by = res_m))

grid_sf <- st_as_sf(grid_df, coords = c("x","y"), crs = epsg_kriging)

# mask grid to AOI
inside_mat <- st_within(grid_sf, aoi_utm, sparse = FALSE)
if (ncol(inside_mat) == 0 || nrow(inside_mat) == 0) {
  stop("Prediction grid or AOI produced an empty containment matrix. 
       Check AOI and resolution.")}
inside <- apply(inside_mat, 1, any)
if (!any(inside)) stop("No prediction points fall inside AOI — adjust resolution or check AOI.")
pred_pts_sp <- as(grid_sf[inside, ], "Spatial")

# 4. Variogram Fitting ---------------------------------------------------------
compute_vario_params <- function(aoi_sf) 
  { bb <- st_bbox(aoi_sf)
  diag_len <- sqrt((bb["xmax"] - bb["xmin"])^2 + (bb["ymax"] - bb["ymin"])^2)
  list(cutoff = diag_len / 2, width = diag_len / 30)}
vp <- compute_vario_params(aoi_utm)

fit_variogram <- function(spdf, z) 
  { fmla <- as.formula(paste0(z, " ~ 1"))
  v <- try(variogram(fmla, locations = spdf, cutoff = vp$cutoff, 
                     width = vp$width), silent = TRUE)
if (inherits(v, "try-error") || nrow(v) == 0) stop("Variogram calculation failed for ", z)
  
# robust initial guesses: nugget = min gamma, sill = max gamma (approx)
nugget_guess <- suppressWarnings(min(v$gamma, na.rm = TRUE))
sill_guess   <- suppressWarnings(max(v$gamma, na.rm = TRUE))
if (!is.finite(nugget_guess)) nugget_guess <- 0
if (!is.finite(sill_guess) || sill_guess <= nugget_guess) sill_guess <- var(spdf[[z]], na.rm = TRUE)
psill_guess  <- max(0, sill_guess - nugget_guess)
  range_guess  <- max(v$dist, na.rm = TRUE) / 3
init <- vgm(psill = psill_guess, model = "Exp", range = range_guess, nugget = nugget_guess)
  
fit <- try(fit.variogram(v, model = init), silent = TRUE)
if (inherits(fit, "try-error") || any(is.na(fit$psill))) {
fit <- try(fit.variogram(v, model = vgm(psill = psill_guess, 
  "Sph", range_guess, nugget_guess)), silent = TRUE)  }

if (inherits(fit, "try-error") || any(is.na(fit$psill))) 
  { stop("Variogram model fitting failed for ", z, 
         ". Inspect data or adjust initial guesses.")}
  
list(v = v, model = fit)}

# fit 3 variograms
var_sand <- fit_variogram(soil_sp, "SAND")
var_silt <- fit_variogram(soil_sp, "SILT")
var_clay <- fit_variogram(soil_sp, "CLAY")

# function to save variogram figure
save_variogram <- function(vario, name) {
  png(file.path(out_fig_dir, paste0("variogram_", name, ".png")),
      width = 1200, height = 900, res = 150)
  
p <- plot(vario$v, vario$model, main = paste("Variogram -", toupper(name)))
  print(p)
  
  dev.off()}

# save all three
save_variogram(var_sand, "SAND")
save_variogram(var_silt, "SILT")
save_variogram(var_clay, "CLAY")

# Convert variogram and model into a ggplot object
vario_to_gg <- function(vario, title_prefix, var_name) {
  
  vdf <- vario$v
  mdf <- variogramLine(vario$model, maxdist = max(vdf$dist))
  
  ggplot(vdf, aes(x = dist, y = gamma)) +
    geom_point() +
    geom_line(data = mdf, aes(x = dist, y = gamma)) +
    theme_bw(base_size = 14, base_family = "serif") +
    labs( x = "Distance", y = "Semivariance",
      title = paste0(title_prefix, " - ", var_name) )}

pA <- vario_to_gg(var_sand, "A", "SAND")
pB <- vario_to_gg(var_silt, "B", "SILT")
pC <- vario_to_gg(var_clay, "C", "CLAY")

combined_plot <- pA / pB / pC  # vertical layout: A above B above C

ggsave(filename = file.path(out_fig_dir, "variogram_ABC_patchwork.png"),
  plot = combined_plot, width = 8, height = 12, dpi = 300)

# 5. Ordinary Kriging ----------------------------------------------------------
krige_one <- function(z, model) 
  { fmla <- as.formula(paste0(z, " ~ 1"))
  gstat::krige(fmla, locations = soil_sp, newdata = pred_pts_sp, model = model)}

kr_sand <- krige_one("SAND", var_sand$model)
kr_silt <- krige_one("SILT", var_silt$model)
kr_clay <- krige_one("CLAY", var_clay$model)

# ensure predictions are aligned (use coordinates from kr_sand)
coords <- coordinates(kr_sand)

kr_df <- data.frame(x_utm = coords[,1], y_utm = coords[,2],
  SAND  = if ("var1.pred" %in% names(kr_sand)) kr_sand$var1.pred else kr_sand@data$var1.pred,
  SILT  = if ("var1.pred" %in% names(kr_silt)) kr_silt$var1.pred else kr_silt@data$var1.pred,
  CLAY  = if ("var1.pred" %in% names(kr_clay)) kr_clay$var1.pred else kr_clay@data$var1.pred)

# 6. Normalize SSC to 100% -----------------------------------------------------
ssc <- kr_df[c("CLAY","SILT","SAND")]
ssc[ssc < 0] <- 0                           # remove negatives

row_sum <- rowSums(ssc)
# avoid division by zero and NAs: if sum == 0, set to 1 (keeps components 0)
row_sum[row_sum == 0] <- 1
ssc_norm <- sweep(ssc, 1, row_sum, "/") * 100
kr_df[c("CLAY","SILT","SAND")] <- ssc_norm

# 7. USDA Classification -------------------------------------------------------
usda <- TT.points.in.classes(kr_df[c("CLAY","SILT","SAND")], class.sys = "USDA.TT")

if (is.factor(usda)) {  kr_df$USDA_texture <- as.character(usda)} else {
  usda_mat <- as.matrix(usda)
  if (ncol(usda_mat) == 0) stop("USDA classification returned empty matrix")
  max_idx <- max.col(usda_mat, ties.method = "first")
  kr_df$USDA_texture <- colnames(usda_mat)[max_idx]}

# 8. Reproject Predictions to WGS84 --------------------------------------------
pred_sf_utm <- st_as_sf(kr_df, coords = c("x_utm", "y_utm"), crs = epsg_kriging, remove = FALSE)

pred_sf_wgs84 <- st_transform(pred_sf_utm, epsg_map)

coords_ll <- st_coordinates(pred_sf_wgs84)
pred_sf_wgs84$lon <- coords_ll[,1]
pred_sf_wgs84$lat <- coords_ll[,2]

# 9. Map Palette ---------------------------------------------------------------
palette_usda <- c(
  "Sand"             = "#F4D06F",
  "Loamy sand"       = "#F7E1A0",
  "Sandy loam"       = "#FFE699",
  "Loam"             = "#8BC34A",
  "Silt loam"        = "#AEDFF7",
  "Silt"             = "#6EC1E4",
  "Sandy clay loam"  = "#F5A65B",
  "Clay loam"        = "#C58C5E",
  "Silty clay loam"  = "#A38BD4",
  "Sandy clay"       = "#E07B39",
  "Silty clay"       = "#7A4E9A",
  "Clay"             = "#D9534F")

used <- sort(unique(pred_sf_wgs84$USDA_texture))
pred_sf_wgs84$USDA_texture <- factor(pred_sf_wgs84$USDA_texture, levels = used)

# match known palette names to used classes
palette_used <- palette_usda[used]
missing <- is.na(palette_used)
if (any(missing)) { palette_used[missing] <- hcl.colors(sum(missing), "Set 3")}
names(palette_used) <- used

# 10. Plot Map -----------------------------------------------------------------
aoi_ll <- st_transform(aoi_sf, epsg_map)

mean_lat <- mean(pred_sf_wgs84$lat, na.rm = TRUE)
lon_step <- res_m / (111320 * cos(mean_lat * pi/180))
lat_step <- res_m / 110540

p_kriging <- ggplot() +
  geom_sf( data = aoi_ll, fill = NA, color = "white", linewidth = 1.5) +
  geom_tile(data = st_drop_geometry(pred_sf_wgs84),
            aes(x = lon, y = lat, fill = USDA_texture),
            width = lon_step, height = lat_step) +
  scale_fill_manual(values = palette_used, drop = FALSE) +
  coord_sf(expand = T) +  labs(x = NULL, y = NULL) +
  theme_bw(base_family = "serif") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  annotation_north_arrow(location = "tr", style = north_arrow_fancy_orienteering(text_family = "serif")) +
  annotation_scale(location = "br", text_family = "serif") +
  ggtitle("USDA Soil Texture Classes (Normalized & Kriged)") + 
  theme(plot.title = element_text(face = "bold", size = 16), axis.text.y = element_text(angle = 90, hjust = 0.5))

print(p_kriging)

ggsave(file.path(out_fig_dir, "USDA_soil_texture_classes_Kriging_map.png"),
       plot = p_kriging, width = 10, height = 8, dpi = 300)


###############################################################################
# ---------------- Soil Texture IDW Interpolation (EXTENSION) ----------------#
###############################################################################

# 1. IDW function (uses soil_sp & pred_pts_sp from kriging section) ------------
run_idw <- function(var, idp = 2, nmax = 20) 
  {fmla <- as.formula(paste0(var, " ~ 1"))
  pred <- gstat::idw(
    formula   = fmla,
    locations = soil_sp,
    newdata   = pred_pts_sp,
    idp       = idp,
    nmax      = nmax)
  pred$var1.pred}

# 2. Run IDW for SAND–SILT–CLAY ------------------------------------------------
idw_df <- data.frame(
  x_utm = pred_pts_sp@coords[,1],
  y_utm = pred_pts_sp@coords[,2],
  SAND  = run_idw("SAND"),
  SILT  = run_idw("SILT"),
  CLAY  = run_idw("CLAY"))

# 3. Normalize SSC to 100% -----------------------------------------------------
ssc <- idw_df[, c("CLAY","SILT","SAND")]
ssc[ssc < 0] <- 0
row_sum <- rowSums(ssc)
row_sum[row_sum == 0] <- 1
ssc_norm <- sweep(ssc, 1, row_sum, "/") * 100
idw_df[, c("CLAY","SILT","SAND")] <- ssc_norm

# 4. USDA classification (same system as kriging) ------------------------------
usda <- TT.points.in.classes(idw_df[c("CLAY","SILT","SAND")], class.sys = "USDA.TT")
if (is.factor(usda)) {idw_df$USDA_texture <- as.character(usda)} else {
  m <- as.matrix(usda)
  idx <- max.col(m, ties.method = "first")
  idw_df$USDA_texture <- colnames(m)[idx]}

# 5. Convert to sf and project to WGS84 (same as kriging) ----------------------
idw_sf_utm <- st_as_sf(idw_df, coords = c("x_utm", "y_utm"),
                       crs = epsg_kriging, remove = FALSE)
idw_sf_wgs84 <- st_transform(idw_sf_utm, epsg_map)
coords_idw <- st_coordinates(idw_sf_wgs84)
idw_sf_wgs84$lon <- coords_idw[,1]
idw_sf_wgs84$lat <- coords_idw[,2]

# 6. Build IDW-specific legend (critical fix) ----------------------------------
used_idw <- sort(unique(idw_sf_wgs84$USDA_texture))

idw_sf_wgs84$USDA_texture <- factor(idw_sf_wgs84$USDA_texture, levels = used_idw)

palette_used_idw <- palette_usda[used_idw]

missing_idw <- is.na(palette_used_idw)
if (any(missing_idw)) {palette_used_idw[missing_idw] <- 
    hcl.colors(sum(missing_idw), "Set 3")}
names(palette_used_idw) <- used_idw

# 7. Plot IDW map (identical style to kriging map) -----------------------------
p_idw <- ggplot() +
  geom_sf(data = aoi_ll, fill = NA, color = "white", linewidth = 1.5) +
  geom_tile(data = st_drop_geometry(idw_sf_wgs84),
    aes(x = lon, y = lat, fill = USDA_texture),
    width = lon_step, height = lat_step) +
  scale_fill_manual(values = palette_used_idw, drop = FALSE) +
  coord_sf(expand = TRUE) + labs(x = NULL, y = NULL) +
  theme_bw(base_family = "serif") +
  theme(legend.position = "bottom", legend.title = element_blank(),
    panel.grid = element_blank()) +
  annotation_north_arrow(location = "tr",
    style = north_arrow_fancy_orienteering(text_family = "serif")) +
  annotation_scale(location = "br", text_family = "serif") +
  ggtitle("USDA Soil Texture Classes (Normalized & IDW)") +
  theme(plot.title = element_text(face = "bold", size = 16),
    axis.text.y = element_text(angle = 90, hjust = 0.5))

print(p_idw)

ggsave(file.path(out_fig_dir, "USDA_soil_texture_classes_IDW_map.png"),
       plot = p_idw, width = 10, height = 8, dpi = 300)


###############################################################################
# ----------- SoilGrids.org Soil Texture Classification (EXTENSION) -----------#
###############################################################################

# 1. Load study area (ensure WGS84) --------------------------------------------
study_area <- st_read("A:/Class_Assignments/SCaM/Data/Raw/study_area.shp")
study_area_sf <- st_transform(study_area, 4326)
study_area <- vect(study_area_sf)  # sf -> terra

# 2. Download SoilGrids raster datasets (5–15 cm layer = depth = 15) -----------
path <- "A:/Class_Assignments/SCaM/Data/Raw/soilgrids_data"

sand <- soil_world("sand", depth = 15, stat = "mean", path = path)
silt <- soil_world("silt", depth = 15, stat = "mean", path = path)
clay <- soil_world("clay", depth = 15, stat = "mean", path = path)

# 3. Crop + mask study area -> View -> Save .tif -------------------------------

sand_c <- mask(crop(sand, study_area), study_area)
silt_c <- mask(crop(silt, study_area), study_area)
clay_c <- mask(crop(clay, study_area), study_area)

plot(sand_c, main = "Sand (Clipped)")
plot(silt_c, main = "Silt (Clipped)")
plot(clay_c, main = "Clay (Clipped)")

writeRaster(sand_c, "A:/Class_Assignments/SCaM/Data/Processed/sand_clipped.tif")
writeRaster(silt_c, "A:/Class_Assignments/SCaM/Data/Processed/silt_clipped.tif")
writeRaster(clay_c, "A:/Class_Assignments/SCaM/Data/Processed/clay_clipped.tif")

# 4. Read and project rasters --------------------------------------------------
SAND.raw <- raster("Data/Processed/clay_clipped.tif")
CLAY.raw <- raster("Data/Processed/clay_clipped.tif")

# Read projection system
crs(CLAY.raw)
projection(CLAY.raw)

# Create projection system (correct syntax)
utm <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Assign/Project to new projection system
CLAY <- projectRaster(from = CLAY.raw, crs = utm, method = "ngb")
crs(CLAY)
projection(CLAY)

SAND <- projectRaster(from = SAND.raw, crs = utm, method = "ngb")

# 5. Creat & join point from rasters -------------------------------------------
point1 <- rasterToPoints(CLAY, fun = function(x) {x>0}, spatial = TRUE)
plot(CLAY)
names(point1)

# Extract sand values
x <- extract(SAND, point1)

# Join the extracted values to point shape file
point2 <- cbind(point1, x)
names(point2)
names(point2) <- c("CLAY", "SAND")
names(point2)

# 6. Convert to data frame -----------------------------------------------------
point3 <- as.data.frame(point2)
names(point3)
head(point3)

# 7. Calculate silt value and create soil texture classes ----------------------
point3$SILT <- 100-point3$CLAY-point3$SAND
head(point3)

# Create soil texture class
point3$Texture <- TT.points.in.classes(tri.data = point3, class.sys = "USDA.TT",
                                       PiC.type = "t")
head(point3)
str(point3)

point3$TextureClass = factor(point3$Texture)
point3$i.texclass = as.numeric(point3$TextureClass)
r = rasterFromXYZ(point3[, c("x", "y", "i.texclass")], crs = utm)
print(r)

r[] = as.numeric(factor(levels(point3$TextureClass))[r[]])

# Define the custom color palette for soil texture classes
manual_palette <- c(
  "#b2182b", "#ffd8b1", "#f4a582", "#92c5de",
  "#fddbc7", "#1b7837", "#2166ac", "#67a9cf",
  "#d1e5f0", "#f7f7f7", "#F4D06F", "#F5A65B"
)

# Now plot the raster with the manual color palette
plot(r, main = "Soil Texture Classes (Rasterized)", col = manual_palette)

# Get the resolution (cell size) of the raster
lon_step <- res(r)[1]  # Longitude resolution
lat_step <- res(r)[2]  # Latitude resolution

# 8. Build ggplot map with USDA Soil Texture classes ---------------------------
p_grid <- ggplot() +
  geom_sf(data = study_area_sf, fill = NA, color = "white", linewidth = 1.5) +  
  geom_tile(data = point3, aes(x = x, y = y, fill = Texture),  
            width = lon_step, height = lat_step) +
  scale_fill_manual(values = manual_palette, drop = FALSE) +  
  coord_sf(expand = TRUE) +  
  labs(x = NULL, y = NULL) +  
  theme_bw(base_family = "serif") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        panel.grid = element_blank()) +
  annotation_north_arrow(location = "tr", style = north_arrow_fancy_orienteering
                         (text_family = "serif")) +
  annotation_scale(location = "br", text_family = "serif") +
  ggtitle("SoilGrids Data (250 m) USDA Texture Classes") +
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(angle = 90, hjust = 0.5))

print(p_grid)

# 9. Save the map as PNG -------------------------------------------------------

out_fig_dir <- "A:/Class_Assignments/SCaM/Outputs/Figures"
ggsave(file.path(out_fig_dir, "SoilGrid_Data_Texture_Classification_map.png"),
       plot = p_grid, width = 10, height = 8, dpi = 300)

# 10. Save the table of x, y, and Texture --------------------------------------
soilgrids_tab <- st_drop_geometry(point3)[, c ("x", "y", "Texture")]
colnames(soilgrids_tab)[3] <- "Texture"
print(soilgrids_tab)

write.csv(soilgrids_tab, file.path(out_tbl_dir, "SoilGridss_table.csv"),
          row.names = FALSE)

#--------------------------------Patch Work-------------------------------------
A <- p_kriging + ggtitle("A: USDA Soil Texture Map (Kriging)")
B <- p_idw     + ggtitle("B: USDA Soil Texture Map (IDW)")
C <- p_grid    + ggtitle("C: USDA Soil Texture Map (SoilGrids)")

combined_plot <- (A | B | C)
plot(combined_plot)

ggsave(file.path(out_fig_dir, "Kriging_IDW_SoilGrids_comparison.png"), 
       plot = combined_plot, width = 12, height = 6, dpi = 300)

#-------------------------table (Kriging vs IDW)--------------------------------

# 1. Extract Kriging table
krig_tab <- st_drop_geometry(pred_sf_wgs84)[, c("lon", "lat", "USDA_texture")]
colnames(krig_tab)[3] <- "USDA_texture_Krig"

# 2. Extract IDW table
idw_tab <- st_drop_geometry(idw_sf_wgs84)[, c("lon", "lat", "USDA_texture")]
colnames(idw_tab)[3] <- "USDA_texture_IDW"

# 3. Merge into one comparison table
comparison_table <- merge(krig_tab, idw_tab, by = c("lon", "lat"), all = TRUE)

# 4. Save to Outputs/Tables
write.csv(comparison_table, file.path(out_tbl_dir, "comparison_table.csv"),
          row.names = FALSE)

head(comparison_table) # View


