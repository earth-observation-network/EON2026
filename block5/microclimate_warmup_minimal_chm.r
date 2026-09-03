# =============================================================================
# Microclimate warmup: same points, six spatial statements
# =============================================================================
#
# Teaching objective:
# Eight station measurements are converted into spatial surfaces. The surfaces
# differ because the models use different transfer assumptions and different
# environmental information. A raster prediction is therefore a modelled
# statement, not an additional observation.
#
# Required files:
#   data/MobilePolsterhaus_spatial_hourly.rds
#   data/dem1_d.tif
#   data/dom1_d.tif
#
# Model comparison:
#   Voronoi       nearest observed station
#   IDW           inverse-distance weighted neighbouring stations
#   LM DEM        one linear temperature–terrain-height relationship
#   LM CHM        one linear temperature–canopy-height relationship
#   RF DEM        position and terrain elevation
#   RF DEM + CHM  position, terrain elevation and canopy height

library(sf)
library(terra)
library(gstat)
library(randomForest)


# -----------------------------------------------------------------------------
# 1. Data and the time slot with the largest spatial contrast
# -----------------------------------------------------------------------------

m <- readRDS("data/MobilePolsterhaus_spatial_hourly.rds")
dem <- rast("data/dem1_d.tif")
dom <- rast("data/dom1_d.tif")
names(dem) <- "altitude"

# DEM and DOM must occupy the same cells before they can be subtracted.
# Bilinear resampling is used because both rasters contain continuous heights.
dom <- terra::resample(dom, dem, method = "bilinear")
names(dom) <- "surface_height"

# The canopy height model is not read as an independent input. It is derived:
# surface elevation (DOM) minus terrain elevation (DEM) = object/canopy height.
chm <- dom - dem
names(chm) <- "chm"

# This hour was selected because it contains the largest contrast between the
# eight stations in the available time series. All interpretations therefore
# refer first of all to this timestamp and its particular weather regime.
target_time <- as.POSIXct(
  "2026-09-01 14:00:00",
  tz = "Europe/Berlin"
)

m <- m[m$time == target_time, ]

# All distance calculations, raster extractions and overlays require a common
# coordinate reference system. The DEM supplies the target CRS.
m <- st_transform(m, terra::crs(dem))
m$temp <- m$air_temp

# Extract both environmental predictors at the eight station locations. These
# values are the only DEM/CHM combinations observed together with temperature.
m$altitude <- terra::extract(dem, terra::vect(m))$altitude
m$chm <- terra::extract(chm, terra::vect(m))$chm

pts <- m[, c("air_channel", "temp", "altitude", "chm")]


# -----------------------------------------------------------------------------
# 2. Prediction area: only 10 m around the station network
# -----------------------------------------------------------------------------

area <- terra::buffer(
  terra::convHull(terra::vect(pts)),
  width = 10
)

# The hull plus buffer defines where predictions are displayed. It is a
# reporting domain, not proof that every included cell is equally well supported.
dem <- terra::mask(terra::crop(dem, area), area)
dom <- terra::mask(terra::crop(dom, area), area)
chm <- terra::mask(terra::crop(chm, area), area)

predictors <- c(dem, chm)

# One table row is created for every raster cell. It contains cell number,
# coordinates, terrain elevation and canopy height. All models predict to these
# exact same locations.
grid <- as.data.frame(predictors, xy = TRUE, cells = TRUE, na.rm = TRUE)

grid_sf <- st_as_sf(
  grid,
  coords = c("x", "y"),
  crs = st_crs(pts),
  remove = FALSE
)

make_map <- function(prediction, map_name) {
  # Model predictions are returned as vectors. This helper writes each value
  # back into its original DEM cell so every result has identical geometry.
  r <- dem
  values(r) <- NA
  values(r)[grid$cell] <- prediction
  names(r) <- map_name
  r
}


# -----------------------------------------------------------------------------
# 3. Six spatial assumptions
# -----------------------------------------------------------------------------

# A: nearest station wins. No environmental raster enters the prediction.
vor <- idw(
  temp ~ 1,
  locations = pts,
  newdata = grid_sf,
  nmax = 1
)
map_vor <- make_map(vor$var1.pred, "Voronoi")


# B: the four nearest stations contribute with inverse squared distance weights.
idw_result <- idw(
  temp ~ 1,
  locations = pts,
  newdata = grid_sf,
  nmax = 4,
  idp = 2
)
map_idw <- make_map(idw_result$var1.pred, "IDW")


# C: one global linear relationship transfers temperature through terrain height.
fit_lm_altitude <- lm(
  temp ~ altitude,
  data = st_drop_geometry(pts)
)
map_lm_altitude <- predict(dem, fit_lm_altitude)
names(map_lm_altitude) <- "LM_altitude"


# D: the same model form is used, but the substantive predictor is CHM height.
# Comparing C and D isolates the effect of environmental-variable choice.
fit_lm_chm <- lm(
  temp ~ chm,
  data = st_drop_geometry(pts)
)
map_lm_chm <- predict(chm, fit_lm_chm)
names(map_lm_chm) <- "LM_CHM"


# E: Random Forest learns nonlinear partitions from x, y and terrain elevation.
xy <- st_coordinates(pts)
pts$x <- xy[, 1]
pts$y <- xy[, 2]

fit_rf_dem <- randomForest(
  temp ~ x + y + altitude,
  data = st_drop_geometry(pts),
  ntree = 500
)

rf_dem_prediction <- predict(
  fit_rf_dem,
  newdata = grid[, c("x", "y", "altitude")]
)
map_rf_dem <- make_map(rf_dem_prediction, "RF_position_altitude")


# F: the same RF additionally receives CHM. Comparing E and F asks whether CHM
# adds useful predictive information within an otherwise identical algorithm.
fit_rf_chm <- randomForest(
  temp ~ x + y + altitude + chm,
  data = st_drop_geometry(pts),
  ntree = 500
)

rf_chm_prediction <- predict(
  fit_rf_chm,
  newdata = grid[, c("x", "y", "altitude", "chm")]
)
map_rf_chm <- make_map(rf_chm_prediction, "RF_position_altitude_CHM")


# -----------------------------------------------------------------------------
# 4. Leave-one-out comparison
# -----------------------------------------------------------------------------

rmse <- function(observed, predicted) {
  sqrt(mean((observed - predicted)^2))
}

# In every iteration one station is excluded, the model is fitted to the other
# seven stations, and the excluded temperature is predicted. This produces one
# genuinely held-out prediction for each station.

lm_altitude_cv <- rep(NA, nrow(pts))
lm_chm_cv <- rep(NA, nrow(pts))
rf_dem_cv <- rep(NA, nrow(pts))
rf_chm_cv <- rep(NA, nrow(pts))

for (i in seq_len(nrow(pts))) {
  train <- pts[-i, ]
  test <- pts[i, ]
  
  lm_altitude_i <- lm(
    temp ~ altitude,
    data = st_drop_geometry(train)
  )
  lm_altitude_cv[i] <- predict(
    lm_altitude_i,
    newdata = st_drop_geometry(test)
  )
  
  lm_chm_i <- lm(
    temp ~ chm,
    data = st_drop_geometry(train)
  )
  lm_chm_cv[i] <- predict(
    lm_chm_i,
    newdata = st_drop_geometry(test)
  )
  
  rf_dem_i <- randomForest(
    temp ~ x + y + altitude,
    data = st_drop_geometry(train),
    ntree = 500
  )
  rf_dem_cv[i] <- predict(
    rf_dem_i,
    newdata = st_drop_geometry(test)
  )
  
  rf_chm_i <- randomForest(
    temp ~ x + y + altitude + chm,
    data = st_drop_geometry(train),
    ntree = 500
  )
  rf_chm_cv[i] <- predict(
    rf_chm_i,
    newdata = st_drop_geometry(test)
  )
}

vor_model <- gstat(
  formula = temp ~ 1,
  locations = pts,
  nmax = 1,
  set = list(idp = 2)
)
vor_cv <- gstat.cv(vor_model, nfold = nrow(pts))

idw_model <- gstat(
  formula = temp ~ 1,
  locations = pts,
  nmax = 4,
  set = list(idp = 2)
)
idw_cv <- gstat.cv(idw_model, nfold = nrow(pts))

model_error <- data.frame(
  model = c(
    "Voronoi", "IDW", "LM terrain elevation", "LM canopy height",
    "RF position + terrain", "RF position + terrain + CHM"
  ),
  MAE = c(
    mean(abs(vor_cv$residual)),
    mean(abs(idw_cv$residual)),
    mean(abs(pts$temp - lm_altitude_cv)),
    mean(abs(pts$temp - lm_chm_cv)),
    mean(abs(pts$temp - rf_dem_cv)),
    mean(abs(pts$temp - rf_chm_cv))
  ),
  RMSE = c(
    sqrt(mean(vor_cv$residual^2)),
    sqrt(mean(idw_cv$residual^2)),
    rmse(pts$temp, lm_altitude_cv),
    rmse(pts$temp, lm_chm_cv),
    rmse(pts$temp, rf_dem_cv),
    rmse(pts$temp, rf_chm_cv)
  )
)

print(model_error)


# =============================================================================
# 5. VISUALISATION
# =============================================================================
#
# The output is deliberately split into two figures.
#
# Figure 1 documents the spatial input:
#   - DEM: terrain elevation
#   - CHM: vegetation height, calculated as DOM minus DEM
#
# Figure 2 compares the six interpolation/model assumptions.
# All temperature maps use exactly the same colours and value range. Therefore,
# one temperature legend is sufficient and differences between panels are real
# model differences rather than artefacts of independently scaled legends.

maps <- c(
  map_vor,
  map_idw,
  map_lm_altitude,
  map_lm_chm,
  map_rf_dem,
  map_rf_chm
)

# The observed minimum and maximum define the common temperature scale.
# No model is allowed to obtain its own visually favourable colour range.
temperature_range <- range(pts$temp)
temperature_colours <- hcl.colors(50, "RdYlBu", rev = TRUE)
xy <- st_coordinates(pts)


# -----------------------------------------------------------------------------
# FIGURE 1: stations and environmental raster information
# -----------------------------------------------------------------------------
#
# The same stations are shown on both rasters. This makes clear where the
# environmental values used by the models originate. The labels contain both
# channel number and measured air temperature; the raster colours represent
# elevation or vegetation height, not temperature.

par(
  mfrow = c(1, 2),
  mar = c(3, 3, 4, 4),
  oma = c(0, 0, 3, 0)
)

environmental_rasters <- c(dem, chm)
environmental_titles <- c(
  "DEM | Terrain elevation (m)",
  "CHM = DSM − DEM | Canopy height (m)"
)

for (i in seq_len(nlyr(environmental_rasters))) {
  plot(
    environmental_rasters[[i]],
    col = hcl.colors(50, "Terrain"),
    main = environmental_titles[i]
  )
  
  # White station symbols remain legible on both raster backgrounds.
  points(pts, pch = 21, bg = "white", cex = 1.2)
  
  text(
    xy,
    labels = paste0(
      "CH", pts$air_channel, "\n",
      round(pts$temp, 1), " °C"
    ),
    pos = 3,
    cex = 0.65
  )
}

title(
  paste0(
    "Panel 1 | Environmental predictors and measurement stations — ",
    format(target_time, "%d.%m.%Y %H:%M")
  ),
  outer = TRUE,
  cex.main = 1.3
)


# -----------------------------------------------------------------------------
# FIGURE 2: six spatial statements from the same eight measurements
# -----------------------------------------------------------------------------
#
# The panels isolate three different questions:
#
#   1. How should observations be transferred by geographic proximity?
#      -> Voronoi and IDW
#
#   2. Which environmental variable is substantively meaningful?
#      -> linear model with terrain elevation versus linear model with CHM
#
#   3. Does CHM add useful information to a flexible model?
#      -> RF with position + DEM versus RF with position + DEM + CHM
#
# MAE and RMSE are leave-one-out errors:
#   MAE  = typical absolute prediction error at an omitted station.
#   RMSE = gives additional weight to individual large errors.
#
# Because only eight stations are available, these values compare the models
# for this network and timestamp. They do not prove universal model validity.

par(
  mfrow = c(2, 3),
  mar = c(3, 3, 4, 3),
  oma = c(3, 0, 3, 0)
)

map_titles <- c(
  "1 | Voronoi: nearest station",
  "2 | IDW: geographic proximity",
  "3 | LM: terrain elevation",
  "4 | LM: canopy height",
  "5 | RF: position + terrain elevation",
  "6 | RF: position + terrain elevation + CHM"
)

for (i in seq_len(nlyr(maps))) {
  plot(
    maps[[i]],
    range = temperature_range,
    col = temperature_colours,
    
    # Only the final map receives a legend. Since every map uses the same
    # range and palette, this is the common legend for all six panels.
    legend = i == nlyr(maps),
    
    main = paste0(
      map_titles[i],
      "\nLOOCV: MAE ", round(model_error$MAE[i], 2),
      " °C | RMSE ", round(model_error$RMSE[i], 2), " °C"
    )
  )
  
  # Station symbols show where observations actually support the surface.
  points(pts, pch = 21, bg = "white", cex = 1)
  text(
    xy,
    labels = paste0("CH", pts$air_channel),
    pos = 3,
    cex = 0.65
  )
}

title(
  "Panel 2 | The same 8 measurements — six different spatial statements",
  outer = TRUE,
  cex.main = 1.3
)

# This footer makes the abbreviated validation measures readable without
# requiring a separate explanatory slide.
mtext(
  paste0(
    "LOOCV: leave out one station, fit with the other seven, and predict the omitted station.  ",
    "MAE: mean absolute error.  RMSE: large errors receive more weight."
  ),
  side = 1,
  outer = TRUE,
  line = 1,
  cex = 0.8
)

# Return to the standard one-panel graphics layout.
par(mfrow = c(1, 1))


# =============================================================================
# Reading the result
# =============================================================================
# Voronoi: which station is nearest?
# IDW:     how does distance transfer the measurements?
# LM altitude: can one terrain-elevation relation explain the pattern?
# LM CHM:      can one canopy-height relation explain the pattern?
# RF DEM:      learns from position (x, y) and terrain elevation.
# RF DEM+CHM:  additionally learns from canopy height.
#
# The maps are not four measurements. They are four different statements made
# from the same eight measurements. CHM adds information, but not new temperature
# observations. The lowest RMSE is not automatically the most defensible map.
# =============================================================================
