# =============================================================================
# Aktuelle Mikroklimadaten: räumliche Auswertung und Interpolation
# =============================================================================
#
# Eingang:
#   MobilePolsterhaus_spatial_hourly.rds
#
# Das RDS ist ein sf-Objekt mit einer Zeile je Messpunkt und Stunde. Verwendet
# werden insbesondere:
#   time          Zeitpunkt
#   air_channel   Nummer des Luftkanals 1 bis 8
#   air_temp      stündlich aggregierte Lufttemperatur in °C
#   geometry      Lage des Messpunkts
#
# Das Skript:
#   1. wählt einen gemeinsamen Messzeitpunkt,
#   2. beschreibt die räumlichen Temperaturunterschiede,
#   3. interpoliert Voronoi und IDW,
#   4. prüft beide Verfahren durch Leave-one-out-Cross-Validation,
#   5. zeigt die räumliche Stützung der Vorhersagen,
#   6. verwendet bei vorhandenem DEM zusätzlich LM Höhe und Random Forest.
#
# Wichtige Einschränkung:
#   Acht Stationen liefern eine räumliche Momentaufnahme, aber keine belastbare
#   hochauflösende Rekonstruktion des gesamten Mikroklimafeldes. Rasterzellen
#   sind Recheneinheiten und keine zusätzliche Beobachtungsinformation.
# =============================================================================


# =============================================================================
# 1. EINSTELLUNGEN
# =============================================================================

# Aktuelles räumliches Stunden-RDS.
data_file <- "data/MobilePolsterhaus_spatial_hourly.rds"

# Optionales Höhenmodell. Wenn diese Datei nicht vorhanden ist oder NULL
# eingetragen wird, werden nur Voronoi und IDW berechnet.
# Beispiel: dem_file <- "DEM1.tif"
dem_file <- "data/DEM1_d.tif"

# Gewünschter Zeitpunkt:
#   NULL  = automatisch jüngster Zeitpunkt mit der maximal verfügbaren Zahl
#           gültiger Luftkanäle
#   oder  = konkreter Zeitstempel, beispielsweise
#           as.POSIXct("2026-09-03 14:00:00", tz = "Europe/Berlin")
target_time <- NULL

# Breite des Puffers um die konvexe Hülle der Stationen. Außerhalb dieses
# kleinen Berichtsraums werden keine Kartenwerte ausgegeben.
buffer_m <- 50

# Rasterweite der Ergebnisdarstellung. Eine feinere Rasterweite erzeugt nur
# glattere Bilder, aber keine zusätzliche Messinformation.
grid_resolution_m <- 1

# Anzahl der nächsten Nachbarn für IDW.
idw_nmax <- 4L

# Potenz der Distanzgewichtung. 2 ist die übliche Standardeinstellung.
idw_power <- 2

# Ergebnisordner.
output_dir <- "interpolation_aktuelle_daten"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# 2. PAKETE
# =============================================================================
# Falls erforderlich:
# install.packages(c("sf", "terra", "gstat", "ggplot2", "FNN"))
#
# Nur wenn ein DEM zusammen mit dem RF-Vergleich verwendet werden soll:
# install.packages("randomForest")

library(sf)
library(terra)
library(gstat)
library(ggplot2)


# =============================================================================
# 3. HILFSFUNKTIONEN
# =============================================================================

# Root Mean Squared Error: typische Größe der Vorhersagefehler, wobei größere
# Fehler stärker gewichtet werden.
rmse_fun <- function(observed, predicted) {
  sqrt(mean((observed - predicted)^2, na.rm = TRUE))
}

# Mean Absolute Error: mittlere absolute Abweichung zwischen Beobachtung und
# Vorhersage.
mae_fun <- function(observed, predicted) {
  mean(abs(observed - predicted), na.rm = TRUE)
}

# Erzeugt aus Vorhersagen in der Reihenfolge des Grids wieder ein SpatRaster.
make_map <- function(prediction, map_name, template, grid_table) {
  result <- template
  terra::values(result) <- NA_real_
  terra::values(result)[grid_table$cell] <- prediction
  names(result) <- map_name
  result
}

# Automatische metrische UTM-Projektion, falls kein DEM das Ziel-KBS vorgibt.
# Das GeoPackage liegt in EPSG:4326; Distanzen dürfen nicht in Grad berechnet
# werden. Für das kleine Untersuchungsgebiet reicht eine lokale UTM-Zone.
local_utm_crs <- function(points) {
  center <- sf::st_coordinates(sf::st_centroid(sf::st_union(points)))
  zone <- floor((center[1] + 180) / 6) + 1
  epsg <- if (center[2] >= 0) 32600 + zone else 32700 + zone
  sf::st_crs(epsg)
}

# Schreibt eine kompakte Matrix als CSV.
write_matrix <- function(x, filename) {
  utils::write.csv(x, file.path(output_dir, filename), row.names = TRUE)
}


# =============================================================================
# 4. AKTUELLE DATEN EINLESEN UND PRÜFEN
# =============================================================================

if (!file.exists(data_file)) {
  stop("RDS-Datei nicht gefunden: ", data_file)
}

all_data <- readRDS(data_file)

if (!inherits(all_data, "sf")) {
  stop("Das RDS muss ein räumliches sf-Objekt enthalten.")
}

required_columns <- c("time", "air_channel", "air_temp")
missing_columns <- setdiff(required_columns, names(all_data))

if (length(missing_columns) > 0) {
  stop(
    "Im aktuellen RDS fehlen folgende Spalten: ",
    paste(missing_columns, collapse = ", ")
  )
}

all_data$time <- as.POSIXct(all_data$time, tz = "Europe/Berlin")
all_data$air_channel <- as.integer(all_data$air_channel)
all_data$air_temp <- as.numeric(all_data$air_temp)


# =============================================================================
# 5. GEMEINSAMEN ZEITPUNKT AUSWÄHLEN
# =============================================================================
# Für eine räumliche Interpolation müssen die Stationswerte denselben Zeitpunkt
# repräsentieren. Bei target_time = NULL wird der jüngste Zeitpunkt gewählt,
# an dem die maximal im Datensatz vorhandene Zahl gültiger Kanäle vorliegt.

valid <- !is.na(all_data$time) &
  !is.na(all_data$air_channel) &
  !is.na(all_data$air_temp)

time_counts <- aggregate(
  all_data$air_channel[valid],
  by = list(time = all_data$time[valid]),
  FUN = function(x) length(unique(x))
)
names(time_counts)[2] <- "n_channels"

if (nrow(time_counts) == 0) {
  stop("Im RDS wurden keine gültigen Lufttemperaturen gefunden.")
}

if (is.null(target_time)) {
  best_times <- time_counts$time[
    time_counts$n_channels == max(time_counts$n_channels)
  ]
  target_time <- max(best_times)
} else {
  target_time <- as.POSIXct(target_time, tz = "Europe/Berlin")
}

points_time <- all_data[
  all_data$time == target_time & !is.na(all_data$air_temp),
]

# Sicherheit gegen doppelte Messwerte desselben Kanals am selben Zeitpunkt.
# Falls Duplikate existieren, wird je Kanal der erste räumliche Datensatz
# verwendet und der Temperaturmittelwert eingesetzt.
if (anyDuplicated(points_time$air_channel)) {
  split_points <- split(points_time, points_time$air_channel)
  points_time <- do.call(
    rbind,
    lapply(split_points, function(x) {
      x$air_temp[1] <- mean(x$air_temp, na.rm = TRUE)
      x[1, ]
    })
  )
}

if (nrow(points_time) < 5L) {
  stop(
    "Am gewählten Zeitpunkt liegen nur ", nrow(points_time),
    " gültige Stationen vor. Mindestens fünf werden benötigt."
  )
}

target_label <- format(target_time, "%Y-%m-%d %H:%M %Z")
target_file_label <- format(target_time, "%Y%m%d_%H%M")

message(
  "Gewählter Zeitpunkt: ", target_label,
  " | gültige Kanäle: ", nrow(points_time)
)


# =============================================================================
# 6. ZIEL-KBS UND OPTIONAL DAS DEM VORBEREITEN
# =============================================================================

use_dem <- !is.null(dem_file) && file.exists(dem_file)

if (!is.null(dem_file) && !file.exists(dem_file)) {
  warning(
    "Das angegebene DEM wurde nicht gefunden. ",
    "Es werden nur Voronoi und IDW berechnet: ", dem_file
  )
}

if (use_dem) {
  dem_original <- terra::rast(dem_file)
  names(dem_original) <- "altitude"
  target_crs <- sf::st_crs(terra::crs(dem_original))

  if (is.na(target_crs)) {
    stop("Das DEM besitzt kein gültiges Koordinatenreferenzsystem.")
  }
} else {
  target_crs <- local_utm_crs(points_time)
}

pts <- sf::st_transform(points_time, target_crs)

if (sf::st_is_longlat(pts)) {
  stop("Das Ziel-KBS muss metrisch/projiziert sein; aktuell ist es geografisch.")
}

pts$temp <- pts$air_temp
pts$channel_label <- paste0("CH", pts$air_channel)


# =============================================================================
# 7. KLEINEN MODELL- UND BERICHTSRAUM FESTLEGEN
# =============================================================================
# Die konvexe Hülle verbindet die äußersten Stationen. Der kleine Puffer macht
# den Rand sichtbar, ohne großflächig außerhalb des Netzes zu extrapolieren.

area <- sf::st_sf(
  geometry = sf::st_buffer(
    sf::st_convex_hull(sf::st_union(sf::st_geometry(pts))),
    dist = buffer_m
  )
)


# =============================================================================
# 8. VORHERSAGERASTER ERZEUGEN
# =============================================================================
# Mit DEM wird dessen Rasterstruktur verwendet. Ohne DEM wird ein neutrales
# Raster mit der eingestellten Zellweite erzeugt. Voronoi und IDW benötigen
# keinen Höhenprädiktor.

if (use_dem) {
  dem <- terra::crop(dem_original, terra::vect(area))
  dem <- terra::mask(dem, terra::vect(area))
  names(dem) <- "altitude"
  template <- dem
} else {
  template <- terra::rast(
    terra::vect(area),
    resolution = grid_resolution_m,
    crs = sf::st_crs(pts)$wkt
  )
  terra::values(template) <- 1
  template <- terra::mask(template, terra::vect(area))
  names(template) <- "domain"
}

grid <- as.data.frame(
  template,
  xy = TRUE,
  cells = TRUE,
  na.rm = TRUE
)

# Die erste vier Spalten sind cell, x, y und Rasterwert.
names(grid)[4] <- if (use_dem) "altitude" else "domain"

grid_sf <- sf::st_as_sf(
  grid,
  coords = c("x", "y"),
  crs = sf::st_crs(pts),
  remove = FALSE
)


# =============================================================================
# 9. DESKRIPTIVE RÄUMLICHE MOMENTAUFNAHME
# =============================================================================
# Die Abweichung bezieht sich nur auf das Mittel der acht Stationen am
# ausgewählten Zeitpunkt. Sie ist keine klimatologische Langzeitanomalie.

network_mean <- mean(pts$temp, na.rm = TRUE)
pts$deviation_from_network_mean <- pts$temp - network_mean

snapshot_summary <- sf::st_drop_geometry(pts)[
  , c(
    "air_channel",
    "channel_label",
    "temp",
    "deviation_from_network_mean"
  )
]
snapshot_summary <- snapshot_summary[
  order(snapshot_summary$air_channel),
]

utils::write.csv(
  snapshot_summary,
  file.path(
    output_dir,
    paste0("stationswerte_", target_file_label, ".csv")
  ),
  row.names = FALSE
)


# =============================================================================
# 10. MODELL 1: VORONOI / NÄCHSTE STATION
# =============================================================================
# Jede Rasterzelle erhält den Wert der nächstgelegenen Station. Das erzeugt
# harte Stationsgebiete und enthält keine Glättung.

vor_prediction <- gstat::idw(
  temp ~ 1,
  locations = pts,
  newdata = grid_sf,
  nmax = 1,
  idp = idw_power
)

map_voronoi <- make_map(
  vor_prediction$var1.pred,
  "Voronoi",
  template,
  grid
)


# =============================================================================
# 11. MODELL 2: IDW
# =============================================================================
# Nahe Stationen erhalten ein größeres Gewicht als weiter entfernte Stationen.
# nmax begrenzt die Berechnung auf die nächsten vier Stationen.

idw_prediction <- gstat::idw(
  temp ~ 1,
  locations = pts,
  newdata = grid_sf,
  nmax = min(idw_nmax, nrow(pts)),
  idp = idw_power
)

map_idw <- make_map(
  idw_prediction$var1.pred,
  "IDW",
  template,
  grid
)


# =============================================================================
# 12. OPTIONALE HÖHENBASIERTE MODELLE
# =============================================================================
# LM und Random Forest werden nur berechnet, wenn ein DEM vorhanden ist.
# Bei acht Stationen ist besonders Random Forest lediglich ein Warn- und
# Vergleichsmodell, keine belastbare Empfehlung.

map_lm <- NULL
map_rf <- NULL
lm_cv <- rep(NA_real_, nrow(pts))
rf_cv <- rep(NA_real_, nrow(pts))

if (use_dem) {
  pts$altitude <- terra::extract(
    dem_original,
    terra::vect(pts)
  )$altitude

  if (anyNA(pts$altitude)) {
    stop("Für mindestens eine Station konnte keine DEM-Höhe extrahiert werden.")
  }

  fit_lm <- stats::lm(
    temp ~ altitude,
    data = sf::st_drop_geometry(pts)
  )

  map_lm <- terra::predict(dem, fit_lm)
  names(map_lm) <- "LM_altitude"

  if (!requireNamespace("randomForest", quietly = TRUE)) {
    warning(
      "Paket 'randomForest' fehlt. Das RF-Vergleichsmodell wird übersprungen."
    )
  } else {
    point_xy <- sf::st_coordinates(pts)
    pts$x <- point_xy[, 1]
    pts$y <- point_xy[, 2]

    fit_rf <- randomForest::randomForest(
      temp ~ x + y + altitude,
      data = sf::st_drop_geometry(pts),
      ntree = 500,
      importance = TRUE
    )

    rf_prediction <- stats::predict(
      fit_rf,
      newdata = grid[, c("x", "y", "altitude")]
    )

    map_rf <- make_map(
      rf_prediction,
      "RF_warning",
      template,
      grid
    )
  }
}


# =============================================================================
# 13. LEAVE-ONE-OUT-CROSS-VALIDATION FÜR VORONOI UND IDW
# =============================================================================
# Jeweils eine Station wird zurückgehalten und aus den übrigen Stationen
# vorhergesagt. Das prüft die Rückvorhersage an Stationsstandorten, nicht die
# flächenhafte Wahrheit zwischen den Stationen.

vor_model <- gstat::gstat(
  formula = temp ~ 1,
  locations = pts,
  nmax = 1,
  set = list(idp = idw_power)
)

vor_cv_result <- gstat::gstat.cv(
  vor_model,
  nfold = nrow(pts),
  verbose = FALSE
)

idw_model <- gstat::gstat(
  formula = temp ~ 1,
  locations = pts,
  nmax = min(idw_nmax, nrow(pts) - 1L),
  set = list(idp = idw_power)
)

idw_cv_result <- gstat::gstat.cv(
  idw_model,
  nfold = nrow(pts),
  verbose = FALSE
)

vor_cv <- vor_cv_result$observed - vor_cv_result$residual
idw_cv <- idw_cv_result$observed - idw_cv_result$residual


# LOOCV der optionalen Modelle.
if (use_dem) {
  for (i in seq_len(nrow(pts))) {
    train <- pts[-i, ]
    test <- pts[i, ]

    fit_lm_i <- stats::lm(
      temp ~ altitude,
      data = sf::st_drop_geometry(train)
    )
    lm_cv[i] <- stats::predict(
      fit_lm_i,
      newdata = sf::st_drop_geometry(test)
    )

    if (requireNamespace("randomForest", quietly = TRUE)) {
      fit_rf_i <- randomForest::randomForest(
        temp ~ x + y + altitude,
        data = sf::st_drop_geometry(train),
        ntree = 500
      )
      rf_cv[i] <- stats::predict(
        fit_rf_i,
        newdata = sf::st_drop_geometry(test)
      )
    }
  }
}


# =============================================================================
# 14. FEHLERKENNWERTE UND STATIONSEINZELFEHLER
# =============================================================================

model_names <- c("Voronoi", "IDW")
predictions <- list(vor_cv, idw_cv)

if (use_dem) {
  model_names <- c(model_names, "LM Höhe")
  predictions <- c(predictions, list(lm_cv))

  if (!all(is.na(rf_cv))) {
    model_names <- c(model_names, "RF Warnmodell")
    predictions <- c(predictions, list(rf_cv))
  }
}

validation_summary <- data.frame(
  model = model_names,
  RMSE_C = vapply(
    predictions,
    function(x) rmse_fun(pts$temp, x),
    numeric(1)
  ),
  MAE_C = vapply(
    predictions,
    function(x) mae_fun(pts$temp, x),
    numeric(1)
  )
)

validation_details <- data.frame(
  channel = pts$channel_label,
  observed_C = pts$temp,
  predicted_Voronoi_C = vor_cv,
  error_Voronoi_C = pts$temp - vor_cv,
  predicted_IDW_C = idw_cv,
  error_IDW_C = pts$temp - idw_cv
)

if (use_dem) {
  validation_details$predicted_LM_C <- lm_cv
  validation_details$error_LM_C <- pts$temp - lm_cv

  if (!all(is.na(rf_cv))) {
    validation_details$predicted_RF_C <- rf_cv
    validation_details$error_RF_C <- pts$temp - rf_cv
  }
}

utils::write.csv(
  validation_summary,
  file.path(
    output_dir,
    paste0("modellguete_", target_file_label, ".csv")
  ),
  row.names = FALSE
)

utils::write.csv(
  validation_details,
  file.path(
    output_dir,
    paste0("loocv_einzelfehler_", target_file_label, ".csv")
  ),
  row.names = FALSE
)


# =============================================================================
# 15. GEOGRAFISCHE STÜTZUNG DER INTERPOLATION
# =============================================================================
# Für jede Rasterzelle werden die Distanz zur nächsten und zur viertnächsten
# Station berechnet. Als Referenz dienen die Distanzen, die beim Zurückhalten
# einer Station in der LOOCV tatsächlich auftreten.
#
# Verhältnis <= 1:
#   Distanz liegt innerhalb des 95-%-Bereichs der CV-Situationen.
# Verhältnis > 1:
#   räumlich weiter reichende Übertragung als in fast allen CV-Situationen.

if (!requireNamespace("FNN", quietly = TRUE)) {
  stop("Für die Stützungsdiagnostik fehlt das Paket 'FNN'.")
}

station_xy <- sf::st_coordinates(pts)
grid_xy <- as.matrix(grid[, c("x", "y")])

station_distances <- as.matrix(sf::st_distance(pts))
diag(station_distances) <- Inf
sorted_station_distances <- t(apply(station_distances, 1, sort))

cv_d1 <- sorted_station_distances[, 1]
k_idw <- min(idw_nmax, nrow(pts) - 1L)
cv_dk <- sorted_station_distances[, k_idw]

grid_neighbors <- FNN::get.knnx(
  data = station_xy,
  query = grid_xy,
  k = k_idw
)

grid_d1 <- grid_neighbors$nn.dist[, 1]
grid_dk <- grid_neighbors$nn.dist[, k_idw]

threshold_d1 <- as.numeric(stats::quantile(cv_d1, 0.95, na.rm = TRUE))
threshold_dk <- as.numeric(stats::quantile(cv_dk, 0.95, na.rm = TRUE))

map_support_voronoi <- make_map(
  grid_d1 / threshold_d1,
  "Voronoi_support_ratio",
  template,
  grid
)

map_support_idw <- make_map(
  grid_dk / threshold_dk,
  "IDW_support_ratio",
  template,
  grid
)

support_summary <- data.frame(
  method = c("Voronoi", paste0("IDW, k = ", k_idw)),
  distance_reference = c("nächste Station", paste0(k_idw, ". Station")),
  cv_q95_distance_m = c(threshold_d1, threshold_dk),
  supported_area_percent = 100 * c(
    mean(grid_d1 <= threshold_d1),
    mean(grid_dk <= threshold_dk)
  )
)

utils::write.csv(
  support_summary,
  file.path(
    output_dir,
    paste0("raeumliche_stuetzung_", target_file_label, ".csv")
  ),
  row.names = FALSE
)


# =============================================================================
# 16. EMPIRISCHES VARIOGRAMM ALS DESKRIPTIVE DIAGNOSE
# =============================================================================
# Bei nur acht Stationen wird bewusst kein theoretisches Variogrammmodell
# erzwungen und kein Kriging durchgeführt. Das empirische Variogramm zeigt nur,
# wie die beobachteten Temperaturunterschiede mit der Distanz zusammenhängen.

empirical_variogram <- gstat::variogram(temp ~ 1, data = pts)


# =============================================================================
# 17. ERGEBNISRASTER SPEICHERN
# =============================================================================

terra::writeRaster(
  map_voronoi,
  file.path(
    output_dir,
    paste0("voronoi_", target_file_label, ".tif")
  ),
  overwrite = TRUE
)

terra::writeRaster(
  map_idw,
  file.path(
    output_dir,
    paste0("idw_", target_file_label, ".tif")
  ),
  overwrite = TRUE
)

terra::writeRaster(
  map_support_voronoi,
  file.path(
    output_dir,
    paste0("support_voronoi_", target_file_label, ".tif")
  ),
  overwrite = TRUE
)

terra::writeRaster(
  map_support_idw,
  file.path(
    output_dir,
    paste0("support_idw_", target_file_label, ".tif")
  ),
  overwrite = TRUE
)

if (!is.null(map_lm)) {
  terra::writeRaster(
    map_lm,
    file.path(
      output_dir,
      paste0("lm_hoehe_", target_file_label, ".tif")
    ),
    overwrite = TRUE
  )
}

if (!is.null(map_rf)) {
  terra::writeRaster(
    map_rf,
    file.path(
      output_dir,
      paste0("rf_warnmodell_", target_file_label, ".tif")
    ),
    overwrite = TRUE
  )
}


# =============================================================================
# 18. BESCHRIFTETE GRAFIK: BEOBACHTETE STATIONSWERTE
# =============================================================================

points_plot <- data.frame(
  x = sf::st_coordinates(pts)[, 1],
  y = sf::st_coordinates(pts)[, 2],
  channel = pts$channel_label,
  temperature = pts$temp,
  deviation = pts$deviation_from_network_mean
)

p_points <- ggplot(points_plot, aes(x = x, y = y)) +
  geom_point(
    aes(colour = temperature),
    size = 5
  ) +
  geom_text(
    aes(label = paste0(channel, "\n", sprintf("%.1f °C", temperature))),
    nudge_y = grid_resolution_m * 4,
    size = 3.5,
    fontface = "bold"
  ) +
  scale_colour_gradient2(
    low = "#2166AC",
    mid = "#FFFFBF",
    high = "#B2182B",
    midpoint = network_mean,
    name = "Temperatur (°C)"
  ) +
  coord_equal() +
  labs(
    title = "Beobachtete Lufttemperatur an den acht Messpunkten",
    subtitle = paste0(
      target_label,
      " | Stationsmittel = ", sprintf("%.1f °C", network_mean)
    ),
    x = "Easting (m)",
    y = "Northing (m)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(p_points)

ggsave(
  file.path(
    output_dir,
    paste0("01_stationswerte_", target_file_label, ".png")
  ),
  p_points,
  width = 8,
  height = 7,
  dpi = 300
)


# =============================================================================
# 19. BESCHRIFTETE GRAFIK: INTERPOLATIONSERGEBNISSE
# =============================================================================
# Alle Karten erhalten dieselbe Temperaturskala. Dadurch können Unterschiede
# zwischen den Verfahren nicht durch voneinander abweichende Farbskalen
# kaschiert oder übertrieben werden.

maps <- c(map_voronoi, map_idw)

if (!is.null(map_lm)) maps <- c(maps, map_lm)
if (!is.null(map_rf)) maps <- c(maps, map_rf)

map_names <- names(maps)
common_range <- range(
  c(pts$temp, terra::minmax(maps)),
  na.rm = TRUE
)

png(
  file.path(
    output_dir,
    paste0("02_interpolationsvergleich_", target_file_label, ".png")
  ),
  width = 2200,
  height = if (terra::nlyr(maps) <= 2) 1100 else 2000,
  res = 220
)

par(
  mfrow = if (terra::nlyr(maps) <= 2) c(1, 2) else c(2, 2),
  mar = c(3, 3, 4, 5)
)

for (i in seq_len(terra::nlyr(maps))) {
  plot(
    maps[[i]],
    range = common_range,
    main = paste0(map_names[i], "\n", target_label),
    axes = TRUE
  )
  points(terra::vect(pts), pch = 21, bg = "white", col = "black", cex = 1.1)
  text(
    sf::st_coordinates(pts),
    labels = paste0(pts$channel_label, ": ", sprintf("%.1f", pts$temp)),
    pos = 3,
    cex = 0.75
  )
}

par(mfrow = c(1, 1))
dev.off()


# =============================================================================
# 20. BESCHRIFTETE GRAFIK: MODELLGÜTE
# =============================================================================

validation_long <- reshape(
  validation_summary,
  varying = c("RMSE_C", "MAE_C"),
  v.names = "error_C",
  timevar = "metric",
  times = c("RMSE", "MAE"),
  direction = "long"
)

p_validation <- ggplot(
  validation_long,
  aes(x = model, y = error_C, fill = metric)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = sprintf("%.2f", error_C)),
    position = position_dodge(width = 0.8),
    vjust = -0.4,
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(RMSE = "#B2182B", MAE = "#2166AC"),
    name = "Fehlermaß"
  ) +
  expand_limits(y = max(validation_long$error_C, na.rm = TRUE) * 1.18) +
  labs(
    title = "Leave-one-out-Fehler der Interpolationsverfahren",
    subtitle = paste0(
      target_label,
      " | kleinere Werte bedeuten bessere Stations-Rückvorhersage"
    ),
    x = "Verfahren",
    y = "Fehler (°C)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(p_validation)

ggsave(
  file.path(
    output_dir,
    paste0("03_modellguete_", target_file_label, ".png")
  ),
  p_validation,
  width = 8,
  height = 5.5,
  dpi = 300
)


# =============================================================================
# 21. BESCHRIFTETE GRAFIKEN: STÜTZUNG UND VARIOGRAMM
# =============================================================================

png(
  file.path(
    output_dir,
    paste0("04_raeumliche_stuetzung_", target_file_label, ".png")
  ),
  width = 2200,
  height = 1100,
  res = 220
)

par(mfrow = c(1, 2), mar = c(3, 3, 4, 5))

plot(
  map_support_voronoi,
  main = "Voronoi-Stützung\nDistanz zur nächsten Station / CV-q95"
)
contour(
  map_support_voronoi,
  levels = 1,
  add = TRUE,
  drawlabels = FALSE,
  lwd = 2
)
points(terra::vect(pts), pch = 21, bg = "white", col = "black")

plot(
  map_support_idw,
  main = paste0(
    "IDW-Stützung\nDistanz zur ", k_idw,
    ". Station / CV-q95"
  )
)
contour(
  map_support_idw,
  levels = 1,
  add = TRUE,
  drawlabels = FALSE,
  lwd = 2
)
points(terra::vect(pts), pch = 21, bg = "white", col = "black")

par(mfrow = c(1, 1))
dev.off()


png(
  file.path(
    output_dir,
    paste0("05_empirisches_variogramm_", target_file_label, ".png")
  ),
  width = 1500,
  height = 1100,
  res = 220
)

plot(
  empirical_variogram$dist,
  empirical_variogram$gamma,
  type = "b",
  pch = 19,
  xlab = "Distanz (m)",
  ylab = "Semivarianz (°C²)",
  main = paste0("Empirisches Variogramm\n", target_label)
)

dev.off()


# =============================================================================
# 22. ERGEBNISSE IN DER KONSOLE
# =============================================================================

cat("\n============================================================\n")
cat("RÄUMLICHE AUSWERTUNG ABGESCHLOSSEN\n")
cat("============================================================\n")
cat("Zeitpunkt: ", target_label, "\n", sep = "")
cat("Stationen: ", nrow(pts), "\n", sep = "")
cat("Stationsmittel: ", sprintf("%.2f °C", network_mean), "\n", sep = "")
cat("Temperaturspanne: ", sprintf(
  "%.2f °C",
  max(pts$temp) - min(pts$temp)
), "\n", sep = "")
cat("DEM verwendet: ", ifelse(use_dem, "ja", "nein"), "\n", sep = "")
cat("Ergebnisordner: ", normalizePath(output_dir), "\n\n", sep = "")

cat("Stationswerte und Abweichungen:\n")
print(snapshot_summary)

cat("\nModellgüte aus Leave-one-out-Cross-Validation:\n")
print(validation_summary)

cat("\nRäumliche Stützung:\n")
print(support_summary)

cat(
  "\nInterpretation:\n",
  "- Voronoi zeigt harte Zuständigkeitsbereiche der nächsten Station.\n",
  "- IDW erzeugt eine geglättete, distanzgewichtete Fläche.\n",
  "- RMSE und MAE bewerten nur Rückvorhersagen an Stationsstandorten.\n",
  "- Ein Stützungsverhältnis über 1 kennzeichnet weiter reichende ",
  "Übertragung als in 95 % der LOOCV-Distanzsituationen.\n",
  "- Bei acht Stationen wird kein belastbares Variogrammmodell und kein ",
  "Kriging erzwungen.\n",
  sep = ""
)


# =============================================================================
# ENDE
# =============================================================================
