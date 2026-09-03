#' Excel-Export einer Ecowitt/Froggit-Station aggregieren und als RDS speichern
#'
#' Eingelesen werden:
#' - Outdoor Temperature
#' - Solar (fc)
#' - Temperature und Humidity der Kanäle 1 bis 8
#' - Temp CH1 bis Temp CH8
#'
#' @param input_file Pfad zur XLSX-Datei.
#' @param output_file Pfad der zu erzeugenden RDS-Datei. Standardmäßig wird
#'   neben der XLSX-Datei eine Datei mit der Endung "_aggregated.rds" angelegt.
#' @param interval Aggregationsintervall, z. B. "1 hour" (Standard),
#'   "30 mins", "10 minutes" oder "1 day". Alternativ Sekunden als Zahl.
#' @param agg_fun Aggregationsfunktion als Name oder Funktion; Standard "mean".
#' @param tz Zeitzone der Zeitstempel; Standard "Europe/Berlin".
#' @param sheet Tabellenblatt; Standard "result_list".
#' @param gpkg_file Optionaler Pfad zum GeoPackage mit den Messpunkten. Wenn
#'   angegeben, werden Klima- und Standortdaten räumlich zusammengeführt.
#' @param gpkg_layer Layername im GeoPackage; Standard "Messpunkte".
#' @param channel_col Spalte mit der Zuordnung zu Temp and Humidity CH1–CH8;
#'   Standard "channel-pos".
#'
#' @return Unsichtbar ein data.frame mit den aggregierten Daten. Zusätzlich
#'   wird dieses data.frame mit saveRDS() gespeichert.
#' @export
ecowitt_excel_to_rds <- function(
    input_file,
    output_file = NULL,
    interval = "1 hour",
    agg_fun = "mean",
    tz = "Europe/Berlin",
    sheet = "result_list",
    gpkg_file = NULL,
    gpkg_layer = "Messpunkte",
    channel_col = "channel-pos") {

  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Das Paket 'readxl' fehlt. Installation: install.packages('readxl')")
  }
  if (!file.exists(input_file)) {
    stop("Eingabedatei nicht gefunden: ", input_file)
  }

  interval_seconds <- parse_interval_seconds(interval)
  fun <- match.fun(agg_fun)

  # Der Stations-Export besitzt zwei Kopfzeilen. Alles wird zunächst als Text
  # gelesen, damit gemischte Kopf- und Datenzeilen keine Spaltentypen verfälschen.
  raw <- readxl::read_excel(
    input_file,
    sheet = sheet,
    col_names = FALSE,
    col_types = "text",
    .name_repair = "minimal"
  )
  raw <- as.data.frame(raw, stringsAsFactors = FALSE, check.names = FALSE)

  if (nrow(raw) < 3L) {
    stop("Das Tabellenblatt enthält keine auswertbaren Daten.")
  }

  group_header <- trimws(as.character(unlist(raw[1, ], use.names = FALSE)))
  sub_header <- trimws(as.character(unlist(raw[2, ], use.names = FALSE)))

  # Leere Zellen der ersten Kopfzeile gehören jeweils zur vorherigen Gruppe.
  for (i in seq_along(group_header)) {
    if ((is.na(group_header[i]) || group_header[i] == "") && i > 1L) {
      group_header[i] <- group_header[i - 1L]
    }
  }

  find_one <- function(group, sub = NULL) {
    hit <- which(group_header == group)
    if (!is.null(sub)) hit <- hit[sub_header[hit] == sub]
    if (length(hit) != 1L) {
      stop(
        "Spalte nicht eindeutig gefunden: ", group,
        if (!is.null(sub)) paste0(" / ", sub) else ""
      )
    }
    hit
  }

  time_col <- find_one("Time")
  selected <- c(
    outdoor_temp = find_one("Outdoor", "Temperature(℃)"),
    solar_fc = find_one("Solar and UVI", "Solar(fc)")
  )

  for (ch in 1:8) {
    group <- paste0("Temp and Humidity CH", ch)
    selected[paste0("th_ch", ch, "_temp")] <- find_one(group, "Temperature(℃)")
    selected[paste0("th_ch", ch, "_humidity")] <- find_one(group, "Humidity(%)")
  }
  for (ch in 1:8) {
    selected[paste0("temp_ch", ch)] <- find_one(
      paste0("Temp CH", ch), "Temperature(℃)"
    )
  }

  dat <- raw[-c(1L, 2L), c(time_col, unname(selected)), drop = FALSE]
  names(dat) <- c("time", names(selected))

  dat$time <- parse_station_time(dat$time, tz = tz)
  for (nm in setdiff(names(dat), "time")) {
    dat[[nm]] <- suppressWarnings(as.numeric(gsub(",", ".", dat[[nm]], fixed = TRUE)))
  }
  dat <- dat[!is.na(dat$time), , drop = FALSE]
  if (!nrow(dat)) stop("Keine gültigen Zeitstempel gefunden.")

  # Zeitklassen beginnen jeweils auf der unteren Intervallgrenze.
  epoch <- as.numeric(dat$time)
  dat$time <- as.POSIXct(
    floor(epoch / interval_seconds) * interval_seconds,
    origin = "1970-01-01",
    tz = tz
  )

  apply_fun <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    value <- tryCatch(
      fun(x, na.rm = TRUE),
      error = function(e) fun(x[!is.na(x)])
    )
    value <- as.numeric(value)[1L]
    if (is.nan(value)) NA_real_ else value
  }

  value_names <- setdiff(names(dat), "time")
  aggregated <- stats::aggregate(
    dat[value_names],
    by = list(time = dat$time),
    FUN = apply_fun
  )
  aggregated <- aggregated[order(aggregated$time), , drop = FALSE]
  rownames(aggregated) <- NULL

  result <- aggregated
  if (!is.null(gpkg_file)) {
    result <- merge_climate_with_locations(
      climate = aggregated,
      gpkg_file = gpkg_file,
      gpkg_layer = gpkg_layer,
      channel_col = channel_col
    )
  }

  if (is.null(output_file)) {
    stem <- tools::file_path_sans_ext(basename(input_file))
    output_file <- file.path(dirname(input_file), paste0(stem, "_aggregated.rds"))
  }
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  saveRDS(result, output_file)

  message(
    "Gespeichert: ", normalizePath(output_file, mustWork = FALSE),
    " (", nrow(result), " Zeilen)"
  )
  invisible(result)
}


merge_climate_with_locations <- function(
    climate,
    gpkg_file,
    gpkg_layer = "Messpunkte",
    channel_col = "channel-pos") {

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Für das GeoPackage wird 'sf' benötigt: install.packages('sf')")
  }
  if (!file.exists(gpkg_file)) {
    stop("GeoPackage nicht gefunden: ", gpkg_file)
  }

  locations <- sf::st_read(gpkg_file, layer = gpkg_layer, quiet = TRUE)
  if (!channel_col %in% names(locations)) {
    repaired <- make.names(channel_col)
    if (repaired %in% names(locations)) {
      channel_col <- repaired
    } else {
      stop("Kanalspalte im GeoPackage nicht gefunden: ", channel_col)
    }
  }
  if (!"Description" %in% names(locations)) {
    stop("Im GeoPackage fehlt die Spalte 'Description'.")
  }

  channel <- suppressWarnings(as.integer(locations[[channel_col]]))
  if (anyNA(channel) || any(!channel %in% 1:8)) {
    stop("Die Kanalpositionen müssen ganzzahlig zwischen 1 und 8 liegen.")
  }

  # IDs der separaten Bodentemperatursensoren laut Sensor-Management.
  soil_id_to_channel <- c(
    B2A2 = 1L,
    B292 = 2L,
    B2B1 = 3L,
    B29A = 4L,
    B318 = 5L,
    B33E = 6L,
    B2EE = 7L,
    B280 = 8L
  )

  extract_id <- function(description, depth) {
    pattern <- paste0("\\b([[:xdigit:]]{4})\\s+", depth, "\\b")
    match <- regexec(pattern, description, ignore.case = TRUE)
    parts <- regmatches(description, match)
    vapply(parts, function(x) {
      if (length(x) >= 2L) toupper(x[2L]) else NA_character_
    }, character(1))
  }

  high_id <- extract_id(locations$Description, "hoch")
  deep_id <- extract_id(locations$Description, "tief")
  unknown_ids <- setdiff(
    unique(stats::na.omit(c(high_id, deep_id))),
    names(soil_id_to_channel)
  )
  if (length(unknown_ids)) {
    stop("Unbekannte Bodensensor-ID(s): ", paste(unknown_ids, collapse = ", "))
  }

  locations$air_channel <- channel
  locations$soil_high_id <- high_id
  locations$soil_deep_id <- deep_id
  locations$soil_high_channel <- unname(soil_id_to_channel[high_id])
  locations$soil_deep_channel <- unname(soil_id_to_channel[deep_id])

  n_time <- nrow(climate)
  expanded <- locations[rep(seq_len(nrow(locations)), each = n_time), ]
  climate_rows <- climate[rep(seq_len(n_time), times = nrow(locations)), ]

  expanded$time <- climate_rows$time
  expanded$outdoor_temp <- climate_rows$outdoor_temp
  expanded$solar_fc <- climate_rows$solar_fc

  get_channel_value <- function(prefix, channels) {
    out <- rep(NA_real_, length(channels))
    for (ch in 1:8) {
      use <- !is.na(channels) & channels == ch
      out[use] <- climate_rows[[paste0(prefix, ch)]][use]
    }
    out
  }

  # channel-pos bezeichnet die kombinierten Luftsensoren.
  expanded$air_temp <- rep(NA_real_, nrow(expanded))
  expanded$air_humidity <- rep(NA_real_, nrow(expanded))
  for (ch in 1:8) {
    use <- expanded$air_channel == ch
    expanded$air_temp[use] <- climate_rows[[paste0("th_ch", ch, "_temp")]][use]
    expanded$air_humidity[use] <- climate_rows[[paste0("th_ch", ch, "_humidity")]][use]
  }

  # IDs aus Description bestimmen die separaten Bodentemperaturkanäle.
  expanded$soil_high_temp <- get_channel_value(
    "temp_ch", expanded$soil_high_channel
  )
  expanded$soil_deep_temp <- get_channel_value(
    "temp_ch", expanded$soil_deep_channel
  )

  expanded
}


parse_interval_seconds <- function(interval) {
  if (is.numeric(interval) && length(interval) == 1L && interval > 0) {
    return(as.numeric(interval))
  }
  if (!is.character(interval) || length(interval) != 1L) {
    stop("'interval' muss z. B. '1 hour', '30 mins' oder Sekunden enthalten.")
  }

  x <- tolower(trimws(interval))
  parts <- regmatches(x, regexec(
    "^([0-9]+(?:[.][0-9]+)?)\\s*(s|sec|secs|second|seconds|min|mins|minute|minutes|h|hour|hours|d|day|days)$",
    x,
    perl = TRUE
  ))[[1]]
  if (!length(parts)) {
    stop("Unbekanntes Intervall: ", interval)
  }

  amount <- as.numeric(parts[2])
  unit <- parts[3]
  factor <- if (unit %in% c("s", "sec", "secs", "second", "seconds")) {
    1
  } else if (unit %in% c("min", "mins", "minute", "minutes")) {
    60
  } else if (unit %in% c("h", "hour", "hours")) {
    3600
  } else {
    86400
  }
  amount * factor
}


parse_station_time <- function(x, tz) {
  x <- trimws(as.character(x))
  out <- as.POSIXct(rep(NA_character_, length(x)), tz = tz)

  # Übliche Zeitformate des Ecowitt/Froggit-Exports.
  formats <- c(
    "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M",
    "%d.%m.%Y %H:%M:%S", "%d.%m.%Y %H:%M"
  )
  for (fmt in formats) {
    missing <- is.na(out)
    out[missing] <- as.POSIXct(x[missing], format = fmt, tz = tz)
  }

  # Falls Excel Zeitstempel als Seriennummern liefert.
  numeric_x <- suppressWarnings(as.numeric(x))
  excel <- is.na(out) & !is.na(numeric_x)
  out[excel] <- as.POSIXct(
    numeric_x[excel] * 86400,
    origin = "1899-12-30",
    tz = tz
  )
  out
}


# Beispiel:
# hourly <- ecowitt_excel_to_rds(
#   "all_MobilePolsterhaus_202609020000-202609022049.xlsx",
#   output_file = "MobilePolsterhaus_hourly.rds"
# )
#
# Für 30-Minuten-Mittelwerte:
# half_hourly <- ecowitt_excel_to_rds(
#   "all_MobilePolsterhaus_202609020000-202609022049.xlsx",
#   output_file = "MobilePolsterhaus_30min.rds",
#   interval = "30 mins"
# )
#
# Klima- und Messpunktdaten räumlich zusammenführen:
# spatial_hourly <- ecowitt_excel_to_rds(
#   "all_MobilePolsterhaus_202609020000-202609022049.xlsx",
#   output_file = "MobilePolsterhaus_spatial_hourly.rds",
#   interval = "1 hour",
#   gpkg_file = "Mikroklima_Stationen_EON25.gpkg"
# )
