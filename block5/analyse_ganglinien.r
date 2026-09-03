# ============================================================================
# Explorative Auswertung der Mikroklima-Ganglinien
# ============================================================================
# Ziel:
#   1. Zeigen, ob die Lufttemperatur-Ganglinien ähnlich verlaufen.
#   2. Zeigen, ob einzelne Kanäle systematisch wärmer oder kälter sind.
#   3. Die mittlere Größe der Unterschiede zwischen den Kanälen quantifizieren.
#   4. Oberflächennahe und tiefere Bodentemperatur vergleichen.
#
# Bewusst NICHT enthalten:
#   - Solarstrahlung
#   - Signifikanztests
#   - Aussagen über langfristige klimatologische Standortunterschiede
#
# Die Daten umfassen nur ungefähr 48 Stunden. Die Ergebnisse sind deshalb
# eine deskriptive Momentaufnahme der beobachteten Ganglinien.
#
# Benötigte Pakete, falls noch nicht installiert:
# install.packages(c("sf", "dplyr", "tidyr", "ggplot2"))
# ============================================================================

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)


# ----------------------------------------------------------------------------
# 1. Einstellungen
# ----------------------------------------------------------------------------

# Pfad zum bereits aggregierten und räumlich verknüpften Datensatz.
input_file <- "data/MobilePolsterhaus_spatial_hourly.rds"

# Alle Tabellen und Abbildungen werden in diesem Ordner gespeichert.
output_dir <- "auswertung_ganglinien"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ----------------------------------------------------------------------------
# 2. Daten einlesen und prüfen
# ----------------------------------------------------------------------------

dat_sf <- readRDS(input_file)

# Für die zeitlichen und statistischen Auswertungen wird die Geometrie nicht
# benötigt. Das ursprüngliche sf-Objekt dat_sf bleibt dabei unverändert.
dat <- dat_sf |>
  st_drop_geometry()

required_columns <- c(
  "time",
  "air_channel",
  "air_temp",
  "air_humidity",
  "soil_high_temp",
  "soil_deep_temp",
  "soil_high_id",
  "soil_deep_id"
)

missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0) {
  stop(
    "Im RDS fehlen folgende benötigte Spalten: ",
    paste(missing_columns, collapse = ", ")
  )
}

dat <- dat |>
  mutate(
    time = as.POSIXct(time),
    air_channel = as.integer(air_channel),
    channel_label = factor(
      paste0("CH", air_channel),
      levels = paste0("CH", 1:8)
    )
  ) |>
  arrange(time, air_channel)


# ----------------------------------------------------------------------------
# 3. Einheitliches Aussehen der Abbildungen
# ----------------------------------------------------------------------------

plot_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )


# ----------------------------------------------------------------------------
# 4. Anomalien gegenüber dem stündlichen Mittel aller Luftkanäle
# ----------------------------------------------------------------------------
# Eine Anomalie ist hier KEINE Abweichung von einem langjährigen Klimamittel.
# Sie ist ausschließlich die Abweichung eines Kanals vom Mittel aller acht
# Kanäle zum selben Zeitpunkt:
#
#   Anomalie(i,t) = Temperatur(i,t) - Mitteltemperatur aller Kanäle(t)
#
# Dadurch wird der gemeinsame Tagesgang entfernt. Positive Werte bedeuten:
# Der Kanal war zu diesem Zeitpunkt wärmer als der Mittelwert aller Kanäle.

air_anomaly <- dat |>
  group_by(time) |>
  mutate(
    spatial_mean_temp = mean(air_temp, na.rm = TRUE),
    temp_anomaly = air_temp - spatial_mean_temp
  ) |>
  ungroup()

p_anomaly <- ggplot(
  air_anomaly,
  aes(
    x = time,
    y = temp_anomaly,
    colour = channel_label,
    group = channel_label
  )
) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.2) +
  labs(
    title = "Abweichung der Lufttemperatur vom stündlichen Stationsmittel",
    subtitle = paste(
      "Positive Werte: wärmer als das Mittel aller Kanäle;",
      "negative Werte: kälter"
    ),
    x = NULL,
    y = "Temperaturanomalie (°C)",
    colour = "Luftkanal"
  ) +
  plot_theme

print(p_anomaly)

ggsave(
  file.path(output_dir, "01_temperaturanomalien.png"),
  p_anomaly,
  width = 10,
  height = 6,
  dpi = 300
)


# ----------------------------------------------------------------------------
# 5. Deskriptive Kennwerte je Luftkanal
# ----------------------------------------------------------------------------
# mean_anomaly zeigt den durchschnittlichen Niveauunterschied eines Kanals zum
# zeitgleichen Mittel aller Kanäle. mean_abs_anomaly ignoriert das Vorzeichen
# und zeigt die typische Größe der Abweichung.

air_summary <- air_anomaly |>
  group_by(air_channel, channel_label) |>
  summarise(
    n = sum(!is.na(air_temp)),
    mean_temp = mean(air_temp, na.rm = TRUE),
    min_temp = min(air_temp, na.rm = TRUE),
    max_temp = max(air_temp, na.rm = TRUE),
    temperature_range = max_temp - min_temp,
    standard_deviation = sd(air_temp, na.rm = TRUE),
    mean_anomaly = mean(temp_anomaly, na.rm = TRUE),
    mean_abs_anomaly = mean(abs(temp_anomaly), na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  air_summary,
  file.path(output_dir, "kennwerte_lufttemperatur.csv"),
  row.names = FALSE
)

print(air_summary)


# ----------------------------------------------------------------------------
# 6. Lufttemperaturen in ein breites Format bringen
# ----------------------------------------------------------------------------
# Für Korrelations-, Bias- und MAE-Matrizen steht danach jeder Kanal in einer
# eigenen Spalte. Falls versehentlich mehrere Werte je Zeitpunkt und Kanal
# vorhanden sind, wird hier ihr Mittelwert verwendet.

air_wide <- dat |>
  group_by(time, channel_label) |>
  summarise(air_temp = mean(air_temp, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(
    names_from = channel_label,
    values_from = air_temp
  ) |>
  arrange(time)

channel_names <- intersect(paste0("CH", 1:8), names(air_wide))

if (length(channel_names) < 2) {
  stop("Für einen Kanalvergleich werden mindestens zwei Luftkanäle benötigt.")
}

temperature_matrix <- as.matrix(air_wide[channel_names])


# Hilfsfunktion: Eine Matrix für ggplot in eine lange Tabelle umwandeln.
matrix_to_long <- function(x, value_name) {
  out <- as.data.frame(as.table(x), stringsAsFactors = FALSE)
  names(out) <- c("channel_y", "channel_x", value_name)
  out$channel_x <- factor(out$channel_x, levels = channel_names)
  out$channel_y <- factor(out$channel_y, levels = rev(channel_names))
  out
}


# ----------------------------------------------------------------------------
# 7. Korrelationsmatrix der absoluten Lufttemperaturen
# ----------------------------------------------------------------------------
# Pearson-r beschreibt, ob zwei Ganglinien gemeinsam steigen und fallen.
# Die Matrix sagt NICHT, ob die Temperaturen auf demselben Niveau liegen.
# Zwei parallel verlaufende Linien können stark korrelieren und trotzdem einen
# deutlichen konstanten Temperaturunterschied besitzen.

cor_raw_matrix <- cor(
  temperature_matrix,
  use = "pairwise.complete.obs",
  method = "pearson"
)

cor_raw_long <- matrix_to_long(cor_raw_matrix, "correlation")

p_cor_raw <- ggplot(
  cor_raw_long,
  aes(x = channel_x, y = channel_y, fill = correlation)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 3.5) +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#2166AC",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson-r"
  ) +
  coord_equal() +
  labs(
    title = "Korrelation der absoluten Lufttemperatur-Ganglinien",
    subtitle = "Hohe Werte zeigen einen ähnlichen zeitlichen Verlauf",
    x = "Luftkanal",
    y = "Luftkanal"
  ) +
  plot_theme

print(p_cor_raw)

ggsave(
  file.path(output_dir, "02_korrelation_rohwerte.png"),
  p_cor_raw,
  width = 8,
  height = 7,
  dpi = 300
)


# ----------------------------------------------------------------------------
# 8. Korrelationsmatrix der stündlichen Temperaturänderungen
# ----------------------------------------------------------------------------
# Die absoluten Temperaturen besitzen einen gemeinsamen Tagesgang. Deshalb
# wird zusätzlich verglichen, ob die Kanäle von einer Stunde zur nächsten in
# ähnlicher Richtung und Stärke reagieren:
#
#   Änderung(t) = Temperatur(t) - Temperatur(t-1)
#
# Diese Matrix ist bei der Frage nach der Ähnlichkeit der kurzfristigen Dynamik
# oft aussagekräftiger als die reine Rohwertkorrelation.

change_matrix <- apply(temperature_matrix, 2, diff)
colnames(change_matrix) <- channel_names

cor_change_matrix <- cor(
  change_matrix,
  use = "pairwise.complete.obs",
  method = "pearson"
)

cor_change_long <- matrix_to_long(cor_change_matrix, "correlation")

p_cor_change <- ggplot(
  cor_change_long,
  aes(x = channel_x, y = channel_y, fill = correlation)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 3.5) +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#2166AC",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson-r"
  ) +
  coord_equal() +
  labs(
    title = "Korrelation der stündlichen Temperaturänderungen",
    subtitle = "Ähnlichkeit der kurzfristigen Erwärmung und Abkühlung",
    x = "Luftkanal",
    y = "Luftkanal"
  ) +
  plot_theme

print(p_cor_change)

ggsave(
  file.path(output_dir, "03_korrelation_aenderungen.png"),
  p_cor_change,
  width = 8,
  height = 7,
  dpi = 300
)


# ----------------------------------------------------------------------------
# 9. Paarweiser Bias
# ----------------------------------------------------------------------------
# Definition der Zellen:
#
#   Bias(Zeile, Spalte) = Mittelwert(Temperatur Zeile - Temperatur Spalte)
#
# Beispiel: +1.2 in der Zelle CH3/CH1 bedeutet, dass CH3 durchschnittlich
# 1.2 °C wärmer als CH1 war. Die Matrix ist deshalb vorzeichenabhängig und
# an der Hauptdiagonalen gleich null.

n_channels <- length(channel_names)
bias_matrix <- matrix(
  NA_real_,
  nrow = n_channels,
  ncol = n_channels,
  dimnames = list(channel_names, channel_names)
)

mae_matrix <- bias_matrix

for (i in seq_len(n_channels)) {
  for (j in seq_len(n_channels)) {
    difference <- temperature_matrix[, i] - temperature_matrix[, j]
    bias_matrix[i, j] <- mean(difference, na.rm = TRUE)
    mae_matrix[i, j] <- mean(abs(difference), na.rm = TRUE)
  }
}

bias_long <- matrix_to_long(bias_matrix, "bias")
max_abs_bias <- max(abs(bias_long$bias), na.rm = TRUE)

p_bias <- ggplot(
  bias_long,
  aes(x = channel_x, y = channel_y, fill = bias)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%+.2f", bias)), size = 3.3) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs_bias, max_abs_bias),
    name = "Bias (°C)"
  ) +
  coord_equal() +
  labs(
    title = "Mittlere Temperaturdifferenz zwischen den Luftkanälen",
    subtitle = "Zellenwert = Kanal der Zeile minus Kanal der Spalte",
    x = "Vergleichskanal (Spalte)",
    y = "Ausgangskanal (Zeile)"
  ) +
  plot_theme

print(p_bias)

ggsave(
  file.path(output_dir, "04_bias_lufttemperatur.png"),
  p_bias,
  width = 8,
  height = 7,
  dpi = 300
)


# ----------------------------------------------------------------------------
# 10. Paarweiser mittlerer absoluter Fehler (MAE)
# ----------------------------------------------------------------------------
# Der MAE zeigt die typische Größe des Unterschieds ohne Vorzeichen:
#
#   MAE = Mittelwert(|Temperatur Kanal i - Temperatur Kanal j|)
#
# Ein kleiner MAE bedeutet, dass zwei Ganglinien nicht nur ähnlich geformt,
# sondern auch auf einem ähnlichen Temperaturniveau liegen.

mae_long <- matrix_to_long(mae_matrix, "mae")

p_mae <- ggplot(
  mae_long,
  aes(x = channel_x, y = channel_y, fill = mae)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", mae)), size = 3.3) +
  scale_fill_gradient(
    low = "white",
    high = "#B2182B",
    limits = c(0, max(mae_long$mae, na.rm = TRUE)),
    name = "MAE (°C)"
  ) +
  coord_equal() +
  labs(
    title = "Mittlere absolute Unterschiede der Lufttemperatur",
    subtitle = "Kleine Werte zeigen ein ähnliches Temperaturniveau",
    x = "Luftkanal",
    y = "Luftkanal"
  ) +
  plot_theme

print(p_mae)

ggsave(
  file.path(output_dir, "05_mae_lufttemperatur.png"),
  p_mae,
  width = 8,
  height = 7,
  dpi = 300
)


# Matrizen zusätzlich als CSV speichern, damit die exakten Werte außerhalb der
# Abbildungen weiterverwendet werden können.
write.csv(
  cor_raw_matrix,
  file.path(output_dir, "matrix_korrelation_rohwerte.csv")
)
write.csv(
  cor_change_matrix,
  file.path(output_dir, "matrix_korrelation_aenderungen.csv")
)
write.csv(
  bias_matrix,
  file.path(output_dir, "matrix_bias.csv")
)
write.csv(
  mae_matrix,
  file.path(output_dir, "matrix_mae.csv")
)


# ----------------------------------------------------------------------------
# 11. Auswertung der Bodentemperaturen
# ----------------------------------------------------------------------------
# Es werden nur Standorte berücksichtigt, an denen mindestens einer der beiden
# Bodensensoren vorhanden ist. Erwartung:
#   - ca. 2 cm: stärkere und schnellere Temperaturschwankung
#   - ca. 25 cm: gedämpfterer Temperaturverlauf

soil <- dat |>
  filter(!is.na(soil_high_temp) | !is.na(soil_deep_temp)) |>
  arrange(air_channel, time) |>
  mutate(
    soil_difference = soil_high_temp - soil_deep_temp
  )

if (nrow(soil) > 0) {

  # Kennwerte je Standort einschließlich Korrelation und Dämpfungsfaktor.
  # Der Dämpfungsfaktor ist die Spannweite in 25 cm geteilt durch die
  # Spannweite in 2 cm. Werte unter 1 bedeuten eine geringere Schwankung in der
  # tieferen Bodenschicht.
  soil_summary <- soil |>
    group_by(air_channel, channel_label, soil_high_id, soil_deep_id) |>
    summarise(
      n_pairs = sum(!is.na(soil_high_temp) & !is.na(soil_deep_temp)),
      mean_high = mean(soil_high_temp, na.rm = TRUE),
      mean_deep = mean(soil_deep_temp, na.rm = TRUE),
      range_high = max(soil_high_temp, na.rm = TRUE) -
        min(soil_high_temp, na.rm = TRUE),
      range_deep = max(soil_deep_temp, na.rm = TRUE) -
        min(soil_deep_temp, na.rm = TRUE),
      mean_difference_high_minus_deep = mean(
        soil_difference,
        na.rm = TRUE
      ),
      mean_absolute_difference = mean(
        abs(soil_difference),
        na.rm = TRUE
      ),
      correlation_high_deep = cor(
        soil_high_temp,
        soil_deep_temp,
        use = "complete.obs"
      ),
      damping_factor = ifelse(
        range_high > 0,
        range_deep / range_high,
        NA_real_
      ),
      .groups = "drop"
    )

  write.csv(
    soil_summary,
    file.path(output_dir, "kennwerte_bodentemperatur.csv"),
    row.names = FALSE
  )

  print(soil_summary)


  # --------------------------------------------------------------------------
  # 11a. Ganglinien beider Bodentiefen
  # --------------------------------------------------------------------------

  soil_long <- soil |>
    select(time, channel_label, soil_high_temp, soil_deep_temp) |>
    pivot_longer(
      cols = c(soil_high_temp, soil_deep_temp),
      names_to = "depth",
      values_to = "temperature"
    ) |>
    filter(!is.na(temperature)) |>
    mutate(
      depth = recode(
        depth,
        soil_high_temp = "ca. 2 cm",
        soil_deep_temp = "ca. 25 cm"
      ),
      depth = factor(depth, levels = c("ca. 2 cm", "ca. 25 cm"))
    )

  p_soil_lines <- ggplot(
    soil_long,
    aes(
      x = time,
      y = temperature,
      colour = depth,
      linetype = depth,
      group = depth
    )
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.1) +
    facet_wrap(~ channel_label, ncol = 2) +
    labs(
      title = "Ganglinien der Bodentemperatur in zwei Tiefen",
      subtitle = "Vergleich der oberflächennahen und tieferen Bodenschicht",
      x = NULL,
      y = "Bodentemperatur (°C)",
      colour = "Messtiefe",
      linetype = "Messtiefe"
    ) +
    plot_theme

  print(p_soil_lines)

  ggsave(
    file.path(output_dir, "06_bodenganglinien.png"),
    p_soil_lines,
    width = 10,
    height = 7,
    dpi = 300
  )


  # --------------------------------------------------------------------------
  # 11b. Zeitlicher Unterschied zwischen oberem und tiefem Bodensensor
  # --------------------------------------------------------------------------
  # Positive Werte bedeuten: Der Boden in ca. 2 cm war wärmer als in ca. 25 cm.
  # Negative Werte bedeuten: Der tiefere Boden war wärmer.

  p_soil_difference <- ggplot(
    soil,
    aes(x = time, y = soil_difference, group = channel_label)
  ) +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.5) +
    geom_line(colour = "#7B3294", linewidth = 0.9) +
    geom_point(colour = "#7B3294", size = 1.1) +
    facet_wrap(~ channel_label, ncol = 2) +
    labs(
      title = "Temperaturunterschied zwischen den Bodentiefen",
      subtitle = "Differenz = ca. 2 cm minus ca. 25 cm",
      x = NULL,
      y = "Temperaturdifferenz (°C)"
    ) +
    plot_theme

  print(p_soil_difference)

  ggsave(
    file.path(output_dir, "07_bodendifferenz.png"),
    p_soil_difference,
    width = 10,
    height = 7,
    dpi = 300
  )


  # --------------------------------------------------------------------------
  # 11c. Dämpfung der Temperaturschwankung mit der Bodentiefe
  # --------------------------------------------------------------------------

  p_damping <- ggplot(
    soil_summary,
    aes(x = channel_label, y = damping_factor)
  ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "grey35"
    ) +
    geom_col(fill = "#008837", width = 0.65) +
    geom_text(
      aes(label = sprintf("%.2f", damping_factor)),
      vjust = -0.4,
      size = 3.8
    ) +
    expand_limits(y = max(1.05, soil_summary$damping_factor, na.rm = TRUE) * 1.1) +
    labs(
      title = "Dämpfung der Temperaturschwankung in ca. 25 cm Tiefe",
      subtitle = paste(
        "Dämpfungsfaktor = Spannweite 25 cm / Spannweite 2 cm;",
        "Werte unter 1 zeigen eine Dämpfung"
      ),
      x = "Messpunkt beziehungsweise Luftkanal",
      y = "Dämpfungsfaktor"
    ) +
    plot_theme

  print(p_damping)

  ggsave(
    file.path(output_dir, "08_daempfungsfaktor_boden.png"),
    p_damping,
    width = 8,
    height = 5,
    dpi = 300
  )

} else {
  message("Keine Bodentemperaturwerte gefunden; Bodenanalyse übersprungen.")
}


# ----------------------------------------------------------------------------
# 12. Kompakte Interpretationshilfe in der R-Konsole
# ----------------------------------------------------------------------------

cat(
  "\nAuswertung abgeschlossen.\n",
  "Ergebnisse gespeichert unter: ", normalizePath(output_dir), "\n\n",
  "Interpretation:\n",
  "- Hohe Korrelation + kleiner MAE: ähnliche Form und ähnliches Niveau.\n",
  "- Hohe Korrelation + großer Bias: ähnliche Form, aber anderes Niveau.\n",
  "- Hohe Rohwert-, aber geringere Änderungskorrelation: gemeinsamer ",
  "Tagesgang, aber unterschiedliche kurzfristige Reaktion.\n",
  "- Wechselnde Anomalien: Unterschiede verändern sich im Tagesverlauf.\n",
  "- Dämpfungsfaktor unter 1: tiefere Bodenschicht schwankt weniger.\n",
  sep = ""
)
