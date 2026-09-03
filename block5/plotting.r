# ============================================================
# Ganglinien aus dem bereits erzeugten räumlichen RDS darstellen
# ============================================================
# Dieses Skript übernimmt KEINE Aufbereitung, Aggregation oder Verknüpfung.
# Es liest nur das fertige RDS ein und erstellt daraus einfache Diagramme.
#
# Benötigte Pakete, falls noch nicht installiert:
# install.packages(c("sf", "dplyr", "tidyr", "ggplot2"))

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)


# ------------------------------------------------------------
# 1. Fertigen Datensatz einlesen
# ------------------------------------------------------------

dat <- readRDS("block5/data/MobilePolsterhaus_spatial_hourly.rds")

# Für Ganglinien wird die räumliche Geometriespalte nicht benötigt.
# Die übrigen Attribute des Datensatzes bleiben erhalten.
dat_plot <- dat |>
  st_drop_geometry() |>
  arrange(air_channel, time)


# Einheitliches, schlichtes Aussehen für alle Abbildungen.
plot_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title.position = "plot"
  )


# ------------------------------------------------------------
# 2. Ganglinien der Lufttemperatur
# ------------------------------------------------------------
# Jeder Farbton entspricht einem Messpunkt bzw. Luftkanal.

p_air_temperature <- ggplot(
  dat_plot,
  aes(
    x = time,
    y = air_temp,
    colour = factor(air_channel),
    group = air_channel
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.3) +
  labs(
    title = "Ganglinien der Lufttemperatur",
    x = NULL,
    y = "Lufttemperatur (°C)",
    colour = "Kanal"
  ) +
  plot_theme

print(p_air_temperature)


# ------------------------------------------------------------
# 3. Ganglinien der relativen Luftfeuchte
# ------------------------------------------------------------

p_air_humidity <- ggplot(
  dat_plot,
  aes(
    x = time,
    y = air_humidity,
    colour = factor(air_channel),
    group = air_channel
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.3) +
  labs(
    title = "Ganglinien der relativen Luftfeuchte",
    x = NULL,
    y = "Relative Luftfeuchte (%)",
    colour = "Kanal"
  ) +
  plot_theme

print(p_air_humidity)


# ------------------------------------------------------------
# 4. Bodentemperatur in beiden Tiefen
# ------------------------------------------------------------
# Die beiden breiten Temperaturspalten werden nur für den Plot in eine
# lange Tabelle umgeformt. Der ursprüngliche Datensatz bleibt unverändert.

soil_plot <- dat_plot |>
  select(time, air_channel, soil_high_temp, soil_deep_temp) |>
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
    )
  )

p_soil_temperature <- ggplot(
  soil_plot,
  aes(
    x = time,
    y = temperature,
    colour = depth,
    linetype = depth,
    group = depth
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.2) +
  facet_wrap(~ air_channel, ncol = 2) +
  labs(
    title = "Ganglinien der Bodentemperatur",
    subtitle = "Nur Messpunkte mit Bodensensoren",
    x = NULL,
    y = "Bodentemperatur (°C)",
    colour = "Messtiefe",
    linetype = "Messtiefe"
  ) +
  plot_theme

print(p_soil_temperature)


# ------------------------------------------------------------
# 5. Direkter Vergleich: Luft und Boden
# ------------------------------------------------------------
# Dieser Plot stellt an jedem Messpunkt Lufttemperatur und – sofern vorhanden –
# die beiden Bodentemperaturen gemeinsam dar.

temperature_plot <- dat_plot |>
  select(
    time,
    air_channel,
    air_temp,
    soil_high_temp,
    soil_deep_temp
  ) |>
  pivot_longer(
    cols = c(air_temp, soil_high_temp, soil_deep_temp),
    names_to = "sensor",
    values_to = "temperature"
  ) |>
  filter(!is.na(temperature)) |>
  mutate(
    sensor = recode(
      sensor,
      air_temp = "Luft",
      soil_high_temp = "Boden ca. 2 cm",
      soil_deep_temp = "Boden ca. 25 cm"
    ),
    sensor = factor(
      sensor,
      levels = c("Luft", "Boden ca. 2 cm", "Boden ca. 25 cm")
    )
  )

p_temperature_comparison <- ggplot(
  temperature_plot,
  aes(
    x = time,
    y = temperature,
    colour = sensor,
    linetype = sensor,
    group = sensor
  )
) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ air_channel, ncol = 2) +
  labs(
    title = "Temperaturverlauf in Luft und Boden",
    x = NULL,
    y = "Temperatur (°C)",
    colour = "Sensor",
    linetype = "Sensor"
  ) +
  plot_theme

print(p_temperature_comparison)


# ------------------------------------------------------------
# 6. Abbildungen optional als PNG speichern
# ------------------------------------------------------------
# Diesen Abschnitt löschen oder mit # auskommentieren, wenn die Diagramme
# ausschließlich in RStudio angezeigt werden sollen.

dir.create("plots_ganglinien", showWarnings = FALSE)

ggsave(
  "plots_ganglinien/01_lufttemperatur.png",
  p_air_temperature,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  "plots_ganglinien/02_luftfeuchte.png",
  p_air_humidity,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  "plots_ganglinien/03_bodentemperatur.png",
  p_soil_temperature,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  "plots_ganglinien/04_luft_boden_vergleich.png",
  p_temperature_comparison,
  width = 11,
  height = 9,
  dpi = 300
)
