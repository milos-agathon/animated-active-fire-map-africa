#############################################
# Animated active fires map of Africa
# Milos Popovic 2022/06/04
#############################################
# libraries we need
libs <- c(
  "tidyverse", "sf", "giscoR",
  "data.table", "classInt",
  "gganimate", "gifski"
)

# install missing libraries
installed_libs <- libs %in% rownames(installed.packages())
if (any(installed_libs == F)) {
  install.packages(libs[!installed_libs])
}

# load libraries
invisible(lapply(libs, library, character.only = T))

# 1. GET THE AFRICA SHP
#----------------------

get_africa_sf <- function(africa) {
  africa <- giscoR::gisco_get_countries(
    year = "2020", epsg = "4326",
    resolution = "10", region = "Africa"
  )

  return(africa)
}

africa <- get_africa_sf()

africa_coords <- st_bbox(africa) %>%
  as.matrix() %>%
  as.vector()

# 2. MAKE A HEXBIN MAP OF AFRICA
#-------------------------------

crsLONGLAT <- "+proj=longlat +datum=WGS84 +no_defs"

get_africa_hex <- function(africa_transfomed, africa_hex, africa_final) {
  # Make grid of circa 100,000 m2

  africa_transfomed <- africa %>%
    st_transform(3575)

  africa_hex <- st_make_grid(africa_transfomed,
    cellsize = (3 * sqrt(3) * 196.18875^2) / 2,
    what = "polygons",
    square = F
  ) %>%
    st_intersection(africa_transfomed) %>%
    st_sf() %>%
    mutate(id = row_number()) %>%
    filter(
      st_geometry_type(.) %in% c("POLYGON", "MULTIPOLYGON")
    ) %>%
    st_cast("MULTIPOLYGON") # transform to multipolygons

  africa_final <- sf::st_transform(africa_hex, crs = crsLONGLAT)

  return(africa_final)
}

africa_final <- get_africa_hex()


# 3. GET THE FIRE DATA
#---------------------

get_fire_data <- function(main_url, map_key, source,
                          area_coords, day_range, date) {
  url <- paste(main_url, map_key, source, area_coords,
    day_range, date,
    sep = "/"
  )

  fire_data <- data.table::fread(url)

  return(fire_data)
}

main_url <- "https://firms.modaps.eosdis.nasa.gov/api/area/csv"
map_key <- "*******************************"
source <- "MODIS_NRT"
area_coords <- paste(africa_coords, sep = "", collapse = ",")
day_range <- 10
date <- Sys.Date() - 1

fire_data <- get_fire_data(
  main_url = main_url, map_key = map_key, source = source,
  area_coords = area_coords, day_range = day_range, date = date
)


# 4. FILTER AND AGGREGATE AFRICA'S FIRE EVENTS
#---------------------------------------------
get_africa_fires <- function(fire_coords, fire_africa) {
  fire_coords <- st_as_sf(fire_data, coords = c(2, 1)) %>%
    st_set_crs(4326) %>%
    st_transform(crsLONGLAT)

  fire_africa <- sf::st_join(fire_coords, africa_final, join = sf::st_within)

  return(fire_africa)
}

fire_africa <- get_africa_fires()


get_aggregated_africa_fires <- function() {
  fire_africa$weighted_fire <- fire_africa$confidence / 100

  fire_africa_sum <- st_drop_geometry(fire_africa) %>%
    group_by(id, acq_date) %>%
    summarise_at(
      vars(weighted_fire),
      list(sum_fire = sum)
    )

  fire_africa_sf <- dplyr::left_join(africa_final, fire_africa_sum, by = "id")

  fire_africa_sf$sum_fire <- round(fire_africa_sf$sum_fire, 0)
  fire_africa_sf$sum_fire[fire_africa_sf$sum_fire == 0] <- NA
  fire_africa_sf$acq_date <- as.Date(fire_africa_sf$acq_date)

  return(fire_africa_sf)
}

fire_africa_sf <- get_aggregated_africa_fires()

# 5. ANIMATE
#-----------

get_intervals <- function(ni, labels, df) {
  df <- drop_na(fire_africa_sf)
  ni <- classIntervals(df$sum_fire,
    n = 5,
    style = "jenks"
  )$brks
  # create categories
  labels <- c()
  for (i in 1:length(ni)) {
    labels <- c(labels, paste0(
      round(ni[i], 0),
      "–",
      round(ni[i + 1], 0)
    ))
  }
  labels <- labels[1:length(labels) - 1]

  df$cat <- cut(df$sum_fire,
    breaks = ni,
    labels = labels,
    include.lowest = T
  )

  return(df)
}

df <- get_intervals()

# animate the map

# FOR WINDOWS USERS
# sysfonts::font_add_google("Montserrat", "Montserrat")
# showtext::showtext_auto()

get_africa_map <- function(africa_map) {
  africa_map <- ggplot(na.omit(df)) +
    geom_sf(
      mapping = aes(
        fill = cat,
        group = interaction(cat, acq_date)
      ),
      color = NA, size = 0
    ) +
    geom_sf(
      data = africa,
      fill = NA,
      color = "white",
      size = 0.05
    ) +
    scale_fill_manual(
      name = "",
      values = rev(c(
        "#3d0965", "#8a335c", "#c26958",
        "#e3aa61", "#f1ef75"
      )),
      drop = F
    ) +
    coord_sf(crs = crsLONGLAT) +
    guides(fill = guide_legend(
      direction = "horizontal",
      keyheight = unit(1.5, units = "mm"),
      keywidth = unit(12.5, units = "mm"),
      title.position = "top",
      title.hjust = .5,
      label.hjust = .5,
      nrow = 1,
      byrow = T,
      reverse = F,
      label.position = "bottom"
    )) +
    theme_minimal() +
    theme(
     # text = element_text(family = "Montserrat"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = c(.5, .125),
      legend.text = element_text(size = 20, color = "white"),
      panel.grid.major = element_line(color = "grey20", size = .2),
      panel.grid.minor = element_blank(),
      plot.title = element_text(
        face = "bold", size = 40,
        color = "white", hjust = .5, vjust = -7
      ),
      plot.caption = element_text(
        size = 14, color = "white",
        hjust = .5, vjust = 16
      ),
      plot.subtitle = element_text(
        size = 60, color = "#c43c4e",
        hjust = .5, vjust = -2
      ),
      plot.margin = unit(c(t = -4, r = -4, b = -4, l = -4), "lines"),
      plot.background = element_rect(fill = "grey20", color = NA),
      panel.background = element_rect(fill = "grey20", color = NA),
      legend.background = element_rect(fill = "grey20", color = NA),
      panel.border = element_blank()
    ) +
    labs(
      x = "",
      y = "",
      title = "Active fires in Africa adjusted by confidence level",
      subtitle = "{as.Date(frame_time)}",
      caption = "©2022 Milos Popovic (https://milospopovic.net)
        MODIS Collection 61 NRT Hotspot/Active Fire Detections MCD14DL
        from NASA FIRMS"
    )

  return(africa_map)
}

africa_map <- get_africa_map()

africa_map <- africa_map +
  transition_time(acq_date) +
  enter_fade() +
  exit_shrink() +
  ease_aes("sine-in-out", interval = .25)

africa_anim <- animate(africa_map,
  nframes = 80,
  duration = 25,
  start_pause = 3,
  end_pause = 30,
  height = 6,
  width = 7.15,
  res = 300,
  units = "in"
)

anim_save("active_fires_africa2.gif", africa_anim)
