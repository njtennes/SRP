
library(ggrepel)
library(maps)

# ---- Site coordinates (approx) ----
sites <- data.frame(
  Site = c("Wilcox, AZ", "Deming, NM", "Kingman, AZ",
           "Casa Grande, AZ", "St. Johns, AZ", "Grand Canyon Junction, AZ", 
           "Medicine Bow, WY", "Encino, NM", "Silver City, NM"),
  Lat  = c(32.2, 32.4, 35.2, 32.87, 34.57, 35.65, 41.9, 34.65, 32.74),
  Lon  = c(-109.8, -108.7, -114.15, -111.74, -109.36, -112.14, -106.2, -105.5, -108.28),
  Type = c("Solar", "Solar", "Solar", "Solar", "Solar", "Both", "Wind", "Wind", "Wind")
)

write.csv(sites, "coordinates.csv", row.names = FALSE)

# map outline
us_states <- map_data("state")
az_nm     <- subset(us_states, region %in% c("arizona", "new mexico", "wyoming", "colorado", "utah"))

#plot limits
x_limits <- c(-115.5, -105.0)
y_limits <- c(31.0,   42.5)


#map
p <- ggplot() +
  geom_polygon(
    data = az_nm,
    aes(x = long, y = lat, group = group),
    fill = NA, color = "grey30", linewidth = 0.4
  ) +
  geom_point(
    data = sites,
    aes(x = Lon, y = Lat, color = Type),  # map Type -> color
    size = 2.8
  ) +
  scale_color_manual(
    values = c(
      "Solar" = "salmon2",
      "Wind"  = "lightblue2",
      "Both"  = "green4"
    ),
    name = "Site Type"
  ) +
  geom_text_repel(
    data = sites,
    aes(x = Lon, y = Lat, label = Site),
    size = 3.4, min.segment.length = 0, seed = 42,
    box.padding = 0.4, point.padding = 0.25
  ) +
  coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE) +
  labs(
    x = "Longitude (°)", y = "Latitude (°)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
print(p)
