# Packages
library(ggplot2)
library(ggrepel)
library(dplyr)
library(maps)

# -----------------------------
# Site coordinates
# -----------------------------
sites <- data.frame(
  Site = c("Willcox, AZ", "Deming, NM", "Kingman, AZ",
           "Casa Grande, AZ", "St. Johns, AZ", "Grand Canyon Junction, AZ"),
  Lat  = c(32.2, 32.4, 35.2, 32.87, 34.57, 35.65),
  Lon  = c(-109.8, -108.7, -114.15, -111.74, -109.36, -112.14)
)

# -----------------------------
# State outlines
# -----------------------------
us_states <- map_data("state")

southwest_states <- us_states %>%
  filter(region %in% c("arizona", "new mexico", "utah", "colorado", "wyoming"))

# -----------------------------
# Function to create ellipse polygons
# -----------------------------
make_ellipse <- function(x0, y0, a, b, angle = 0, n = 200) {
  t <- seq(0, 2 * pi, length.out = n)
  x <- a * cos(t)
  y <- b * sin(t)
  
  theta <- angle * pi / 180
  x_rot <- x * cos(theta) - y * sin(theta)
  y_rot <- x * sin(theta) + y * cos(theta)
  
  data.frame(
    x = x0 + x_rot,
    y = y0 + y_rot
  )
}

# -----------------------------
# Rough schematic desert regions
# -----------------------------
mojave <- make_ellipse(-114.2, 36.0, 1.9, 1.5)
sonoran <- make_ellipse(-112.8, 32.6, 2.9, 2.5, angle = 8)
chihuahuan <- make_ellipse(-108.0, 32.0, 1.9, 1.6)

mojave$desert <- "Mojave"
sonoran$desert <- "Sonoran"
chihuahuan$desert <- "Chihuahuan"

deserts_poly <- rbind(mojave, sonoran, chihuahuan)

# -----------------------------
# Plot limits
# -----------------------------
x_limits <- c(-115.8, -106.0)
y_limits <- c(31.0, 37.4)

# -----------------------------
# Plot
# -----------------------------
p <- ggplot() +
  
  # base map first
  geom_polygon(
    data = southwest_states,
    aes(x = long, y = lat, group = group),
    fill = "grey97",
    color = "grey50",
    linewidth = 0.4
  ) +
  
  # desert shading on top of map
  geom_polygon(
    data = deserts_poly,
    aes(x = x, y = y, group = desert, fill = desert),
    alpha = 0.25,
    color = NA
  ) +
  
  scale_fill_manual(
    values = c(
      "Mojave" = "#C98E8E",
      "Sonoran" = "#D8A06A",
      "Chihuahuan" = "#A8A27A"
    )
  ) +
  
  # site points
  geom_point(
    data = sites,
    aes(x = Lon, y = Lat),
    size = 2.8,
    color = "black"
  ) +
  
  # site labels
  geom_text_repel(
    data = sites,
    aes(x = Lon, y = Lat, label = Site),
    size = 3.5,
    color = "black",
    seed = 42,
    box.padding = 0.35,
    point.padding = 0.2,
    min.segment.length = 0,
    segment.color = "grey35",
    segment.linewidth = 0.35
  ) +
  
  # desert labels
  annotate(
    "text",
    x = -114.2, y = 36.15,
    label = "Mojave\n(Winter Rainfall)",
    fontface = "bold",
    size = 4,
    color = "#102A5E"
  ) +
  annotate(
    "text",
    x = -112.8, y = 32.55,
    label = "Sonoran\n(Bimodal Rainfall)",
    fontface = "bold",
    size = 4,
    color = "#102A5E"
  ) +
  annotate(
    "text",
    x = -108.0, y = 31.95,
    label = "Chihuahuan\n(Summer Monsoon)",
    fontface = "bold",
    size = 4,
    color = "#102A5E"
  ) +
  
  coord_fixed(
    ratio = 1.18,
    xlim = x_limits,
    ylim = y_limits,
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude (°)",
    y = "Latitude (°)",
    caption = "Shaded regions are approximate representations of generalized desert climate zones."
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "grey20"),
    plot.caption = element_text(color = "grey35", hjust = 0)
  )

p

# Optional export
ggsave("site_map_deserts.png", p, width = 8.5, height = 5.5, dpi = 400, bg = "white")
p
