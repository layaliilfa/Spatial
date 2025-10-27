# ==== Step 1b: Plot lokasi GPS di peta SF Bay Area ====
required <- c("ggplot2", "sf", "rnaturalearth", "rnaturalearthdata")
to_install <- setdiff(required, rownames(installed.packages()))
if(length(to_install)) install.packages(to_install)
library(ggplot2); library(sf); library(rnaturalearth); library(rnaturalearthdata)

# Konversi ke sf (longitude–latitude)
pts <- st_as_sf(dat, coords = c("Longitude", "Latitude"), crs = 4326)

# Ambil shapefile world (termasuk California)
world <- ne_countries(scale = "medium", returnclass = "sf")

# Batasi area sekitar San Francisco Bay
# ==== versi zoom lebih jauh ====
bbox_sf <- st_bbox(c(
  xmin = -124.0,  # lebih ke barat
  xmax = -120.5,  # lebih ke timur
  ymin = 36.0,    # lebih ke selatan
  ymax = 39.8     # lebih ke utara
), crs = st_crs(pts))


ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70") +
  geom_sf(data = pts, color = "red", size = 1.5, alpha = 0.9) +
  coord_sf(xlim = c(bbox_sf["xmin"], bbox_sf["xmax"]),
           ylim = c(bbox_sf["ymin"], bbox_sf["ymax"])) +
  annotate("segment", x = -121.2, xend = -121.2, y = 36.7, yend = 37.2,
           arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  annotate("text", x = -121.25, y = 37.3, label = "N", size = 4, hjust = 0.5) +
  labs(title = "San Francisco Bay Area GPS Stations",
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 12)

# ==== Step 3b: Visual 3-panel per k ====
library(ggplot2)
library(patchwork)
library(sf)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)

# base map
world <- ne_countries(scale = "medium", returnclass = "sf")

# fungsi untuk buat 3 panel sekaligus
make_tripanel <- function(k){
  set.seed(42)
  km <- kmeans(scale(dat[,c("Velocity (E)","Velocity (N)")]), centers=k, nstart=25)
  dat$cluster <- factor(km$cluster)
  pts_sf <- st_as_sf(dat, coords=c("Longitude","Latitude"), crs=4326)
  
  p1 <- ggplot(dat, aes(`Velocity (E)`,`Velocity (N)`, color=cluster))+
    geom_point(size=2)+
    labs(title=paste("Velocity Space (k=",k,")"), x="E [mm/yr]", y="N [mm/yr]")+
    theme_minimal(base_size=11)+theme(legend.position="none")
  
  p2 <- ggplot()+
    geom_sf(data=world, fill="grey95", color="grey70")+
    geom_sf(data=pts_sf, aes(color=cluster), size=1.5)+
    coord_sf(xlim=c(-124,-120.5), ylim=c(36,39.8))+
    labs(title="Map View", x="Longitude", y="Latitude")+
    theme_minimal(base_size=11)+theme(legend.position="none")
  
  p3 <- ggplot(dat, aes(Longitude, Latitude))+
    geom_segment(aes(xend=Longitude+`Velocity (E)`/500,
                     yend=Latitude+`Velocity (N)`/500,
                     color=cluster),
                 arrow=arrow(length=unit(0.12,"cm")))+
    coord_sf(xlim=c(-124,-120.5), ylim=c(36,39.8))+
    labs(title="Velocity Field")+
    theme_minimal(base_size=11)+theme(legend.position="none")
  
  (p1 | p2 | p3) + plot_annotation(title=paste("K-means clustering (k =",k,")"))
}

# Contoh: pilih 3 nilai k representatif
make_tripanel(3)
make_tripanel(5)
make_tripanel(7)

# ==== FINAL clean version (tanpa fault line) ====
library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)
library(patchwork)

world <- ne_countries(scale = "medium", returnclass = "sf")

# hasil clustering k = 4
set.seed(42)
vel_scaled <- scale(dat[,c("Velocity (E)","Velocity (N)")])
km4 <- kmeans(vel_scaled, centers = 4, nstart = 50)
dat$cluster <- factor(km4$cluster)

# warna & simbol
cols   <- c("#d62728","#ff7f0e","#2ca02c","#9467bd")
shapes <- c(15,17,18,16)
names(cols) <- names(shapes) <- levels(dat$cluster)

# === Panel kiri: velocity space (1 garis per blok) ===
centroids <- dat %>%
  group_by(cluster) %>%
  summarise(ve_mean = mean(`Velocity (E)`),
            vn_mean = mean(`Velocity (N)`))

p1 <- ggplot(dat, aes(`Velocity (E)`, `Velocity (N)`, color=cluster, shape=cluster)) +
  geom_point(size=2) +
  geom_segment(data=centroids,
               aes(x=0, y=0, xend=ve_mean, yend=vn_mean),
               color="black", linewidth=0.8, arrow=arrow(length=unit(0.15,"cm"))) +
  scale_color_manual(values=cols) +
  scale_shape_manual(values=shapes) +
  coord_equal() +
  theme_minimal(base_size=12) +
  theme(legend.position="none") +
  labs(title="GPS Horizontal Velocities (k = 4)",
       x="Velocity East (mm/yr)", y="Velocity North (mm/yr)") +
  annotate("text", x=-27, y=33, label="Pacific Block (PB)", fontface="bold", size=3.8) +
  annotate("text", x=-19, y=25, label="Bay Block (BB)", fontface="bold", size=3.8) +
  annotate("text", x=-14, y=17, label="East Bay Block (EBB)", fontface="bold", size=3.8) +
  annotate("text", x=-10, y=9,  label="Sierra Nevada - Great Valley Block (SNGVB)", fontface="bold", size=3.3) +
  annotate("text", x=-27, y=1, label="Reduced χ² = 20.40   RMS = 1.79", size=3, hjust=0, fontface="italic")

# === Panel kanan: map view tanpa fault line ===
pts_sf <- st_as_sf(dat, coords=c("Longitude","Latitude"), crs=4326)

p2 <- ggplot() +
  geom_sf(data=world, fill="grey95", color="grey80") +
  geom_sf(data=pts_sf, aes(color=cluster, shape=cluster), size=2) +
  scale_color_manual(values=cols) +
  scale_shape_manual(values=shapes) +
  coord_sf(xlim=c(-124,-120.5), ylim=c(36,39.8)) +
  labs(title="Map of Clusters (k = 4)", x="Longitude", y="Latitude") +
  theme_minimal(base_size=12) +
  theme(legend.position="none") +
  annotate("text", x=-123.3, y=37.2, label="PB", size=4, fontface="bold") +
  annotate("text", x=-122.4, y=37.9, label="BB", size=4, fontface="bold") +
  annotate("text", x=-122.1, y=38.7, label="EBB", size=4, fontface="bold") +
  annotate("text", x=-122.8, y=39.6, label="SNGVB", size=4, fontface="bold")

# gabung dua panel
(p1 | p2) +
  plot_annotation(
    title="K-Means Clustering of GPS Velocities (San Francisco Bay Area)",
    subtitle="k = 4 — PB, BB, EBB, SNGVB (after Simpson et al., 2012)",
    theme = theme(plot.title=element_text(face="bold", size=13))
  )
