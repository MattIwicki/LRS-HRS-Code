source("scripts/optimal_seg_functions.R")

library(terra)
library(sf)
library(rgeoda)
library(terra)
library(farver)
library(supercells)
library(regional)
library(ggplot2)
library(segmetric)

library(beepr)

raster = rast("data/Poznan/ready/orto_190831.tif")
raster_lab = raster_rgb2lab(raster)

slic = supercells(raster_lab, step = 14, compactness = .27)
plotRGB(raster)
plot(slic["geometry"], add = T, border = "red")

slic = lab2rgb(slic)

w = queen_weights(slic)
smst = optimal_segmentation(slic[, c("r", "g", "b")], w, algo_fun = redcap, max_clusters = 150)

# --- defining optimal segments number ---

# plotting lv (mean local variance of all superpixels) to find where the graph
# becomes smooth - let's say it does become smooth with 30 segments and more

ggplot(smst, aes(x = seg_num, y = lv)) + geom_line() + geom_vline(xintercept = 40, col = "red") + xlab("Number of clusters") + ylab("LV")

ggplot(smst, aes(x = seg_num, y = roclv)) + geom_line() + geom_vline(xintercept = 40, col = "red") + geom_vline(xintercept = 49, col = "grey") +
  xlab("Number of clusters") + ylab("ROC-LV")

clustered = cluster_superpixels(slic[, c("r", "g", "b")], 49, w)
plot(raster)
plot(clustered$geometry, add = T, lwd = 1, border = "red")

ref = read_sf("data/orto4_osm_cleanup.gpkg")

ref = st_transform(ref, st_crs(clustered))

metric_obj = sm_read(ref_sf = ref, seg_sf = clustered)

plot(metric_obj)

temp_summary = sm_compute(metric_obj, "AFI") %>% summary()

temp_summary

metric_summaries = c(metric_summaries, toString(temp_summary))

# 49
#        AFI
#  -14.3351

