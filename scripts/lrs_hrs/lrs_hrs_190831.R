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
slic_2 = supercells(raster_lab, step = 17, compactness = .7)

slic = supercells(raster_lab, step = 14, compactness = .27)
plotRGB(raster)
plot(slic["geometry"], add = T, border = "red")
draw_rect_on_raster(0.42, 0.8, 0.52, 0.97, raster)
draw_rect_on_raster(0.09, 0.35, .35, .48, raster)
draw_rect_on_raster(0.71, 0.58, .78, .9, raster)

w = queen_weights(slic)
smst = optimal_segmentation(slic[, c("L", "A", "B")], w, algo_fun = redcap, max_clusters = 150)

# --- defining optimal segments number ---

# plotting lv (mean local variance of all superpixels) to find where the graph
# becomes smooth - let's say it does become smooth with 30 segments and more

ggplot(smst, aes(x = seg_num, y = lv)) + geom_line() + geom_vline(xintercept = 50, col = "red") + xlab("Number of clusters") + ylab("LV")
# finding peaks

smst = smst[1:148,]

spikes = get_peaks(smst, 50)

ggplot(smst, aes(x = seg_num, y = roclv)) + geom_line() + geom_vline(xintercept = spikes, col = "gray") +
  xlab("Number of clusters") + ylab("ROC-LV")

# calculating mean inhomogeneity for peaks' segnum
mean_inh = peaks_mean_inh(spikes, slic[,c("L", "A", "B")], raster_lab, w, 0.05)

beep()

ggplot(mean_inh, aes(x = index, y = mean_inh)) + geom_line() + geom_point() + 
  geom_vline(xintercept = 10, col = "gray") +
  geom_vline(xintercept = 87, col = "gray") + 
  labs(x = "Number of segments", y = "Mean inhomogeneity")
ggplot(mean_inh, aes(x = index, y = mean_isol)) + geom_line() + geom_point() +
  geom_vline(xintercept = 10, col = "gray") +
  geom_vline(xintercept = 87, col = "gray") +
  labs(x = "Number of segments", y = "Mean isolation")

ggplot(mean_inh, aes(x = index, y = mean_inh)) + geom_line() + geom_point() +
  geom_line(aes(y = mean_isol), color = "red") + geom_point(aes(y = mean_isol), color = "red") +
  scale_y_continuous(sec.axis = sec_axis(~., name="mean_isol"))



mean_inh[8,]
mean_inh[20,]
mean_inh[21,]

clustered = cluster_superpixels(slic[, c("L", "A", "B")], 87, w)
plotRGB(raster)
plot(clustered$geometry, add = T, lwd = 1, border = "red")


ref = read_sf("data/orto4_osm_cleanup.gpkg")

ref = st_transform(ref, st_crs(clustered))

metric_obj = sm_read(ref_sf = ref, seg_sf = clustered)

plot(metric_obj)

temp_summary = sm_compute(metric_obj, "AFI") %>% sm_compute("D_index") %>% sm_compute("OMerging") %>%
  sm_compute("UMerging") %>% sm_compute("recall") %>% summary()

temp_summary

metric_summaries = c(metric_summaries, toString(temp_summary))

# 10 
# AFI      D_index     OMerging     UMerging       recall 
# -163.4477751    0.6407899  134.3774249    0.2626128    0.6994615 

# 87
# AFI    D_index   OMerging   UMerging     recall 
# -6.9789672  0.6112813  2.4153782  0.7438823  0.5123075 

metric_summaries
names(metric_summaries) = c("fc, high isol", "mc, lowest inh", "mc, highest isol")
metric_summaries

# The best result is "many clusters, lowest inhomogeneity" - most metrics are most advantagious
# For best boundary recall the best result is "few clusters, high isolation"
# but this is highly undersegmented.
