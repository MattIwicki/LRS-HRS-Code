source("scripts/lrs_hrs_source.R")

library(terra)
library(sf)
library(rgeoda)
library(terra)
library(farver)
library(supercells)
library(regional)
library(ggplot2)


raster = rast("data/rasters/orto_170519.tif")
raster_lab = raster_rgb2lab(raster)
slic_2 = supercells(raster_lab, step = 17, compactness = .7)

slic = supercells(raster_lab, step = 14, compactness = .27)
plotRGB(raster)
plot(slic["geometry"], add = T, border = "red")

w = queen_weights(slic)
smst = optimal_segmentation(slic[, c("L", "A", "B")], w, algo_fun = redcap, max_clusters = 150)

# --- defining optimal segments number ---

# plotting lv (mean local variance of all superpixels) to find where the graph
# becomes smooth - let's say it does become smooth with 30 segments and more

ggplot(smst, aes(x = seg_num, y = lv)) + geom_line() + geom_vline(xintercept = 40, col = "red") + xlab("Number of clusters") + ylab("LV")

# finding peaks

smst = smst[1:148,]

spikes = get_peaks(smst, 40)

ggplot(smst, aes(x = seg_num, y = roclv)) + geom_line() + geom_vline(xintercept = spikes, col = "gray") + 
  xlab("Number of clusters") + ylab("ROC-LV")

# calculating mean inhomogeneity for peaks' segnum
mean_inh = peaks_mean_inh(spikes, slic[,c("L", "A", "B")], raster_lab, w, 0.05)

beep()

spikes[7]

ggplot(mean_inh, aes(x = index, y = mean_inh)) + geom_line() + geom_point() + 
  geom_vline(xintercept = 24, col = "gray") + 
  geom_vline(xintercept = 54, col = "gray") + 
  labs(x = "Number of segments", y = "Mean inhomogeneity")

ggplot(mean_inh, aes(x = index, y = mean_isol)) + geom_line() + geom_point() +
  geom_vline(xintercept = 24, col = "gray") + 
  geom_vline(xintercept = 54, col = "gray") + 
  labs(x = "Number of segments", y = "Mean isolation")


clustered = cluster_superpixels(slic[, c("L", "A", "B")], 54, w)
plot(raster)
plot(clustered$geometry, add = T, lwd = 1, border = "red")

ref = read_sf("data/osm_reference.gpkg")

ref = st_transform(ref, st_crs(clustered))

metric_obj = sm_read(ref_sf = ref, seg_sf = clustered)

plot(metric_obj)

temp_summary = sm_compute(metric_obj, "AFI") %>% sm_compute("D_index") %>% sm_compute("OMerging") %>%
  sm_compute("UMerging") %>% sm_compute("recall") %>% summary()

temp_summary

metric_summaries = c(metric_summaries, toString(temp_summary))

# 24 
#         AFI     D_index    OMerging    UMerging      recall 
# -63.6723885   0.5371534  46.4873257   0.3552743   0.7646601 

# 54
#        AFI    D_index   OMerging   UMerging     recall 
# -5.4000624  0.6226288  1.6523608  0.7931964  0.4731620

#        AFI     D_index    OMerging    UMerging      recall 
# -12.7834440   0.5898102   6.0196550   0.6277903   0.5976134 

metric_summaries
names(metric_summaries) = c("fc, high isol", "mc, lowest inh", "mc, highest isol")
metric_summaries
