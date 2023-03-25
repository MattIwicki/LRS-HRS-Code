source("scripts/lrs_hrs_source.R")

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

raster = rast("data/rasters/orto_180407.tif")
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

ggplot(smst, aes(x = seg_num, y = lv)) + geom_line() + geom_vline(xintercept = 60, col = "red") + xlab("Number of clusters") + ylab("LV")

# finding peaks
smst = smst[1:148,]

spikes = get_peaks(smst, 60)
length(spikes)
ggplot(smst, aes(x = seg_num, y = roclv)) + geom_line() + geom_vline(xintercept = spikes, col = "gray") +
  xlab("Number of segments") + ylab("ROC-LV")

# calculating mean inhomogeneity for peaks' segnum
mean_inh = peaks_mean_inh(spikes, slic[,c("L", "A", "B")], raster_lab, w, 0.05)

beep()

ggplot(mean_inh, aes(x = index, y = mean_inh)) + geom_line() + geom_point() + 
  geom_vline(xintercept = 27, col = "gray") +
  geom_vline(xintercept = 78, col = "gray") + 
  labs(x = "Number of segments", y = "Mean inhomogeneity")
ggplot(mean_inh, aes(x = index, y = mean_isol)) + geom_line() + geom_point() +
  geom_vline(xintercept = 27, col = "gray") +
  geom_vline(xintercept = 78, col = "gray") +
  labs(x = "Number of segments", y = "Mean isolation")


clustered = cluster_superpixels(slic[, c("L", "A", "B")], 78, w)
plotRGB(raster)
plot(clustered$geometry, add = T, lwd = 1, border = "red")


ref = read_sf("data/osm_reference.gpkg")

ref = st_transform(ref, st_crs(clustered))


metric_obj = sm_read(ref_sf = ref, seg_sf = clustered)

plot(metric_obj)

temp_summary = sm_compute(metric_obj, "AFI") %>% sm_compute("D_index") %>% sm_compute("OMerging") %>%
  sm_compute("UMerging") %>% sm_compute("recall") %>% summary()

temp_summary

metric_summaries = c(metric_summaries, toString(temp_summary))

# 27 
#         AFI     D_index    OMerging    UMerging      recall 
# -16.8947829   0.5154926  12.0487016   0.3926926   0.6882677 

# 78
#        AFI    D_index   OMerging   UMerging     recall
# -3.4988736  0.5911763  1.3549356  0.7147376  0.4861250

metric_summaries
names(metric_summaries) = c("fc, high isol", "mc, lowest inh", "mc, highest isol")
metric_summaries

# The best result is "many clusters, lowest inhomogeneity" - most metrics are most advantagious
# For best boundary recall the best result is "few clusters, high isolation"
# but this is highly undersegmented.
