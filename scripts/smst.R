source("scripts/lrs_hrs_source.R")


# choose the image to be processed
# Choose_image is a data frame containing values needed for this script:
# file name, smooth_threshold, lrs and hrs number of segments.
# possible values: "2017-05-19", "2018-04-07",
#                  "2019-08-31", "2020-03-28",
#                  "2020-04-06", "2021-10-01"
# or: 1, 2, 3, 4, 5, 6
current_image = choose_image["2017-05-19",]

# load the RGB raster and convert it to CIELAB
raster = rast(paste0("data/rasters/", current_image$file))
raster_lab = raster_rgb2lab(raster)

# create superpixels and plot them onto the raster
slic = supercells(raster_lab, step = 14, compactness = .27)
plotRGB(raster)
plot(slic["geometry"], add = T, border = "red")

#conver CIELAB superpixel values to RGB
slic = vect_lab2rgb(slic)

# calculate superpixel weights for the REDCAP algorithm
# and calculate lv and roclv
w = queen_weights(slic)
lv_roclv = calculate_lv_roclv(slic[, c("r", "g", "b")], w, algo_fun = redcap, max_clusters = 150)

# plot the lv graph to find the smooth threshold
ggplot(lv_roclv, aes(x = seg_num, y = lv)) + geom_line() +
  geom_vline(xintercept = current_image$smooth_threshold, col = "red") +
  xlab("Number of clusters") + ylab("LV")

# choose peak in roclv
ggplot(lv_roclv, aes(x = seg_num, y = roclv)) + geom_line() +
  geom_vline(xintercept = current_image$smooth_threshold, col = "red") +
  geom_vline(xintercept = current_image$smst_segments, col = "grey") +
  xlab("Number of clusters") + ylab("ROC-LV")

# SMST result
smst = cluster_superpixels(slic[, c("r", "g", "b")], 53, w)
plot(raster)
plot(clustered$geometry, add = T, lwd = 1, border = "red")

# AFI metric
ref = read_sf("data/osm_reference.gpkg")
ref = st_transform(ref, st_crs(smst))

metric_obj = sm_read(ref_sf = ref, seg_sf = smst)
afi = sm_compute(metric_obj, "AFI") %>% summary()
afi

