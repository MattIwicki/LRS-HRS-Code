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

# calculate superpixel weights for the REDCAP algorithm
# and calculate lv and roclv - the target variable is called smst
# to show which part of this process is based on the SMST paper (Wang et al., 2018)
w = queen_weights(slic)
smst = calculate_lv_roclv(slic[, c("L", "A", "B")], w, algo_fun = redcap, max_clusters = 150)

# plot the lv graph to find the smooth threshold
ggplot(smst, aes(x = seg_num, y = lv)) + geom_line() +
  geom_vline(xintercept = current_image$smooth_threshold, col = "red") + xlab("Number of clusters") + ylab("LV")

# find peaks in roclv
# - the peaks variable is called spikes in order not to interfere with
# pracma package's namespace
smst = smst[1:148,]
spikes = get_peaks(smst, current_image$smooth_threshold)

# plot the found peaks
ggplot(smst, aes(x = seg_num, y = roclv)) + geom_line() + geom_vline(xintercept = spikes, col = "gray") + 
  xlab("Number of clusters") + ylab("ROC-LV")

# calculate mean inhomogeneity and isolation
# for the numbers of segments defined by peaks
inh_isol = peaks_mean_inh_isol(spikes, slic[,c("L", "A", "B")], raster_lab, w, 0.05)

# plot the inhomogeneity and isolation graphs
ggplot(inh_isol, aes(x = index, y = mean_inh)) + geom_line() + geom_point() + 
  geom_vline(xintercept = current_image$lrs_segments, col = "gray") + 
  geom_vline(xintercept = current_image$hrs_segments, col = "gray") + 
  labs(x = "Number of segments", y = "Mean inhomogeneity")

ggplot(inh_isol, aes(x = index, y = mean_isol)) + geom_line() + geom_point() +
  geom_vline(xintercept = current_image$lrs_segments, col = "gray") + 
  geom_vline(xintercept = current_image$hrs_segments, col = "gray") + 
  labs(x = "Number of segments", y = "Mean isolation")

# LRS result
lrs = cluster_superpixels(slic[, c("L", "A", "B")], current_image$lrs_segments, w)
plotRGB(raster)
plot(lrs$geometry, add = T, lwd = 1, border = "red")

# HRS result
hrs = cluster_superpixels(slic[, c("L", "A", "B")], current_image$hrs_segments, w)
plotRGB(raster)
plot(hrs$geometry, add = T, lwd = 1, border = "red")

# AFI metric
ref = read_sf("data/osm_reference.gpkg")
ref = st_transform(ref, st_crs(lrs))

afi = calculate_afi_lrs_hrs(lrs, hrs, ref)
afi
