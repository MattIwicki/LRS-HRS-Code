#install.packages(c("terra", "sf", "rgeoda", "supercells", "regional",
#                   "segmetric", "pracma", "ggplot2", "farver", "modEvA", "osmdata"))

library(pracma)
library(sf)
library(rgeoda)
library(terra)
library(modEvA)
library(farver)
library(supercells)
library(regional)
library(ggplot2)
library(osmdata)
library(segmetric)


vect_lab2rgb = function(slic){
  slic_df = st_drop_geometry(slic[,c("L","A","B")])
  
  # LAB 0-1 -> LAB
  slic_df[,1] = slic_df[,1] * 100
  slic_df[,2] = (slic_df[,2] * 256) - 128
  slic_df[,3] = (slic_df[,3] * 256) - 128
  
  # LAB -> RGB
  slic_df_rgb = convert_colour(slic_df, from = "lab", to = "rgb")
  slic$r = slic_df_rgb$r
  slic$g = slic_df_rgb$g
  slic$b = slic_df_rgb$b
  return(slic)
}

raster_rgb2lab = function(s){
  # RGB -> RGB 0-1 (standardized)
  s1 = s/255
  
  # LAB
  library(farver)
  convert_colour2 = function(r, g, b, ...){
    spectrum = cbind(r, g, b)
    convert_colour(spectrum, ...)
  }
  s2 = lapp(s1, convert_colour2, "rgb", "lab")
  
  # LAB 0-1 (standardized)
  s3a = app(s2[[1]], range01)
  s3b = app(s2[[2]], range01)
  s3c = app(s2[[3]], range01)
  
  s3 = c(s3a, s3b, s3c)
  names(s3) = c("L", "A", "B")
  return(s3)
}

cluster_superpixels = function(superpixels, k, w, algo_fun = redcap){
  # perform clustering
  clusters = algo_fun(k, w, superpixels, cpu_threads=1, method = "fullorder-completelinkage")
  
  # assign superpixels to clusters
  superpixels$cluster = clusters$Clusters
  
  # merge superpixels into clusters
  superpixels = aggregate(superpixels, by = list(superpixels$cluster), mean)
  return(superpixels)
}

calculate_lv_roclv = function(superpixels, w, algo_fun = redcap, max_clusters = 70, echo = F){
  # this is a partial implementation of the SMST process (Wang et al., 2018)
  # [!] superpixels must only only contain wanted variables, eg superpixels[,c("L", "A", "B")]
  
  attributes_number = ncol(slic) - 1
  smst = data.frame(seg_num = 1:max_clusters, lv = 1:max_clusters, roclv = 1:max_clusters, time_elapsed = 1:max_clusters)
  
  last_lv = 0
  for(i in max_clusters:1){
    if(echo){
      print(i)
    }
    start_i = Sys.time()
    clusters = algo_fun(i, w, superpixels, cpu_threads=1, method = "fullorder-completelinkage")
    temp = superpixels
    temp$cluster = clusters$Cluster
    temp = st_drop_geometry(temp)
    sum_subtree = 0
    for(subtree in 1:i){
      sum_attributes = 0
      for(attribute in 1:attributes_number){
        mean_attribute = mean(temp[temp$cluster == subtree, attribute])
        sum_attributes = sum_attributes + sum((temp[temp$cluster == subtree, attribute] - mean_attribute)^2)
      }
      sum_subtree = sum_subtree + sqrt(sum_attributes/ nrow(temp[temp$cluster == subtree,]))
    }
    
    smst$lv[i] = sum_subtree / i
    smst$roclv[i] = ((smst$lv[i] - last_lv) / last_lv) * 100
    last_lv = smst$lv[i]
    end_i = Sys.time()
    
    smst$time_elapsed[i] = end_i - start_i
    if(echo){
      print(end_i - start_i)
    }
  }
  return(smst[-70,])
}

get_peaks_old = function(smst, smooth_threshold, diff_threshold){
  
  spikes = findpeaks(smst[, "roclv"])
  #spikes = c(0, 0, diff(diff(smst[, "roclv"])))
  #spikes[1:smooth_threshold] = 1
  #spikes[spikes < diff_threshold] = 1
  #spikes[smooth_threshold] = 0
  return(spikes)
}

get_peaks = function(smst, smooth_threshold, below_threshold = 3, above_threshold = 1.2){
  smst1 = smst[1:(smooth_threshold - 1),"roclv"]
  smst2 = smst[smooth_threshold:length(smst[, "roclv"]),"roclv"]
  
  # peaks below threshold, positive and negative
  spikes_max1 = findpeaks(smst1, threshold = below_threshold)
  spikes_min1 = findpeaks(-smst1, threshold = below_threshold)
  
  # peaks above threshold, positive and negative
  spikes_max2 = findpeaks(smst2, threshold = above_threshold)
  spikes_max2[,2] = spikes_max2[,2] + smooth_threshold
  spikes_min2 = findpeaks(-smst2, threshold = above_threshold)
  spikes_min2[,2] = spikes_min2[,2] + smooth_threshold
  
  # clean-up 
  spikes = c(spikes_max1[,2], spikes_min1[,2], spikes_max2[,2], spikes_min2[,2])
  spikes = unique(spikes)
  spikes = sort(spikes)
  return(spikes)
}
  
plot_roclv = function(smst, spikes){
  ggplot(smst, aes(x = seg_num, y = roclv)) + geom_line() + geom_vline(xintercept = which(spikes != 1), col = "gray") +
    geom_vline(xintercept = min(which(spikes > 1)), col = "green") + geom_vline(xintercept = which(spikes == 0), col = "red")
}

peaks_mean_inh = function(spikes, superpixels, raster, w, sample_size){
  indexes = spikes
  result = data.frame(index = indexes, mean_inh = 1:length(indexes), mean_isol = 1:length(indexes))
  
  for(i in 1:length(indexes)){
    print(i)
    superpixels1 = cluster_superpixels(superpixels, k = i, w = w)
    superpixels1$inh = reg_inhomogeneity(superpixels1, raster, sample_size = sample_size)
    superpixels1$isol = reg_isolation(superpixels1, raster, sample_size = sample_size)
    result[i,"mean_inh"] = mean(superpixels1$inh)
    result[i, "mean_isol"] = mean(superpixels1$isol)
  }
  return(result)
}

peaks_mean_inh2 = function(spikes, superpixels, raster, w, sample_size){
  indexes = which(spikes > 1)
  result = data.frame(index = indexes, mean_inh = 1:length(indexes), mean_isol = 1:length(indexes))
  
  for(i in 1:length(indexes)){
    print(i)
    superpixels1 = cluster_superpixels(superpixels, k = i, w = w)
    superpixels1$inh = reg_inhomogeneity(superpixels1, raster, sample_size = sample_size)
    superpixels1$isol = reg_isolation(superpixels1, raster, sample_size = sample_size)
    result[i,"mean_inh"] = mean(superpixels1$inh)
    result[i, "mean_isol"] = mean(superpixels1$isol)
  }
  return(result)
}

# functions not used in this project

download_osm_from_raster = function(raster) {
  # reprojecting raster to fit the OSM projection
  temp_raster = project(raster, "+proj=longlat +datum=WGS84")
  bbox = st_bbox(temp_raster)
  
  osm = osmdata_sf(opq(bbox = bbox))
  
  polygons = osm$osm_polygons$geometry
  polygons = st_transform(polygons, crs = crs(raster))
  
  #clipping
  raster_bbox = st_bbox(raster)
  clip_polygon = st_polygon(list(cbind(c(raster_bbox[1], raster_bbox[3], raster_bbox[3], raster_bbox[1], raster_bbox[1]),
                                       c(raster_bbox[2], raster_bbox[2], raster_bbox[4], raster_bbox[4], raster_bbox[2]))))
  clip_polygon = st_sfc(clip_polygon, crs = st_crs(polygons))
  clipped_polygons = st_intersection(polygons, clip_polygon)
  
  return(clipped_polygons)
}

draw_rect_on_raster = function(xmin, ymin, xmax, ymax, raster, color = "red", width = 1.5){
  bbox = st_bbox(raster)
  difference_x = bbox[3] - bbox[1]
  difference_y = bbox[4] - bbox[2]
  rect(bbox[1] + (difference_x * xmin),
       bbox[2] + (difference_y * ymin),
       bbox[1] + (difference_x * xmax),
       bbox[2] + (difference_y * ymax),
       border = color, lwd = width)
}
