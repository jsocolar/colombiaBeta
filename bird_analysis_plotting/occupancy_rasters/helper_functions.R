# functions
require(stars); require(dplyr)

# calculate sum & max occupancies across square regions centered on each 
# centroid. Computing the max allows for exclusion of regions that are out of 
# range (i.e. have very low occupancies due to being peripheral to a species'
# range)
calc_regional_occupancies <- function(forest_probs, pasture_probs, centroids, grid_size) {
    # note grid_size supplied as a resolution, so divide to get a diameter
    polys_buff <- st_buffer(centroids, grid_size/2, endCapStyle = "SQUARE")
    forest_sum <- aggregate(forest_probs, polys_buff, FUN = sum, use_gdal = T)
    pasture_sum <- aggregate(pasture_probs, polys_buff, FUN = sum, use_gdal = T)
    
    forest_max <- aggregate(forest_probs, polys_buff, FUN = max, use_gdal = T)
    pasture_max <- aggregate(pasture_probs, polys_buff, FUN = max, use_gdal = T)
    
    list(forest_sum = st_as_sf(forest_sum), 
         forest_max = st_as_sf(forest_max), 
         pasture_sum = st_as_sf(pasture_sum), 
         pasture_max = st_as_sf(pasture_max))
}

# calculate summary measures of the difference between pasture and forest for 
# supplied regions. Note that this is summarising, across the community present 
# in a region, the absolute/relative difference between expected number of 
# points occupied in forest and pasture 
# note: works with probabilities, *not* logodds
calc_regional_summary <- function(regional_occs, centroids, threshold) {
    # do thresholding
    forest_max_mat <- sf_to_mat(regional_occs$forest_max)
    forest_max_mat[is.na(forest_max_mat)] <- 0
    pasture_max_mat <- sf_to_mat(regional_occs$pasture_max)
    pasture_max_mat[is.na(pasture_max_mat)] <- 0
    p_above_threshold <- (forest_max_mat > threshold) | (pasture_max_mat > threshold)
    
    forest_sum_mat <- sf_to_mat(regional_occs$forest_sum) 
    forest_sum_mat <- forest_sum_mat * p_above_threshold
    
    pasture_sum_mat <- sf_to_mat(regional_occs$pasture_sum) 
    pasture_sum_mat <- pasture_sum_mat * p_above_threshold
    
    # calculate ratio/difference
    rel_diff_mat <- forest_sum_mat/pasture_sum_mat
    abs_diff_mat <- forest_sum_mat - pasture_sum_mat
    
    centroids %>%
        mutate(avg_abs_diff = matrixStats::rowMeans2(abs_diff_mat, na.rm=T),
               median_abs_diff = matrixStats::rowMedians(abs_diff_mat, na.rm=T),
               avg_ratio = matrixStats::rowMeans2(rel_diff_mat, na.rm=T),
               avg_logratio = matrixStats::rowMeans2(log(rel_diff_mat), na.rm=T),
               median_logratio = matrixStats::rowMedians(log(rel_diff_mat), na.rm=T))
}

# calculate summary measures of the difference in occupancy between pasture and 
# forest at the point level. Note that this is summarising, across the community 
# present at each point, the absolute/relative difference in occupancy between 
# forest and pasture at that point in space
# note: works with probabilities, *not* logodds
calc_point_summary <- function(forest_points, pasture_points, centroids, threshold) {
    forest_points_mat <- sf_to_mat(forest_points)
    pasture_points_mat <- sf_to_mat(pasture_points)
    
    # do thresholding
    fp_mat2 <- forest_points_mat
    fp_mat2[is.na(fp_mat2)] <- 0
    pp_mat2 <- pasture_points_mat
    pp_mat2[is.na(pp_mat2)] <- 0
    p_above_threshold <- (fp_mat2 > threshold) | (pp_mat2 > threshold)
    
    # calculate ratio
    rel_diff_mat <- forest_points_mat/pasture_points_mat
    abs_diff_mat <- forest_points_mat - pasture_points_mat
    
    centroids %>%
        mutate(avg_abs_diff = matrixStats::rowMeans2(abs_diff_mat, na.rm=T),
               median_abs_diff = matrixStats::rowMedians(abs_diff_mat, na.rm=T),
               avg_ratio = matrixStats::rowMeans2(rel_diff_mat, na.rm=T),
               avg_logratio = matrixStats::rowMeans2(log(rel_diff_mat), na.rm=T),
               median_logratio = matrixStats::rowMedians(log(rel_diff_mat), na.rm=T))
}

# helper function to extract matrix of species x point occupancies from sf object
sf_to_mat <- function(x) {
    as_tibble(x) %>%
        dplyr::select(-geometry) %>%
        as.matrix
}
