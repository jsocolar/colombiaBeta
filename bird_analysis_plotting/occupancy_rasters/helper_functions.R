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
calc_regional_summary <- function(regional_occs, centroids, 
                                  apply_threshold = TRUE, threshold_value = 0.1) {
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
# note:
#   (1) *significantly* faster than using st_apply across x and y
#   (2) works with probabilities, *not* logodds
calc_pt_summary <- function(forest_probs, pasture_probs, 
                            apply_threshold=TRUE, threshold_value = 0.1) {
    if(apply_threshold) {
        above_threshold <- c(forest_probs, pasture_probs) %>%
            setNames(c("p_forest", "p_pasture")) %>%
            mutate(above_threshold = ifelse(p_forest > 0.01 | p_pasture > 0.01, 1, NA)) %>%
            select(above_threshold)
        forest_probs <- forest_probs * above_threshold
        pasture_probs <- pasture_probs * above_threshold
    }
    
    fmat <- st_as_sf(forest_probs, na.rm = FALSE) %>%
        sf_to_mat
    pmat <- st_as_sf(pasture_probs, na.rm = FALSE) %>%
        sf_to_mat
    
    rel_diff <- fmat/pmat
    
    forest_probs %>% 
        slice(1, along="species") %>%
        mutate(avg_ratio = matrixStats::rowMeans2(rel_diff, na.rm=T),
               avg_logratio = matrixStats::rowMeans2(log(rel_diff), na.rm=T),
               median_logratio = matrixStats::rowMedians(log(rel_diff), na.rm=T)) %>%
        select(-1) # note: need to remove attribute here, rather than start
}

# helper function to extract matrix of species x point occupancies from sf object
sf_to_mat <- function(x) {
    as_tibble(x) %>%
        dplyr::select(-geometry) %>%
        as.matrix
}

classify_thresh <- function(DT, threshold=0.1) {
    for (i in names(DT)) {
        DT[is.infinite(get(i)), (i):=NA]
        DT[get(i) < threshold, (i):=0]
        DT[get(i) >= threshold, (i):=1]
    }
}

trim_stars <- function(x) {
    n_attr <- length(x)
    dim1_list <- c()
    dim2_list <- c()
    for(i in 1:n_attr) {
        dim1 = apply(x[[i]], 1, function(x) !all(is.na(x)))
        dim1 = which(dim1)
        dim1 = dim1[1]:dim1[length(dim1)]
        dim2 = apply(x[[i]], 2, function(x) !all(is.na(x)))
        dim2 = which(dim2)
        dim2 = dim2[1]:dim2[length(dim2)]
        dim1_list[[i]] <- dim1
        dim2_list[[i]] <- dim2
    }
    d1u <- unlist(dim1_list)
    d2u <- unlist(dim2_list)
    d1u <- sort(unique(d1u))
    d2u <- sort(unique(d2u))
    
    x = x[, d1u, d2u]
    x = st_normalize(x)
    return(x)
}

calc_regional_summ_v3 <- function(forest_probs, pasture_probs, 
                                  d = c(0, 2, 5, 10, 12, 20), point_spacing=5, 
                                  threshold = .1) {
    # note: 
    # - point spacing is expressed in terms of every nth cell to use as a cell
    # for generating regions from
    # - d is expressed in terms of number of cells to buffer around focal cell
    # cell_res <- 2
    # d = c(0, 2, 5, 10, 12, 20)
    # point_spacing=5
    # threshold = .1
    
    start.time <- Sys.time()
    
    # extract indexing in the stars object 
    col_index_stars <- forest_probs %>%
        slice(1, along="species") %>%
        mutate(id_cell = 1:n(), id_x=NA, id_y = NA) %>%
        select(id_cell, id_x, id_y)
    
    xy_dim <- dim(col_index_stars)
    col_index_stars[["id_x"]] <- matrix(rep(1:xy_dim[2], each=xy_dim[1]), 
                                        nrow=xy_dim[2], ncol=xy_dim[1])
    col_index_stars[["id_y"]] <- matrix(rep(1:xy_dim[1], xy_dim[2]), 
                                        nrow=xy_dim[2], ncol=xy_dim[1])
    
    col_index_sf <- st_as_sf(col_index_stars, as_points = TRUE)
    col_index_dt <- as.data.table(col_index_sf)
    
    # manage memory
    rm(col_index_sf, col_index_stars)
    
    # store x and y ids
    x_ids <- seq(1:max(col_index_dt$id_x))
    y_ids <- seq(1:max(col_index_dt$id_y))
    
    # convert sf object to dt and append cell indexing, then remove all-NA rows
    forest_sf <- st_as_sf(forest_probs, na.rm=F, as_points = T)
    forest_dt <- as.data.table(forest_sf)
    rm(forest_sf)
    forest_dt[,id_cell := col_index_dt$id_cell]
    forest_dt <- forest_dt[rowSums(!is.na(forest_dt[,..species_names])) > 0, ]
    
    pasture_sf <- st_as_sf(pasture_probs, na.rm=F, as_points=T)
    pasture_dt <- as.data.table(pasture_sf)
    rm(pasture_sf)
    pasture_dt[,id_cell := col_index_dt$id_cell]
    pasture_dt <- pasture_dt[forest_dt[,"id_cell"], on="id_cell"]
    
    # update col_index_dt to only have rows with at least one occupied cell
    col_index_dt <- col_index_dt[forest_dt, on="id_cell", .(id_cell, id_x, id_y)]
    
    # specify spacing grid 
    dt_grid <- expand.grid(id_x = as.integer(seq(1, max(col_index_dt$id_x), point_spacing)), 
                           id_y = as.integer(seq(1, max(col_index_dt$id_y), point_spacing))) %>%
        as.data.table
    
    # join in focal cells (i.e. cells to buffer from)
    # joining on x and y position
    dt_focal_cells <- col_index_dt[dt_grid, on = c("id_x", "id_y")]
    # remove cells without data (id_cell inherit from col_index_dt, which has 
    # all-NA cells removed )
    dt_focal_cells <- dt_focal_cells[!is.na(id_cell),]
    
    # iteratively specify offsets and calculate summary across region
    print(Sys.time() - start.time)
    start.time <- Sys.time()
    out_list <- vector("list", length(d))
    for(i in 1:length(d)) {
        d_i <- d[i]
        
        # specify offsets
        dt_offsets <- as.data.table(expand.grid(id_x_offset = -d_i:d_i, 
                                                id_y_offset = -d_i:d_i, 
                                                id_cell = dt_focal_cells$id_cell))
        
        # join based on id_cell to generate dt with cell x and y coordinates,
        # plus the buffer offsets
        # dt_full is the full set of cells that need forest and pasture probs
        # appending
        # update x and y coordinates with the offsets
        dt_offsets[dt_focal_cells, on="id_cell", `:=`(id_x_offset = id_x + id_x_offset, 
                                                      id_y_offset = id_y + id_y_offset)]
        setnames(dt_offsets, "id_x_offset", "id_x")
        setnames(dt_offsets, "id_y_offset", "id_y")
        setnames(dt_offsets, "id_cell", "id_focal_cell")
        
        # join with col_index_dt to identify NA cells in buffer region and
        # remove
        dt_offsets <- col_index_dt[dt_offsets, on=c("id_x", "id_y")]
        dt_offsets <- dt_offsets[!is.na(id_cell),]
        
        # merge in forest & pasture
        forest_matched <- forest_dt[dt_offsets, on = "id_cell"]
        pasture_matched <- pasture_dt[dt_offsets, on = "id_cell"]
        
        # summarise
        forest_summ <- forest_matched[,lapply(.SD, sum, na.rm=T), by="id_focal_cell", .SDcols=species_names]
        forest_max <- suppressWarnings(
            forest_matched[,lapply(.SD, max, na.rm=T), by="id_focal_cell", .SDcols=species_names]
        )
        
        pasture_summ <- pasture_matched[,lapply(.SD, sum, na.rm=T), by="id_focal_cell", .SDcols=species_names]
        pasture_max <- suppressWarnings(
            pasture_matched[,lapply(.SD, max, na.rm=T), by="id_focal_cell", .SDcols=species_names]
        )
        
        # update pasture_max & forest_max     
        classify_thresh(pasture_max, threshold); classify_thresh(forest_max, threshold)
        above_threshold <- copy(pasture_max == 1 | forest_max == 1)
        above_threshold[above_threshold == FALSE] <- NA
        
        # calculate relative diffs
        pasture_summ <- pasture_summ * above_threshold
        forest_summ <- forest_summ * above_threshold
        rel_diff_mat <- as.matrix(forest_summ/pasture_summ)
        rel_diff_mat <- rel_diff_mat[,!colnames(rel_diff_mat) %in% c("id_focal_cell", "n_cell")]
        
        out_dt <- forest_summ[,c("id_focal_cell")]
        setnames(out_dt, "id_focal_cell", "id_cell")
        out_dt[, `:=`(avg_ratio = matrixStats::rowMeans2(rel_diff_mat, na.rm=T),
                      avg_logratio = matrixStats::rowMeans2(log(rel_diff_mat), na.rm=T),
                      median_logratio = matrixStats::rowMedians(log(rel_diff_mat), na.rm=T))]
        
        # bind geometry info back in
        out_dt <- forest_dt[out_dt, on="id_cell", 
                            .(avg_ratio, avg_logratio, median_logratio, geometry)]
        print(Sys.time() - start.time)
        
        out_list[[i]] <- out_dt %>%
            st_as_sf %>%
            st_rasterize(., dx = 2*point_spacing*1e3, dy=2*point_spacing*1e3)
    }
    if(length(d) == 1) {
        return(out_list[[1]])
    } else {
        do.call(c, c(out_list, list(along="region_size"))) %>%
            st_set_dimensions(., 3, values=(2*d + 1))
    }
}

calc_regional_summ_v4 <- function(forest_probs, pasture_probs, 
                                  d = c(0, 2, 5, 10, 12, 20), point_spacing=5, 
                                  threshold = .1) {
    # note: 
    # - point spacing is expressed in terms of every nth cell to use as a cell
    # for generating regions from
    # - d is expressed in terms of number of cells to buffer around focal cell
    start.time <- Sys.time()
    
    # extract indexing in the stars object 
    col_index_stars <- forest_probs %>%
        slice(1, along="species") %>%
        mutate(id_cell = 1:n(), id_x=NA, id_y = NA) %>%
        select(id_cell, id_x, id_y)
    
    xy_dim <- dim(col_index_stars)
    col_index_stars[["id_x"]] <- matrix(rep(1:xy_dim[2], each=xy_dim[1]), 
                                        nrow=xy_dim[2], ncol=xy_dim[1])
    col_index_stars[["id_y"]] <- matrix(rep(1:xy_dim[1], xy_dim[2]), 
                                        nrow=xy_dim[2], ncol=xy_dim[1])
    
    col_index_sf <- st_as_sf(col_index_stars, as_points = TRUE)
    col_index_dt <- as.data.table(col_index_sf)
    
    # manage memory
    rm(col_index_sf, col_index_stars)
    
    # store x and y ids
    x_ids <- seq(1:max(col_index_dt$id_x))
    y_ids <- seq(1:max(col_index_dt$id_y))
    
    # convert sf object to dt and append cell indexing, then remove all-NA rows
    forest_sf <- st_as_sf(forest_probs, na.rm=F, as_points = T)
    forest_dt <- as.data.table(forest_sf)
    # manage memory
    rm(forest_sf)
    # add cell ids
    forest_dt[,id_cell := col_index_dt$id_cell]
    forest_dt <- forest_dt[rowSums(!is.na(forest_dt[,..species_names])) > 0, ]
    
    pasture_sf <- st_as_sf(pasture_probs, na.rm=F, as_points=T)
    pasture_dt <- as.data.table(pasture_sf)
    # manage memory
    rm(pasture_sf)
    # add cell ids
    pasture_dt[,id_cell := col_index_dt$id_cell]
    pasture_dt <- pasture_dt[forest_dt[,"id_cell"], on="id_cell"]
    
    # update col_index_dt to only have rows with at least one occupied cell
    col_index_dt <- col_index_dt[forest_dt, on="id_cell", .(id_cell, id_x, id_y)]
    
    # specify spacing grid 
    dt_grid <- expand.grid(id_x = as.integer(seq(1, max(col_index_dt$id_x), point_spacing)), 
                           id_y = as.integer(seq(1, max(col_index_dt$id_y), point_spacing))) %>%
        as.data.table
    
    # join in focal cells (i.e. cells to buffer from)
    # joining on x and y position
    dt_focal_cells <- col_index_dt[dt_grid, on = c("id_x", "id_y")]
    # remove cells without data (id_cell inherit from col_index_dt, which has 
    # all-NA cells removed )
    dt_focal_cells <- dt_focal_cells[!is.na(id_cell),]
    
    # iteratively specify offsets and calculate summary across region
    print(Sys.time() - start.time)
    start.time <- Sys.time()
    out_list <- vector("list", length(d))
    for(i in 1:length(d)) {
        d_i <- d[i]
        
        # run in chunks if memory requirements become large
        # approx. mem requirements for forest_matched and pasture_matched 
        # (largest objects in pipeline, by a margin)
        rmem <- (nrow(dt_focal_cells) * (d_i*2 + 1)^2) * 8 * (5 + length(species_names)) * 2
        n_blocks <- ceiling(rmem/1e9/20)
        v <- dt_focal_cells$id_cell
        
        if(n_blocks > 1) {
            chunked_index <- split(v, cut(v, n_blocks))
        } else {
            chunked_index <- list(v)
        }
        out_chunk_list <- vector("list", n_blocks)
        
        for(j in 1:n_blocks) {
            
            # specify offsets
            dt_offsets <- as.data.table(expand.grid(id_x_offset = -d_i:d_i, 
                                                    id_y_offset = -d_i:d_i, 
                                                    id_cell = chunked_index[[j]]))
            
            # join based on id_cell to generate dt with cell x and y coordinates,
            # plus the buffer offsets
            # dt_full is the full set of cells that need forest and pasture probs
            # appending
            # update x and y coordinates with the offsets
            dt_offsets[dt_focal_cells, on="id_cell", `:=`(id_x_offset = id_x + id_x_offset, 
                                                          id_y_offset = id_y + id_y_offset)]
            
            setnames(dt_offsets, "id_x_offset", "id_x")
            setnames(dt_offsets, "id_y_offset", "id_y")
            setnames(dt_offsets, "id_cell", "id_focal_cell")
            
            # join with col_index_dt to identify NA cells in buffer region and
            # remove
            dt_offsets <- col_index_dt[dt_offsets, on=c("id_x", "id_y")]
            dt_offsets <- dt_offsets[!is.na(id_cell),]
            
            # merge in forest & pasture
            forest_matched <- forest_dt[dt_offsets, on = "id_cell"]
            pasture_matched <- pasture_dt[dt_offsets, on = "id_cell"]
            
            # summarise
            forest_summ <- forest_matched[,lapply(.SD, sum, na.rm=T), by="id_focal_cell", .SDcols=species_names]
            forest_max <- suppressWarnings(
                forest_matched[,lapply(.SD, max, na.rm=T), by="id_focal_cell", .SDcols=species_names]
            )
            
            pasture_summ <- pasture_matched[,lapply(.SD, sum, na.rm=T), by="id_focal_cell", .SDcols=species_names]
            pasture_max <- suppressWarnings(
                pasture_matched[,lapply(.SD, max, na.rm=T), by="id_focal_cell", .SDcols=species_names]
            )
            
            # update pasture_max & forest_max     
            classify_thresh(pasture_max, threshold); classify_thresh(forest_max, threshold)
            above_threshold <- copy(pasture_max == 1 | forest_max == 1)
            above_threshold[above_threshold == FALSE] <- NA
            
            # calculate relative diffs
            pasture_summ <- pasture_summ * above_threshold
            forest_summ <- forest_summ * above_threshold
            rel_diff_mat <- as.matrix(forest_summ/pasture_summ)
            rel_diff_mat <- rel_diff_mat[,!colnames(rel_diff_mat) %in% c("id_focal_cell", "n_cell")]
            
            out_dt <- forest_summ[,c("id_focal_cell")]
            setnames(out_dt, "id_focal_cell", "id_cell")
            out_dt[, `:=`(avg_ratio = matrixStats::rowMeans2(rel_diff_mat, na.rm=T),
                          avg_logratio = matrixStats::rowMeans2(log(rel_diff_mat), na.rm=T),
                          median_logratio = matrixStats::rowMedians(log(rel_diff_mat), na.rm=T))]
            
            # bind geometry info back in
            out_dt <- forest_dt[out_dt, on="id_cell", 
                                .(avg_ratio, avg_logratio, median_logratio, geometry)]
            
            out_chunk_list[[j]] <- out_dt
        }
        
        print(Sys.time() - start.time)
        
        out_list[[i]] <- rbindlist(out_chunk_list) %>%
            st_as_sf %>%
            st_rasterize(., dx = 2*point_spacing*1e3, dy=2*point_spacing*1e3)
    }
    
    # return
    if(length(d) == 1) {
        return(out_list[[1]])
    } else {
        do.call(c, c(out_list, list(along="region_size"))) %>%
            st_set_dimensions(., 3, values=(2*d + 1))
    }
}




