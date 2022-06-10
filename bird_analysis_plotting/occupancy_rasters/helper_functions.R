# functions

# sf species x point occupancies to matrix ----
# note: probably redundant- remove at a later point
sf_to_mat <- function(x) {
    as_tibble(x) %>%
        dplyr::select(-geometry) %>%
        as.matrix
}

# tidy dt after thresholding ----
# note: probably redundant- remove at a later point
classify_thresh <- function(DT, threshold=0.1) {
    for (i in names(DT)) {
        DT[is.infinite(get(i)), (i):=NA]
        DT[get(i) < threshold, (i):=0]
        DT[get(i) >= threshold, (i):=1]
    }
}

# trim NAs from edge of stars ----
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

# fns to swap cell ids and xy positions ----
# rewrite for arbitrary dims (e.g. if change resolution at a later point)?
cell_to_y_pos <- function(id_cell, xy_dim = c(924, 679)) {
    ceiling(id_cell/xy_dim[2])
    
}

cell_to_x_pos <- function(id_cell, xy_dim = c(924, 679)) {
    x <- id_cell %% xy_dim[2]
    x[x==0] <- xy_dim[2]
    x
}

xy_to_cell <- function(id_x, id_y, xy_dim = c(924, 679)) {
    (id_y - 1) * xy_dim[2] + id_x
}

# calculate regional ratios ----
calc_regional_summ_v7 <- function(pred_info, buffer = c(2, 5, 10, 15), 
                                  point_spacing=5, 
                                  threshold = .1, 
                                  xy_lookup = xy_lookup) {
    # note: 
    # - point spacing is expressed in terms of every nth cell to buffer from to 
    # generate regions (1 would be at a 2km resolution)
    # - buffer is expressed in terms of number of cells to buffer around focal 
    # cell (2 is equivalent to a 10x10km region, i.e. 25 points)
    
    # calc_OR <- function(p_forest, p_pasture) (p_forest/(1-p_forest))/(p_pasture/(1 - p_pasture))
    # calc_RR <- function(p_forest, p_pasture) p_forest/p_pasture
    
    # store x and y ids
    xy_dim <- c(max(xy_lookup$id_x), max(xy_lookup$id_y))
    x_ids <- seq(1:xy_dim[1])
    y_ids <- seq(1:xy_dim[2])
    
    # specify spacing grid 
    dt_grid <- expand.grid(id_x = as.integer(seq(1, xy_dim[1], point_spacing)), 
                           id_y = as.integer(seq(1, xy_dim[2], point_spacing))) %>%
        as.data.table
    
    # join in focal cells (i.e. cells to buffer from)
    # joining on x and y position
    dt_grid[xy_lookup, on=c("id_x", "id_y"), id_cell := id_cell]
    
    # get cells that contain values
    nonNA_cells <- unique(pred_info[,id_cell])
    
    # note: pred_info only contains cells with some data
    dt_focal_cells <- dt_grid[id_cell %in% nonNA_cells, ]
    
    # iteratively specify offsets and calculate summary across region
    start.time <- Sys.time()
    out_list <- vector("list", length(buffer))
    for(i in 1:length(buffer)) {
        print(paste0("buffer:", buffer[i]))
        buffer_i <- buffer[i]
        
        # run in chunks if memory requirements become large
        # note: calculation now out of date, but seems to prevent running out 
        # of RAM on my device. Can relax on HPC. 
        rmem <- (nrow(dt_focal_cells) * (buffer_i*2 + 1)^2) * 8 * (5 + 1614) * 2
        n_blocks <- ceiling(rmem/1e9/10)
        v <- dt_focal_cells$id_cell
        
        if(n_blocks > 1) {
            chunked_index <- split(v, cut(v, n_blocks))
        } else {
            chunked_index <- list(v)
        }
        out_chunk_list <- vector("list", n_blocks)
        for(j in 1:n_blocks) {
            # specify offsets
            dt_offsets <- as.data.table(expand.grid(id_x_offset = -buffer_i:buffer_i, 
                                                    id_y_offset = -buffer_i:buffer_i, 
                                                    id_cell = chunked_index[[j]]))
            
            # join based on id_cell to generate dt with cell x and y coordinates,
            # plus the buffer offsets
            dt_offsets[dt_focal_cells, on="id_cell", `:=`(id_x_offset = id_x + id_x_offset, 
                                                          id_y_offset = id_y + id_y_offset)]
            
            setnames(dt_offsets, "id_x_offset", "id_x")
            setnames(dt_offsets, "id_y_offset", "id_y")
            setnames(dt_offsets, "id_cell", "id_focal_cell")
            
            # exchange x and y for cell info (for merge below)
            dt_offsets[xy_lookup, on=c("id_x", "id_y"), id_cell := id_cell]
            dt_offsets[,`:=`(id_x = NULL, id_y = NULL)]
            
            # remove any NA cells in buffer region
            dt_offsets <- dt_offsets[id_cell %in% nonNA_cells,]
            
            # calculate number of nonNA cells in each regional buffer
            n_cell <- dt_offsets[,list(n_cell = length(unique(id_cell))), by="id_focal_cell"]
            
            # note: code is somewhat unwieldy from here on out. We need to join the
            # species-level forest and pasture probs, and rel_diffs, into the 
            # dt_offsets dt
            dt_offsets <- pred_info[dt_offsets, on = "id_cell", allow.cartesian=T]
            
            # ..then compute the max prob for each species x region
            dt_max <- dt_offsets[, list(
                max_prob = max(c(max(p_forest, na.rm=T), max(p_pasture, na.rm=T)))
            ), by=c("id_focal_cell", "species")]
            
            # ..trim lookup df to only species that exceed the threshold
            dt_max <- dt_max[max_prob > threshold,]
            
            # manage memory
            dt_max[,max_prob := NULL]
            
            # ..remove species not in region (i.e. not above threshold)
            SR_summ <- dt_offsets[dt_max, on = c("id_focal_cell", "species")
            ][, list(sr_pt = .N), by = c("id_cell", "id_focal_cell")]
            
            SR_summ2 <- SR_summ[,list(mean_sr_pt = mean(sr_pt)), by="id_focal_cell"]
            
            # calculate beta
            SR_gamma <- dt_offsets[dt_max, on = c("id_focal_cell", "species")
            ][, list(sr_region = length(unique(species))), 
              by = "id_focal_cell"]
            
            SR_tot <- SR_gamma[SR_summ2, on="id_focal_cell"]
            SR_tot[, beta := sr_region/mean_sr_pt]
            
            # manage memory
            rm(dt_max)
            
            # summarise log ratios, by (1) region, and (2) points (within regions)
            ## 1: calculate the number of occupied points within a region for 
            ## forest and pasture
            dt_summ <- dt_offsets[, list(
                sum_forest = sum(p_forest), 
                sum_pasture = sum(p_pasture)
            ), by=c("id_focal_cell", "species")]
            
            ## .. calculate the median log-ratio across species within each region
            dt_summ <- dt_summ[, list(
                median_log_ratio = median(log(sum_forest/sum_pasture), na.rm=T)
            ), by=c("id_focal_cell")]
            
            ## 2: calculate the median-sensitivity species on each point
            ## note: have to calculate medians here rather than at the outset 
            ## as need to threshold first
            dt_summ_pt <- dt_offsets[
                , list(log_ratio = median(log(p_forest/p_pasture), na.rm=T)), 
                by=c("id_focal_cell", "id_cell")]
            
            ## .. calculate the average (mean) sensitivity across points within a 
            # region
            dt_summ_pt <- dt_summ_pt[, list(mean_pt_log_ratio = mean(log_ratio, na.rm=T)), 
                                     by=c("id_focal_cell")]
            
            # join together n_cell, sr, and loss sensitivity dts
            out_chunk_list[[j]] <- dt_summ[n_cell, on = "id_focal_cell"
            ][SR_tot, on="id_focal_cell"
            ][dt_summ_pt, on = "id_focal_cell"]
            
            # manage memory
            rm(dt_offsets); gc()
        }
        
        print(Sys.time() - start.time)
        
        out_list[[i]] <- rbindlist(out_chunk_list)
    }
    out_list
}

# memory efficient row removal ----
# note: not currently used
remove_rows <- function(dt, threshold) {
    cols <- names(dt)
    above_threshold <- dt$p_forest > threshold | dt$p_pasture > threshold
    dt_subset = data.table(temp_name = dt[[cols[1]]][above_threshold])
    setnames(dt_subset, "temp_name", cols[1])
    for (i in 1:length(cols)){
        col <- cols[1]
        dt_subset[, (col) := dt[[col]][above_threshold]]
        dt[, (col) := NULL] #delete
    }
    return(dt_subset)
}
