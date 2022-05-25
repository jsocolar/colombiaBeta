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
cell_to_y_pos <- function(id_cell) {
    ceiling(id_cell/679)
    
}

cell_to_x_pos <- function(id_cell) {
    x <- id_cell %% 679
    x[x==0] <- 679
    x
}

xy_to_cell <- function(id_x, id_y) {
    (id_x-1) * 679 + id_y
}

# calculate regional ratios ----
calc_regional_summ_v5 <- function(pred_info,
                                  d = c(0, 2, 5, 10), point_spacing=5, 
                                  threshold = .1) {
    # note: 
    # - point spacing is expressed in terms of every nth cell to use as a cell
    # for generating regions from
    # - d is expressed in terms of number of cells to buffer around focal cell
    start.time <- Sys.time()
    # d = c(0, 2, 5, 10, 12, 20)
    # point_spacing=10
    # threshold = .1
    
    # extract indexing in the stars object 
    # xyinfo <- readRDS("outputs/xy_info_lookup.rds")
    xy_dim <- c(679, 925) #dim(col_index_stars)
    
    # store x and y ids
    x_ids <- seq(1:679)
    y_ids <- seq(1:925)
    
    # specify spacing grid 
    dt_grid <- expand.grid(id_x = as.integer(seq(1, xy_dim[2], point_spacing)), 
                           id_y = as.integer(seq(1, xy_dim[1], point_spacing))) %>%
        as.data.table
    
    # join in focal cells (i.e. cells to buffer from)
    # joining on x and y position
    dt_grid[,id_cell := xy_to_cell(id_x, id_y)]
    
    # get cells that contain values
    nonNA_cells <- unique(pred_info[,id_cell])
    
    # note: pred_info only contains cells with some data
    dt_focal_cells <- dt_grid[id_cell %in% nonNA_cells, ]
    
    # iteratively specify offsets and calculate summary across region
    print(Sys.time() - start.time)
    start.time <- Sys.time()
    out_list <- vector("list", length(d))
    for(i in 1:length(d)) {
        d_i <- d[i]
        # d_i <- 2
        # run in chunks if memory requirements become large
        # approx. mem requirements for forest_matched and pasture_matched 
        # (largest objects in pipeline, by a margin)
        rmem <- (nrow(dt_focal_cells) * (d_i*2 + 1)^2) * 8 * (5 + 1614) * 2
        n_blocks <- ceiling(rmem/1e9/20)
        v <- dt_focal_cells$id_cell
        
        if(n_blocks > 1) {
            chunked_index <- split(v, cut(v, n_blocks))
        } else {
            chunked_index <- list(v)
        }
        out_chunk_list <- vector("list", n_blocks)
        for(j in 1:n_blocks) {
            # j <- 1
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
            
            # exchange x and y for cell info (for merge below)
            dt_offsets[,id_cell := xy_to_cell(id_x, id_y)]
            dt_offsets[,`:=`(id_x = NULL, id_y = NULL)]
            
            # remove any NA cells in buffer region
            dt_offsets <- dt_offsets[id_cell %in% nonNA_cells,]
            
            # calculate number of cells in each regional buffer
            n_cell <- dt_offsets[,list(n_cell = length(unique(id_cell))), by="id_focal_cell"]
            
            # note: could briefly save a bit of memory here (~4GB) by overwriting
            # base pred_info, but not obviously worth it for the additional 
            # complications it introduces (e.g. read times)
            # also, important to compute n_cell before thresholding
            # pred_info_local <- pred_info[p_forest > threshold | p_pasture > threshold, ]
            
            # merge in forest & pasture
            dt_offsets <- pred_info[dt_offsets, on = "id_cell", allow.cartesian=TRUE]
            
            dt_max <- dt_offsets[,list(max_forest = max(p_forest, na.rm=T), 
                                       max_pasture = max(p_pasture, na.rm=T)), 
                                 by=c("id_focal_cell", "species")]
            dt_max <- dt_max[max_forest > threshold | max_pasture > threshold,]
            
            # remove species not in region (i.e. not above threshold)
            dt_offsets <- dt_offsets[dt_max[,.(id_focal_cell, species)], 
                                     on = c("id_focal_cell", "species")]
            
            # calculate beta
            SR_summ <- dt_offsets[,list(sr_pt = .N), by = c("id_cell", "id_focal_cell")]
            SR_summ2 <- SR_summ[,list(mean_sr_pt = mean(sr_pt)), by="id_focal_cell"]
            
            SR_gamma <- dt_offsets[,list(sr_region = length(unique(species))), 
                                   by = c("id_focal_cell")]
            
            SR_tot <- SR_gamma[SR_summ2, on="id_focal_cell"]
            SR_tot[, beta := 1 - mean_sr_pt/sr_region]
            
            # manage memory
            dt_offsets[,id_cell := NULL]
            rm(dt_max)
            
            # summarise log ratios
            dt_summ <- dt_offsets[, list(
                median_log_ratio = median(sum(p_forest)/sum(p_pasture), na.rm = T)
            ), by=c("id_focal_cell")]
            
            out_chunk_list[[j]] <- dt_summ[n_cell, on = "id_focal_cell"][SR_tot, on="id_focal_cell"]
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
