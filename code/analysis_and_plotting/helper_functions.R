# functions
# funs 
extract_stars_indexing <- function(stars_object) {
    # drop 3rd dimension, if exists
    if(length(dim(stars_object )) != 2) stars_object <- stars_object[,,,1]
    
    col_index_stars <- stars_object %>%
        mutate(id_cell = 1:n(), id_x=NA, id_y = NA) %>%
        select(id_cell, id_x, id_y)
    
    xy_dim <- dim(col_index_stars)
    col_index_stars[["id_x"]] <- matrix(rep(1:xy_dim[2], each=xy_dim[1]), 
                                        nrow=xy_dim[2], ncol=xy_dim[1])
    col_index_stars[["id_y"]] <- matrix(rep(1:xy_dim[1], xy_dim[2]), 
                                        nrow=xy_dim[2], ncol=xy_dim[1])
    
    col_index_sf <- st_as_sf(col_index_stars, as_points = TRUE)
    as.data.table(col_index_sf)
}

# 
calculate_grid <- function(xy_dim, point_spacing, nonNA_cells, xy_lookup) {
    # specify spacing grid (each cell is a 'focal cell' to buffer from)
    dt_grid <- expand.grid(id_x = as.integer(seq(1, xy_dim[1], point_spacing)), 
                           id_y = as.integer(seq(1, xy_dim[2], point_spacing))) %>%
        as.data.table
    
    # append cell id from each cell's x and y position via join with xy_lookup
    dt_grid[xy_lookup, on=c("id_x", "id_y"), id_cell := id_cell]
    
    # trim dt_grid to nonNA_cells via inner join
    nonNA_cells[dt_grid, on="id_cell", nomatch=0]
}

calculate_offsets <- function(id_cells, dt_focal_cells, buffer,nonNA_cells, 
                              xy_lookup) {
    dt_offsets <- as.data.table(expand.grid(id_x_offset = -buffer:buffer, 
                                            id_y_offset = -buffer:buffer, 
                                            id_cell = id_cells))
    
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
    
    # remove any NA cells in buffer region and return
    dt_offsets[nonNA_cells, on="id_cell", nomatch=0]
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
    -ceiling(id_cell/xy_dim[2])
    
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
                                  threshold_base = .01,
                                  recompute_base = TRUE,
                                  xy_lookup = xy_lookup) {
    # note: 
    # - point spacing is expressed in terms of every nth cell to buffer from to 
    # generate regions (1 would be at a 2km resolution)
    # - buffer is expressed in terms of number of cells to buffer around focal 
    # # cell (2 is equivalent to a 10x10km region, i.e. 25 points)
    # point_spacing = 20
    # threshold = 1
    # buffer = 2
    # i <- 1
    # j <- 1
    
    # calculate expected number of forest points within 2x2km square
    pred_info[, `:=`(n_forest = p_forest * (2000^2/(pi*100^2)), 
                     n_pasture = p_pasture * (2000^2/(pi*100^2)))]
    
    # exclude extremely low probabilities from computation
    pred_info[!(p_forest > threshold_base | p_pasture > threshold_base), 
              `:=`(n_forest = NA, n_pasture=NA)]
    
    # for species that occur on at lesat 
    pred_info[(n_forest > threshold | n_pasture > threshold),`:=`(log_rd = log(n_forest/n_pasture))]
    
    # store x and y ids
    xy_dim <- c(max(xy_lookup$id_x), max(xy_lookup$id_y))
    x_ids <- seq(1:xy_dim[1])
    y_ids <- seq(1:xy_dim[2])
    
    # get cells that contain values (pred_info does not contain NA cells)
    nonNA_cells <- unique(xy_lookup[,list(id_cell)])
    
    # sr_pt <- pred_info[!is.na(log_rd), .(sr_pt = .N), by="id_cell"]
    
    # trim dt_grid to nonNA_cells via inner join
    dt_focal_cells <- calculate_grid(xy_dim, point_spacing, nonNA_cells, xy_lookup)
    
    dt_pt_summ <- pred_info[,list(median_log_rd_pt = median(log_rd, na.rm=T),
                                  prop_winners = sum(log_rd < 0, na.rm=T)/sum(!is.na(log_rd)),
                                  sr_pt = sum(!is.na(log_rd))), by="id_cell"]
    
    # iteratively specify offsets and calculate summary across region
    start.time <- Sys.time()
    out_list <- vector("list", length(buffer))
    
    for(i in 1:length(buffer)) {
        buffer_i <- buffer[i]
        print(paste0("buffer:", buffer_i))
        
        # run in chunks if memory requirements become large
        # note: calculation now out of date, but seems to prevent running out 
        # of RAM on my device. Can relax on HPC. 
        # rmem <- (nrow(dt_focal_cells) * (buffer_i*2 + 1)^2) * 8 * (5 + 1614) * 2
        rmem <- (nrow(dt_focal_cells) * (buffer_i*2 + 1)^2)
        n_blocks <- ceiling(rmem/4e5)
        v <- dt_focal_cells$id_cell
        
        if(n_blocks > 1) {
            chunked_index <- split(v, cut(v, n_blocks))
        } else {
            chunked_index <- list(v)
        }
        out_chunk_list <- vector("list", n_blocks)
        
        for(j in 1:n_blocks) {
            print(paste0("chunk ", j, " of ", n_blocks))
            
            # buffer_i <- 10
            # calculate offsets for cell ids in the jth chunk
            dt_offsets <- calculate_offsets(chunked_index[[j]], dt_focal_cells, 
                                            buffer = buffer_i, nonNA_cells, 
                                            xy_lookup)
            
            # calculate number of nonNA cells in each regional buffer
            n_cell <- dt_offsets[,list(n_cell = length(unique(id_cell))), 
                                 by="id_focal_cell"]
            
            pt_summ <- dt_offsets[dt_pt_summ, on="id_cell", nomatch=0]
            pt_summ <- pt_summ[, .(avg_log_ratio = mean(median_log_rd_pt), 
                                   avg_prop_winners = mean(prop_winners), 
                                   avg_sr_pt = mean(sr_pt)), by="id_focal_cell"]
            
            # Join the species-level forest and pasture probs into dt_offsets
            # (note: nomatch = 0, because we don't need to include rows where a 
            # species has an NA probability (i.e. doesn't have a row in pred_info)
            dt_offsets <- pred_info[dt_offsets, on = "id_cell",
                                    allow.cartesian=T, nomatch=0]
            
            dt_offsets[,mean(log_rd, na.rm=T), by=c("id_cell")]
            
            ## Summarising regional and point-level species richness
            # Compute the max prob for each species x region, to threshold by 
            # in the next step (this is the slow bit)
            dt_max <- dt_offsets[, list(n_forest = sum(n_forest, na.rm=T), 
                                        n_pasture = sum(n_pasture, na.rm=T)
            ), by=c("id_focal_cell", "species")]
            
            # ..trim lookup df to only species that exceed the threshold
            dt_max <- dt_max[n_forest > threshold | n_pasture > threshold,]
            dt_max[,log_rd := log(n_forest/n_pasture)]
            
            dt_region_summ <- dt_max[,list(med_log_ratio_region = median(log_rd), 
                                           avg_log_ratio_region = mean(log_rd), 
                                           prop_winners = sum(log_rd < 0)/.N,
                                           sr_region = .N), by="id_focal_cell"]
            
            #dt_region_summ <- dt_region_summ[mean_n_sp_region, on="id_focal_cell"]
            
            dt_out <- dt_region_summ[pt_summ, on="id_focal_cell",]
            out_chunk_list[[j]] <- dt_out
            
        }
        
        print(Sys.time() - start.time)
        out_list[[i]] <- rbindlist(out_chunk_list)
    }
    out <- rbindlist(out_list, idcol="id")
    labs <- paste0((buffer * 2 + 1)*2, "km")
    out[,labs := factor(labs[id], levels=labs)]
    out
}

# redundant functions ----
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

