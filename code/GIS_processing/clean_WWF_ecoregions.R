# Tidy WWF ecoregions; plot map and save output. 
# Trim WWF ecoregions to the limits of Colombia, select regions to predict over, 
# simplify geometry of Amazonian regions, and clip unforested eastern region of
# llanos. 

# packages
library(dplyr); library(sf); library(ggplot2); library(stars)

AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 
    +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# get Colombia outline
CO_low <- rnaturalearth::ne_countries(country = "Colombia", 
                                      scale = "small", 
                                      returnclass = "sf")
CO <- rnaturalearth::ne_countries(country = "Colombia", 
                                  scale = "large", returnclass = "sf") %>%
    st_crop(., CO_low)

# read WWF ecoregions
WWF_regions <- read_sf("data/official/wwf_terr_ecos.shp") %>%
    st_make_valid() %>%
    st_crop(CO)

# crop WWF ecoregions to Colombia borders
WWF_CO <- WWF_regions %>%
    st_intersection(CO)

# get bird survey points
birds <- readRDS("outputs/birds.RDS")
point_locs <- birds %>%
    select(site, point, lon, lat) %>%
    unique %>%
    st_as_sf(., coords=c("lon", "lat"), crs="wgs84")

# subset ecoregions to those we are going to project over
# sampled_within <- t2 %>%
#     filter(rowSums(st_intersects(., point_locs, sparse = F)) != 0) %>%
#     pull(ECO_NAME) %>%
#     c(., "Caqueta moist forests", "Apure-Villavicencio dry forests", 
#       "Magdalena Valley dry forests")

region_subset <- c("Apure-Villavicencio dry forests", 
                   "Caqueta moist forests", 
                   "Cauca Valley montane forests",
                   "Cordillera Oriental montane forests",
                   "Eastern Cordillera real montane forests",
                   "Llanos",
                   "Magdalena Valley dry forests", 
                   "Magdalena Valley montane forests",
                   "Napo moist forests", 
                   "Northern Andean páramo", 
                   "Northwestern Andean montane forests", 
                   "Santa Marta montane forests", 
                   # possibly need to do something with this last one?
                   "Magdalena-Urabá moist forests")

WWF_subset <- WWF_CO %>%
    filter(ECO_NAME %in% region_subset) %>%
    group_by(ECO_NAME) %>%
    summarise(geometry = st_union(geometry)) %>%
    st_transform(., crs=AEAstring)

# save first draft of ecoregions
WWF_subset %>%
    ggplot() + geom_sf(aes(fill=ECO_NAME), col=NA) +
    geom_sf(data=point_locs) +
    geom_sf(data=CO, fill=NA) +
    # scale_fill_viridis_d() +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"), 
          panel.grid = element_blank())
ggsave("figures/map_WWF_terr_ecoregions_first_draft.png")

# Tidying polygons ----
# simplify geometry of Amazonian ecoregions & trim Llanos to only include 
# western portions
Amazon_regions <- c("Caqueta moist forests", "Napo moist forests")

WWF_subset2 <- WWF_subset %>%
    filter(!(ECO_NAME %in% Amazon_regions))


Amazon_outline <- WWF_subset %>%
    filter(ECO_NAME %in% Amazon_regions) %>%
    st_buffer(., 10) %>%
    st_convex_hull()

catch_region <- c()
for(region in Amazon_regions) {
    Amazon_region_outline <- WWF_subset %>%
        filter(ECO_NAME == region) %>%
        st_buffer(., 10) %>%
        st_convex_hull()
    
    Amazon_region <- WWF_subset %>%
        filter(ECO_NAME == region) 
    
    diff_region <- st_difference(Amazon_region_outline, Amazon_region)
    catch_region[[region]] <- st_cast(diff_region, "POLYGON") %>%
        mutate(area = st_area(geometry)) %>%
        filter(area == max(area)) %>%
        st_difference(Amazon_region_outline, .)
}

WWF_subset3 <- bind_rows(catch_region) %>%
    bind_rows(WWF_subset2, .) %>%
    select(-(ECO_NAME.2:area))

WWF_subset3 %>%
    ggplot() + geom_sf(aes(fill=ECO_NAME), col=NA) +
    geom_sf(data=point_locs) +
    geom_sf(data=CO, fill=NA) +
    # scale_fill_viridis_d() +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"), 
          panel.grid = element_blank())

## clip Llanos ----
# for now just arbitrarily clip at 70 degrees west
clip_coords <- c(xmin = -70.5, ymin=3, xmax = -60, ymax=20)
class(clip_coords) <- "bbox"    
clip <- st_as_sfc(clip_coords) %>%
    st_as_sf(.,  crs="wgs84") %>%
    rename(geometry = x) %>%
    st_transform(., st_crs(WWF_subset3))


Llanos_crop <- st_read("inputs/E_llanos_2.kml") %>%
    st_transform(., st_crs(WWF_subset3))

Llanos_cropped <- WWF_subset3 %>%
    filter(ECO_NAME == "Llanos") %>%
    st_difference(., smoothr::smooth(Llanos_crop)) %>%
    select(-Name, -Description) %>%
    st_zm(., drop=TRUE)
    

WWF_subset4 <- bind_rows(WWF_subset3 %>% filter(ECO_NAME != "Llanos"), 
                         Llanos_cropped)

# check
ggplot(WWF_subset4) + geom_sf()

# plot final map ----
# save first draft of ecoregions
set.seed(1)
WWF_subset4 %>%
    mutate(ECO_NAME = factor(ECO_NAME, levels = sample(ECO_NAME))) %>%
    ggplot() + 
    geom_sf(data=CO, fill="grey90", col="black") +
    geom_sf(aes(fill=ECO_NAME), col="black", alpha=.5) +
    geom_sf(data=point_locs, pch=21, fill="white", size=2) +
    # scale_fill_manual(values = cols) +
    scale_fill_viridis_d() +
    scale_x_continuous(breaks=seq(-80, 60, 5)) +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"), 
          panel.grid = element_blank(), 
          legend.position = "bottom") +
    labs(fill = "") +
    guides(fill = guide_legend(ncol=3))
ggsave("figures/map_WWF_terr_ecoregions.png")

saveRDS(WWF_subset4, "outputs/WWF_terrestrial_ecoregions.rds")

WWF_subset4 %>%
    summarise(geometry = st_union(geometry)) %>%
    saveRDS(., "outputs/WWF_terr_ecoregions_mask.rds")
