This repository contains code for analysis of bird and dung beetle data from Colombia.

The core scripts required for analysis are listed below. **Bold** headings organize the codebase by topic.  Solid bullets correspond to code files.  Hollow sub-bullets indicate each script that produces objects required for the script to run (generally objects saved to disk, except in cases of scripts called via `source()`).  Scripts can be run in order of first appearance in this list without problems.

#### Baseline species lists
* species_lists.R   (baseline species list for Colombia; also creates birdlife range map objects)
* nf_species_list.R   (baseline species list for neotropical forests)
    * species_lists.R

#### Sampling point data
*	GEE_setup.sh   (setup for Earth Engine)
*	coord_processing.R   (import and process coordinates for all points)
    *	GEE_setup.sh
*	points_formatting.R    (Add point-level covariates, including ALOS elevations from Earth Engine via reticulate; *still need to do precipitation and wildlife-friendliness*)
    * GEE_setup.sh
    * coord_processing.R

#### Survey data
* bird_import_and_cleaning.R   (Import and format bird survey data)
    * species_lists.R
    * points_formatting.R

#### Species ranges
*	hydrosheds_extraction.R   (create polygons for biogeographic clipping)
*	ayerbe_maps.R   (read and process Ayerbe range maps)
    * species_lists.R
*	combined_bird_maps.R   (create properly buffered maps--originally these were combined Ayerbe-BirdLife maps, but now they are just biogeographically buffered Ayerbe maps).
    * species_lists.R
    * bird_import_and_cleaning.R
    * hydrosheds_extraction.R
    * ayerbe_maps.R

#### Species traits
*	birdlife_scraper.R   (scrape birdlife website for trait info)
    * nf_species_list.R
*	parker_standardization.R  (update Parker-Stotz-Fitzpatrick databases to BirdLife taxonomy *still requires some additional checking by hand*)
    * nf_species_list.R
    * birdlife_scraper.R
    * species_lists.R
*	elevations_prep_and_exploration.R   (read in elevational range limits file and explore data)
    * species_lists.R
*	species_covariate_formatting.R   (create file with species covariates)
    * species_lists.R
    * parker_standardization.R
    * combined_bird_maps.R
    * elevations_prep_and_exploration.R
    * birdlife_scraper.R
*	migratory_dates.R   (read in file of migratory dates and format)
    * species_lists.R

#### Format for analysis
*	format_for_analysis.R   (merge various covariate files, also create new covariates from certain files and range maps)
    * combined_bird_maps.R
    * bird_import_and_cleaning.R
    * elevations_prep_and_exploration.R
    * points_formatting.R
    * species_covariate_formatting.R
    * migratory_dates.R

#### Analysis
*	occupancyMod_full.stan   (stan code for biogeographically clipped multi-species occupancy model with within-chain parallelization)
*	bird_occupancy_full.R   (fit the occupancy model described in occupancyMod_full.R via cmdstanR on local computer)
    * format_for_analysis.R
