# Script Details ---------------------------------------------------------------
###project: Exploratory HMSC Models for OEPA Dataset
###author(s): MW
###goal(s): 
###date(s): May 2026
###note(s):

# Housekeeping -----------------------------------------------------------------
### load necessary libraries ----
# install.packages("librarian")
librarian::shelf(tidyverse, readxl, scales, broom, purrr, dataRetrieval,
                 mvabund, vegan, cowplot, patchwork, stringr, splitstackshape, 
                 readxl, janitor, skimr, Hmsc, ape, taxize)

### define custom functions ----
nacheck <- function(df) {
      na_count_per_column <- sapply(df, function(x) sum(is.na(x)))
      print(na_count_per_column)
}


# Environmental ----------------------------------------------------------------
### load necessary data --------------------------------------------------------
lu <- read_csv('localdata/landuse.csv') |> 
      janitor::clean_names() |> 
      mutate(across(open_water:wetlands, ~ .x / 100)) |> 
      mutate(barriers = rowSums(across(starts_with("segment_")))) |> 
      select(gage_year, barren_land:wetlands, watershed_area = area_m2, barriers) |> 
      mutate(gage_year = paste0('0', gage_year, sep = ''))
glimpse(lu)
skim(lu)

hy <- read_csv('localdata/hydro.csv') |> 
      janitor::clean_names() |> 
      mutate(log_magnitude = log(magnitude)) |>
      select(gage_year, site_no, lat, lon, log_magnitude, magnitude, frequency, duration, timing, cv) |> 
      group_by(site_no) |> 
      mutate(across(log_magnitude:cv, ~as.numeric(scale(.x)), .names = "scaled_site_{.col}")) |> 
      ungroup() |> 
      mutate(across(log_magnitude:cv, ~as.numeric(scale(.x)), .names = "scaled_across_{.col}"))
glimpse(hy)
skim(hy)

env <- hy |> left_join(lu)

# Biotic -----------------------------------------------------------------------
### load necessary data --------------------------------------------------------
fi <- read_csv('localdata/fish.csv') |> 
      janitor::clean_names() |> 
      
      # add leading zero to gage_year, lost along the way somehow
      mutate(gage_year = paste0('0', gage_year, sep = '')) |> 
      
      # remove nonstandard and longline data (lost 6.4% of data)
      filter(geartype %in% c('BOAT', 'WADE')) |> 
      
      # remove samples based on CAP filter 1 (lost 0% of remaining data)
      filter(!(dist < 100 & da < 200)) |> 
      
      # remove samples based on CAP filter 2 (lost 1.2% of remaining data)
      filter(!(dist < 300 & da > 200 & gear == "A")) |> 
      
      # remove hybrid species (lost 2.9% of remaining data)
      filter(!grepl(" x ", fishname)) |> 
      
      # retain desired columns
      select(
            gage_year, 
            site = storet, 
            sheet,
            station = station_id, 
            river, 
            rivercode, 
            huc, 
            lat = latdd, 
            lon = longdd,
            fishname, 
            counts, 
            weight,
            dist
      ) |> 
      
      # create some new columns
      mutate(
            sample   = paste0(site, sheet, sep = ''),
            cpue     = counts/dist,
            presence = case_when(
                  counts >= 1 ~ 1,
                  TRUE ~ 0
            )
      ) |> 
      
      # zero-fill the dataset
      complete(
            nesting(gage_year, site, sample, lat, lon, dist),
            fishname,
            fill = list(counts   = 0,
                        presence = 0, 
                        cpue     = 0)
      ) |> 
      
      # keep desired columns for moving forward
      dplyr::select(gage_year, site, sample, lat, lon, dist, common_name = fishname, cpue, counts, presence)
skim(fi)

spp_list <- fi |> select(common_name) |> distinct() |> arrange() |> mutate(common_name = common_name |> str_trim())
skim(spp_list)
glimpse(spp_list)

### grab taxonomic information all species in dataset --------------------------
sci_lookup <- comm2sci(
      com = spp_list$common_name,
      db = 'itis',
      simplify = FALSE
)
tsns <- get_tsn(spp_list$common_name, searchtype = "common")

# Traits -----------------------------------------------------------------------
### load necessary data --------------------------------------------------------
tr <- read_csv('localdata/traits.csv') |> 
      janitor::clean_names()
skim(tr)
