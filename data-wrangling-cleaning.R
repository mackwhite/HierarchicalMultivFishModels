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
                 readxl, janitor, skimr, Hmsc, ape, taxize, rfishbase, 
                 GGally, corrplot, rotl, Rphylopars, ape, fishtree)

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

spp_list <- fi |> select(common_name) |> distinct() |> arrange() |> mutate(common_name = common_name |> str_trim()) |> mutate(type = 'og')
skim(spp_list)
glimpse(spp_list)

### grab taxonomic information all species in dataset --------------------------

#### clean up taxonomic information --------------------------------------------
sci_names <- rfishbase::common_to_sci(spp_list$common_name) |> 
      janitor::clean_names() |> 
      rename(
            scientific_name = species,
            common_name     = com_name
      ) |> 
      select(scientific_name, common_name) |> 
      distinct()

dups <- sci_clean |>
      group_by(common_name) |>
      mutate(n_sci_per_common = n()) |>
      group_by(scientific_name) |>
      mutate(n_common_per_sci = n()) |>
      ungroup() |>
      mutate(is_dup = n_sci_per_common > 1 | n_common_per_sci > 1) |> 
      filter(is_dup) |>
      arrange(common_name, scientific_name)

sci_clean <- sci_names |> 
      mutate(
            scientific_name = case_when(
                  common_name == 'Gizzard shad'        ~ 'Dorosoma cepedianum',
                  common_name == 'Rock bass'           ~ 'Ambloplites rupestris',
                  common_name == 'Brown trout'         ~ 'Salmo trutta',
                  common_name == 'Grass carp'          ~ 'Ctenopharyngodon idella',
                  common_name == 'Silver chub'         ~ 'Macrhybopsis storeriana',
                  common_name == 'Spotted bass'        ~ 'Micropterus punctulatus',
                  common_name == 'White perch'         ~ 'Morone americana',
                  common_name == 'Bowfin'              ~ 'Amia calva',
                  common_name == 'Central stoneroller' ~ 'Campostoma anomalum',
                  common_name == 'Creek chub'          ~ 'Semotilus atromaculatus',
                  common_name == 'Goldfish'            ~ 'Carassius auratus',
                  common_name == 'Paddlefish'          ~ 'Polyodon spathula',
                  common_name == 'Sea lamprey'         ~ 'Petromyzon marinus',
                  common_name == 'Silver shiner'       ~ 'Notropis photogenis',
                  TRUE ~ scientific_name
            )
      ) |> 
      
      # capitalize all common names 
      mutate(common_name = str_to_title(common_name)) |> 
      distinct()

final_spp_list <- spp_list |>
      left_join(sci_clean) |> 
      mutate(common_name = case_when(
            common_name == 'Warmouth Sunfish'    ~ 'Warmouth',
            common_name == 'North Brook Lamprey' ~ 'Northern Brook Lamprey',
            TRUE ~ common_name
      )) |> 
      mutate(scientific_name = case_when(
            common_name == 'Warmouth'               ~ 'Lepomis gulosus',
            common_name == "Northern Brook Lamprey" ~ 'Ichthyomyzon fossor',
            TRUE ~ scientific_name
      )) |> 
      distinct() |> 
      filter(scientific_name != 'Clupanodon thrissa')
glimpse(final_spp_list)

#### fit check -----------------------------------------------------------------
fi_tax0 <- fi |> 
      mutate(common_name = case_when(
            common_name == 'Warmouth Sunfish'    ~ 'Warmouth',
            common_name == 'North Brook Lamprey' ~ 'Northern Brook Lamprey',
            TRUE ~ common_name
      )) |> 
      left_join(final_spp_list) |> 
      mutate(common_name = case_when(
            common_name == 'Stonecat Madtom'      ~ 'Stonecat',
            common_name == 'South. Redbelly Dace' ~ 'Southern Redbelly Dace',
            TRUE ~ common_name
      )) |> 
      mutate(
            scientific_name = case_when(
                  common_name == 'Eastern Banded Killifish' ~ 'Fundulus diaphanus',
                  common_name == 'Stonecat'                 ~ 'Noturus flavus',
                  common_name == 'Southern Redbelly Dace'   ~ 'Chrosomus erythrogaster',
                  TRUE ~ scientific_name
            )
      ) |> 
      dplyr::select(-type)
skim(fi_tax0)

spp_list <- fi_tax0 |> select(common_name, scientific_name) |> distinct()
glimpse(spp_list)

#### pull remaining taxonomic classifications ----------------------------------
class_list <- classification(
      spp_list$scientific_name,
      db = "itis"
)

class_long <- class_list |>
      keep(is.data.frame) |>                                  
      imap_dfr(\(df, sci) tibble(scientific_name = sci,
                                 rank = df$rank,
                                 name = df$name,
                                 id   = df$id)) |> 
      filter(rank %in% c("kingdom", "phylum", "class", "order",
                         "family", "genus", "species")) |>
      select(-id) |>                                          
      pivot_wider(names_from = rank, values_from = name)

taxonomy_dirty <- fi_tax0 |> 
      select(common_name, scientific_name) |> 
      distinct() |> 
      left_join(class_long)
write_csv(taxonomy_dirty, 'localdata/taxonomy_dirty.csv')

taxonomy_clean <- read_csv('localdata/taxonomy_clean.csv')
skim(taxonomy_clean)

# Traits -----------------------------------------------------------------------
### load necessary data --------------------------------------------------------
tr <- read_csv('localdata/traits.csv') |> 
      janitor::clean_names() |> 
      select(
            scientific_name = sciname,
            max_body_length = max_body_l,
            length_mat      = leng_mat,
            longevity       = longevi,
            age_maturity    = age_mat,
            fecundity       = fecund,
            egg_size        = egg_size,
            parental_care   = parent_c
      )
skim(tr)

oepa_traits <- taxonomy_clean |> 
      left_join(tr) |> 
      distinct()
glimpse(oepa_traits)
skim(oepa_traits)

oepa_traits_dirty <- oepa_traits |>
      filter(if_any(max_body_length:parental_care, is.na))

fb_names <- validate_names(oepa_traits_dirty$scientific_name)

name_xwalk <- tibble(
      scientific_name = oepa_traits_dirty$scientific_name,
      fb_name = fb_names
)

name_xwalk |> filter(is.na(fb_name))
sp_to_query <- name_xwalk$fb_name |> na.omit() |> unique()

fb_species   <- species(sp_to_query)
fb_maturity  <- maturity(sp_to_query)
fb_fecund    <- fecundity(sp_to_query)
fb_repro     <- reproduction(sp_to_query)

sp_summary <- fb_species |>
      filter(LTypeMaxM == "TL" | is.na(LTypeMaxM)) |>   # only TL records
      select(Species,
             fb_max_body_length = Length,                # cm TL
             fb_longevity       = LongevityWild) |>
      mutate(fb_max_body_length = fb_max_body_length * 10)  # cm → mm

mat_summary <- fb_maturity |>
      filter(is.na(Type1) | Type1 == "TL") |>           # TL only
      group_by(Species) |>
      summarise(
            fb_length_mat   = median(LengthMatMin, na.rm = TRUE) * 10,  # cm → mm
            fb_age_maturity = median(AgeMatMin,    na.rm = TRUE),
            .groups = "drop"
      ) |>
      mutate(across(starts_with("fb_"), \(x) ifelse(is.nan(x), NA, x)))

fec_summary <- fb_fecund |>
      group_by(Species) |>
      summarise(
            fb_fecundity = (median(FecundityMax, na.rm = TRUE) + median(FecundityMin, na.rm = TRUE))/2,
            
            .groups = "drop"
      ) |>
      mutate(fb_fecundity = ifelse(is.nan(fb_fecundity), NA, fb_fecundity))

fb_combined <- name_xwalk |>
      left_join(sp_summary,  by = c("fb_name" = "Species")) |>
      left_join(mat_summary, by = c("fb_name" = "Species")) |>
      left_join(fec_summary, by = c("fb_name" = "Species"))
fb_combined

oepa_traits_filled <- oepa_traits |>
      left_join(fb_combined |> select(scientific_name, starts_with("fb_")),
                by = "scientific_name") |>
      mutate(
            max_body_length = coalesce(max_body_length, fb_max_body_length),
            length_mat      = coalesce(length_mat,      fb_length_mat),
            longevity       = coalesce(longevity,       fb_longevity),
            age_maturity    = coalesce(age_maturity,    fb_age_maturity),
            fecundity       = coalesce(fecundity,       fb_fecundity)
      ) |>
      select(-starts_with("fb_"))
glimpse(oepa_traits_filled)

##### visualize relationship among traits --------------------------------------
traits_only <- oepa_traits_filled |>
      select(max_body_length:parental_care)

cor_mat <- cor(traits_only, use = "pairwise.complete.obs", method = "spearman")

corrplot(cor_mat,
         method = "color",
         type = "upper",
         addCoef.col = "black", 
         tl.col = "black",
         tl.srt = 45,
         diag = FALSE)

ggpairs(traits_only,
        upper = list(continuous = wrap("cor", method = "spearman", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.5, size = 1)),
        diag  = list(continuous = wrap("densityDiag", alpha = 0.5))) +
      theme_bw()

##### use taxonomic imputation for remaining NAs -------------------------------
