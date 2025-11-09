################################################################################
# # 02_build_city_ez_flags.R ── Enrich 1990‐Tract Flags with Population & Spatial Dimension --
################################################################################
# Purpose:
#   PART 1: Pull, clean & merge NHGIS population (incl. under-18) to 1990 tracts
#   PART 2: Map 1990 tracts → 1990 places (centroid-within), aggregate to cities
#
# Inputs:
#   data_clean/1_intermediate/tracts90_flags.rds
#   data_raw/nhgis0004_csv/nhgis0004_ds120_1990_tract.csv          (NHGIS counts)
#   data_raw/nhgis0005_shape-nhgis0005_shapefile_tl2000_us_place_1990/US_place_1990.shp (1990 Place Polygons)
#
# Outputs:
#   data_clean/1_intermediate/tracts90_w_pop.rds
#   data_clean/1_intermediate/city_treatment_flags.rds
################################################################################

## 0. Total Script Setup: Libraries & Paths ------------------------------------------------
source(here::here("code", "00_setup.R")) #Load packages and anchor project root

library(readr)
library(dplyr)
library(tidyr)
library(haven)
library(stringr)
library(purrr)
library(sf)
library(tigris)
options(tigris_use_cache = TRUE)

# Paths: 

#a. tract flags from Script 01:
tracts90_flags_path <- here("data_clean", "1_intermediate", "tracts90_flags.rds")

#b. raw NHGIS extract -1990 population data:
nhgis_raw_path      <- here("data_raw", "nhgis0004_csv", "nhgis0004_ds120_1990_tract.csv")

#c. raw NHGIS extract -1990 place shape data:
place90_raw_path      <- here("data_raw", "nhgis0005_shape-nhgis0005_shapefile_tl2000_us_place_1990", "US_place_1990.shp")

# ------------------------------------------------------------------------------
### PART 1: PULL, CLEAN, & MERGE NHGIS TRACT-LEVEL POPULATION DATA TO TRACT90
# ------------------------------------------------------------------------------

## 1. Load Existing 1990-Tract Flags -------------------------------------------
  #tracts90 has: tract1990, ez1, app, future, treat_yr (and placeholder pop fields)
tracts90 <- read_rds(tracts90_flags_path) #from script 1

## 2. Load & Clean NHGIS 1990 Tract-Level Population Data ----------------------

nhgis_1990_raw <- read_csv(nhgis_raw_path, show_col_types = FALSE)

# keep only needed fields and rename:
nhgis_1990 <- nhgis_1990_raw %>%
  transmute(
    STATEA    = STATEA,  #2-digit state
    COUNTYA   = COUNTYA, #3-digit county
    TRACTA    = TRACTA,  #4-6 digits tract (4 digits unless tract has two digit suffix) 
    
    # totals & age buckets:
    pop_total_1990 = as.integer(ET1001),   #total population in the tract
    age_under1   = as.integer(ET3001),
    age_1_2      = as.integer(ET3002),
    age_3_4      = as.integer(ET3003),
    age_5        = as.integer(ET3004),
    age_6        = as.integer(ET3005),
    age_7_9      = as.integer(ET3006),
    age_10_11    = as.integer(ET3007),
    age_12_13    = as.integer(ET3008),
    age_14       = as.integer(ET3009),
    age_15       = as.integer(ET3010),
    age_16       = as.integer(ET3011),
    age_17       = as.integer(ET3012)
  ) %>%
  
  mutate(
    #TRACTA is 4 digits when no suffix; pad to 6 by appending '00'
    TRACTA_2 = if_else(
      str_length(TRACTA) == 4,
      str_c(str_pad(TRACTA, width = 4, side = "left", pad = "0"), "00"),
      TRACTA
    ),
    tract1990 = str_c(STATEA, COUNTYA, TRACTA_2),
    # compute under-16:
    pop_under16_1990 = age_under1 + age_1_2 + age_3_4 + age_5 +
                       age_6 + age_7_9 + age_10_11 + age_12_13 +
                       age_14 + age_15,
    # compute under-18:
    pop_under18_1990 = pop_under16_1990 + age_16 + age_17
  ) %>%
    # keep only new columns
  select(
    tract1990,
    pop_total_1990,
    pop_under16_1990,
    pop_under18_1990
  )

# check: how many flagged tracts missing in NHGIS?
cat("Flagged tracts missing from NHGIS:",
    tracts90 %>% anti_join(nhgis_1990, by = "tract1990") %>% nrow(), "\n")

## 3. Merge NHGIS Covariates into 1990-Tract Flags -----------------------------
tracts90_w_pop <- nhgis_1990 %>%
  left_join(
    tracts90 %>% select(tract1990, ez1, app, future, treat_yr),
    by = "tract1990"
  ) %>%
  mutate(
    pop1990        = pop_total_1990,
    u16_pop_1990   = pop_under16_1990,
    u18_pop_1990   = pop_under18_1990,
    u18_share_1990 = if_else(pop1990 > 0, u18_pop_1990 / pop1990, NA_real_),
    u16_share_1990 = if_else(pop1990 > 0, u16_pop_1990 / pop1990, NA_real_),
    across(c(ez1, app, future), ~replace_na(.x, FALSE))
  ) %>%
  select(tract1990, pop1990, u16_pop_1990, u18_pop_1990,
         u16_share_1990, u18_share_1990, ez1, app, future, treat_yr)



# Keep all flagged tracts; otherwise require positive population
tracts90_w_pop <- tracts90_w_pop %>%
  filter(
    (ez1 == TRUE) | (app == TRUE) | (future == TRUE) |
      (!is.na(pop1990) & pop1990 > 0 )
    )

#build state_fips 
state_fips <- tracts90_w_pop$tract1990 %>% substr(1, 2) %>% unique()
valid_fips <- unique(tigris::fips_codes$state_code)
state_fips <- intersect(state_fips, valid_fips)


# Sanity:counts retained
tracts90_w_pop %>%
  summarise(
    keep_ez1    = sum(ez1),
    keep_app    = sum(app),
    keep_future = sum(future),
    total_rows  = n()
  ) %>% 
  print()

# Save 
write_rds(tracts90_w_pop,
          here("data_clean/1_intermediate/tracts90_w_pop.rds"))

# ------------------------------------------------------------------------------
### PART 2: TRANSLATE 1990 TRACTS → 1990 PLACE FIPS (centroid-within)
# ------------------------------------------------------------------------------
## 4. Pull 1990 place polygons (NHGIS TL2000 vintage) --------------------------
place90_raw <- st_read(place90_raw_path, quiet = TRUE)

places90 <- place90_raw %>% 
  filter(ENTITY == "P", STATE %in% state_fips) %>%  # keep ALL place classes
  mutate(
    place_fips = paste0(STATE, stringr::str_pad(FPL90, 5, pad = "0")),         # 7-digit 1990 state+place code
    place_name = stringr::str_squish(stringr::str_to_title(NAME))
  ) %>% 
  group_by(place_fips, place_name) %>%
  summarise(
    geometry = sf::st_union(geometry),
    .groups = "drop"
  ) %>%
  sf::st_make_valid() # align CRS after we build tracts90_sf in Step 5.

## 5. Build 1990 tract polygons (cb=TRUE; we only need to centroid them) -------

tracts90_sf <- purrr::map_dfr(
  state_fips,
  \(st) {
    sf_tmp <- tigris::tracts(state = st, year = 1990, cb = TRUE, progress_bar = FALSE)
    sf_tmp %>%
      mutate(
        tract1990 = paste0(STATEFP, COUNTYFP, TRACTBASE, TRACTSUF)
      ) %>%
      select(tract1990, geometry)
  }
)

# Drop water-only “.99” tracts in BOTH tables
tracts90_sf   <- tracts90_sf   %>% filter(!str_detect(tract1990, "99$"))
tracts90_w_pop <- tracts90_w_pop %>% filter(!str_detect(tract1990, "99$"))

## 6. Merge flags & population to 1990 tract polygons --------------------------
tracts90_sf <- tracts90_sf %>%
  left_join(tracts90_w_pop, by = "tract1990")

# Normalize NAs for non-flagged tracts (prevents NA cascades)
tracts90_sf <- tracts90_sf %>%
  mutate(
    across(c(ez1, app, future), ~replace_na(.x, FALSE)),
    pop1990       = replace_na(pop1990, 0L),
    u16_pop_1990  = replace_na(u16_pop_1990, 0L),
    u18_pop_1990  = replace_na(u18_pop_1990, 0L),
    u16_share_1990 = if_else(pop1990 > 0, u16_pop_1990 / pop1990, 0),
    u18_share_1990 = if_else(pop1990 > 0, u18_pop_1990 / pop1990, 0)
  )

## 7. Ensure valid geometries & common CRS -------------------------------------
tracts90_sf <- tracts90_sf %>% sf::st_make_valid()
places90    <- places90    %>% sf::st_transform(sf::st_crs(tracts90_sf))

## 8. Assign each tract to one city via centroid-within (1990 places) ----------
# Use centroid of largest polygon to handle multipart tracts safely
tract_cent <- tracts90_sf %>% 
  sf::st_make_valid() %>%
  sf::st_centroid(of_largest_polygon = TRUE)

tract_city_main <- sf::st_join(
  tract_cent,
  places90,
  join = sf::st_within,
  left = TRUE
) %>%
  sf::st_drop_geometry() %>%
  select(tract1990, place_fips, place_name) %>%
  distinct(tract1990, place_fips, .keep_all = TRUE)

# Diagnostics
message("Tracts with no 1990-place match: ",
        sum(is.na(tract_city_main$place_fips)))

flagged_missing_place <- tracts90_w_pop %>%
  filter(ez1 | app | future) %>%
  anti_join(tract_city_main, by = "tract1990") %>%
  nrow()
message("Flagged tracts without city assignment: ", flagged_missing_place)


## 9. Merge city codes back into tract table -----------------------------------
tracts90_geo <- tracts90_w_pop %>%
  left_join(tract_city_main, by = "tract1990")%>%
  filter(pop1990 > 0 | ez1 | app | future) %>%
  filter(!is.na(place_fips))     # only city-assigned rows

# Save 
write_rds(tracts90_geo,
          here("data_clean/1_intermediate/tracts90_geo.rds"))

## 10. Aggregate to City Level — baseline (centroid→1990 places; boolean flags)
city_treat <- tracts90_geo %>% 
  group_by(place_fips, place_name) %>% 
  summarise(
    # city totals from tract-level vectors
    pop90_city   = sum(pop1990,      na.rm = TRUE),
    u18_90_city  = sum(u18_pop_1990, na.rm = TRUE),
    vap90_city  = pop90_city - u18_90_city,
    
    # population-weighted shares: use tract-level pop in the numerators
    share_ez     = if_else(pop90_city > 0,
                            sum(pop1990[ez1], na.rm = TRUE) / pop90_city, NA_real_),
    share_app    = if_else(pop90_city > 0,
                            sum(pop1990[app], na.rm = TRUE) / pop90_city, NA_real_),
    share_future = if_else(pop90_city > 0,
                            sum(pop1990[future], na.rm = TRUE) / pop90_city,  NA_real_),
    
    # earliest treatment year among city’s tracts
    treat_yr     = suppressWarnings(min(treat_yr, na.rm = TRUE)),
    .groups = "drop"
  ) %>% 
  mutate(
    treat_yr = ifelse(is.infinite(treat_yr), NA_real_, treat_yr),
    # dummies (BGK-style precedence: treated > app > future)
    treated_dummy  = share_ez >= 0.10,
    nearmiss_dummy = !treated_dummy & (share_app >= 0.10),
    future_dummy   = !treated_dummy & !nearmiss_dummy & (share_future >= 0.10)
  )


# Quick roll-up
city_treat %>%
  summarise(
    n_cities  = n(),
    treated   = sum(treated_dummy),
    nearmiss  = sum(nearmiss_dummy),
    future    = sum(future_dummy),
    med_pop90 = median(pop90_city),
    med_vap90 = median(vap90_city)
  ) %>% print()

## 11. Save city-level treatment file ------------------------------------------
write_rds(city_treat,
          here("data_clean/1_intermediate/city_treat_flags.rds"))

