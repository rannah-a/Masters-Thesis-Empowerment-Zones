################################################################################
# 03_build_city_covars.R —— Add 1990 socio-economic covariates (BGK) ----------
################################################################################
#  Inputs
#    • data_clean/1_intermediate/city_treatment_flags.rds  (Script 02 output)
#    • data_raw/BGK_raw/.../city_level_data.dta            (Busso–Gregory–Kline)
#
#  Output
#    • data_clean/1_intermediate/city_covars_1990.rds      (cross-section)
################################################################################

## 0. Setup -------------------------------------------------------------------
library(here); source(here("code","00_setup.R"))
library(tidyverse)
library(haven)       # read_dta()

## 1. Load city treatment file ------------------------------------------------
city_treat <- read_rds(
  here("data_clean","1_intermediate","city_treat_flags.rds")
)

## 2. Read & filter BGK city dataset -----------------------------------------
bgk_city_path <- here(
  "data_raw","BGK_raw","112607-V1","DataAppendix","Data",
  "PublicData","CityData","city_level_data.dta")

bgk_city_raw  <- read_dta(bgk_city_path)


bgk_1990 <- bgk_city_raw %>% 
  filter(year == 3) %>%                                 # 1990 in BGK coding
  mutate(
    statefp    = substr(tractid, 1, 2),
    placefp5   = str_pad(city, 5, pad = "0"),
    place_id   = paste0(statefp, placefp5),             # 7-digit state+place
    place_name = str_squish(str_to_title(true_city))
  ) %>% 
  group_by(place_id, place_name) %>% 
  # Collapse only numeric columns within each (state, place) using max (or any aggregator you prefer)
  summarise(
    across(
      .cols = where(is.numeric),
      .fns  = ~ suppressWarnings(max(.x, na.rm = TRUE))
    ),
    .groups = "drop"
  ) %>% 
  # Turn any ±Inf (from all-missing groups) into NA
  mutate(
    across(where(is.numeric), ~ ifelse(is.infinite(.x), NA_real_, .x))
  ) %>% 
  # Build the clean set of covariates you want to carry forward
  transmute(
    place_id, place_name,
    bgk_pop1990      = population_city,
    child_share_90   = child_city,              # already a fraction
    black_city_90    = black_city,
    latino_city_90   = latino_city,
    white_city_90    = white_city,
    poverty_rate_90  = poverty_city,            # %
    unemploy_90      = unemployment_city,       # %
    ln_wage_90       = ifelse(wage_city > 0, log(wage_city), NA_real_),
    welfare_90       = welfare_city,            # per 1,000 pop
    crime_90         = crimerate_city,          # per 1,000
    gov_emp_share_90 = pctgov_city              # %
  ) %>% 
  filter(!is.na(place_name) & place_name != "")

## 3. Merge by place_fips -----------------------------------------------------
city_covars90 <- city_treat %>% 
  left_join(bgk_1990, by = c("place_fips" = "place_id"))

## 4. Quick diagnostics---------------------------------------------------------
city_covars90 %>% 
  summarise(
    missing_covars = sum(is.na(child_share_90)),
    treated_cities = sum(treated_dummy),
    nearmiss_cities= sum(nearmiss_dummy),
    future_cities  = sum(future_dummy)
  ) %>% print()

#Match rate with BGK covariates
diag_match <- city_covars90 %>%
  summarise(
    n_cities          = n(),
    matched_bgk       = sum(!is.na(bgk_pop1990)),
    match_rate        = matched_bgk / n_cities,
    treated_cities    = sum(treated_dummy),
    nearmiss_cities   = sum(nearmiss_dummy),
    future_cities     = sum(future_dummy)
  )
print(diag_match)

# Compare our 1990 pop (from NHGIS) vs BGK’s
pop_compare <- city_covars90 %>%
  filter(!is.na(bgk_pop1990)) %>%
  summarise(
    corr_pop       = cor(pop90_city, bgk_pop1990, use = "complete.obs"),
    med_abs_diff   = median(abs(pop90_city - bgk_pop1990), na.rm = TRUE),
    p90_abs_diff   = quantile(abs(pop90_city - bgk_pop1990), 0.90, na.rm = TRUE)
  )
print(pop_compare)


## 5. Save --------------------------------------------------------------------
write_rds(
  city_covars90,
  here("data_clean","1_intermediate","city_treat_covars_1990.rds")
)


