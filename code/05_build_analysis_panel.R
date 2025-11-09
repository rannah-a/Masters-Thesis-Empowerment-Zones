################################################################################
# 05_build_analysis_panel.R ── Build Final Analysis Panel 
###############################################################################
# Inputs:
#   data_clean/1_intermediate/city_panel.rds
# Outputs:
#   data_clean/1_intermediate/city_panel_analysis_mayor.rds
#   data_clean/1_intermediate/city_panel_analysis_council.rds
################################################################################

## 0. Setup---------------------------------------------------------------------
source(here::here("code","00_setup.R"))
library(dplyr)
library(stringr)
library(readr)
library(fixest)
library(lubridate)
library(modelsummary)

## 1. Load combined city panel from script 4-----------------------------------
city_panel <- read_rds(here("data_clean","1_intermediate","city_panel.rds"))


## 2. Keep only Mayor / City Council; build month-aware POST -------------------
city_panel <- city_panel %>%
  filter(office_consolidated %in% c("Mayor","City Council")) %>%
  mutate(
    # Adoption is December 1994 (EZ). Use month-aware post for Mayors at treat_yr.
    post = dplyr::case_when(
      # Mayors: post if after treat_yr OR (in treat_yr AND month_final >= 12)
      treated_dummy == 1L & office_consolidated == "Mayor" ~
        as.integer(year > treat_yr | (year == treat_yr & coalesce(as.integer(month_final), 0L) >= 12L)),
      # Councils: conservative (almost no December city elections)
      treated_dummy == 1L & office_consolidated == "City Council" ~
        as.integer(year > treat_yr),
      TRUE ~ 0L
    )
  )

## 3. Event time & cohort tags (no window trimming here) -----------------------
panel_core <- city_panel %>%
  mutate(
    rel_year = ifelse(!is.na(treat_yr), year - treat_yr, NA_real_),
    cohort   = ifelse(treated_dummy == 1L, treat_yr, NA_real_),
    group_t  = treated_dummy,
    group_nm = nearmiss_dummy,
    group_fu = future_dummy
  )

## 4. Measurement hygiene (baseline rules) -------------------------------------
  # - Mayor turnout: cap >120% => NA (keep row for FE)
  # - Council per-seat turnout (votes_per_seat_pct_vap): cap >200 => NA
  # - Keep controls only if city VAP >= 1,000 (treated/near-miss/future kept regardless)
panel_core <- panel_core %>%
  mutate(
    turnout_pct = ifelse(
      office_consolidated=="Mayor" & !is.na(turnout_pct) & turnout_pct > 120,
      NA_real_, turnout_pct
    ),
    votes_per_seat_pct_vap = ifelse(
      office_consolidated=="City Council" & !is.na(votes_per_seat_pct_vap) & votes_per_seat_pct_vap > 200,
      NA_real_, votes_per_seat_pct_vap
    )
  ) %>%
  mutate(vap_ok = vap90_city >= 1000) %>%
  filter(vap_ok | treated_dummy | nearmiss_dummy | future_dummy) %>%
  select(-vap_ok)

## 5. Coverage rules for FE (don’t trim event window here) ---------------------
  # - Treated cities: need ≥1 pre & ≥1 post observation
  # - Controls: require ≥2 total observations so FE estimable
coverage_by_city <- panel_core %>%
  group_by(place_fips, office_consolidated) %>%
  summarise(
    n_pre   = sum(!is.na(rel_year) & rel_year < 0,  na.rm = TRUE),
    n_post  = sum(!is.na(rel_year) & rel_year >= 0, na.rm = TRUE),
    n_total = dplyr::n(),
    .groups = "drop"
  )

keep_rules <- coverage_by_city %>%
  mutate(
    keep_treated = (n_pre >= 1 & n_post >= 1),
    keep_control = (n_total >= 2)
  )

panel_keep <- panel_core %>%
  left_join(keep_rules, by = c("place_fips","office_consolidated")) %>%
  filter(ifelse(treated_dummy == 1L, keep_treated, keep_control)) %>%
  select(-keep_treated, -keep_control)

  


## 6. Save split analysis panels ----------------------------------------------------
panel_mayor   <- panel_keep %>% filter(office_consolidated == "Mayor")
panel_council <- panel_keep %>% filter(office_consolidated == "City Council")

write_rds(panel_mayor,   here("data_clean","1_intermediate", "city_panel_analysis_mayor.rds"))
write_rds(panel_council, here("data_clean","1_intermediate", "city_panel_analysis_council.rds"))


## 7. Diagnostics
cat("\n== Dimensions (unique cities / years / rows) ==\n")
panel_mayor   %>% summarise(n_cities=n_distinct(place_fips), n_years=n_distinct(year), n_rows=n())   %>% print()
panel_council %>% summarise(n_cities=n_distinct(place_fips), n_years=n_distinct(year), n_rows=n())   %>% print()

cat("\n== Key uniqueness (should be zero rows) ==\n")
panel_mayor   %>% count(place_fips, year)   %>% filter(n>1) %>% print()
panel_council %>% count(place_fips, year)   %>% filter(n>1) %>% print()

cat("\n== Join health: 1990 covariates present (should be 0) ==\n")
panel_mayor   %>% summarise(na_covars = sum(is.na(pop90_city)))   %>% print()
panel_council %>% summarise(na_covars = sum(is.na(pop90_city)))   %>% print()

cat("\n== Group sizes (treated / near-miss / future / untreated rows) ==\n")
panel_mayor %>% summarise(
  treated   = sum(treated_dummy,  na.rm=TRUE),
  nearmiss  = sum(nearmiss_dummy, na.rm=TRUE),
  future    = sum(future_dummy,   na.rm=TRUE),
  untreated = sum(!treated_dummy & !nearmiss_dummy & !future_dummy)
) %>% print()

panel_council %>% summarise(
  treated   = sum(treated_dummy,  na.rm=TRUE),
  nearmiss  = sum(nearmiss_dummy, na.rm=TRUE),
  future    = sum(future_dummy,   na.rm=TRUE),
  untreated = sum(!treated_dummy & !nearmiss_dummy & !future_dummy)
) %>% print()

cat("\n== Treated support in event time (rows) ==\n")
rel_support <- function(df){
  df %>%
    filter(treated_dummy == 1L, !is.na(treat_yr), !is.na(year)) %>%
    mutate(rel_year = year - treat_yr) %>%
    summarise(
      min_rel = min(rel_year, na.rm=TRUE),
      p10     = quantile(rel_year, .10, na.rm=TRUE),
      p50     = quantile(rel_year, .50, na.rm=TRUE),
      p90     = quantile(rel_year, .90, na.rm=TRUE),
      max_rel = max(rel_year, na.rm=TRUE),
      n_pre   = sum(rel_year < 0, na.rm=TRUE),
      n_post  = sum(rel_year >= 0, na.rm=TRUE)
    )
}
rel_support(panel_mayor)   %>% print()
rel_support(panel_council) %>% print()

cat("\n== Share of treated cities with both pre & post ==\n")
prepost <- function(df){
  df %>%
    filter(treated_dummy == 1L, !is.na(treat_yr)) %>%
    mutate(rel_year = year - treat_yr) %>%
    group_by(place_fips) %>%
    summarise(has_pre = any(rel_year < 0), has_post = any(rel_year >= 0), .groups="drop") %>%
    summarise(n_treated = n(), both_pre_post = sum(has_pre & has_post), share_both = both_pre_post/n_treated)
}
prepost(panel_mayor)   %>% print()
prepost(panel_council) %>% print()

cat("\n== Outcome coverage & scales ==\n")
diag_outcome <- function(df, y){
  stopifnot(y %in% names(df))
  df %>%
    summarise(
      nonmiss_rows   = sum(!is.na(.data[[y]])),
      share_nonmiss  = mean(!is.na(.data[[y]])),
      p10            = quantile(.data[[y]], .1, na.rm=TRUE),
      p50            = median(.data[[y]], na.rm=TRUE),
      p90            = quantile(.data[[y]], .9, na.rm=TRUE),
      min_val        = suppressWarnings(min(.data[[y]], na.rm=TRUE)),
      max_val        = suppressWarnings(max(.data[[y]], na.rm=TRUE))
    )
}

cat("\n-- Mayor: dem_share / turnout_pct / inc_win --\n")
diag_outcome(panel_mayor, "dem_share")     %>% print()
diag_outcome(panel_mayor, "turnout_pct")   %>% print()
diag_outcome(panel_mayor, "inc_win")       %>% print()

cat("\n-- Council: dem_share / votes_per_seat_pct_vap / inc_win_share --\n")
diag_outcome(panel_council, "dem_share")                 %>% print()
diag_outcome(panel_council, "votes_per_seat_pct_vap")    %>% print()
diag_outcome(panel_council, "inc_win_share")             %>% print()

cat("\n== Range checks for shares in [0,1] (should be empty) ==\n")
bad_mayor <- panel_mayor %>%
  filter(dem_share < 0 | dem_share > 1 |
           female_cand_share < 0 | female_cand_share > 1 |
           female_win_share  < 0 | female_win_share  > 1)
bad_council <- panel_council %>%
  filter(dem_share < 0 | dem_share > 1 |
           inc_win_share < 0 | inc_win_share > 1 |
           female_cand_share < 0 | female_cand_share > 1 |
           female_win_share  < 0 | female_win_share  > 1)
cat("mayor bad rows:", nrow(bad_mayor), "\n")
cat("council bad rows:", nrow(bad_council), "\n")

cat("\n== 1994 treated MAYOR rows: month-aware POST sanity ==\n")
panel_mayor %>%
  filter(treated_dummy == 1L, year == treat_yr) %>%
  count(month_final, post) %>%
  arrange(month_final) %>%
  print(n = Inf)

cat("\n== Covariate standardized differences (treated vs all controls) ==\n")
std_diff <- function(df, var){
  with(df, {
    m1 <- mean(get(var)[treated_dummy == 1L], na.rm=TRUE)
    m0 <- mean(get(var)[treated_dummy != 1L],  na.rm=TRUE)
    s  <- sqrt((var(get(var)[treated_dummy == 1L], na.rm=TRUE) +
                  var(get(var)[treated_dummy != 1L],  na.rm=TRUE))/2)
    (m1 - m0) / s
  })
}
covs <- c("ln_wage_90","poverty_rate_90","unemploy_90",
          "black_city_90","latino_city_90","white_city_90",
          "child_share_90","crime_90","gov_emp_share_90")
tibble::tibble(
  cov = covs,
  sdiff_mayor   = sapply(covs, \(v) std_diff(panel_mayor,   v)),
  sdiff_council = sapply(covs, \(v) std_diff(panel_council, v))
) %>% print(n = Inf)

cat("\nDone. Estimation starts in Script 6 (incl. quick unbalanced baseline + balanced EB models).\n")

###==================================================================================================
### PART TWO:5B) Patch policy arms (EZ1 / SEZ1) and rebuild analysis panels
#      — no upstream edits required
# ================================

library(dplyr); library(stringr)

# Helper: find a city's FIPS safely by name + state
.fips_for <- function(df, name_pat, state) {
  cand <- df %>%
    mutate(name_lc = str_to_lower(place_name_treat)) %>%
    filter(str_detect(name_lc, name_pat), state_abb == state) %>%
    distinct(place_fips, place_name_treat, state_abb)
  if (nrow(cand) == 0) {
    warning(sprintf("No match for %s, %s", name_pat, state))
    return(character(0))
  }
  if (nrow(cand) > 1) {
    message("Multiple matches for ", name_pat, " (", state, "); taking first:\n",
            capture.output(print(cand)))
  }
  cand$place_fips[1]
}

# Work off the raw city_panel (already loaded earlier in Script 5)
fips_xwalk <- city_panel %>% distinct(place_fips, place_name_treat, state_abb)

# Canonical cities (BGK Round I EZ-1) + SEZ-1
fips_atl <- .fips_for(fips_xwalk, "atlanta",       "GA")
fips_bal <- .fips_for(fips_xwalk, "baltimore",     "MD")
fips_chi <- .fips_for(fips_xwalk, "chicago",       "IL")
fips_det <- .fips_for(fips_xwalk, "detroit",       "MI")
fips_nyc <- .fips_for(fips_xwalk, "new york",      "NY")
fips_phi <- .fips_for(fips_xwalk, "philadelphia",  "PA")
fips_cam <- .fips_for(fips_xwalk, "camden",        "NJ")

fips_la  <- .fips_for(fips_xwalk, "los angeles",   "CA")  # SEZ-1
fips_cle <- .fips_for(fips_xwalk, "cleveland",     "OH")  # SEZ-1

ez1_set  <- c(fips_atl, fips_bal, fips_chi, fips_det, fips_nyc, fips_phi, fips_cam) %>% unique() %>% setdiff(NA_character_)
sez1_set <- c(fips_la, fips_cle) %>% unique() %>% setdiff(NA_character_)

# --- Patch baseline treatment only for EZ1; store SEZ1 separately ---

patch_policy <- function(df) {
  df %>%
    mutate(
      policy_arm = case_when(
        place_fips %in% ez1_set  ~ "EZ1",
        place_fips %in% sez1_set ~ "SEZ1",
        TRUE                     ~ "None"
      ),
      # Baseline treated set = EZ1 only
      treated_dummy = if_else(policy_arm == "EZ1", 1L, 0L),
      # Treat year only for EZ1 in the baseline panel
      treat_yr      = if_else(policy_arm == "EZ1", 1994, NA_real_),
      
      # Keep SEZ-1 timing for later robustness, but don't affect baseline filters
      sez1_year     = if_else(policy_arm == "SEZ1", 1999, NA_real_),
      
      # Recompute post/rel_year for EZ1 baseline
      post          = if_else(policy_arm == "EZ1" & !is.na(year), year >= 1994, FALSE),
      rel_year      = if_else(policy_arm == "EZ1", year - 1994, NA_real_),
      
      # Any city newly treated should not be flagged as "future" any more
      future_dummy  = if_else(policy_arm == "EZ1", 0L, future_dummy)
    )
}

# Recreate the "panel_core" (measurement trims identical to your Step 2)
panel_core2 <- city_panel %>%
  # keep only Mayor / Council, rebuild post based on (potentially) updated treat_yr in patch step
  filter(office_consolidated %in% c("Mayor", "City Council")) %>%
  patch_policy() %>%
  mutate(
    # Conservative measurement trims (same rules you had)
    turnout_pct = if_else(
      office_consolidated == "Mayor" & turnout_pct > 120, NA_real_, turnout_pct
    ),
    votes_per_seat_pct_vap = if_else(
      office_consolidated == "City Council" & votes_per_seat_pct_vap > 200,
      NA_real_, votes_per_seat_pct_vap
    ),
    vap_ok = vap90_city >= 1000
  ) %>%
  # Don’t drop SEZ-1 rows here — only drop on VAP for tiny cities unless in treated/near-miss/future
  filter(vap_ok | treated_dummy == 1L | nearmiss_dummy == 1L | future_dummy == 1L) %>%
  select(-vap_ok)

# Coverage window (same as your Step 3; applied to treated only)
panel_window2 <- panel_core2 %>%
  filter(is.na(rel_year) | dplyr::between(rel_year, -10, 10))

# Keep treated with ≥1 pre AND ≥1 post; keep all controls
coverage_by_city2 <- panel_core2 %>%
  group_by(place_fips, office_consolidated) %>%
  summarise(
    n_pre  = sum(!is.na(rel_year) & rel_year < 0,  na.rm = TRUE),
    n_post = sum(!is.na(rel_year) & rel_year >= 0, na.rm = TRUE),
    .groups = "drop"
  )

keep_rules2 <- coverage_by_city2 %>%
  mutate(keep_treated = (n_pre >= 1 & n_post >= 1),
         keep_control = TRUE)

panel_keep2 <- panel_window2 %>%
  left_join(keep_rules2, by = c("place_fips","office_consolidated")) %>%
  filter(if_else(treated_dummy == 1L, keep_treated, keep_control)) %>%
  select(-keep_treated, -keep_control)

# Split & overwrite analysis panels
panel_mayor   <- panel_keep2 %>% filter(office_consolidated == "Mayor")
panel_council <- panel_keep2 %>% filter(office_consolidated == "City Council")

write_rds(panel_mayor,   here::here("data_clean","1_intermediate","city_panel_analysis_mayor.rds"))
write_rds(panel_council, here::here("data_clean","1_intermediate","city_panel_analysis_council.rds"))

# --- Checks: did we patch as intended? ---

cat("\n== Policy arm counts ==\n")
panel_keep2 %>% count(office_consolidated, policy_arm) %>% print(n = Inf)

cat("\n== EZ1 treated city counts (should include Baltimore, Chicago, NYC, Philly, plus ATL, DET, Camden) ==\n")
panel_keep2 %>%
  filter(policy_arm == "EZ1") %>%
  summarise(
    treated_cities = n_distinct(place_fips),
    treated_rows   = sum(treated_dummy == 1L),
    min_year       = min(year, na.rm = TRUE),
    max_year       = max(year, na.rm = TRUE)
  ) %>% print()

cat("\n== Newly treated (by name) — pre/post support ==\n")
panel_keep2 %>%
  filter(policy_arm == "EZ1",
         place_fips %in% c(fips_bal, fips_chi, fips_nyc, fips_phi)) %>%
  mutate(rel_year = year - 1994) %>%
  group_by(place_fips, place_name_treat, office_consolidated) %>%
  summarise(n_pre = sum(rel_year < 0, na.rm=TRUE),
            n_post = sum(rel_year >= 0, na.rm=TRUE),
            years = paste(sort(unique(year)), collapse = ","),
            .groups = "drop") %>%
  arrange(place_name_treat, office_consolidated) %>% print(n = Inf)

cat("\n== SEZ1 cities should be present but NOT treated in baseline ==\n")
panel_keep2 %>%
  filter(place_fips %in% sez1_set) %>%
  distinct(place_fips, place_name_treat, office_consolidated, treated_dummy, sez1_year) %>%
  arrange(place_name_treat, office_consolidated) %>% print(n = Inf)

cat("\n== Final treated city counts per panel ==\n")
panel_mayor   %>% summarise(n_cities = n_distinct(place_fips[treated_dummy == 1L]),
                            treated_rows = sum(treated_dummy == 1L)) %>% print()
panel_council %>% summarise(n_cities = n_distinct(place_fips[treated_dummy == 1L]),
                            treated_rows = sum(treated_dummy == 1L)) %>% print()

cat("\nPatched policy arms applied. Proceed to Script 6 (unbalanced baseline + EB-balanced models).\n")


##2. Fix pre treat NYC city council 
# Target: New York City (place_fips 3651000), City Council, treated, 1994
nyc_cc_1994 <- with(panel_council,
                    place_fips == "3651000" &
                      office_consolidated == "City Council" &
                      treated_dummy == 1L &
                      year == 1994)

# Sanity before
panel_council %>%
  filter(nyc_cc_1994) %>%
  select(place_name_treat, office_consolidated, year, treated_dummy, post) %>%
  print(n = Inf)

# Flip that row to PRE (post = FALSE)
panel_council <- panel_council %>%
  mutate(post = if_else(nyc_cc_1994, FALSE, post))

# If the panel happens to carry a rel_year column already, set it to -1 for that row
if ("rel_year" %in% names(panel_council)) {
  panel_council <- panel_council %>%
    mutate(rel_year = if_else(nyc_cc_1994, -1, rel_year))
}

# Sanity after
panel_council %>%
  filter(nyc_cc_1994) %>%
  select(place_name_treat, office_consolidated, year, treated_dummy, post, dplyr::any_of("rel_year")) %>%
  print(n = Inf)

# (Re)save if you want this persisted
write_rds(panel_council, here::here("data_clean","1_intermediate","city_panel_analysis_council.rds"))

