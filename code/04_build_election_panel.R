################################################################################
# 04_build_election_panel.R —— Assemble City-level Election Outcomes  ----------
################################################################################
#  Purpose
#   • Bring in city EZ-treatment flags (Script 02 output + Script 03 covariates)
#   • Build *separate* panels for Mayor and City Council
#   • Fix incumbent-win logic
#   • Avoid double counting (keep the contest month with the most votes each year)
#   • Save mayor_panel.rds, council_panel.rds, and (optionally) city_panel_long.rds
#
#  Inputs
#    data_clean/1_intermediate/city_treat_covars_1990.rds   (from Script 03)
#    data_raw/ALEGD/Replication/ledb_candidatelevel.rds (ALEGD)                  # raw election files
#    
#
#  Outputs
#   data_clean/2_intermediate/mayor_panel.rds
#   data_clean/2_intermediate/council_panel.rds
#   data_clean/2_intermediate/city_panel_long.rds   (optional convenience file)
################################################################################

## 0.  Set-up ------------------------------------------------------------------
library(here)        # project-root paths
library(tidyverse)   # dplyr, tidyr, readr, stringr, …
library(janitor)     # clean_names(), tabyl(), …
library(sf)          # only if any spatial joins remain
library(arrow)       # if your ALEGD files are in Parquet/Feather

source(here("code", "00_setup.R"))     # ensures packages & here::i_am()


##  1. Paths--------------------------------------------------------------------
covars90_path  <- here("data_clean", "1_intermediate",
                       "city_treat_covars_1990.rds")   # flags + 1990 covars

ledb_raw_path  <- here("data_raw", "ALEGD_raw",
                       "Replication", "ledb_candidatelevel.rds")

##  2. Load city flags + 1990 covariates (script 3 output) ---------------------

city_covars90 <- read_rds(covars90_path)

# Recommended: rename columns 
city_covars90 <- city_covars90 %>%
  dplyr::rename(place_name_treat = place_name.x) %>%
  dplyr::rename(place_name_bgk = place_name.y)

#Keep  only needed columns 
base_cities <- city_covars90 %>%
  select(place_fips, place_name_treat, pop90_city, u18_90_city, vap90_city, share_ez, share_app, share_future,
         treat_yr, treated_dummy, nearmiss_dummy, future_dummy,
         # BGK covariates:
         place_name_bgk, bgk_pop1990, child_share_90, black_city_90, latino_city_90, white_city_90,
         poverty_rate_90, unemploy_90, ln_wage_90, welfare_90, crime_90, gov_emp_share_90) 


## 3. Load ALEGD candidate-leve & subset to only Mayor / Council----------------
ledb_cand <- read_rds(here(ledb_raw_path)) %>%
  select(fips, geo_name, state_abb, office_consolidated,
         year, month, district, votes, n_winners, winner,
         incumbent, pid_est, gender_est, race_est, ledb_candid) %>% 
  filter(office_consolidated %in% c("Mayor","City Council")) %>%
  mutate(
    month = as.integer(month),
    district = tidyr::replace_na(district, "_ALL")
  )

# Helpers: 
is_female <- function(x) x == "F"
is_dem    <- function(x) x %in% c("D")
is_poc    <- function(x) x %in% c("black","hispanic","asian","other")

## 4. Contest-level summaries (city × year × month × office) -------------------

# =========================== MAYOR PANEL ======================================
  # Choose the "final" month per city-year:
  # rule = month with max total votes; tie -> latest month wins.
  # Also return diagnostics about how many months occurred in that year.
  pick_final_month_mayor <- function(mayor_month_df){
    # one row per fips-year-month; uses 'total_votes_raw' already aggregated
    base <- mayor_month_df %>%
      dplyr::select(fips, year, month, total_votes_raw)
    
    choice <- base %>%
      dplyr::group_by(fips, year) %>%
      dplyr::arrange(dplyr::desc(total_votes_raw), dplyr::desc(month), .by_group = TRUE) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::mutate(month_final = month) %>%
      dplyr::select(fips, year, month_final)
    
    diag <- base %>%
      dplyr::group_by(fips, year) %>%
      dplyr::summarise(
        n_months = dplyr::n(),
        months_list = paste(sort(unique(month)), collapse = ","),
        max_votes = max(total_votes_raw, na.rm = TRUE),
        .groups = "drop"
      )
    
    dplyr::left_join(choice, diag, by = c("fips","year"))
  }

# A) contest-month summaries (one row per city-year-month)
mayor_month <- ledb_cand %>%
  dplyr::filter(office_consolidated == "Mayor") %>%
  dplyr::group_by(fips, geo_name, state_abb, year, month) %>%
  dplyr::summarise(
    total_votes_raw = sum(votes, na.rm = TRUE),
    dem_votes_raw   = sum(dplyr::if_else(is_dem(pid_est), votes, 0), na.rm = TRUE),
    
    # set dem_share = NA in nonpartisan/unknown cases (no party labels at all)
    any_party_label = any(!is.na(pid_est) & pid_est %in% c("D","R")),
    dem_share       = dplyr::if_else(
      any_party_label & total_votes_raw > 0,
      dem_votes_raw / total_votes_raw, NA_real_),
    
    n_candidates = dplyr::n_distinct(ledb_candid),
    
    # incumbent winner? (any winning row with incumbent==1)
    inc_win = as.integer(any(winner == "win" & incumbent == 1, na.rm = TRUE)),
    
    # candidate-side composition (shares among candidates on the ballot)
    female_cand_share = mean(gender_est == "F", na.rm = TRUE),
    poc_cand_share    = mean(is_poc(race_est),  na.rm = TRUE),
    
    # winner composition (robust to 1+ seats, though mayor should be 1)
    wins_total      = sum(winner == "win", na.rm = TRUE),
    female_wins     = sum(winner == "win" & gender_est == "F", na.rm = TRUE),
    poc_wins        = sum(winner == "win" & is_poc(race_est),  na.rm = TRUE),
    female_win_share = dplyr::if_else(wins_total > 0, female_wins / wins_total, NA_real_),
    poc_win_share    = dplyr::if_else(wins_total > 0, poc_wins    / wins_total, NA_real_),
    
    # winner id by raw votes in THIS month (convenient for QC)
    winner_id = ledb_candid[which.max(votes)],
    
    .groups = "drop"
  )

# B) pick the final month per city-year
mayor_final_info <- pick_final_month_mayor(mayor_month)

# & keep only the final month row, but attach diagnostics
mayor_by_cityyear <- mayor_month %>%
  dplyr::inner_join(
    mayor_final_info,
    by = c("fips","year","month" = "month_final")
  ) %>%
  dplyr::mutate(
    had_multiple_months = n_months > 1L,
    # simple “general-ish” sanity flag: final not in Oct/Nov is unusual
    final_not_oct_or_nov = !(month %in% c(10L, 11L))
  ) %>%
  dplyr::transmute(
    place_fips = fips, geo_name, state_abb, year,
    office_consolidated = "Mayor",
    month_final = month, n_months, months_list, had_multiple_months, final_not_oct_or_nov,
    total_votes = total_votes_raw,
    dem_votes   = dem_votes_raw,
    dem_share,
    n_candidates, inc_win, winner_id,
    female_cand_share, poc_cand_share,
    female_win_share,  poc_win_share
  )


# C) Merge covars/flags and build turnout (keep specials, just flagged)
mayor_panel <- base_cities %>%
  dplyr::left_join(mayor_by_cityyear, by = "place_fips") %>%
  dplyr::filter(!is.na(year)) %>%                         # drop city-years that never had a mayor row
  dplyr::mutate(
    turnout_pct = dplyr::if_else(
      vap90_city > 0 & !is.na(total_votes),
      100 * total_votes / vap90_city, NA_real_
    ),
    turnout_pct = dplyr::if_else(turnout_pct > 120, NA_real_, turnout_pct)
  )

# D) Sanity
  # How often do we see multiple months for mayors?
  mayor_panel %>%
    dplyr::count(had_multiple_months, final_not_oct_or_nov) %>%
    print(n = Inf)
  
  # Any “final” months outside Oct/Nov? (not errors—just flagging)
  mayor_panel %>%
    dplyr::filter(final_not_oct_or_nov) %>%
    dplyr::count(state_abb, geo_name, year, month_final) %>%
    head(20) %>% print(n = Inf)

# ========================= COUNCIL PANEL ======================================
  # Helper: pick the “final” month per (fips, year, district)
  # Rule = month with the largest seat-normalized total votes; ties → later month
  pick_final_month_district <- function(df){
    df %>%
      dplyr::group_by(fips, year, district, month) %>%
      dplyr::summarise(votes_m = sum(total_votes_ps, na.rm = TRUE), .groups = "drop_last") %>%
      dplyr::group_by(fips, year, district) %>%
      dplyr::arrange(dplyr::desc(votes_m), dplyr::desc(month), .by_group = TRUE) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::select(fips, year, district, month_final = month)
  }

# A) contest-month summaries (one row per fips × year × district × month)
council_month <- ledb_cand %>%
  dplyr::filter(office_consolidated == "City Council") %>%
  dplyr::group_by(fips, geo_name, state_abb, year, district, month) %>%
  dplyr::summarise(
    # seats in this district-month (robust to the small declared/marked mismatch)
    n_seats_decl    = suppressWarnings(max(n_winners, na.rm = TRUE)),
    n_seats_marked  = sum(winner == "win", na.rm = TRUE),
    n_seats         = dplyr::coalesce(pmax(n_seats_decl, n_seats_marked), 0),

    # raw totals
    total_votes_raw = sum(votes, na.rm = TRUE),
    dem_votes_raw   = sum(dplyr::if_else(is_dem(pid_est), votes, 0), na.rm = TRUE),

    # seat-normalized vote totals (so at-large 3-seat races are comparable to SMDs)
    total_votes_ps  = total_votes_raw / pmax(n_seats, 1),
    dem_votes_ps    = dem_votes_raw   / pmax(n_seats, 1),
    dem_share       = dplyr::if_else(total_votes_ps > 0, dem_votes_ps / total_votes_ps, NA_real_),

    # incumbency success in this district-month (for later seat-weighted aggregation)
    inc_wins        = sum(as.integer(winner == "win" & incumbent == 1), na.rm = TRUE),
    wins_total      = sum(winner == "win", na.rm = TRUE),

    # candidate composition (district-month)
    n_candidates_ct   = dplyr::n_distinct(ledb_candid),
    female_cand_share = mean(gender_est == "F", na.rm = TRUE),
    poc_cand_share    = mean(is_poc(race_est),  na.rm = TRUE),

    # winner composition (district-month, seat-weighted)
    female_wins     = sum(as.integer(winner == "win" & gender_est == "F"), na.rm = TRUE),
    poc_wins        = sum(as.integer(winner == "win" & is_poc(race_est)), na.rm = TRUE),

    .groups = "drop"
  )

# B) pick the final month per (fips, year, district)
council_final_month <- pick_final_month_district(council_month)

# C) keep only that month for each district, then roll up districts → city-year
council_by_cityyear <- council_month %>%
  dplyr::inner_join(council_final_month,
                    by = c("fips","year","district","month" = "month_final")) %>%
  dplyr::group_by(fips, geo_name, state_abb, year) %>%
  dplyr::summarise(
    # per-seat vote intensity (sum across districts)
    total_votes_ps = sum(total_votes_ps, na.rm = TRUE),
    dem_votes_ps   = sum(dem_votes_ps,   na.rm = TRUE),
    dem_share      = dplyr::if_else(total_votes_ps > 0, dem_votes_ps / total_votes_ps, NA_real_),

    # seats & incumbents (sum of district seats elected this year)
    n_seats_total  = sum(pmax(n_seats, 0), na.rm = TRUE),
    inc_win_share  = dplyr::if_else(n_seats_total > 0,
                                    sum(inc_wins, na.rm = TRUE) / n_seats_total, NA_real_),

    # contest structure
    n_contests     = dplyr::n(),                                  # districts (final)
    n_candidates   = sum(n_candidates_ct, na.rm = TRUE),

    # candidate composition (weight by #candidates in each district)
    female_cand_share = dplyr::if_else(n_candidates > 0,
      sum(female_cand_share * n_candidates_ct, na.rm = TRUE) / n_candidates, NA_real_),
    poc_cand_share = dplyr::if_else(n_candidates > 0,
      sum(poc_cand_share * n_candidates_ct, na.rm = TRUE) / n_candidates, NA_real_),

    # winner composition (seat-weighted)
    female_win_share = dplyr::if_else(n_seats_total > 0,
      sum(female_wins, na.rm = TRUE) / n_seats_total, NA_real_),
    poc_win_share    = dplyr::if_else(n_seats_total > 0,
      sum(poc_wins,    na.rm = TRUE) / n_seats_total, NA_real_),

    .groups = "drop"
  ) %>%
  dplyr::mutate(office_consolidated = "City Council",
                place_fips = fips) %>%
  dplyr::select(place_fips, geo_name, state_abb, year, office_consolidated,
                total_votes_ps, dem_votes_ps, dem_share,
                n_seats_total, inc_win_share, n_contests, n_candidates,
                female_cand_share, poc_cand_share, female_win_share, poc_win_share)

# D) merge covars/flags and build votes-per-seat intensity (% of city VAP)
council_panel <- base_cities %>%
  dplyr::left_join(council_by_cityyear, by = "place_fips") %>%
  dplyr::filter(!is.na(year)) %>%
  dplyr::mutate(
    votes_per_seat_pct_vap = dplyr::if_else(
      vap90_city > 0 & !is.na(total_votes_ps),
      100 * total_votes_ps / vap90_city, NA_real_
    )
  )

# E) Sanity
cat("\n== Council diagnostics ==\n")

  # How many months appear per city-year in the raw council contests?
  c_mo_any <- ledb_cand %>%
    dplyr::filter(office_consolidated == "City Council") %>%
    dplyr::group_by(fips, year) %>%
    dplyr::summarise(n_months_any = dplyr::n_distinct(month),
                     months_any   = paste(sort(unique(month)), collapse=","),
                     .groups="drop")
  
  # How many distinct final months across districts in a city-year?
  c_mo_final <- council_final_month %>%
    dplyr::group_by(fips, year) %>%
    dplyr::summarise(n_months_final = dplyr::n_distinct(month_final),
                     months_final   = paste(sort(unique(month_final)), collapse=","),
                     .groups="drop")
  
  c_diag <- c_mo_any %>%
    dplyr::left_join(c_mo_final, by = c("fips","year")) %>%
    dplyr::mutate(had_multiple_months_any   = n_months_any   > 1,
                  had_multiple_months_final = n_months_final > 1)
  
  print(head(c_diag, 20))
  
  # Seat totals sanity (should be reasonable integers)
  cat("\nSample of (city-year) seats elected and vote intensity:\n")
  council_panel %>%
    dplyr::select(place_fips, geo_name, year, n_seats_total, total_votes_ps, votes_per_seat_pct_vap) %>%
    dplyr::arrange(place_fips, year) %>%
    dplyr::sample_n(20) %>%
    print()
  



## 5. Combined EDA panel-------------------------------------------------------- 
city_panel <- bind_rows(
  mayor_panel   %>% mutate(panel = "mayor"),
  council_panel %>% mutate(panel = "council")
)


## 6. Save ---------------------------------------------------------------------
write_rds(mayor_panel,   here("data_clean","1_intermediate","mayor_panel.rds"))
write_rds(council_panel, here("data_clean","1_intermediate","council_panel.rds"))
write_rds(city_panel,    here("data_clean","1_intermediate","city_panel.rds"))


## 7. More Diagnostics ---------------------------------------------------------
mayor_panel %>% summarise(n_cities=n_distinct(place_fips), n_rows=n())
council_panel %>% summarise(n_cities=n_distinct(place_fips), n_rows=n())
council_panel %>% summarise(mean_inc_win_share = mean(inc_win_share, na.rm=TRUE))

# Core Integretity (both panels)
  library(tidyverse)
  
  # 0) dimensions
  mayor_panel   %>% summarise(n_cities=n_distinct(place_fips), n_years=n_distinct(year), n_rows=n())
  council_panel %>% summarise(n_cities=n_distinct(place_fips), n_years=n_distinct(year), n_rows=n())
  
  # 1) uniqueness of keys
  mayor_panel   %>% count(place_fips, year)   %>% filter(n>1)   # should be empty
  council_panel %>% count(place_fips, year)   %>% filter(n>1)   # should be empty
  
  # 2) join health (covars/flags should be present)
  mayor_panel   %>% summarise(na_covars = sum(is.na(pop90_city)))
  council_panel %>% summarise(na_covars = sum(is.na(pop90_city)))
  
  # 3) treatment group sizes (same logic as before)
  for(df in list(mayor_panel=mayor_panel, council_panel=council_panel)){
    print(names(df))
  }
  mayor_panel %>% summarise(
    treated   = sum(treated_dummy,  na.rm=TRUE),
    nearmiss  = sum(nearmiss_dummy, na.rm=TRUE),
    future    = sum(future_dummy,   na.rm=TRUE),
    untreated = sum(!treated_dummy & !nearmiss_dummy & !future_dummy)
  )
  council_panel %>% summarise(
    treated   = sum(treated_dummy,  na.rm=TRUE),
    nearmiss  = sum(nearmiss_dummy, na.rm=TRUE),
    future    = sum(future_dummy,   na.rm=TRUE),
    untreated = sum(!treated_dummy & !nearmiss_dummy & !future_dummy)
  )
  
  # 4) range checks for shares (should all be in [0,1])
  rng_check <- function(df, cols){
    df %>% summarise(across(all_of(cols),
                            ~ c(min=min(.x, na.rm=TRUE), max=max(.x, na.rm=TRUE)) %>% list()))
  }
  rng_check(mayor_panel,
            c("dem_share","female_cand_share","poc_cand_share","female_win_share","poc_win_share"))
  rng_check(council_panel,
            c("dem_share","inc_win_share","female_cand_share","poc_cand_share","female_win_share","poc_win_share"))
  
  # 5) name/FIPS consistency (one name per FIPS)
  mayor_panel   %>% group_by(place_fips)   %>% summarise(n_names=n_distinct(place_name_treat))   %>% filter(n_names>1)
  council_panel %>% group_by(place_fips)   %>% summarise(n_names=n_distinct(place_name_treat))   %>% filter(n_names>1)
  

# Timing & DiD Readiness (both panels)
  
  # rel_year for treated cities only; how wide is support?
  rel_support <- function(df){
    df %>%
      filter(treated_dummy, !is.na(treat_yr), !is.na(year)) %>%
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
  rel_support(mayor_panel)
  rel_support(council_panel)
  
  # share of treated cities with at least one pre- and one post- observation
  prepost <- function(df){
    df %>%
      filter(treated_dummy, !is.na(treat_yr), !is.na(year)) %>%
      mutate(rel_year = year - treat_yr) %>%
      group_by(place_fips) %>%
      summarise(has_pre = any(rel_year < 0), has_post = any(rel_year >= 0), .groups="drop") %>%
      summarise(
        n_treated = n(),
        both_pre_post = sum(has_pre & has_post),
        share_both = both_pre_post / n_treated
      )
  }
  prepost(mayor_panel)
  prepost(council_panel)
  

# Basic Prior Stats Replicated
  
  # obs per city
  mayor_panel %>% count(place_fips) %>% summarise(min_obs=min(n), median=median(n), max_obs=max(n), mean_obs=mean(n))
  council_panel %>% count(place_fips) %>% summarise(min_obs=min(n), median=median(n), max_obs=max(n), mean_obs=mean(n))
  
  # years covered
  mayor_panel %>% summarise(first_year=min(year, na.rm=TRUE), last_year=max(year, na.rm=TRUE))
  council_panel %>% summarise(first_year=min(year, na.rm=TRUE), last_year=max(year, na.rm=TRUE))
  
  # EZ share vs Dem share correlation (crude)
  mayor_panel   %>% filter(!is.na(share_ez), !is.na(dem_share))   %>% summarise(cor_ez_dem = cor(share_ez, dem_share))
  council_panel %>% filter(!is.na(share_ez), !is.na(dem_share))   %>% summarise(cor_ez_dem = cor(share_ez, dem_share))

# Mayor-Panel Specific 
  # 1) one row per city-year (after final-month pick)
  mayor_panel %>% count(place_fips, year) %>% filter(n>1)
  
  # 2) turnout sanity and outliers
  mayor_panel %>% summarise(
    over_100 = mean(turnout_pct > 100, na.rm=TRUE),
    p50      = median(turnout_pct, na.rm=TRUE),
    p90      = quantile(turnout_pct, .9, na.rm=TRUE),
    max      = max(turnout_pct, na.rm=TRUE)
  )
  
  # 3) incumbent win rate
  mayor_panel %>% summarise(inc_rate = mean(inc_win==1, na.rm=TRUE))
  
  # 4) contest structure
  mayor_panel %>% summarise(
    median_candidates = median(n_candidates, na.rm=TRUE),
    p90_candidates    = quantile(n_candidates, .9, na.rm=TRUE),
    max_candidates    = max(n_candidates, na.rm=TRUE)
  )
  
  # 5) composition ranges
  mayor_panel %>% summarise(
    female_cand_p50 = median(female_cand_share, na.rm=TRUE),
    poc_cand_p50    = median(poc_cand_share,    na.rm=TRUE),
    female_win_p50  = median(female_win_share,  na.rm=TRUE),
    poc_win_p50     = median(poc_win_share,     na.rm=TRUE)
  )
  
  # 6) pre/post descriptive shift (treated only; informal)
  mayor_panel %>%
    filter(treated_dummy, !is.na(treat_yr)) %>%
    mutate(post = year >= treat_yr) %>%
    group_by(post) %>%
    summarise(
      n = n(),
      mean_dem = mean(dem_share, na.rm=TRUE),
      mean_turnout = mean(turnout_pct, na.rm=TRUE),
      mean_inc = mean(inc_win, na.rm=TRUE)
    )

  # For mayors: verify the chosen (city,year) is the month with max votes
  mayor_check <- ledb_cand %>%
    filter(office_consolidated=="Mayor") %>%
    group_by(fips, year, month) %>%
    summarise(votes_m=sum(votes,na.rm=TRUE), .groups="drop_last") %>%
    group_by(fips, year) %>%
    mutate(rank = rank(-votes_m, ties.method="min"),
           is_latest = month==max(month)) %>%
    ungroup()
  
  # How often does our kept month have rank==1 OR (tie on votes & latest month)?
  # (Assumes you carried month_final in your build; if not, re-derive it with the helper.)
  
# Council-Panel Specific 
  # 1) one row per city-year
  council_panel %>% count(place_fips, year) %>% filter(n>1)
  
  # 2) seat & contest sanity
  council_panel %>% summarise(
    med_seats  = median(n_seats_total, na.rm=TRUE),
    max_seats  = max(n_seats_total, na.rm=TRUE),
    med_cont   = median(n_contests, na.rm=TRUE),
    p90_cont   = quantile(n_contests, .9, na.rm=TRUE),
    max_cont   = max(n_contests, na.rm=TRUE)
  )
  
  # 3) per-seat vote intensity sanity (our “turnout” analog)
  council_panel %>% summarise(
    p50_vps = median(votes_per_seat_pct_vap, na.rm=TRUE),
    p90_vps = quantile(votes_per_seat_pct_vap, .9, na.rm=TRUE),
    max_vps = max(votes_per_seat_pct_vap, na.rm=TRUE)
  )
  
  # 4) incumbent seat share distribution
  council_panel %>% summarise(
    inc_p50 = median(inc_win_share, na.rm=TRUE),
    inc_p90 = quantile(inc_win_share, .9, na.rm=TRUE)
  )
  
  # 5) composition metrics (candidate vs winner)
  council_panel %>% summarise(
    female_cand_p50 = median(female_cand_share, na.rm=TRUE),
    poc_cand_p50    = median(poc_cand_share,    na.rm=TRUE),
    female_win_p50  = median(female_win_share,  na.rm=TRUE),
    poc_win_p50     = median(poc_win_share,     na.rm=TRUE)
  )
  
  # 6) pre/post descriptive shift (treated only; informal)
  council_panel %>%
    filter(treated_dummy, !is.na(treat_yr)) %>%
    mutate(post = year >= treat_yr) %>%
    group_by(post) %>%
    summarise(
      n = n(),
      mean_dem   = mean(dem_share, na.rm=TRUE),
      mean_vps   = mean(votes_per_seat_pct_vap, na.rm=TRUE),
      mean_inc   = mean(inc_win_share, na.rm=TRUE)
    )
  
# Cross-Panel Alignment
  # cities that appear in both panels (any year)
  inner_join(
    mayor_panel   %>% distinct(place_fips) %>% mutate(in_mayor=1),
    council_panel %>% distinct(place_fips) %>% mutate(in_council=1),
    by="place_fips"
  ) %>% summarise(n_both = n())
  
  # treated city overlap
  inner_join(
    mayor_panel   %>% filter(treated_dummy)   %>% distinct(place_fips),
    council_panel %>% filter(treated_dummy)   %>% distinct(place_fips),
    by="place_fips"
  ) %>% summarise(treated_in_both = n())
  
#Outlier & Hygiene Checks 
  # extreme values for quick manual review
  mayor_panel %>% arrange(desc(turnout_pct)) %>% select(place_name_treat, state_abb, year, turnout_pct) %>% slice_head(n=10)
  council_panel %>% arrange(desc(votes_per_seat_pct_vap)) %>% select(place_name_treat, state_abb, year, votes_per_seat_pct_vap) %>% slice_head(n=10)
  
  # implausible shares outside [0,1] or NA denom effects
  bad_mayor <- mayor_panel %>%
    filter(dem_share < 0 | dem_share > 1 |
             female_cand_share < 0 | female_cand_share > 1 |
             female_win_share < 0 | female_win_share > 1)
  bad_council <- council_panel %>%
    filter(dem_share < 0 | dem_share > 1 |
             inc_win_share < 0 | inc_win_share > 1 |
             female_cand_share < 0 | female_cand_share > 1 |
             female_win_share < 0 | female_win_share > 1)
  
  nrow(bad_mayor); nrow(bad_council)

# Richer Summaries 
  #skim each panel 
  if (!requireNamespace("skimr", quietly=TRUE)) install.packages("skimr")
  skimr::skim(mayor_panel)
  skimr::skim(council_panel)
  
  #yearly trends: 
  library(ggplot2)
  
  mayor_panel %>%
    group_by(year) %>%
    summarise(avg_candidates = mean(n_candidates, na.rm=TRUE)) %>%
    ggplot(aes(year, avg_candidates)) + geom_line() + labs(title="Mayor: Avg # Candidates")
  
  council_panel %>%
    group_by(year) %>%
    summarise(avg_vps = mean(votes_per_seat_pct_vap, na.rm=TRUE)) %>%
    ggplot(aes(year, avg_vps)) + geom_line() + labs(title="Council: Votes-per-Seat % of VAP")
  

