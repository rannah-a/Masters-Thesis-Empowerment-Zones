###############################################################################
# 07_robustness_A4R_EB.R —— Robustness & sensitivity checks
###############################################################################

###############################################################################
# ROBUSTNESS CHECK 
# #1 Reuse baseline A4R Division→Region + EB (exact→same kernel fallback),
              #   but refit weights under donor-pool restrictions:
              #   - "near_miss_only"
              #   - "near_miss_plus_future"
              #   - "future_only"
###############################################################################

here::i_am("code/07_robustness_models.R")

## Packages -----------------------------------------------------------------
suppressPackageStartupMessages({
  for (p in c("here","tidyverse","fixest","ebal")) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  }
  library(here); library(tidyverse); library(fixest); library(ebal)
})

theme_set(theme_minimal())

## Load the panels you used in Script 6 -------------------------------------
# (We need raw panels to rebuild city_table & donors; weights must be refit)
panel_mayor   <- readr::read_rds(here("data_clean","1_intermediate","city_panel_analysis_mayor.rds"))
panel_council <- readr::read_rds(here("data_clean","1_intermediate","city_panel_analysis_council.rds"))

## Prep helpers (respect stored month-aware post) --------------------------
  prep_unbalanced <- function(df) {
    df %>%
      mutate(
        # Use existing post if present (month-aware). If not, fall back safely.
        post2 = case_when(
          "post" %in% names(df)                    ~ as.integer(!!as.name("post")),
          "rel_year" %in% names(df) & treated_dummy == 1L ~ as.integer(rel_year >= 0),
          treated_dummy == 1L & !is.na(treat_yr)  ~ as.integer(year >= treat_yr),
          TRUE                                     ~ 0L
        ),
        # sunab cohort g: treat year for treated, 0 = never-treated
        g = if_else(treated_dummy == 1L & !is.na(treat_yr),
                    as.integer(treat_yr), 0L),
        year = as.integer(year),
        treated_dummy = as.integer(treated_dummy)
      )
  }

pm <- prep_unbalanced(panel_mayor)
pc <- prep_unbalanced(panel_council)

stopifnot(all(c("place_fips","year","treated_dummy","g","post2") %in% names(pm)))
stopifnot(all(c("place_fips","year","treated_dummy","g","post2") %in% names(pc)))

# --- robust sdiff + weight diagnostics (unchanged) ----------------------------
.wvar_unbiased <- function(x, w = NULL){
  x <- as.numeric(x); if (is.null(w)) w <- rep(1, length(x))
  ok <- is.finite(x) & is.finite(w) & w >= 0; x <- x[ok]; w <- w[ok]
  if (length(x) < 2 || sum(w) <= 0) return(NA_real_)
  w <- w / sum(w); mu <- sum(w * x)
  eff_n <- (sum(w)^2) / sum(w^2); if (!is.finite(eff_n) || eff_n <= 1) return(NA_real_)
  sum(w * (x - mu)^2) * (eff_n / (eff_n - 1))
}
sdiff <- function(x, t, w = NULL){
  t <- as.integer(t == 1L)
  if (is.null(w)) {
    m1 <- mean(x[t==1], na.rm=TRUE); v1 <- var(x[t==1], na.rm=TRUE)
    m0 <- mean(x[t==0], na.rm=TRUE); v0 <- var(x[t==0], na.rm=TRUE)
  } else {
    w <- as.numeric(w); w[!is.finite(w) | w < 0] <- 0
    mw <- mean(w, na.rm = TRUE); if (is.finite(mw) && mw > 0) w <- w / mw
    m1 <- tryCatch(weighted.mean(x[t==1], w[t==1], na.rm=TRUE), error = function(e) NA_real_)
    m0 <- tryCatch(weighted.mean(x[t==0], w[t==0], na.rm=TRUE), error = function(e) NA_real_)
    v1 <- .wvar_unbiased(x[t==1], w[t==1]); v0 <- .wvar_unbiased(x[t==0], w[t==0])
  }
  s  <- sqrt((v1 + v0)/2); if (!is.finite(s) || s == 0) return(NA_real_)
  (m1 - m0) / s
}
weight_diag <- function(w_ctrl){
  w <- w_ctrl[is.finite(w_ctrl) & w_ctrl > 0]
  if (!length(w)) return(tibble(donors_kept=0, sum_w=0, ess=0, ess_share=0,
                                max_w=NA_real_, p99_w=NA_real_, p95_w=NA_real_))
  ess <- (sum(w)^2) / sum(w^2)
  tibble(donors_kept=length(w), sum_w=sum(w), ess=ess, ess_share=ess/length(w),
         max_w=max(w), p99_w=quantile(w, .99, names=FALSE), p95_w=quantile(w, .95, names=FALSE))
}

# --- Census state → division/region (unchanged) --------------------------------
state_div <- tribble(
  ~state_fips,~division,~region,
  "01",6,3,"02",9,4,"04",8,4,"05",7,3,"06",9,4,"08",8,4,"09",1,1,"10",5,3,"11",5,3,
  "12",5,3,"13",5,3,"15",9,4,"16",8,4,"17",3,2,"18",3,2,"19",4,2,"20",4,2,"21",6,3,
  "22",7,3,"23",1,1,"24",5,3,"25",1,1,"26",3,2,"27",4,2,"28",6,3,"29",4,2,"30",8,4,
  "31",4,2,"32",8,4,"33",1,1,"34",2,1,"35",8,4,"36",2,1,"37",5,3,"38",4,2,"39",3,2,
  "40",7,3,"41",9,4,"42",2,1,"44",1,1,"45",5,3,"46",4,2,"47",6,3,"48",7,3,"49",8,4,
  "50",1,1,"51",5,3,"53",9,4,"54",5,3,"55",3,2,"56",8,4
)


## --- City-level donor flags from row-level dummies ---------------------------
derive_city_donor_flags <- function(panel_df) {
  panel_df %>%
    group_by(place_fips) %>%
    summarise(
      near_miss_ctrl = as.integer(any(nearmiss_dummy %in% c(TRUE, 1L))),
      future_ctrl    = as.integer(any(future_dummy  %in% c(1L, TRUE))),
      .groups = "drop"
    )
}

## Quick sanity (optional; prints once):
pm %>%
  derive_city_donor_flags() %>%
  summarise(n_nearmiss = sum(near_miss_ctrl==1),
            n_future   = sum(future_ctrl==1)) %>% print()

### --- Build city table (always creates donor flags internally) ---------------
build_city_table <- function(panel_df, y, pre_years = 1989:1993){
  
  df_geo <- panel_df %>%
    distinct(place_fips, .keep_all = TRUE) %>%
    mutate(state_fips = substr(place_fips, 1, 2)) %>%
    left_join(state_div, by = "state_fips")
  
  # city-level donor flags (from panel dummies)
  donor_flags <- derive_city_donor_flags(panel_df)
  
  covs <- df_geo %>%
    transmute(
      place_fips,
      ln_pop90       = log(pmax(pop90_city, 1)),
      poverty_rate_90, unemploy_90, black_city_90,
      division = as.integer(division),
      region   = as.integer(region)
    ) %>%
    left_join(donor_flags, by = "place_fips") %>%
    mutate(
      near_miss_ctrl = replace_na(near_miss_ctrl, 0L),
      future_ctrl    = replace_na(future_ctrl,    0L)
    )
  
  pre_y <- panel_df %>%
    filter(year %in% pre_years) %>%
    select(place_fips, year, !!sym(y)) %>%
    group_by(place_fips) %>%
    summarise(mean_pre = mean(.data[[y]], na.rm = TRUE), .groups="drop")
  
  trt <- panel_df %>%
    group_by(place_fips) %>%
    summarise(treated_city = as.integer(any(treated_dummy == 1L)), .groups="drop")
  
  covs %>%
    left_join(pre_y, by="place_fips") %>%
    left_join(trt,  by="place_fips") %>%
    mutate(treated_city = replace_na(treated_city, 0L))
}

## ----------  Small adapter so the wrapper can “find” a builder ----------
# Your runner prints and returns a list with $weights — perfect.
a4r_eb_weights <- function(panel_df, outcome, restriction){
  run_A4R_EB_with_restriction(
    panel_df    = panel_df,
    panel_name  = "Panel",     # just for console labeling
    outcome     = outcome,
    restrict    = restriction
  )
}
# Alias (the wrapper will also try this name):
get_a4r_weights <- a4r_eb_weights

# --- A4R donor pool (unchanged), with an extra "restrict" filter ---------------
# `restrict` ∈ {"near_miss_only","near_miss_plus_future","future_only", "all"}
donor_pool_union_restricted <- function(city_tab, vars,
                                        min_div = 20, min_reg = 40,
                                        restrict = c("all","near_miss_only","near_miss_plus_future","future_only")){
  restrict <- match.arg(restrict)
  donors_union <- character(0)
  trt_rows <- city_tab %>% filter(treated_city == 1L)
  
  # choose the control-eligibility mask once
  ctrl_mask <- case_when(
    restrict == "near_miss_only"       ~ city_tab$near_miss_ctrl == 1L,
    restrict == "near_miss_plus_future"~ (city_tab$near_miss_ctrl == 1L | city_tab$future_ctrl == 1L),
    restrict == "future_only"          ~ city_tab$future_ctrl == 1L,
    TRUE                               ~ city_tab$treated_city == 0L  # "all"
  )
  
  for (i in seq_len(nrow(trt_rows))){
    this_div <- trt_rows$division[i]
    this_reg <- trt_rows$region[i]
    
    cand_div <- city_tab %>%
      filter(treated_city == 0L, division == this_div, ctrl_mask) %>%
      drop_na(all_of(vars)) %>% pull(place_fips) %>% unique()
    
    if (length(cand_div) >= min_div) { donors_union <- union(donors_union, cand_div); next }
    
    cand_reg <- city_tab %>%
      filter(treated_city == 0L, region == this_reg, ctrl_mask) %>%
      drop_na(all_of(vars)) %>% pull(place_fips) %>% unique()
    
    if (length(cand_reg) >= min_reg) { donors_union <- union(donors_union, cand_reg); next }
    
    cand_nat <- city_tab %>%
      filter(treated_city == 0L, ctrl_mask) %>%
      drop_na(all_of(vars)) %>% pull(place_fips) %>% unique()
    
    donors_union <- union(donors_union, cand_nat)
  }
  donors_union
}

# --- EB with the SAME fallback you used in Script 6 (unchanged) ---------------
eb_fit_with_fallback <- function(city_tab, vars, donors_union){
  Tset <- city_tab %>% filter(treated_city == 1L) %>% drop_na(all_of(vars))
  Cset <- city_tab %>% filter(treated_city == 0L, place_fips %in% donors_union) %>% drop_na(all_of(vars))
  
  target <- colMeans(as.matrix(Tset[, vars, drop=FALSE]), na.rm=TRUE)
  Xc     <- as.matrix(Cset[, vars, drop=FALSE])
  
  eb_try <- try(ebal::ebalance(target.margins = target, X = Xc,
                               base.weights = rep(1, nrow(Xc)),
                               maxit = 10000, constraint.tolerance = 1e-8),
                silent = TRUE)
  
  if (inherits(eb_try, "try-error") || is.null(eb_try$w)) {
    message("EB exact failed; switching to baseline kernel fallback (same as Script 6).")
    
    xc_s  <- scale(Xc)
    tgt_s <- as.numeric(scale(matrix(target, nrow=1),
                              center = attr(xc_s, "scaled:center"),
                              scale  = attr(xc_s, "scaled:scale")))
    j_black <- match("black_city_90", vars)
    lambda_black <- 1.05
    
    diff_mat <- xc_s - matrix(tgt_s, nrow = nrow(xc_s), ncol = ncol(xc_s), byrow = TRUE)
    d2 <- rowSums(diff_mat^2)
    if (!is.na(j_black)) d2 <- d2 + (lambda_black^2 - 1) * (diff_mat[, j_black]^2)
    
    h  <- median(d2, na.rm=TRUE); if (!is.finite(h) || h <= 1e-8) h <- mean(d2, na.rm=TRUE) + 1
    w  <- exp(- d2 / h)
    if (mean(w, na.rm=TRUE) > 0) w <- w / mean(w, na.rm=TRUE)
  } else {
    w <- as.numeric(eb_try$w)
    if (mean(w, na.rm=TRUE) > 0) w <- w / mean(w, na.rm=TRUE)
  }
  
  ctrl_weights <- tibble(place_fips = Cset$place_fips, w_eb = as.numeric(w))
  city_tab %>%
    left_join(ctrl_weights, by="place_fips") %>%
    mutate(w_eb = case_when(
      treated_city == 1L            ~ 1,
      is.finite(w_eb) & w_eb > 0    ~ w_eb,
      TRUE                          ~ 0
    )) %>%
    select(place_fips, treated_city, all_of(vars), w_eb)
}

# --- One runner that only changes the donor restriction -----------------------
run_A4R_EB_with_restriction <- function(panel_df, panel_name, outcome,
                                        restrict = c("near_miss_only","near_miss_plus_future","future_only","all")){
  restrict <- match.arg(restrict)
  cat("\n\n=== ", panel_name, " | ", outcome, " — A4R+EB with donors = ", restrict, " ===\n", sep="")
  city_tab <- build_city_table(panel_df, y = outcome, pre_years = 1989:1993)
  vars <- c("mean_pre", "ln_pop90", "poverty_rate_90", "unemploy_90", "black_city_90")
  
  print(city_tab %>% summarise(treated_cities = sum(treated_city==1L),
                               donor_cities   = sum(treated_city==0L)))
  
  donors_union <- donor_pool_union_restricted(city_tab, vars,
                                              min_div = 20, min_reg = 40,
                                              restrict = restrict)
  
  wtab <- eb_fit_with_fallback(city_tab, vars, donors_union)
  
  bal <- tibble(var = vars) %>%
    mutate(
      sdiff_unw = map_dbl(var, ~ sdiff(wtab[[.x]], wtab$treated_city, NULL)),
      sdiff_w   = map_dbl(var, ~ sdiff(wtab[[.x]], wtab$treated_city,
                                       w = ifelse(wtab$treated_city==1, 1, wtab$w_eb)))
    )
  wd  <- weight_diag(wtab$w_eb[wtab$treated_city==0])
  
  cat("\n-- Donor pool size kept (after restriction & A4R) --\n")
  print(tibble(donors_kept = sum(wtab$treated_city==0 & wtab$w_eb>0),
               total_controls = sum(wtab$treated_city==0)))
  
  cat("\n-- Std. diffs: before vs after weighting --\n")
  print(bal %>% mutate(across(starts_with("sdiff"), \(x) round(x, 3))))
  
  cat("\n-- Weight diagnostics (controls) --\n")
  print(wd %>% mutate(across(where(is.numeric), \(x) round(x, 3))))
  
  list(weights = wtab %>% select(place_fips, w_eb), balance = bal, wdiag = wd)
}

###  Attach A4R+EB weights, run models, save-------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr)
  library(fixest); library(readr)
  library(here)
})

#--- Output dirs (separate from Script 6 by prefixing 07_) ---------------------
out_root   <- here("data_clean", "2_output")
out_models <- file.path(out_root, "models")
out_tables <- file.path(out_root, "tables")
dir.create(out_models, showWarnings = FALSE, recursive = TRUE)
dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)

#--- Helper: find the weight builder from Sections 0–2 -------------------------
# Your Sec.0–2 wrapper should exist already. We try common function names:
get_a4r_weights_wrapper <- function(panel_df, outcome, restriction) {
  # preferred name
  if (exists("a4r_eb_weights")) {
    return(a4r_eb_weights(panel_df, outcome, restriction))
  }
  # fallback alias
  if (exists("get_a4r_weights")) {
    return(get_a4r_weights(panel_df, outcome, restriction))
  }
  stop("Could not find an A4R weight function. Expected `a4r_eb_weights()` or `get_a4r_weights()`.")
}

#--- Helper: standardize weights tibble to (place_fips, w_a4r) ----------------
standardize_w <- function(w_any) {
  # accept list(returned_by_wrapper) or bare tibble
  w_tbl <- if (is.list(w_any) && "weights" %in% names(w_any)) w_any$weights else w_any
  stopifnot(is.data.frame(w_tbl), "place_fips" %in% names(w_tbl))
  wcol <- setdiff(names(w_tbl), "place_fips")
  stopifnot(length(wcol) >= 1)
  w_tbl %>%
    select(place_fips, !!sym(wcol[1])) %>%
    rename(w_a4r = !!sym(wcol[1])) %>%
    mutate(w_a4r = as.numeric(w_a4r))
}

#--- Helper: attach weights (treated=1, non-donors=0) --------------------------
attach_a4r_weights <- function(panel_df, w_tbl) {
  panel_df %>%
    left_join(w_tbl, by = "place_fips") %>%
    mutate(
      w_a4r = dplyr::case_when(
        treated_dummy == 1L           ~ 1,                       # treated get unit weight
        is.finite(w_a4r) & w_a4r > 0  ~ w_a4r,                   # donors keep their weight
        TRUE                          ~ 0                        # controls outside donors = 0
      )
    )
}

#--- Estimators (weighted) -----------------------------------------------------
est_att_w <- function(df, y, wcol = "w_a4r") {
  d <- df %>% filter(!is.na(.data[[y]]), .data[[wcol]] > 0)
  feols(
    as.formula(paste0(y, " ~ post2:treated_dummy | place_fips + year")),
    data = d, weights = as.formula(paste0("~", wcol)), cluster = ~ place_fips
  )
}

est_es_w <- function(df, y, wcol = "w_a4r", ref = -1) {
  d <- df %>% filter(!is.na(.data[[y]]), .data[[wcol]] > 0)
  feols(
    as.formula(paste0(y, " ~ sunab(g, year, ref.p = ", ref, ") | place_fips + year")),
    data = d, weights = as.formula(paste0("~", wcol)), cluster = ~ place_fips
  )
}

#--- Small extractor for a compact ATT table -----------------------------------
att_row <- function(model, coef_name = "post2:treated_dummy") {
  ct <- fixest::coeftable(model)
  if (!coef_name %in% rownames(ct)) return(tibble(estimate = NA_real_, se = NA_real_, p_value = NA_real_))
  tibble(
    estimate = unname(ct[coef_name, 1]),
    se       = unname(ct[coef_name, 2]),
    p_value  = unname(ct[coef_name, 4])
  )
}

#--- Panels & outcomes ---------------------------------------------------------
PANELS <- list(
  Mayor = list(df = pm,  outcomes = c("dem_share", "turnout_pct")),
  Council = list(df = pc, outcomes = c("dem_share", "votes_per_seat_pct_vap"))
)

RESTRICTIONS <- c("near_miss_only", "near_miss_plus_future", "future_only")

#--- Loop: build weights -> attach -> estimate -> save -------------------------
att_collect <- list()  # will bind_rows into a CSV at end

for (rstr in RESTRICTIONS) {
  message("\n===== Running A4R-EB robustness for restriction: ", rstr, " =====")
  
  for (pname in names(PANELS)) {
    df_panel <- PANELS[[pname]]$df
    for (y in PANELS[[pname]]$outcomes) {
      
      # 1) Build city-level weights for this restriction × outcome
      w_any  <- get_a4r_weights_wrapper(df_panel, y, rstr)
      w_tbl  <- standardize_w(w_any)  # (place_fips, w_a4r)
      
      # 2) Attach to panel
      df_w <- attach_a4r_weights(df_panel, w_tbl)
      
      # 3) Weighted TWFE ATT and SUNAB ES
      m_att <- est_att_w(df_w, y, wcol = "w_a4r")
      m_es  <- est_es_w(df_w, y, wcol = "w_a4r", ref = -1)
      
      # 4) Save models with clean 07_ stubs
      stub <- paste0("07_A4R_", rstr, "_", str_to_lower(pname), "_", y)
      saveRDS(m_att, file = file.path(out_models, paste0(stub, "_ATT.rds")))
      saveRDS(m_es,  file = file.path(out_models, paste0(stub, "_ES.rds")))
      
      # 5) Add to compact ATT table
      att_collect[[length(att_collect) + 1]] <- dplyr::bind_cols(
        restriction = rstr,
        panel       = pname,
        outcome     = y,
        N           = nrow(df_w %>% filter(!is.na(.data[[y]]), w_a4r > 0)),
        cities      = dplyr::n_distinct(df_w$place_fips[df_w$w_a4r > 0]),
        att_row(m_att)
      )
      
      # 6) Console breadcrumb (nice to have)
      cat(sprintf("Saved: %s_[ATT|ES].rds | panel=%s outcome=%s N=%s cities=%s\n",
                  stub, pname, y,
                  nrow(df_w %>% filter(!is.na(.data[[y]]), w_a4r > 0)),
                  dplyr::n_distinct(df_w$place_fips[df_w$w_a4r > 0])))
    }
  }
}

#--- Write combined ATT CSV for all restrictions/panels/outcomes ---------------
att_tbl <- dplyr::bind_rows(att_collect) %>%
  relocate(restriction, panel, outcome, N, cities, estimate, se, p_value)

readr::write_csv(att_tbl, file.path(out_tables, "07_A4R_ATT_all.csv"))

cat("\nDone. Models saved under 2_output/models with 07_A4R_* stubs; ATT table at 2_output/tables/07_A4R_ATT_all.csv\n")




###############################################################################
# Robustness A-1b: Balance pre-paths + city trends
    # - Expands EB targets: mean_pre + pre_slope + several pre-lags
    # - Adds city-specific linear trends in ATT & SUNAB
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(stringr)
  library(fixest); library(readr); library(ggplot2); library(ggplotify); library(gridExtra)
})

theme_set(theme_minimal())

# ---------- Small utilities (reuse from earlier parts if already defined) -----

# unbiased weighted variance (already in your file; redefined here safely)
.wvar_unbiased <- function(x, w = NULL){
  x <- as.numeric(x); if (is.null(w)) w <- rep(1, length(x))
  ok <- is.finite(x) & is.finite(w) & w >= 0; x <- x[ok]; w <- w[ok]
  if (length(x) < 2 || sum(w) <= 0) return(NA_real_)
  w <- w / sum(w); mu <- sum(w * x)
  eff_n <- (sum(w)^2) / sum(w^2); if (!is.finite(eff_n) || eff_n <= 1) return(NA_real_)
  sum(w * (x - mu)^2) * (eff_n / (eff_n - 1))
}

sdiff <- function(x, t, w = NULL){
  t <- as.integer(t == 1L)
  if (is.null(w)) {
    m1 <- mean(x[t==1], na.rm=TRUE); v1 <- var(x[t==1], na.rm=TRUE)
    m0 <- mean(x[t==0], na.rm=TRUE); v0 <- var(x[t==0], na.rm=TRUE)
  } else {
    w <- as.numeric(w); w[!is.finite(w) | w < 0] <- 0
    mw <- mean(w, na.rm = TRUE); if (is.finite(mw) && mw > 0) w <- w / mw
    m1 <- tryCatch(weighted.mean(x[t==1], w[t==1], na.rm=TRUE), error = function(e) NA_real_)
    m0 <- tryCatch(weighted.mean(x[t==0], w[t==0], na.rm=TRUE), error = function(e) NA_real_)
    v1 <- .wvar_unbiased(x[t==1], w[t==1]); v0 <- .wvar_unbiased(x[t==0], w[t==0])
  }
  s  <- sqrt((v1 + v0)/2); if (!is.finite(s) || s == 0) return(NA_real_)
  (m1 - m0) / s
}

weight_diag <- function(w_ctrl){
  w <- w_ctrl[is.finite(w_ctrl) & w_ctrl > 0]
  if (!length(w)) return(tibble(donors_kept=0, sum_w=0, ess=0, ess_share=0,
                                max_w=NA_real_, p99_w=NA_real_, p95_w=NA_real_))
  ess <- (sum(w)^2) / sum(w^2)
  tibble(donors_kept=length(w), sum_w=sum(w), ess=ess, ess_share=ess/length(w),
         max_w=max(w), p99_w=quantile(w, .99, names=FALSE), p95_w=quantile(w, .95, names=FALSE))
}

# ---------- Donor pool rule (A4R) and EB fallback (reuse your logic) ----------
donor_pool_union_restricted <- function(city_tab, vars,
                                        min_div = 20, min_reg = 40,
                                        restrict = c("all","near_miss_only","near_miss_plus_future","future_only")){
  restrict <- match.arg(restrict)
  donors_union <- character(0)
  trt_rows <- city_tab %>% filter(treated_city == 1L)
  
  ctrl_mask <- dplyr::case_when(
    restrict == "near_miss_only"        ~ city_tab$near_miss_ctrl == 1L,
    restrict == "near_miss_plus_future" ~ (city_tab$near_miss_ctrl == 1L | city_tab$future_ctrl == 1L),
    restrict == "future_only"           ~ city_tab$future_ctrl == 1L,
    TRUE                                ~ city_tab$treated_city == 0L
  )
  
  for (i in seq_len(nrow(trt_rows))){
    this_div <- trt_rows$division[i]; this_reg <- trt_rows$region[i]
    
    cand_div <- city_tab %>% filter(treated_city==0L, division==this_div, ctrl_mask) %>%
      tidyr::drop_na(dplyr::all_of(vars)) %>% dplyr::pull(place_fips) %>% unique()
    if (length(cand_div) >= min_div) { donors_union <- union(donors_union, cand_div); next }
    
    cand_reg <- city_tab %>% filter(treated_city==0L, region==this_reg, ctrl_mask) %>%
      tidyr::drop_na(dplyr::all_of(vars)) %>% dplyr::pull(place_fips) %>% unique()
    if (length(cand_reg) >= min_reg) { donors_union <- union(donors_union, cand_reg); next }
    
    cand_nat <- city_tab %>% filter(treated_city==0L, ctrl_mask) %>%
      tidyr::drop_na(dplyr::all_of(vars)) %>% dplyr::pull(place_fips) %>% unique()
    donors_union <- union(donors_union, cand_nat)
  }
  donors_union
}

eb_fit_with_fallback <- function(city_tab, vars, donors_union){
  Tset <- city_tab %>% filter(treated_city == 1L) %>% tidyr::drop_na(dplyr::all_of(vars))
  Cset <- city_tab %>% filter(treated_city == 0L, place_fips %in% donors_union) %>%
    tidyr::drop_na(dplyr::all_of(vars))
  
  if (nrow(Cset) == 0) {
    stop("EB: No control donors with complete values for targets: ",
         paste(vars, collapse = ", "),
         ". Check the [DIAG] table above and consider using fewer pre-lag targets.")
  }
  
  target <- colMeans(as.matrix(Tset[, vars, drop=FALSE]), na.rm=TRUE)
  Xc     <- as.matrix(Cset[, vars, drop=FALSE])
  
  # ---- robust manual standardization (avoid NA from base::scale) ----
  center <- colMeans(Xc, na.rm = TRUE)
  scalev <- apply(Xc, 2, sd)
  scalev[!is.finite(scalev) | scalev < 1e-8] <- 1  # guard: zero/NA sd -> 1
  
  xc_s  <- sweep(Xc,    2, center, FUN = "-")
  xc_s  <- sweep(xc_s,  2, scalev, FUN = "/")
  tgt_s <- (target - center) / scalev
  
  j_black <- match("black_city_90", vars)
  lambda_black <- 1.05
  
    diff_mat <- xc_s - matrix(tgt_s, nrow = nrow(xc_s), ncol = ncol(xc_s), byrow = TRUE)
    d2 <- rowSums(diff_mat^2)
    if (!is.na(j_black)) d2 <- d2 + (lambda_black^2 - 1) * (diff_mat[, j_black]^2)
    
    # ---- sharp bandwidth and sticky kernel (no uniform escape) ----
    h <- stats::median(d2[is.finite(d2)], na.rm = TRUE) / 5
    if (!is.finite(h) || h <= 1e-8) h <- max(stats::median(d2[is.finite(d2)], na.rm = TRUE), 1)
    
    w <- exp(- d2 / h)
    w[!is.finite(w)] <- 0
    w <- pmax(w, 1e-12)
    
    ess_dbg <- (sum(w)^2) / sum(w^2)
    message(sprintf("[KERN] h=%.3f  sum=%.6f  min=%.6g  max=%.6g  ESS=%.2f/%d",
                    h, sum(w), min(w), max(w), ess_dbg, length(w)))
    
  # robust normalization
  m <- mean(w, na.rm = TRUE)
  if (is.finite(m) && m > 0) {
    w <- w / m
  } else {
    stop("Fallback/EB returned non-finite weights after normalization. ",
         "Likely too few donors with complete targets. See the [DIAG] counts printed above.")
  } 
  
  ctrl_weights <- tibble(place_fips = Cset$place_fips, w_eb = as.numeric(w))
  city_tab %>%
    left_join(ctrl_weights, by="place_fips") %>%
    mutate(w_eb = dplyr::case_when(
      treated_city == 1L         ~ 1,
      is.finite(w_eb) & w_eb > 0 ~ w_eb,
      TRUE                       ~ 0
    )) %>%
    select(place_fips, treated_city, dplyr::all_of(vars), w_eb)
}


# ---------- PRETREND builder: adds pre-lags and pre-slope to EB targets --------
# Choose the specific pre-years you want to match (keep it sparse to avoid drop-off)
.PRE_YEARS  <- 1989:1993
.PRE_LAGS   <- c(1989, 1991, 1993)   # can change to 1988/1990/1992 if you prefer

build_city_table_pretrend <- function(panel_df, y, pre_years = .PRE_YEARS, pre_lags = .PRE_LAGS){
  
  df_geo <- panel_df %>%
    distinct(place_fips, .keep_all = TRUE) %>%
    mutate(state_fips = substr(place_fips, 1, 2)) %>%
    left_join(state_div, by = "state_fips")
  
  donor_flags <- derive_city_donor_flags(panel_df)
  
  covs <- df_geo %>%
    transmute(
      place_fips,
      ln_pop90       = log(pmax(pop90_city, 1)),
      poverty_rate_90, unemploy_90, black_city_90,
      division = as.integer(division),
      region   = as.integer(region)
    ) %>%
    left_join(donor_flags, by="place_fips") %>%
    mutate(
      near_miss_ctrl = replace_na(near_miss_ctrl, 0L),
      future_ctrl    = replace_na(future_ctrl,    0L)
    )
  
  # mean_pre and pre_slope (on pre_years)
  pre_path <- panel_df %>%
    filter(year %in% pre_years) %>%
    select(place_fips, year, !!sym(y)) %>%
    group_by(place_fips) %>%
    summarize(
      mean_pre  = mean(.data[[y]], na.rm = TRUE),
      pre_slope = {
        yy <- .data[[y]]; tt <- year
        if (sum(is.finite(yy)) >= 2) {
          as.numeric(coef(lm(yy ~ tt))[2])
        } else NA_real_
      },
      .groups = "drop"
    )
  
  # selected pre-lag levels (wide)
  pre_lag_wide <- panel_df %>%
    filter(year %in% pre_lags) %>%
    select(place_fips, year, !!sym(y)) %>%
    mutate(var = paste0("y_", year)) %>%
    select(-year) %>%
    tidyr::pivot_wider(names_from = var, values_from = !!sym(y))
  
  trt <- panel_df %>%
    group_by(place_fips) %>%
    summarise(treated_city = as.integer(any(treated_dummy == 1L)), .groups = "drop")
  
  covs %>%
    left_join(pre_path,    by = "place_fips") %>%
    left_join(pre_lag_wide,by = "place_fips") %>%
    left_join(trt,         by = "place_fips") %>%
    mutate(treated_city = replace_na(treated_city, 0L))
}

# ---------- ATT & ES with city-specific linear trends --------------------------
est_att_w_trend <- function(df, y, wcol = "w_a4r"){
  d <- df %>% filter(!is.na(.data[[y]]), .data[[wcol]] > 0) %>% mutate(trend = year)
  feols(
    as.formula(paste0(y, " ~ post2:treated_dummy + i(place_fips, trend) | place_fips + year")),
    data = d, weights = as.formula(paste0("~", wcol)), cluster = ~ place_fips
  )
}

est_es_w_trend <- function(df, y, wcol = "w_a4r", ref = -1){
  d <- df %>% filter(!is.na(.data[[y]]), .data[[wcol]] > 0) %>% mutate(trend = year)
  feols(
    as.formula(paste0(y, " ~ sunab(g, year, ref.p = ", ref, ") + i(place_fips, trend) | place_fips + year")),
    data = d, weights = as.formula(paste0("~", wcol)), cluster = ~ place_fips
  )
}

# ---------- Plot helpers for ES (same style as Script 6) -----------------------
es_tidy <- function(es_model) {
  ct <- as.data.frame(fixest::coeftable(es_model)); ct$term <- rownames(ct)
  keep <- grepl("^year::", ct$term)
  tibble::tibble(
    rel  = as.integer(sub("^year::", "", ct$term[keep])),
    beta = ct$Estimate[keep],
    se   = ct$`Std. Error`[keep]
  ) %>%
    dplyr::mutate(lo = beta - 1.96 * se, hi = beta + 1.96 * se) %>%
    dplyr::arrange(rel)
}

plot_es <- function(es_model, title, ylab, file_stub){
  df <- es_tidy(es_model)
  p <- ggplot(df, aes(x = rel, y = beta)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
    geom_line() + geom_point(size = 1.6) +
    labs(title = title, x = "Event time (years from treatment; ref = -1)", y = ylab)
  ggsave(file.path(out_figs, paste0(file_stub, ".png")), plot = p, width = 7, height = 4.5, dpi = 300)
  p
}

# ---------- Runner: build targets + EB; attach; estimate; save; print ---------
run_pretrend_repair <- function(panel_df, panel_name, outcome,
                                restrict = "all",
                                pre_years = .PRE_YEARS, pre_lags = .PRE_LAGS){
  
  cat("\n\n=== PRE-TREND REPAIR — ", panel_name, " | ", outcome,
      " — targets = mean_pre + pre_slope + {", paste(pre_lags, collapse=", "),
      "} ; donors = ", restrict, " ===\n", sep="")
  
  # 1) City table w/ pre-path targets
  city_tab <- build_city_table_pretrend(panel_df, y = outcome,
                                        pre_years = pre_years, pre_lags = pre_lags)
  
  # ===== DIAGNOSTIC: availability of EB targets among controls =====
  lag_names <- paste0("y_", pre_lags)  # pre_lags is defined earlier in the function
  vars_all  <- c(
    "mean_pre", "pre_slope",
    lag_names,
    "ln_pop90", "poverty_rate_90", "unemploy_90", "black_city_90"
  )
  
  cat("\n[DIAG] Non-missing counts among CONTROLS for EB targets:\n")
  print(
    city_tab %>%
      dplyr::filter(treated_city == 0L) %>%
      dplyr::summarise(dplyr::across(all_of(vars_all), ~sum(is.finite(.)))) %>%
      tidyr::pivot_longer(everything(), names_to = "var", values_to = "n_ctrl") %>%
      dplyr::arrange(var)
  )
  cat("\n")
  # ===== end diagnostic =====
  
  # EB variables (order matters only for printing)
  lag_names <- paste0("y_", pre_lags)
  vars <- c("mean_pre", "pre_slope", lag_names,
            "ln_pop90","poverty_rate_90","unemploy_90","black_city_90")
  
  # quick counts
  print(city_tab %>% summarise(treated_cities = sum(treated_city==1L),
                               donor_cities   = sum(treated_city==0L)))
  
  # 2) Donor pool (A4R with optional restriction)
  donors_union <- donor_pool_union_restricted(city_tab, vars,
                                              min_div = 20, min_reg = 40,
                                              restrict = restrict)
  
  # 3) EB with fallback
  wtab <- eb_fit_with_fallback(city_tab, vars, donors_union)
  
  # 4) Diagnostics: std diffs before/after & weight diag
  bal <- tibble(var = vars) %>%
    mutate(
      sdiff_unw = map_dbl(var, ~ sdiff(wtab[[.x]], wtab$treated_city, NULL)),
      sdiff_w   = map_dbl(var, ~ sdiff(wtab[[.x]], wtab$treated_city,
                                       w = ifelse(wtab$treated_city==1, 1, wtab$w_eb)))
    )
  wd  <- weight_diag(wtab$w_eb[wtab$treated_city==0])
  
  cat("\n-- Donor pool size kept (after A4R) --\n")
  print(tibble(donors_kept = sum(wtab$treated_city==0 & wtab$w_eb>0),
               total_controls = sum(wtab$treated_city==0)))
  
  cat("\n-- Std. diffs: before vs after weighting --\n")
  print(bal %>% mutate(across(starts_with("sdiff"), \(x) round(x, 3))))
  
  cat("\n-- Weight diagnostics (controls) --\n")
  print(wd %>% mutate(across(where(is.numeric), \(x) round(x, 3))))
  
  # Save balance CSV
  bal_csv <- file.path(out_tables, paste0("07_PRETREND_balance_", tolower(panel_name), "_", outcome, ".csv"))
  readr::write_csv(bal, bal_csv)
  cat("\nSaved balance diagnostics to:\n  - ", bal_csv, "\n", sep="")
  
  # 5) Attach weights to panel and run models with city linear trends
  df_w <- panel_df %>%
    left_join(wtab %>% select(place_fips, w_eb), by="place_fips") %>%
    mutate(w_a4r = if_else(treated_dummy==1L, 1, replace_na(w_eb, 0)))
  
  # ATT with trends
  m_att <- est_att_w_trend(df_w, outcome, wcol = "w_a4r")
  # SUNAB with trends
  m_es  <- est_es_w_trend(df_w, outcome, wcol = "w_a4r", ref = -1)
  
  # Wald joint test of pre-treatment leads (-5..-1) = 0
  ct_names <- rownames(coeftable(m_es))
  pre_terms <- ct_names[grepl("^year::-(5|4|3|2|1)$", ct_names)]
  wald_res  <- if (length(pre_terms)) fixest::wald(m_es, pre_terms) else NULL
  
  # Save models & tiny ATT row
  stub <- paste0("07_PRETREND_", tolower(panel_name), "_", outcome)
  saveRDS(m_att, file = file.path(out_models, paste0(stub, "*ATT_wTrend.rds")))
  saveRDS(m_es,  file = file.path(out_models, paste0(stub, "*ES_wTrend.rds")))
  if (!is.null(wald_res)) saveRDS(wald_res, file = file.path(out_models, paste0(stub, "_ES_preWALD.rds")))
  
  att_row <- {
    ct <- fixest::coeftable(m_att)
    if ("post2:treated_dummy" %in% rownames(ct)) {
      tibble(estimate = unname(ct["post2:treated_dummy",1]),
             se       = unname(ct["post2:treated_dummy",2]),
             p_value  = unname(ct["post2:treated_dummy",4]))
    } else tibble(estimate = NA_real_, se = NA_real_, p_value = NA_real_)
  }
  att_csv <- file.path(out_tables, paste0("07_PRETREND_ATTrow_", tolower(panel_name), "_", outcome, ".csv"))
  readr::write_csv(att_row, att_csv)
  
  cat("\n-- Weighted TWFE ATT with city trends --\n"); print(summary(m_att))
  if (!is.null(wald_res)) { cat("\nJoint test (pre-leads = 0):\n"); print(wald_res) }
  
  # Plot & save ES
  p <- plot_es(m_es,
               title = paste0(panel_name, " — Event study (weighted; pre-trend adjusted)"),
               ylab  = paste0("Effect on ", outcome),
               file_stub = stub)
  cat("\nSaved ES plot to:\n  - ", file.path(out_figs, paste0(stub, ".png")), "\n", sep="")
  
  invisible(list(weights = wtab %>% select(place_fips, w_eb),
                 balance = bal, wdiag = wd,
                 att = m_att, es = m_es, wald = wald_res))
}

# ===================== RUN PRE-TREND REPAIR FOR ALL OUTCOMES ===================

PANELS <- list(
  Mayor   = list(df = pm, outcomes = c("dem_share","turnout_pct")),
  Council = list(df = pc, outcomes = c("dem_share","votes_per_seat_pct_vap"))
)

pret_results <- list()
for (pname in names(PANELS)){
  df_panel <- PANELS[[pname]]$df
  for (y in PANELS[[pname]]$outcomes){
    pret_results[[paste(pname,y,sep="::")]] <-
      run_pretrend_repair(df_panel, panel_name = pname, outcome = y, restrict = "all")
  }
}

cat("\nDone: PRE-TREND robustness finished for all four outcomes.\n")



###############################################################################
# Robustness A-1b: City Trends ONLY 

###############################################################################
# --- CITY-TREND-ONLY spec (A-1a) ---------------------------------------------
run_citytrend_only <- function(panel_df, panel_name, outcome,
                               restrict = "all",
                               target_vars = c("mean_pre",
                                               "ln_pop90","poverty_rate_90",
                                               "unemploy_90","black_city_90"),
                               pre_years = 1989:1993){
  
  cat("\n\n=== CITY-TREND ONLY — ", panel_name, " | ", outcome,
      " — targets = {", paste(target_vars, collapse=", "),
      "} ; donors = ", restrict, " ===\n", sep="")
  
  # 1) City table with mean_pre and covariates (no pre-lag levels)
  city_tab <- build_city_table_pretrend(panel_df, y = outcome,
                                        pre_years = pre_years,
                                        pre_lags  = integer(0))
  
  # 2) Donor pool & EB weights with fallback using simple targets
  donors_union <- donor_pool_union_restricted(city_tab, target_vars,
                                              min_div = 20, min_reg = 40,
                                              restrict = restrict)
  wtab <- eb_fit_with_fallback(city_tab, target_vars, donors_union)
  
  # Diagnostics
  bal <- tibble::tibble(var = target_vars) |>
    dplyr::mutate(
      sdiff_unw = purrr::map_dbl(var, ~ sdiff(wtab[[.x]], wtab$treated_city, NULL)),
      sdiff_w   = purrr::map_dbl(var, ~ sdiff(wtab[[.x]], wtab$treated_city,
                                              w = ifelse(wtab$treated_city==1, 1, wtab$w_eb)))
    )
  wd <- weight_diag(wtab$w_eb[wtab$treated_city==0])
  print(dplyr::mutate(bal, across(starts_with("sdiff"), ~round(.,3))))
  print(dplyr::mutate(wd, across(where(is.numeric), ~round(.,3))))
  
  # 3) Attach weights and estimate ATT/ES with city linear trends
  df_w <- panel_df |>
    dplyr::left_join(wtab |> dplyr::select(place_fips, w_eb), by="place_fips") |>
    dplyr::mutate(w_a4r = dplyr::if_else(treated_dummy==1L, 1, tidyr::replace_na(w_eb, 0)))
  
  m_att <- est_att_w_trend(df_w, outcome, wcol = "w_a4r")  # has i(place_fips, trend)
  m_es  <- est_es_w_trend(df_w, outcome, wcol = "w_a4r", ref = -1)
  
  # Save
  stub <- paste0("07A_CITYTREND_", tolower(panel_name), "_", outcome)
  readr::write_csv(bal, file.path(out_tables, paste0(stub, "_balance.csv")))
  saveRDS(m_att, file = file.path(out_models, paste0(stub, "_ATT.rds")))
  saveRDS(m_es,  file = file.path(out_models, paste0(stub, "_ES.rds")))
  
  # Wald test of pre-leads (full window; you can trim later)
  ct_names <- rownames(fixest::coeftable(m_es))
  pre_terms <- ct_names[grepl("^year::-(5|4|3|2|1)$", ct_names)]
  wald_res  <- if (length(pre_terms)) fixest::wald(m_es, pre_terms) else NULL
  if (!is.null(wald_res)) print(wald_res)
  
  # Plot ES (same helper you already have)
  p <- plot_es(m_es,
               title = paste0(panel_name, " — Event study (weighted; city-trend adj)"),
               ylab  = paste0("Effect on ", outcome),
               file_stub = stub)
  print(p)
  
  invisible(list(weights = wtab, balance = bal, wdiag = wd, att = m_att, es = m_es, wald = wald_res))
}

run_citytrend_only(pm, "Mayor", "dem_share", restrict="all")
run_citytrend_only(pm, "Mayor", "turnout_pct", restrict="all")
run_citytrend_only(pc, "Council", "dem_share", restrict="all")
run_citytrend_only(pc, "Council", "votes_per_seat_pct_vap", restrict="all")


###############################################################################
# Robustness G: Shocks 
# Division×Year FE (with city trends) 
# + 2-way clustering + Division Consistent 
###############################################################################
# --- ONE STEP: Division×Year FE (with city trends) + 2-way clustering ---

#  Ensure these helpers exist once in your session
make_shock_keys <- function(df){
  df %>%
    mutate(
      # state_fips already exists in your build_* helpers; if not, add substr(place_fips,1,2)
      div_year   = interaction(division, year, drop = TRUE),   # Division × Year FE key
      state_year = interaction(state_fips, year, drop = TRUE)  # for clustering
    )
}

est_att_divyear <- function(df, outcome, wcol = "w_a4r"){
  d <- df %>%
    filter(is.finite(.data[[outcome]]), .data[[wcol]] > 0) %>%
    make_shock_keys() %>%
    mutate(trend = year)
  
  feols(
    as.formula(paste0(outcome, " ~ post2:treated_dummy + i(place_fips, trend) | place_fips + div_year")),
    data    = d,
    weights = as.formula(paste0("~", wcol)),
    cluster = ~ place_fips + state_year   # two-way clustering
  )
}

est_es_divyear <- function(df, outcome, wcol = "w_a4r", ref = -1){
  d <- df %>%
    filter(is.finite(.data[[outcome]]), .data[[wcol]] > 0) %>%
    make_shock_keys() %>%
    mutate(trend = year)
  
  feols(
    as.formula(paste0(outcome, " ~ sunab(g, year, ref.p = ", ref, ") + i(place_fips, trend) | place_fips + div_year")),
    data    = d,
    weights = as.formula(paste0("~", wcol)),
    cluster = ~ place_fips + state_year   # two-way clustering
  )
}

# Re-run and STORE the return object
res_mayor_dem <- run_citytrend_only(pm, "Mayor", "dem_share", restrict = "all")

# Extract the weights table (has columns: place_fips, w_eb)
w_mayor_dem <- res_mayor_dem$weights

# Attach weights
# 0) Bring division onto the panel (join via the state FIPS)
#    Assumes `state_div` exists with columns: state_fips, division, region
df_w <- pm %>%
  dplyr::mutate(state_fips = substr(place_fips, 1, 2)) %>%
  dplyr::left_join(state_div, by = "state_fips") %>%   # adds `division` (and `region` if present)
  dplyr::left_join(w_mayor_dem, by = "place_fips") %>%
  dplyr::mutate(
    w_a4r   = dplyr::if_else(treated_dummy == 1L, 1, tidyr::replace_na(w_eb, 0)),
    division = as.factor(division),                    # make sure it's a factor
    trend   = year
  )

# quick sanity check
df_w %>% dplyr::summarise(n_div = dplyr::n_distinct(division), share_na_div = mean(is.na(division)))

# 1) ATT with Division×Year FE + city trends + 2-way clustering
m_att_divyr <- est_att_divyear(df_w, outcome = "dem_share", wcol = "w_a4r")
cat("\n[ATT | Division×Year FE + City trends]\n"); print(summary(m_att_divyr))

# 2) ES (SUNAB) with Division×Year FE + city trends + 2-way clustering
m_es_divyr  <- est_es_divyear(df_w, outcome = "dem_share", wcol = "w_a4r", ref = -1)

# pre-leads joint test (-5..-1)
ctn <- rownames(coeftable(m_es_divyr))
pre_terms <- ctn[grepl("^year::-(5|4|3|2|1)$", ctn)]
if (length(pre_terms)) {
  cat("\n[Pre-leads joint test | Division×Year FE]\n")
  print(fixest::wald(m_es_divyr, pre_terms))
}

###############################################################################
## ================== DIVISION-CONSISTENT WEIGHTS (City-trend targets) ==================
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(stringr)
})

# ---- utilities (safe to redeclare) ----
.wvar_unbiased <- function(x, w = NULL){
  x <- as.numeric(x); if (is.null(w)) w <- rep(1, length(x))
  ok <- is.finite(x) & is.finite(w) & w >= 0; x <- x[ok]; w <- w[ok]
  if (length(x) < 2 || sum(w) <= 0) return(NA_real_)
  w <- w / sum(w); mu <- sum(w * x)
  eff_n <- (sum(w)^2) / sum(w^2); if (!is.finite(eff_n) || eff_n <= 1) return(NA_real_)
  sum(w * (x - mu)^2) * (eff_n / (eff_n - 1))
}

sdiff <- function(x, t, w = NULL){
  t <- as.integer(t == 1L)
  if (is.null(w)) {
    m1 <- mean(x[t==1], na.rm=TRUE); v1 <- var(x[t==1], na.rm=TRUE)
    m0 <- mean(x[t==0], na.rm=TRUE); v0 <- var(x[t==0], na.rm=TRUE)
  } else {
    w <- as.numeric(w); w[!is.finite(w) | w < 0] <- 0
    mw <- mean(w, na.rm = TRUE); if (is.finite(mw) && mw > 0) w <- w / mw
    m1 <- tryCatch(weighted.mean(x[t==1], w[t==1], na.rm=TRUE), error = function(e) NA_real_)
    m0 <- tryCatch(weighted.mean(x[t==0], w[t==0], na.rm=TRUE), error = function(e) NA_real_)
    v1 <- .wvar_unbiased(x[t==1], w[t==1]); v0 <- .wvar_unbiased(x[t==0], w[t==0])
  }
  s <- sqrt((v1 + v0)/2); if (!is.finite(s) || s == 0) return(NA_real_)
  (m1 - m0) / s
}

weight_diag <- function(w_ctrl){
  w <- w_ctrl[is.finite(w_ctrl) & w_ctrl > 0]
  if (!length(w)) return(tibble(donors_kept=0, sum_w=0, ess=0, ess_share=0,
                                max_w=NA_real_, p99_w=NA_real_, p95_w=NA_real_))
  ess <- (sum(w)^2) / sum(w^2)
  tibble(donors_kept=length(w), sum_w=sum(w), ess=ess, ess_share=ess/length(w),
         max_w=max(w), p99_w=quantile(w, .99, names=FALSE), p95_w=quantile(w, .95, names=FALSE))
}

# ---- build city table for CITY-TREND targets (mean_pre + 1990 covariates + division) ----
build_city_table_citytrend <- function(panel_df, y, pre_years = 1989:1993){
  if (!all(c("place_fips","year","treated_dummy") %in% names(panel_df)))
    stop("Panel must have place_fips, year, treated_dummy.")
  if (!exists("state_div")) stop("state_div (state_fips, division) not found in env.")
  df_geo <- panel_df %>%
    distinct(place_fips, .keep_all = TRUE) %>%
    mutate(state_fips = substr(place_fips, 1, 2)) %>%
    left_join(state_div, by = "state_fips")
  pre_path <- panel_df %>%
    filter(year %in% pre_years) %>%
    select(place_fips, year, !!sym(y)) %>%
    group_by(place_fips) %>%
    summarise(mean_pre = mean(.data[[y]], na.rm = TRUE), .groups = "drop")
  trt <- panel_df %>%
    group_by(place_fips) %>%
    summarise(treated_city = as.integer(any(treated_dummy == 1L)), .groups = "drop")
  df_geo %>%
    transmute(place_fips,
              division = as.integer(division),
              ln_pop90 = log(pmax(pop90_city, 1)),
              poverty_rate_90, unemploy_90, black_city_90) %>%
    left_join(pre_path, by = "place_fips") %>%
    left_join(trt,      by = "place_fips") %>%
    mutate(treated_city = replace_na(treated_city, 0L))
}

# ---- EB within-division, with safe kernel fallback, normalized WITHIN each division ----
division_consistent_weights <- function(city_tab, vars){
  if (!"division" %in% names(city_tab)) stop("division not in city_tab.")
  div_list <- sort(unique(city_tab$division[city_tab$treated_city == 1L]))
  all_w <- list()
  
  for (dv in div_list){
    cat(sprintf("\n[DIV %s] -- EB on same-division donors\n", dv))
    
    # Treated in this division (targets)
    Tset <- city_tab %>%
      dplyr::filter(treated_city == 1L, division == dv) %>%
      dplyr::select(dplyr::all_of(vars))
    
    # Controls in this division (complete cases only)
    Ckeep <- city_tab %>%
      dplyr::filter(treated_city == 0L, division == dv) %>%
      dplyr::select(place_fips, dplyr::all_of(vars)) %>%
      tidyr::drop_na()
    
    if (nrow(Ckeep) == 0) {
      warning(sprintf("[DIV %s] No complete-case donors; skipping (weights=0).", dv))
      next
    }
    
    target <- colMeans(as.matrix(Tset), na.rm = TRUE)
    Xc     <- as.matrix(Ckeep[, vars, drop = FALSE])
    
    # Try EB first
    eb_try <- try(ebal::ebalance(target.margins = target, X = Xc,
                                 base.weights = rep(1, nrow(Xc)),
                                 maxit = 10000, constraint.tolerance = 1e-8),
                  silent = TRUE)
    
    if (!inherits(eb_try, "try-error") && !is.null(eb_try$w)) {
      w <- as.numeric(eb_try$w)
      msg <- "[DIV %s] EB exact converged."
    } else {
      # Safe kernel fallback (robust standardization)
      center <- colMeans(Xc, na.rm = TRUE)
      scalev <- apply(Xc, 2, sd); scalev[!is.finite(scalev) | scalev < 1e-8] <- 1
      xc_s  <- sweep(Xc, 2, center, "-"); xc_s <- sweep(xc_s, 2, scalev, "/")
      tgt_s <- (target - center) / scalev
      diff_mat <- xc_s - matrix(tgt_s, nrow = nrow(xc_s), ncol = ncol(xc_s), byrow = TRUE)
      d2 <- rowSums(diff_mat^2)
      h  <- stats::quantile(d2[is.finite(d2)], 0.25, na.rm = TRUE)
      if (!is.finite(h) || h <= 1e-8) h <- max(stats::median(d2[is.finite(d2)], na.rm = TRUE), 1)
      w <- exp(- d2 / h); w[!is.finite(w)] <- 0
      msg <- sprintf("[DIV %s] EB failed; using kernel fallback (h=%.3f).", dv, h)
    }
    
    # Normalize within division
    m <- mean(w, na.rm = TRUE)
    if (!is.finite(m) || m <= 0) {
      warning(sprintf("[DIV %s] Non-finite mean weight; setting uniform.", dv))
      w <- rep(1, length(w)); m <- 1
    }
    w <- w / m
    
    # Weights aligned with the complete-case donor table
    cw <- tibble::tibble(place_fips = Ckeep$place_fips, w_eb = as.numeric(w))
    
    # Sub-division table for diagnostics and balance
    sub_tab <- city_tab %>%
      dplyr::filter(division == dv) %>%
      dplyr::left_join(cw, by = "place_fips") %>%
      dplyr::mutate(w_use = dplyr::if_else(treated_city == 1L, 1, tidyr::replace_na(w_eb, 0)))
    
    # Diagnostics (within this division)
    donors_kept <- sum(sub_tab$treated_city == 0L & sub_tab$w_use > 0)
    ess <- (sum(sub_tab$w_use[sub_tab$treated_city == 0L])^2) /
      sum(sub_tab$w_use[sub_tab$treated_city == 0L]^2)
    cat(sprintf(msg, dv), " donors_kept=", donors_kept,
        " ESS=", round(ess, 2), "/", sum(sub_tab$treated_city == 0L), "\n", sep = "")
    
    bal <- tibble::tibble(var = vars) %>%
      dplyr::mutate(
        sdiff_unw = purrr::map_dbl(var, ~ sdiff(sub_tab[[.x]], sub_tab$treated_city, NULL)),
        sdiff_w   = purrr::map_dbl(var, ~ sdiff(sub_tab[[.x]], sub_tab$treated_city, w = sub_tab$w_use))
      )
    print(dplyr::mutate(bal, dplyr::across(dplyr::starts_with("sdiff"), ~round(., 3))))
    
    all_w[[as.character(dv)]] <- cw %>% dplyr::mutate(division = dv)
  }
  
  if (!length(all_w)) {
    warning("No divisions produced weights; returning zeros.")
    return(city_tab %>% dplyr::transmute(place_fips, w_eb = 0))
  }
  dplyr::bind_rows(all_w) %>% dplyr::select(place_fips, w_eb) %>% dplyr::distinct()
}

# ---- runner: compute division-consistent weights for a panel/outcome ----
run_citytrend_divisionconsistent <- function(panel_df, panel_name, outcome){
  cat("\n=== CITY-TREND (Division-consistent weights) — ", panel_name, " | ", outcome, " ===\n", sep="")
  city_tab <- build_city_table_citytrend(panel_df, y = outcome, pre_years = 1989:1993)
  vars <- c("mean_pre","ln_pop90","poverty_rate_90","unemploy_90","black_city_90")
  # DIAG: treated divisions
  tdiv <- sort(unique(city_tab$division[city_tab$treated_city==1L]))
  cat("Treated divisions:", paste(tdiv, collapse=", "), "\n")
  # compute weights
  wtab <- division_consistent_weights(city_tab, vars)
  # return weights + a small summary
  list(
    weights = wtab,
    treated_divisions = tdiv
  )
}

# -------------------- USAGE EXAMPLE --------------------
res_div <- run_citytrend_divisionconsistent(pm, "Mayor", "dem_share")

# ================== RUN: Mayor | dem_share  (Div-consistent weights + Div×Year FE + city trends) ==================

suppressPackageStartupMessages({
  library(dplyr); library(fixest)
})

# 0) Ensure we have division-consistent weights (place_fips, w_eb)
if (!exists("res_div") || !"weights" %in% names(res_div)) {
  # assumes run_citytrend_divisionconsistent() is already defined in your session
  res_div <- run_citytrend_divisionconsistent(pm, "Mayor", "dem_share")
}
w_div <- res_div$weights %>% dplyr::select(place_fips, w_eb)

# 1) Prepare weighted panel with Division and state-year keys
#    Assumes `state_div` exists with columns: state_fips, division, (region optional)
df_w <- pm %>%
  dplyr::mutate(state_fips = substr(place_fips, 1, 2)) %>%
  dplyr::left_join(state_div, by = "state_fips") %>%
  dplyr::left_join(w_div, by = "place_fips") %>%
  dplyr::mutate(
    w_a4r     = dplyr::if_else(treated_dummy == 1L, 1, tidyr::replace_na(w_eb, 0)),
    division  = as.factor(division),
    div_year  = interaction(division, year, drop = TRUE),
    trend     = year,
    state_year = interaction(state_fips, year, drop = TRUE)
  )

# 2) ATT with Division×Year FE + city-specific linear trends + 2-way clustering
m_att_divyr <- feols(
  dem_share ~ post2:treated_dummy + i(place_fips, trend) | place_fips + div_year,
  data = df_w,
  weights = ~ w_a4r,
  cluster = ~ place_fips + state_year
)

# 3) ES (SUNAB) with Division×Year FE + city-specific linear trends + 2-way clustering
m_es_divyr <- feols(
  dem_share ~ sunab(g, year, ref.p = -1) + i(place_fips, trend) | place_fips + div_year,
  data = df_w,
  weights = ~ w_a4r,
  cluster = ~ place_fips + state_year
)

# 4) Compact printout for ATT
cat("\n================== ATT (Mayor | dem_share) — Div×Year FE + city trends ==================\n")
att_ct <- fixest::coeftable(m_att_divyr)
if ("post2:treated_dummy" %in% rownames(att_ct)) {
  att_est <- unname(att_ct["post2:treated_dummy","Estimate"])
  att_se  <- unname(att_ct["post2:treated_dummy","Std. Error"])
  att_p   <- unname(att_ct["post2:treated_dummy","Pr(>|t|)"])
  cat(sprintf("ATT (post×treated):  %.3f   [SE=%.3f]   p=%.3f\n", att_est, att_se, att_p))
} else {
  cat("ATT (post×treated) coefficient not found in model.\n")
}
cat(sprintf("Obs: %d | Cities: %d | Div×Year FE: %d\n",
            nobs(m_att_divyr),
            length(fixef(m_att_divyr)$place_fips),
            length(fixef(m_att_divyr)$div_year)))

# 5) Joint test of pre-treatment leads (-5..-1) in the ES
cat("\n================== Pre-leads joint test (ES) ==================\n")
ctn <- rownames(coeftable(m_es_divyr))
pre_terms <- ctn[grepl("^year::-(5|4|3|2|1):cohort::", ctn)]
if (length(pre_terms) > 0) {
  jw <- fixest::wald(m_es_divyr, pre_terms)
  cat(sprintf("Wald stat = %.4f | p = %.6f | df1 = %s | df2 = %s | VCOV = %s\n",
              jw$stat, jw$p, jw$df1, jw$df2, jw$vcov))
} else {
  cat("No pre-lead terms detected at lags -5..-1 (check your ES window / ref period).\n")
}


# 6) (Optional) Print a one-line event-time summary (middle few coefficients)
# You can comment this block out if you want it even leaner.
es_ct <- as.data.frame(coeftable(m_es_divyr))
es_ct$term <- rownames(es_ct)
et <- es_ct[grepl("^year::", es_ct$term), c("term","Estimate","Std. Error","Pr(>|t|)")]
if (nrow(et) > 0) {
  keep <- et[order(et$term), ]
  cat("\nEvent-time (subset):\n")
  print(utils::head(keep, 10), row.names = FALSE)
}

###############################################################################
# Mayor | turnout_pct
# Step A: City-trend-only robustness (EB over all untreated controls)
# Step B: Division-consistent weights + Division×Year FE + city trends
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(fixest); library(stringr)
})

# ---------- Small helpers ----------
.wvar_unbiased <- function(x, w = NULL){
  x <- as.numeric(x); if (is.null(w)) w <- rep(1, length(x))
  ok <- is.finite(x) & is.finite(w) & w >= 0; x <- x[ok]; w <- w[ok]
  if (length(x) < 2 || sum(w) <= 0) return(NA_real_)
  w <- w / sum(w); mu <- sum(w * x)
  eff_n <- (sum(w)^2) / sum(w^2); if (!is.finite(eff_n) || eff_n <= 1) return(NA_real_)
  sum(w * (x - mu)^2) * (eff_n / (eff_n - 1))
}
sdiff <- function(x, t, w = NULL){
  t <- as.integer(t == 1L)
  if (is.null(w)) {
    m1 <- mean(x[t==1], na.rm=TRUE); v1 <- var(x[t==1], na.rm=TRUE)
    m0 <- mean(x[t==0], na.rm=TRUE); v0 <- var(x[t==0], na.rm=TRUE)
  } else {
    w <- as.numeric(w); w[!is.finite(w) | w < 0] <- 0
    mw <- mean(w, na.rm = TRUE); if (is.finite(mw) && mw > 0) w <- w / mw
    m1 <- tryCatch(weighted.mean(x[t==1], w[t==1], na.rm=TRUE), error = function(e) NA_real_)
    m0 <- tryCatch(weighted.mean(x[t==0], w[t==0], na.rm=TRUE), error = function(e) NA_real_)
    v1 <- .wvar_unbiased(x[t==1], w[t==1]); v0 <- .wvar_unbiased(x[t==0], w[t==0])
  }
  s  <- sqrt((v1 + v0)/2); if (!is.finite(s) || s == 0) return(NA_real_)
  (m1 - m0) / s
}

# ---------- Builders ----------
# City-trend targets: mean_pre + (1990 covariates). Includes division for later joins.
build_city_table_citytrend <- function(panel_df, y, pre_years = 1989:1993){
  if (!all(c("place_fips","year","treated_dummy") %in% names(panel_df)))
    stop("Panel must have place_fips, year, treated_dummy.")
  if (!exists("state_div")) stop("state_div (state_fips, division) not found in env.")
  
  df_geo <- panel_df %>%
    distinct(place_fips, .keep_all = TRUE) %>%
    mutate(state_fips = substr(place_fips, 1, 2)) %>%
    left_join(state_div, by = "state_fips")
  
  pre_path <- panel_df %>%
    filter(year %in% pre_years) %>%
    select(place_fips, year, !!sym(y)) %>%
    group_by(place_fips) %>%
    summarise(mean_pre = mean(.data[[y]], na.rm = TRUE), .groups = "drop")
  
  trt <- panel_df %>%
    group_by(place_fips) %>%
    summarise(treated_city = as.integer(any(treated_dummy == 1L)), .groups = "drop")
  
  df_geo %>%
    transmute(place_fips,
              division = as.integer(division),
              ln_pop90 = log(pmax(pop90_city, 1)),
              poverty_rate_90, unemploy_90, black_city_90) %>%
    left_join(pre_path, by = "place_fips") %>%
    left_join(trt,      by = "place_fips") %>%
    mutate(treated_city = tidyr::replace_na(treated_city, 0L))
}

# ---------- A) CITY-TREND ONLY (all controls) ----------
run_citytrend_only <- function(panel_df, panel_name, outcome){
  cat("\n=== CITY-TREND ONLY — ", panel_name, " | ", outcome,
      " — targets = {mean_pre, ln_pop90, poverty_rate_90, unemploy_90, black_city_90} ; donors = all ===\n",
      sep="")
  city_tab <- build_city_table_citytrend(panel_df, y = outcome, pre_years = 1989:1993)
  vars <- c("mean_pre","ln_pop90","poverty_rate_90","unemploy_90","black_city_90")
  
  # Build EB over all controls; kernel fallback with sharp-ish bandwidth if EB fails
  Tset <- city_tab %>% filter(treated_city==1L) %>% select(dplyr::all_of(vars))
  Ckeep <- city_tab %>% filter(treated_city==0L) %>% select(place_fips, dplyr::all_of(vars)) %>% tidyr::drop_na()
  target <- colMeans(as.matrix(Tset), na.rm = TRUE)
  Xc     <- as.matrix(Ckeep[, vars, drop=FALSE])
  
  eb_try <- try(ebal::ebalance(target.margins = target, X = Xc,
                               base.weights = rep(1, nrow(Xc)),
                               maxit = 10000, constraint.tolerance = 1e-8),
                silent = TRUE)
  
  if (!inherits(eb_try, "try-error") && !is.null(eb_try$w)) {
    w <- as.numeric(eb_try$w)
    h_dbg <- NA_real_
  } else {
    center <- colMeans(Xc, na.rm = TRUE)
    scalev <- apply(Xc, 2, sd); scalev[!is.finite(scalev) | scalev < 1e-8] <- 1
    xc_s  <- sweep(Xc, 2, center, "-"); xc_s <- sweep(xc_s, 2, scalev, "/")
    tgt_s <- (target - center) / scalev
    diff_mat <- xc_s - matrix(tgt_s, nrow = nrow(xc_s), ncol = ncol(xc_s), byrow = TRUE)
    d2 <- rowSums(diff_mat^2)
    h_dbg <- stats::quantile(d2[is.finite(d2)], 0.25, na.rm = TRUE)
    if (!is.finite(h_dbg) || h_dbg <= 1e-8) h_dbg <- max(stats::median(d2[is.finite(d2)], na.rm = TRUE), 1)
    w <- exp(- d2 / h_dbg); w[!is.finite(w)] <- 0
  }
  m <- mean(w, na.rm=TRUE); if (!is.finite(m) || m <= 0) { w <- rep(1, length(w)); m <- 1 }
  w <- w / m
  ess <- (sum(w)^2) / sum(w^2)
  
  cat(sprintf("[KERN] h=%s  sum=%.6f  min=%.6g  max=%.6g  ESS=%.2f/%d\n",
              ifelse(is.na(h_dbg), "EB", format(h_dbg, digits=4)),
              sum(w), min(w), max(w), ess, length(w)))
  
  wtab <- tibble::tibble(place_fips = Ckeep$place_fips, w_eb = as.numeric(w))
  
  # balance table
  sub_tab <- city_tab %>% left_join(wtab, by="place_fips") %>%
    mutate(w_use = if_else(treated_city==1L, 1, tidyr::replace_na(w_eb, 0)))
  bal <- tibble::tibble(var = vars) %>%
    mutate(sdiff_unw = map_dbl(var, ~ sdiff(sub_tab[[.x]], sub_tab$treated_city, NULL)),
           sdiff_w   = map_dbl(var, ~ sdiff(sub_tab[[.x]], sub_tab$treated_city, w=sub_tab$w_use)))
  print(mutate(bal, across(starts_with("sdiff"), ~round(.,3))))
  
  # attach weights to panel and run ATT + ES with city trends
  df_w <- panel_df %>%
    left_join(wtab, by = "place_fips") %>%
    mutate(w_a4r = if_else(treated_dummy==1L, 1, tidyr::replace_na(w_eb, 0)),
           trend = year)
  
  m_att <- feols(
    as.formula(paste0(outcome, " ~ post2:treated_dummy + i(place_fips, trend) | place_fips + year")),
    data = df_w, weights = ~ w_a4r, cluster = ~ place_fips
  )
  cat("\n[ATT | City trends only]\n"); print(summary(m_att))
  
  m_es  <- feols(
    as.formula(paste0(outcome, " ~ sunab(g, year, ref.p = -1) + i(place_fips, trend) | place_fips + year")),
    data = df_w, weights = ~ w_a4r, cluster = ~ place_fips
  )
  
  # joint test pre-leads (-5..-1), flexible to label variants
  ctn <- rownames(coeftable(m_es))
  pre_terms <- ctn[grepl("^year::-(5|4|3|2|1)", ctn)]
  if (length(pre_terms)) {
    cat("\n[Pre-leads joint test | City trends only]\n"); print(fixest::wald(m_es, pre_terms))
  } else {
    cat("\n[Pre-leads joint test | City trends only] No pre-lead terms detected in labels.\n")
  }
  
  invisible(list(weights = wtab, att = m_att, es = m_es, balance = bal))
}

# ---------- B) Division-consistent weights ----------
division_consistent_weights <- function(city_tab, vars){
  if (!"division" %in% names(city_tab)) stop("division not in city_tab.")
  div_list <- sort(unique(city_tab$division[city_tab$treated_city == 1L]))
  all_w <- list()
  for (dv in div_list){
    cat(sprintf("\n[DIV %s] -- EB on same-division donors\n", dv))
    
    Tset <- city_tab %>%
      dplyr::filter(treated_city == 1L, division == dv) %>%
      dplyr::select(dplyr::all_of(vars))
    
    Ckeep <- city_tab %>%
      dplyr::filter(treated_city == 0L, division == dv) %>%
      dplyr::select(place_fips, dplyr::all_of(vars)) %>%
      tidyr::drop_na()
    
    if (nrow(Ckeep) == 0) {
      warning(sprintf("[DIV %s] No complete-case donors; skipping (weights=0).", dv))
      next
    }
    
    target <- colMeans(as.matrix(Tset), na.rm = TRUE)
    Xc     <- as.matrix(Ckeep[, vars, drop = FALSE])
    
    eb_try <- try(ebal::ebalance(target.margins = target, X = Xc,
                                 base.weights = rep(1, nrow(Xc)),
                                 maxit = 10000, constraint.tolerance = 1e-8),
                  silent = TRUE)
    
    if (!inherits(eb_try, "try-error") && !is.null(eb_try$w)) {
      w <- as.numeric(eb_try$w)
      msg <- sprintf("[DIV %s] EB exact converged.", dv)
    } else {
      center <- colMeans(Xc, na.rm = TRUE)
      scalev <- apply(Xc, 2, sd); scalev[!is.finite(scalev) | scalev < 1e-8] <- 1
      xc_s  <- sweep(Xc, 2, center, "-"); xc_s <- sweep(xc_s, 2, scalev, "/")
      tgt_s <- (target - center) / scalev
      diff_mat <- xc_s - matrix(tgt_s, nrow = nrow(xc_s), ncol = ncol(xc_s), byrow = TRUE)
      d2 <- rowSums(diff_mat^2)
      h  <- stats::quantile(d2[is.finite(d2)], 0.25, na.rm = TRUE)
      if (!is.finite(h) || h <= 1e-8) h <- max(stats::median(d2[is.finite(d2)], na.rm = TRUE), 1)
      w <- exp(- d2 / h); w[!is.finite(w)] <- 0
      msg <- sprintf("[DIV %s] EB failed; using kernel fallback (h=%.3f).", dv, h)
    }
    
    m <- mean(w, na.rm = TRUE)
    if (!is.finite(m) || m <= 0) { warning(sprintf("[DIV %s] Non-finite mean; uniform.", dv)); w <- rep(1, length(w)); m <- 1 }
    w <- w / m
    
    cw <- tibble::tibble(place_fips = Ckeep$place_fips, w_eb = as.numeric(w))
    
    sub_tab <- city_tab %>%
      dplyr::filter(division == dv) %>%
      dplyr::left_join(cw, by = "place_fips") %>%
      dplyr::mutate(w_use = dplyr::if_else(treated_city == 1L, 1, tidyr::replace_na(w_eb, 0)))
    
    donors_kept <- sum(sub_tab$treated_city == 0L & sub_tab$w_use > 0)
    ess <- (sum(sub_tab$w_use[sub_tab$treated_city == 0L])^2) /
      sum(sub_tab$w_use[sub_tab$treated_city == 0L]^2)
    cat(msg, " donors_kept=", donors_kept, " ESS=", round(ess,2), "/", sum(sub_tab$treated_city == 0L), "\n", sep = "")
    
    vars_chk <- vars
    bal <- tibble::tibble(var = vars_chk) %>%
      dplyr::mutate(
        sdiff_unw = purrr::map_dbl(var, ~ sdiff(sub_tab[[.x]], sub_tab$treated_city, NULL)),
        sdiff_w   = purrr::map_dbl(var, ~ sdiff(sub_tab[[.x]], sub_tab$treated_city, w = sub_tab$w_use))
      )
    print(dplyr::mutate(bal, dplyr::across(dplyr::starts_with("sdiff"), ~round(., 3))))
    
    all_w[[as.character(dv)]] <- cw %>% dplyr::mutate(division = dv)
  }
  if (!length(all_w)) {
    warning("No divisions produced weights; returning zeros.")
    return(city_tab %>% dplyr::transmute(place_fips, w_eb = 0))
  }
  dplyr::bind_rows(all_w) %>% dplyr::select(place_fips, w_eb) %>% dplyr::distinct()
}

run_citytrend_divisionconsistent <- function(panel_df, panel_name, outcome){
  cat("\n=== CITY-TREND (Division-consistent weights) — ", panel_name, " | ", outcome, " ===\n", sep="")
  city_tab <- build_city_table_citytrend(panel_df, y = outcome, pre_years = 1989:1993)
  vars <- c("mean_pre","ln_pop90","poverty_rate_90","unemploy_90","black_city_90")
  tdiv <- sort(unique(city_tab$division[city_tab$treated_city==1L]))
  cat("Treated divisions:", paste(tdiv, collapse=", "), "\n")
  wtab <- division_consistent_weights(city_tab, vars)
  list(weights = wtab, treated_divisions = tdiv)
}

# ========================= RUNS FOR Mayor | turnout_pct =========================

# A) City-trend-only
res_turnout_citytrend <- run_citytrend_only(pm, "Mayor", "turnout_pct")

# B) Division-consistent weights + Division×Year FE + city trends + 2-way clustering
res_div_turnout <- run_citytrend_divisionconsistent(pm, "Mayor", "turnout_pct")
w_div <- res_div_turnout$weights %>% dplyr::select(place_fips, w_eb)

# Prepare weighted panel with division & state-year keys
df_w_turnout <- pm %>%
  dplyr::mutate(state_fips = substr(place_fips, 1, 2)) %>%
  dplyr::left_join(state_div, by = "state_fips") %>%
  dplyr::left_join(w_div, by = "place_fips") %>%
  dplyr::mutate(
    w_a4r     = dplyr::if_else(treated_dummy == 1L, 1, tidyr::replace_na(w_eb, 0)),
    division  = as.factor(division),
    div_year  = interaction(division, year, drop = TRUE),
    trend     = year,
    state_year = interaction(state_fips, year, drop = TRUE)
  )

# ATT with Div×Year FE + city trends + 2-way clustering
m_att_divyr_turnout <- feols(
  turnout_pct ~ post2:treated_dummy + i(place_fips, trend) | place_fips + div_year,
  data = df_w_turnout,
  weights = ~ w_a4r,
  cluster = ~ place_fips + state_year
)

cat("\n================== ATT (Mayor | turnout_pct) — Div×Year FE + city trends ==================\n")
att_ct <- fixest::coeftable(m_att_divyr_turnout)
if ("post2:treated_dummy" %in% rownames(att_ct)) {
  cat(sprintf("ATT (post×treated):  %.3f   [SE=%.3f]   p=%.3f\n",
              unname(att_ct["post2:treated_dummy","Estimate"]),
              unname(att_ct["post2:treated_dummy","Std. Error"]),
              unname(att_ct["post2:treated_dummy","Pr(>|t|)"])))
} else {
  cat("ATT (post×treated) coefficient not found in model.\n")
}
cat(sprintf("Obs: %d | Cities: %d | Div×Year FE: %d\n",
            nobs(m_att_divyr_turnout),
            length(fixef(m_att_divyr_turnout)$place_fips),
            length(fixef(m_att_divyr_turnout)$div_year)))

# ES with Div×Year FE + city trends + 2-way clustering
m_es_divyr_turnout <- feols(
  turnout_pct ~ sunab(g, year, ref.p = -1) + i(place_fips, trend) | place_fips + div_year,
  data = df_w_turnout,
  weights = ~ w_a4r,
  cluster = ~ place_fips + state_year
)

cat("\n================== Pre-leads joint test (ES) — turnout_pct ==================\n")
ctn <- rownames(coeftable(m_es_divyr_turnout))
pre_terms <- ctn[grepl("^year::-(5|4|3|2|1)", ctn)]
if (length(pre_terms) > 0) {
  jw <- fixest::wald(m_es_divyr_turnout, pre_terms)
  cat(sprintf("Wald stat = %.4f | p = %.6f | df1 = %s | df2 = %s | VCOV = %s\n",
              jw$stat, jw$p, jw$df1, jw$df2, jw$vcov))
} else {
  cat("No pre-lead terms detected at lags -5..-1 (check ES window / ref period).\n")
}

# (Optional) quick peek at first few event-time coefficients
es_ct <- as.data.frame(coeftable(m_es_divyr_turnout)); es_ct$term <- rownames(es_ct)
et <- es_ct[grepl("^year::", es_ct$term), c("term","Estimate","Std. Error","Pr(>|t|)")]
if (nrow(et) > 0) {
  keep <- et[order(et$term), ]
  cat("\nEvent-time (subset):\n")
  print(utils::head(keep, 10), row.names = FALSE)
}

###############################################################################
# Council | dem_share
# Step A: City-trend-only robustness (EB over all untreated controls)
# Step B: Division-consistent weights + Division×Year FE + city trends
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(fixest); library(stringr)
})

# ---------- Step A: City-trend-only ----------
res_council_citytrend <- run_citytrend_only(pc, "Council", "dem_share")
# res_council_citytrend$weights  # (if you want the A-step weights/table)
# res_council_citytrend$att
# res_council_citytrend$es

# ---------- Step B: Division-consistent weights ----------
res_div_council <- run_citytrend_divisionconsistent(pc, "Council", "dem_share")
w_div <- res_div_council$weights %>% dplyr::select(place_fips, w_eb)

# Prepare weighted panel with Division and state-year keys
df_w_council <- pc %>%
  dplyr::mutate(state_fips = substr(place_fips, 1, 2)) %>%
  dplyr::left_join(state_div, by = "state_fips") %>%
  dplyr::left_join(w_div, by = "place_fips") %>%
  dplyr::mutate(
    w_a4r      = dplyr::if_else(treated_dummy == 1L, 1, tidyr::replace_na(w_eb, 0)),
    division   = as.factor(division),
    div_year   = interaction(division, year, drop = TRUE),
    trend      = year,
    state_year = interaction(state_fips, year, drop = TRUE)
  )

# ATT with Div×Year FE + city trends + 2-way clustering
m_att_divyr_council <- feols(
  dem_share ~ post2:treated_dummy + i(place_fips, trend) | place_fips + div_year,
  data    = df_w_council,
  weights = ~ w_a4r,
  cluster = ~ place_fips + state_year
)

cat("\n================== ATT (Council | dem_share) — Div×Year FE + city trends ==================\n")
att_ct <- fixest::coeftable(m_att_divyr_council)
if ("post2:treated_dummy" %in% rownames(att_ct)) {
  cat(sprintf("ATT (post×treated):  %.3f   [SE=%.3f]   p=%.3f\n",
              unname(att_ct["post2:treated_dummy","Estimate"]),
              unname(att_ct["post2:treated_dummy","Std. Error"]),
              unname(att_ct["post2:treated_dummy","Pr(>|t|)"])))
} else {
  cat("ATT (post×treated) coefficient not found in model.\n")
}
cat(sprintf("Obs: %d | Cities: %d | Div×Year FE: %d\n",
            nobs(m_att_divyr_council),
            length(fixef(m_att_divyr_council)$place_fips),
            length(fixef(m_att_divyr_council)$div_year)))

# ES (SUNAB) with Div×Year FE + city trends + 2-way clustering
m_es_divyr_council <- feols(
  dem_share ~ sunab(g, year, ref.p = -1) + i(place_fips, trend) | place_fips + div_year,
  data    = df_w_council,
  weights = ~ w_a4r,
  cluster = ~ place_fips + state_year
)

cat("\n================== Pre-leads joint test (ES) — Council dem_share ==================\n")
ctn <- rownames(coeftable(m_es_divyr_council))
pre_terms <- ctn[grepl("^year::-(5|4|3|2|1)", ctn)]
if (length(pre_terms) > 0) {
  jw <- fixest::wald(m_es_divyr_council, pre_terms)
  cat(sprintf("Wald stat = %.4f | p = %.6f | df1 = %s | df2 = %s | VCOV = %s\n",
              jw$stat, jw$p, jw$df1, jw$df2, jw$vcov))
} else {
  cat("No pre-lead terms detected at lags -5..-1 (check ES window / ref period).\n")
}

# (Optional) peek at a few event-time coefficients
es_ct <- as.data.frame(coeftable(m_es_divyr_council)); es_ct$term <- rownames(es_ct)
et <- es_ct[grepl("^year::", es_ct$term), c("term","Estimate","Std. Error","Pr(>|t|)")]
if (nrow(et) > 0) {
  keep <- et[order(et$term), ]
  cat("\nEvent-time (subset):\n")
  print(utils::head(keep, 10), row.names = FALSE)
}

  # ===== G-shock (Council | dem_share): national EB weights + Div×Year FE + city trends =====
suppressPackageStartupMessages({ library(dplyr); library(fixest) })

# 0) Get national EB weights from your city-trend-only step (re-use if already computed)
if (!exists("res_council_citytrend") || !"weights" %in% names(res_council_citytrend)) {
  res_council_citytrend <- run_citytrend_only(pc, "Council", "dem_share")  # defined earlier
}
w_nat <- res_council_citytrend$weights %>% dplyr::select(place_fips, w_eb)

# 1) Join weights + add Division×Year and state×year keys
df_w <- pc %>%
  dplyr::mutate(state_fips = substr(place_fips, 1, 2)) %>%
  dplyr::left_join(state_div, by = "state_fips") %>%     # needs state_div: (state_fips, division)
  dplyr::left_join(w_nat, by = "place_fips") %>%
  dplyr::mutate(
    w_a4r      = dplyr::if_else(treated_dummy == 1L, 1, tidyr::replace_na(w_eb, 0)),
    division   = as.factor(division),
    div_year   = interaction(division, year, drop = TRUE),
    trend      = year,
    state_year = interaction(state_fips, year, drop = TRUE)
  )

# 2) ATT with Div×Year FE + city-specific linear trends + 2-way clustering
m_att <- feols(
  dem_share ~ post2:treated_dummy + i(place_fips, trend) | place_fips + div_year,
  data = df_w, weights = ~ w_a4r, cluster = ~ place_fips + state_year
)

cat("\n================== ATT (Council | dem_share) — national EB + Div×Year FE + city trends ==================\n")
ct <- fixest::coeftable(m_att)
if ("post2:treated_dummy" %in% rownames(ct)) {
  cat(sprintf("ATT (post×treated):  %.3f   [SE=%.3f]   p=%.3f\n",
              unname(ct["post2:treated_dummy","Estimate"]),
              unname(ct["post2:treated_dummy","Std. Error"]),
              unname(ct["post2:treated_dummy","Pr(>|t|)"])))
} else {
  cat("ATT (post×treated) coefficient not found in model.\n")
}
cat(sprintf("Obs: %d | Cities: %d | Div×Year FE: %d\n",
            nobs(m_att), length(fixef(m_att)$place_fips), length(fixef(m_att)$div_year)))

# 3) Event-study (SUNAB) with same FE, trends, weights, and clustering
m_es <- feols(
  dem_share ~ sunab(g, year, ref.p = -1) + i(place_fips, trend) | place_fips + div_year,
  data = df_w, weights = ~ w_a4r, cluster = ~ place_fips + state_year
)

# 4) Joint pre-lead test (-5..-1)
cat("\n================== Pre-leads joint test (ES) — Council dem_share ==================\n")
ctn <- rownames(coeftable(m_es))
pre_terms <- ctn[grepl("^year::-(5|4|3|2|1):cohort::", ctn)]
if (length(pre_terms) > 0) {
  jw <- fixest::wald(m_es, pre_terms)
  cat(sprintf("Wald stat = %.4f | p = %.6f | df1 = %s | df2 = %s | VCOV = %s\n",
              jw$stat, jw$p, jw$df1, jw$df2, jw$vcov))
} else {
  cat("No pre-lead terms detected at lags -5..-1 (check ES window / ref period).\n")
}

###############################################################################
# Council | votes_per_seat_pct_vap  (turnout proxy)
# A) City-trend-only (pretrend repair check, EB over all controls)
# B) Division-consistent weights + Division×Year FE + city trends (+ 2-way cluster)
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(fixest)
})

# ---------- A) City-trend-only ----------
res_council_vps_citytrend <- run_citytrend_only(pc, "Council", "votes_per_seat_pct_vap")
# You can inspect: res_council_vps_citytrend$weights / $att / $es

# ---------- B) Division-consistent weights ----------
res_div_council_vps <- run_citytrend_divisionconsistent(pc, "Council", "votes_per_seat_pct_vap")
w_div <- res_div_council_vps$weights %>% dplyr::select(place_fips, w_eb)

# Build weighted panel with Division×Year FE key and state×year cluster key
df_w_council_vps <- pc %>%
  dplyr::mutate(state_fips = substr(place_fips, 1, 2)) %>%
  dplyr::left_join(state_div, by = "state_fips") %>%
  dplyr::left_join(w_div, by = "place_fips") %>%
  dplyr::mutate(
    w_a4r      = dplyr::if_else(treated_dummy == 1L, 1, tidyr::replace_na(w_eb, 0)),
    division   = as.factor(division),
    div_year   = interaction(division, year, drop = TRUE),
    trend      = year,
    state_year = interaction(state_fips, year, drop = TRUE)
  )

# ATT with Div×Year FE + city-specific linear trends + two-way clustering
m_att_divyr_council_vps <- feols(
  votes_per_seat_pct_vap ~ post2:treated_dummy + i(place_fips, trend) | place_fips + div_year,
  data = df_w_council_vps,
  weights = ~ w_a4r,
  cluster = ~ place_fips + state_year
)

cat("\n================== ATT (Council | votes_per_seat_pct_vap) — Div×Year FE + city trends ==================\n")
ct <- fixest::coeftable(m_att_divyr_council_vps)
if ("post2:treated_dummy" %in% rownames(ct)) {
  cat(sprintf("ATT (post×treated):  %.3f   [SE=%.3f]   p=%.3f\n",
              unname(ct["post2:treated_dummy","Estimate"]),
              unname(ct["post2:treated_dummy","Std. Error"]),
              unname(ct["post2:treated_dummy","Pr(>|t|)"])))
} else {
  cat("ATT (post×treated) coefficient not found in model.\n")
}
cat(sprintf("Obs: %d | Cities: %d | Div×Year FE cells: %d\n",
            nobs(m_att_divyr_council_vps),
            length(fixef(m_att_divyr_council_vps)$place_fips),
            length(fixef(m_att_divyr_council_vps)$div_year)))

# Event-study (SUNAB) with same FE, trends, weights, clustering
m_es_divyr_council_vps <- feols(
  votes_per_seat_pct_vap ~ sunab(g, year, ref.p = -1) + i(place_fips, trend) | place_fips + div_year,
  data = df_w_council_vps,
  weights = ~ w_a4r,
  cluster = ~ place_fips + state_year
)

cat("\n================== Pre-leads joint test (ES) — Council votes_per_seat_pct_vap ==================\n")
ctn <- rownames(coeftable(m_es_divyr_council_vps))
pre_terms <- ctn[grepl("^year::-(5|4|3|2|1)", ctn)]
if (length(pre_terms) > 0) {
  jw <- fixest::wald(m_es_divyr_council_vps, pre_terms)
  cat(sprintf("Wald stat = %.4f | p = %.6f | df1 = %s | df2 = %s | VCOV = %s\n",
              jw$stat, jw$p, jw$df1, jw$df2, jw$vcov))
} else {
  cat("No pre-lead terms detected at lags -5..-1 (check ES window / ref period).\n")
}

# (Optional) peek at first few event-time coefficients
es_ct <- as.data.frame(coeftable(m_es_divyr_council_vps)); es_ct$term <- rownames(es_ct)
et <- es_ct[grepl("^year::", es_ct$term), c("term","Estimate","Std. Error","Pr(>|t|)")]
if (nrow(et) > 0) {
  keep <- et[order(et$term), ]
  cat("\nEvent-time (subset):\n"); print(utils::head(keep, 10), row.names = FALSE)
}


###############################################################################
# COUNCIL — Light checks + Preferred spec w/ wild bootstrap + G-shock variant
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(fixest)
})

# --- Safety checks ---
stopifnot(exists("pc"), exists("state_div"))
if (!exists("run_citytrend_only")) stop("run_citytrend_only() not found in session.")

# --- Quick diagnostic helper (light checks) ---
council_quick_checks <- function(panel_df, outcome){
  cat("\n================== LIGHT CHECKS ::", outcome, "==================\n")
  # 1) Basic coverage by treated status
  cov1 <- panel_df %>%
    summarize(
      cities_total   = n_distinct(place_fips),
      cities_treated = n_distinct(place_fips[treated_dummy == 1]),
      cities_ctrl    = n_distinct(place_fips[treated_dummy == 0]),
      obs_total      = n(),
      obs_nonmiss    = sum(is.finite(.data[[outcome]]), na.rm = TRUE)
    )
  print(cov1)
  
  # 2) Non-missing counts per city (helps spot sparse outcomes)
  cov2 <- panel_df %>%
    group_by(place_fips) %>%
    summarize(n_nonmiss = sum(is.finite(.data[[outcome]])),
              ever_treated = as.integer(any(treated_dummy == 1L)),
              .groups = "drop") %>%
    summarize(
      p25 = quantile(n_nonmiss, .25, na.rm = TRUE),
      p50 = quantile(n_nonmiss, .50, na.rm = TRUE),
      p75 = quantile(n_nonmiss, .75, na.rm = TRUE),
      min = min(n_nonmiss, na.rm = TRUE),
      max = max(n_nonmiss, na.rm = TRUE)
    )
  cat("Non-missing outcome years per city (distribution):\n"); print(cov2)
  
  # 3) Simple pre/post counts within treated cities (not a test; just context)
  cov3 <- panel_df %>%
    filter(treated_dummy == 1L) %>%
    mutate(pre = as.integer(year < min(year[g >= 0], na.rm = TRUE))) %>%
    group_by(place_fips) %>%
    summarize(pre_obs = sum(pre == 1 & is.finite(.data[[outcome]])),
              post_obs = sum(pre == 0 & is.finite(.data[[outcome]])),
              .groups = "drop") %>%
    summarize(
      pre_med  = median(pre_obs,  na.rm = TRUE),
      post_med = median(post_obs, na.rm = TRUE)
    )
  cat("Median pre/post usable obs in treated cities:\n"); print(cov3)
  invisible(NULL)
}

# --- Wild cluster bootstrap helper (Rademacher; city clusters) ---
wb_att <- function(model, param = "post2:treated_dummy", B = 999){
  if (!requireNamespace("fwildclusterboot", quietly = TRUE)) {
    return(list(ok = FALSE, msg = "fwildclusterboot not installed; skipping wild bootstrap."))
  }
  out <- try(
    fwildclusterboot::boottest(
      model,
      clustid = "place_fips",
      param = param,
      B = B,
      bootcluster = "place_fips",     # wild cluster at city level
      type = "rademacher",
      impose_null = TRUE
    ),
    silent = TRUE
  )
  if (inherits(out, "try-error")) {
    list(ok = FALSE, msg = "boottest() failed; skipping wild bootstrap.")
  } else {
    list(ok = TRUE, pval = out$p_val, statistic = out$t_stat)
  }
}

# --- G-shock runner (national weights + Div×Year FE + city trends) ---
run_gshock_nat_divyear <- function(panel_df, outcome, res_citytrend_only_weights){
  w_nat <- res_citytrend_only_weights %>% select(place_fips, w_eb)
  
  df_w <- panel_df %>%
    mutate(state_fips = substr(place_fips, 1, 2)) %>%
    left_join(state_div, by = "state_fips") %>%
    left_join(w_nat, by = "place_fips") %>%
    mutate(
      w_a4r      = if_else(treated_dummy == 1L, 1, tidyr::replace_na(w_eb, 0)),
      division   = as.factor(division),
      div_year   = interaction(division, year, drop = TRUE),
      trend      = year,
      state_year = interaction(state_fips, year, drop = TRUE)
    )
  
  # ATT
  m_att <- feols(
    as.formula(paste0(outcome, " ~ post2:treated_dummy + i(place_fips, trend) | place_fips + div_year")),
    data = df_w, weights = ~ w_a4r, cluster = ~ place_fips + state_year
  )
  
  ct <- fixest::coeftable(m_att)
  if ("post2:treated_dummy" %in% rownames(ct)) {
    att_est <- unname(ct["post2:treated_dummy","Estimate"])
    att_se  <- unname(ct["post2:treated_dummy","Std. Error"])
    att_p   <- unname(ct["post2:treated_dummy","Pr(>|t|)"])
    cat(sprintf("ATT (post×treated):  %.3f   [SE=%.3f]   p=%.3f\n", att_est, att_se, att_p))
    wb <- wb_att(m_att, "post2:treated_dummy", B = 999)
    if (isTRUE(wb$ok)) cat(sprintf("   Wild bootstrap p (place_fips): %.3f\n", wb$pval)) else cat("   ", wb$msg, "\n", sep = "")
  } else {
    cat("ATT (post×treated) coefficient not found in model.\n")
  }
  cat(sprintf("Obs: %d | Cities: %d | Div×Year FE cells: %d\n",
              nobs(m_att), length(fixef(m_att)$place_fips), length(fixef(m_att)$div_year)))
  
  # ES + prelead test
  m_es <- feols(
    as.formula(paste0(outcome, " ~ sunab(g, year, ref.p = -1) + i(place_fips, trend) | place_fips + div_year")),
    data = df_w, weights = ~ w_a4r, cluster = ~ place_fips + state_year
  )
  ctn <- rownames(coeftable(m_es))
  pre_terms <- ctn[grepl("^year::-(5|4|3|2|1)", ctn)]
  if (length(pre_terms) > 0) {
    jw <- fixest::wald(m_es, pre_terms)
    cat(sprintf("Pre-leads joint test (ES): Wald = %.4f | p = %.6f | df1 = %s | df2 = %s | VCOV = %s\n",
                jw$stat, jw$p, jw$df1, jw$df2, jw$vcov))
  } else {
    cat("Pre-leads joint test (ES): no pre-lead terms detected in labels.\n")
  }
  
  
  invisible(list(att = m_att, es = m_es))
}

# =========================== RUNS (Council) ===========================
for (outcome in c("dem_share", "votes_per_seat_pct_vap")) {
  
  # (0) Light checks
  council_quick_checks(pc, outcome)
  
  # (A) Preferred spec: City-trend-only (national donors), with wild bootstrap
  cat("\n================== City-trend-only (national donors) ::", outcome, "==================\n")
  res_city <- run_citytrend_only(pc, "Council", outcome)   # prints balance + pre-lead test
  m_att_A  <- res_city$att
  
  # Wild bootstrap p-value on ATT
  ctA <- fixest::coeftable(m_att_A)
  if ("post2:treated_dummy" %in% rownames(ctA)) {
    cat("\n[ATT summary | City-trend-only]\n")
    cat(sprintf("Estimate = %.3f | SE = %.3f | p = %.3f\n",
                unname(ctA["post2:treated_dummy","Estimate"]),
                unname(ctA["post2:treated_dummy","Std. Error"]),
                unname(ctA["post2:treated_dummy","Pr(>|t|)"])))
    wbA <- wb_att(m_att_A, "post2:treated_dummy", B = 999)
    if (isTRUE(wbA$ok)) cat(sprintf("Wild bootstrap p (place_fips): %.3f\n", wbA$pval)) else cat(wbA$msg, "\n")
  }
  
  # (B) Peace-of-mind G-shock: national EB weights + Div×Year FE + city trends
  cat("\n================== G-SHOCK :: national weights + Div×Year FE + city trends ::",
      outcome, "==================\n")
  invisible(run_gshock_nat_divyear(pc, outcome, res_city$weights))
}

cat("\n---- DONE (Council: dem_share & votes_per_seat_pct_vap) ----\n")

