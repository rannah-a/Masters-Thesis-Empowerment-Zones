###############################################################################
# 06_baseline_estimation.R —— Main results (event study + DiD ATT)
###############################################################################
here::i_am("code/06_main_estimation.R")


### PART 1: QUICK UNBALANCED BASELINE (diagnostics + ATT)=====================

## 0) Packages ----------------------------------------------------------------
pkgs <- c("here","dplyr","stringr","fixest", "ggplotify", "gridExtra")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only = TRUE))
library(ggplot2)
theme_set(ggplot2::theme_minimal())

## 1) Load analysis panels (from Script 5) ------------------------------------
panel_mayor   <- readr::read_rds(here::here("data_clean","1_intermediate","city_panel_analysis_mayor.rds"))
panel_council <- readr::read_rds(here::here("data_clean","1_intermediate","city_panel_analysis_council.rds"))

## 2) Prep helpers (respect stored month-aware post) --------------------------
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

## 3) Quick diagnostics (sample sizes & support) -------------------------------
quick_diag <- function(df, y, panel_label){
  message("\n== ", panel_label, " :: ", y, " ==")
  d0 <- df %>% filter(!is.na(.data[[y]]))
  
  # Treated vs control counts; treated pre/post cells
  tt <- with(d0, table(treated_dummy, post2, useNA = "ifany"))
  print(tt)
  
  # How many clusters / rows in estimation sample?
  cat("Rows in model sample:", nrow(d0),
      "| Cities:", dplyr::n_distinct(d0$place_fips), "\n")
  
  # Treated event-time support (coarse)
  if ("treat_yr" %in% names(d0)) {
    d0 %>%
      mutate(rel = if ("rel_year" %in% names(d0)) rel_year else ifelse(treated_dummy==1L, year - treat_yr, NA_real_)) %>%
      filter(treated_dummy == 1L, !is.na(rel)) %>%
      group_by(rel) %>% summarise(n_rows = dplyr::n(), n_cities = n_distinct(place_fips), .groups="drop") %>%
      arrange(rel) %>% print(n = 30)
  }
  
  # Scale check treated vs control (medians)
  d0 %>%
    mutate(group = if_else(treated_dummy==1L,"treated","control")) %>%
    group_by(group) %>%
    summarise(p10 = quantile(.data[[y]], .1, na.rm=TRUE),
              p50 = median(.data[[y]], na.rm=TRUE),
              p90 = quantile(.data[[y]], .9, na.rm=TRUE),
              .groups="drop") %>%
    print()
  invisible(d0)
}

cat("\n### QUICK DIAGNOSTICS (unbalanced) ###\n")
pm_dem  <- quick_diag(pm, "dem_share",  "Mayor")
pm_turn <- quick_diag(pm, "turnout_pct","Mayor")
pc_dem  <- quick_diag(pc, "dem_share",  "Council")
pc_vps  <- quick_diag(pc, "votes_per_seat_pct_vap","Council")

## 4) Quick UNBALANCED ATT (TWFE DiD) -----------------------------------------
# Uses stored month-aware post2 and clusters by city. No balancing yet.

est_att <- function(df, y){
  d <- df %>% filter(!is.na(.data[[y]]))
  # interaction post2:treated_dummy gives ATT (treated x post)
  fixest::feols(
    as.formula(paste0(y, " ~ post2:treated_dummy | place_fips + year")),
    data = d, cluster = ~ place_fips
  )
}

cat("\n### QUICK UNBALANCED ATT (TWFE) ###\n")
m_mayor_dem   <- est_att(pm, "dem_share");               print(summary(m_mayor_dem))
m_mayor_turn  <- est_att(pm, "turnout_pct");             print(summary(m_mayor_turn))
m_council_dem <- est_att(pc, "dem_share");               print(summary(m_council_dem))
m_council_vps <- est_att(pc, "votes_per_seat_pct_vap");  print(summary(m_council_vps))

## 5) (Optional) SUNAB event-studies, still unbalanced -------------------------
# These are handy to eyeball pre-trends; we keep them lightweight here.
est_es <- function(df, y, ref = -1){
  d <- df %>% filter(!is.na(.data[[y]]))
  feols(
    as.formula(paste0(y, " ~ sunab(g, year, ref.p = ", ref, ") | place_fips + year")),
    data = d, cluster = ~ place_fips
  )
}

cat("\n### QUICK UNBALANCED ES (SUNAB) — optional ###\n")
es_mayor_dem   <- est_es(pm, "dem_share");               print(summary(es_mayor_dem))
es_council_dem <- est_es(pc, "dem_share");               print(summary(es_council_dem))

## 6) Tiny ATT table (on-screen) -----------------------------------------------
att_row <- function(mod){
  ct <- fixest::coeftable(mod)
  if (!"post2:treated_dummy" %in% rownames(ct)) return(tibble())
  tibble(
    estimate = unname(ct["post2:treated_dummy",1]),
    se       = unname(ct["post2:treated_dummy",2]),
    p_value  = unname(ct["post2:treated_dummy",4])
  )
}

att_quick <- bind_rows(
  tibble(panel="Mayor",   outcome="dem_share")             |> bind_cols(att_row(m_mayor_dem)),
  tibble(panel="Mayor",   outcome="turnout_pct")           |> bind_cols(att_row(m_mayor_turn)),
  tibble(panel="Council", outcome="dem_share")             |> bind_cols(att_row(m_council_dem)),
  tibble(panel="Council", outcome="votes_per_seat_pct_vap")|> bind_cols(att_row(m_council_vps))
)
cat("\n### QUICK ATT SUMMARY (unbalanced) ###\n"); print(att_quick, n = Inf)

# (Balanced / EB versions will follow after this block in Script 6.)


## =============================================================================
## PART 2: Division→Region donor pools + EB (exact→approx) with diagnostics
## Targets: mean_pre(y), ln_pop90, poverty_rate_90, unemploy_90, black_city_90
## City-invariant control weights -> use in weighted TWFE / SUNAB
## =============================================================================
suppressPackageStartupMessages({
  if (!requireNamespace("ebal", quietly = TRUE)) install.packages("ebal")
  library(dplyr); library(tidyr); library(purrr); library(fixest); library(ebal)
})

## --- Robust sdiff + weight diagnostics (no GUI deps) -------------------------
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

## --- State→Division/Region map (Census) --------------------------------------
state_div <- tribble(
  ~state_fips,~division,~region,
  "01",6,3,"02",9,4,"04",8,4,"05",7,3,"06",9,4,"08",8,4,"09",1,1,"10",5,3,"11",5,3,
  "12",5,3,"13",5,3,"15",9,4,"16",8,4,"17",3,2,"18",3,2,"19",4,2,"20",4,2,"21",6,3,
  "22",7,3,"23",1,1,"24",5,3,"25",1,1,"26",3,2,"27",4,2,"28",6,3,"29",4,2,"30",8,4,
  "31",4,2,"32",8,4,"33",1,1,"34",2,1,"35",8,4,"36",2,1,"37",5,3,"38",4,2,"39",3,2,
  "40",7,3,"41",9,4,"42",2,1,"44",1,1,"45",5,3,"46",4,2,"47",6,3,"48",7,3,"49",8,4,
  "50",1,1,"51",5,3,"53",9,4,"54",5,3,"55",3,2,"56",8,4
)
## (Divisions 1..9; Regions: 1=Northeast, 2=Midwest, 3=South, 4=West)

## --- Build city-level covars + pre-outcome and geo tags ----------------------
## CHANGE (A): add ln_pop90 to covariates/targets
build_city_table <- function(panel_df, y, pre_years = 1989:1993){
  df_geo <- panel_df %>%
    distinct(place_fips, .keep_all = TRUE) %>%
    mutate(state_fips = substr(place_fips, 1, 2)) %>%
    left_join(state_div, by = "state_fips")
  
  covs <- df_geo %>%
    transmute(place_fips,
              ln_pop90       = log(pmax(pop90_city, 1)),
              poverty_rate_90, unemploy_90, black_city_90,
              division = as.integer(division), region = as.integer(region))
  
  pre_y <- panel_df %>%
    filter(year %in% pre_years) %>%
    select(place_fips, year, !!sym(y)) %>%
    group_by(place_fips) %>%
    summarise(mean_pre = mean(.data[[y]], na.rm = TRUE), .groups="drop")
  
  trt <- panel_df %>%
    group_by(place_fips) %>%
    summarise(treated_city = as.integer(any(treated_dummy == 1L)), .groups="drop")
  
  covs %>% left_join(pre_y, by="place_fips") %>% left_join(trt, by="place_fips") %>%
    mutate(treated_city = replace_na(treated_city, 0L))
}

## --- Build donor pools: division first, else region, else national ------------
donor_pool_union <- function(city_tab, vars, min_div = 20, min_reg = 40){
  donors_union <- character(0)
  trt_rows <- city_tab %>% filter(treated_city == 1L)
  
  for (i in seq_len(nrow(trt_rows))){
    this_div <- trt_rows$division[i]
    this_reg <- trt_rows$region[i]
    
    cand_div <- city_tab %>%
      filter(treated_city == 0L, division == this_div) %>%
      drop_na(all_of(vars)) %>% pull(place_fips) %>% unique()
    
    if (length(cand_div) >= min_div) {
      donors_union <- union(donors_union, cand_div); next
    }
    
    cand_reg <- city_tab %>%
      filter(treated_city == 0L, region == this_reg) %>%
      drop_na(all_of(vars)) %>% pull(place_fips) %>% unique()
    
    if (length(cand_reg) >= min_reg) {
      donors_union <- union(donors_union, cand_reg); next
    }
    
    # fallback: national donors (complete cases)
    cand_nat <- city_tab %>%
      filter(treated_city == 0L) %>%
      drop_na(all_of(vars)) %>% pull(place_fips) %>% unique()
    
    donors_union <- union(donors_union, cand_nat)
  }
  
  donors_union
}

## --- Exact EB; fallback to approximate EB via Gaussian kernel -----------------
## CHANGE (B): anisotropic kernel with a mild up-weight on black_city_90
eb_fit_with_fallback <- function(city_tab, vars, donors_union){
  # treated set (complete)
  Tset <- city_tab %>% filter(treated_city == 1L) %>% drop_na(all_of(vars))
  Cset <- city_tab %>% filter(treated_city == 0L, place_fips %in% donors_union) %>%
    drop_na(all_of(vars))
  
  # targets: treated means
  target <- colMeans(as.matrix(Tset[, vars, drop=FALSE]), na.rm=TRUE)
  Xc     <- as.matrix(Cset[, vars, drop=FALSE])
  
  # try exact entropy balancing
  eb_try <- try(ebal::ebalance(target.margins = target, X = Xc,
                               base.weights = rep(1, nrow(Xc)),
                               maxit = 10000, constraint.tolerance = 1e-8),
                silent = TRUE)
  
  if (inherits(eb_try, "try-error") || is.null(eb_try$w)) {
    message("EB exact failed; switching to approximate EB (kernel).")
    
    ## --- approximate weights: Gaussian kernel (standardized), mild λ on black
    xc_s  <- scale(Xc)  # donors standardized per-column
    tgt_s <- as.numeric(scale(matrix(target, nrow=1),
                              center = attr(xc_s, "scaled:center"),
                              scale  = attr(xc_s, "scaled:scale")))
    
    # identify the column index for black_city_90
    j_black <- match("black_city_90", vars)
    # mild emphasis (λ in [1.0, 1.1]); default 1.05
    lambda_black <- 1.05
    
    # squared distances with anisotropic weight on the black dimension
    diff_mat <- xc_s - matrix(tgt_s, nrow = nrow(xc_s), ncol = ncol(xc_s), byrow = TRUE)
    d2 <- rowSums(diff_mat^2)
    if (!is.na(j_black)) {
      d2 <- d2 + (lambda_black^2 - 1) * (diff_mat[, j_black]^2)
    }
    
    # bandwidth and weights
    h  <- median(d2, na.rm=TRUE); if (!is.finite(h) || h <= 1e-8) h <- mean(d2, na.rm=TRUE) + 1
    w  <- exp(- d2 / h)
    if (mean(w, na.rm=TRUE) > 0) w <- w / mean(w, na.rm=TRUE)
    
  } else {
    w <- as.numeric(eb_try$w)
    if (mean(w, na.rm=TRUE) > 0) w <- w / mean(w, na.rm=TRUE)
  }
  
  # stitch weights to a city-level table
  ctrl_weights <- tibble(place_fips = Cset$place_fips, w_eb = as.numeric(w))
  out <- city_tab %>%
    left_join(ctrl_weights, by="place_fips") %>%
    mutate(w_eb = case_when(
      treated_city == 1L            ~ 1,
      is.finite(w_eb) & w_eb > 0    ~ w_eb,
      TRUE                          ~ 0
    ))
  out %>% select(place_fips, treated_city, all_of(vars), w_eb)
}

## --- One-stop wrapper: build donors, fit EB, print diagnostics ----------------
run_divregion_eb <- function(panel_df, panel_name, outcome){
  cat("\n\n=== ", panel_name, " | Outcome: ", outcome, " — Division→Region + EB ===\n", sep="")
  city_tab <- build_city_table(panel_df, y = outcome, pre_years = 1989:1993)
  
  # CHANGE (A): include ln_pop90 among targets
  vars <- c("mean_pre", "ln_pop90", "poverty_rate_90", "unemploy_90", "black_city_90")
  
  # who is treated?
  trt_ct <- city_tab %>% summarise(treated_cities = sum(treated_city==1L),
                                   donor_cities   = sum(treated_city==0L))
  print(trt_ct)
  
  # donor pool (union across treated cities): division → region → national
  donors_union <- donor_pool_union(city_tab, vars, min_div = 20, min_reg = 40)
  
  # EB with fallback
  wtab <- eb_fit_with_fallback(city_tab, vars, donors_union)
  
  # diagnostics: balance and weights (controls only)
  bal <- tibble(var = vars) %>%
    mutate(
      sdiff_unw = map_dbl(var, ~ sdiff(wtab[[.x]], wtab$treated_city, NULL)),
      sdiff_w   = map_dbl(var, ~ sdiff(wtab[[.x]], wtab$treated_city,
                                       w = ifelse(wtab$treated_city==1, 1, wtab$w_eb)))
    )
  wd  <- weight_diag(wtab$w_eb[wtab$treated_city==0])
  
  cat("\n-- Donor pool size (union of per-city pools) --\n")
  print(tibble(donors_kept = sum(wtab$treated_city==0 & wtab$w_eb>0),
               total_controls = sum(wtab$treated_city==0)))
  
  cat("\n-- Std. diffs (treated vs controls): before vs after weighting --\n")
  print(bal %>% mutate(across(starts_with("sdiff"), \(x) round(x, 3))))
  
  cat("\n-- Weight diagnostics (controls) --\n")
  print(wd %>% mutate(across(where(is.numeric), \(x) round(x, 3))))
  
  invisible(list(weights = wtab %>% select(place_fips, w_eb),
                 balance = bal, wdiag = wd))
}

## --- Attach weights to a panel and run weighted TWFE quickly ------------------
attach_weights <- function(panel_df, wtab){
  panel_df %>%
    left_join(wtab, by="place_fips") %>%
    mutate(w_eb = if_else(treated_dummy==1L, 1,
                          replace_na(w_eb, 0)))
}
est_att_w <- function(df, y){
  d <- df %>% filter(!is.na(.data[[y]]), w_eb > 0)
  feols(as.formula(paste0(y, " ~ post2:treated_dummy | place_fips + year")),
        data = d, weights = ~ w_eb, cluster = ~ place_fips)
}

## ===================== RUN (example, both panels) ============================
## Mayor
res_divreg_mayor_dem  <- run_divregion_eb(pm, "Mayor",   "dem_share")
pm_w                   <- attach_weights(pm, res_divreg_mayor_dem$weights)
cat("\n-- Weighted TWFE ATT (Mayor, dem_share) --\n"); print(summary(est_att_w(pm_w, "dem_share")))

res_divreg_mayor_turn <- run_divregion_eb(pm, "Mayor",   "turnout_pct")
pm_w2                  <- attach_weights(pm, res_divreg_mayor_turn$weights)
cat("\n-- Weighted TWFE ATT (Mayor, turnout_pct) --\n"); print(summary(est_att_w(pm_w2, "turnout_pct")))

## Council
res_divreg_cc_dem     <- run_divregion_eb(pc, "Council", "dem_share")
pc_w                  <- attach_weights(pc, res_divreg_cc_dem$weights)
cat("\n-- Weighted TWFE ATT (Council, dem_share) --\n"); print(summary(est_att_w(pc_w, "dem_share")))

res_divreg_cc_vps     <- run_divregion_eb(pc, "Council", "votes_per_seat_pct_vap")
pc_w2                 <- attach_weights(pc, res_divreg_cc_vps$weights)
cat("\n-- Weighted TWFE ATT (Council, votes_per_seat_pct_vap) --\n"); print(summary(est_att_w(pc_w2, "votes_per_seat_pct_vap")))


## ===================== PART 3: SAVE MODELS / PLOTS ===========================

# deps for saving outputs and pretty tables (no wild bootstrap)
pkgs_more <- c("readr","gridExtra","ggplotify","glue","scales","rempsyc")
for (p in pkgs_more) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(readr); library(gridExtra); library(ggplotify); library(glue); library(scales); library(rempsyc)

# output dirs
out_root   <- here::here("data_clean", "2_output")
out_models <- file.path(out_root, "models")
out_tables <- file.path(out_root, "tables")
out_figs   <- file.path(out_root, "figures")
dir.create(out_models, showWarnings = FALSE, recursive = TRUE)
dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(out_figs,   showWarnings = FALSE, recursive = TRUE)

# helpers to save small things
save_model <- function(obj, name) saveRDS(obj, file = file.path(out_models, paste0(name, ".rds")))
save_ctable <- function(model, coef_name, file_stub) {
  ct <- fixest::coeftable(model)
  if (!coef_name %in% rownames(ct)) return(invisible(NULL))
  row <- ct[coef_name, , drop = FALSE]
  tibble::tibble(
    term     = coef_name,
    estimate = as.numeric(row[1]),
    std_err  = as.numeric(row[2]),
    t_value  = as.numeric(row[3]),
    p_value  = as.numeric(row[4])
  ) |>
    readr::write_csv(file.path(out_tables, paste0(file_stub, ".csv")))
}

# tidy ES (sunab) and plotting
es_tidy <- function(es_model) {
  ct <- as.data.frame(fixest::coeftable(es_model))
  ct$term <- rownames(ct)
  keep <- grepl("^year::", ct$term)
  tibble::tibble(
    rel  = as.integer(sub("^year::", "", ct$term[keep])),
    beta = ct$Estimate[keep],
    se   = ct$`Std. Error`[keep]
  ) |>
    dplyr::mutate(lo = beta - 1.96 * se, hi = beta + 1.96 * se) |>
    dplyr::arrange(rel)
}
plot_es <- function(es_model, title, ylab, file_stub){
  df <- es_tidy(es_model)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = rel, y = beta)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1.6) +
    ggplot2::labs(title = title, x = "Event time (years from treatment; ref = -1)", y = ylab) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(out_figs, paste0(file_stub, ".png")), plot = p, width = 7, height = 4.5, dpi = 300)
  p
}

# --------------------- weight diagnostics for table ----------------------------
donor_ess_info <- function(df_w, y){
  d <- df_w %>% dplyr::filter(!is.na(.data[[y]]), w_eb > 0)
  donors <- d %>% dplyr::filter(treated_dummy == 0L) %>% dplyr::distinct(place_fips, w_eb)
  n_donors <- nrow(donors)
  if (n_donors == 0) return(list(donors = 0L, ess = NA_real_))
  w <- donors$w_eb
  ess <- (sum(w)^2) / sum(w^2)
  list(donors = n_donors, ess = ess)
}

# ---------------- weighted ATT (no bootstrap; clustered by city) ---------------
att_w <- function(df_w, y, stub){
  d  <- df_w %>% dplyr::filter(!is.na(.data[[y]]), w_eb > 0)
  m  <- fixest::feols(
    as.formula(paste0(y, " ~ post2:treated_dummy | place_fips + year")),
    data = d, weights = ~ w_eb, cluster = ~ place_fips
  )
  save_model(m, paste0(stub, "_TWFE"))
  save_ctable(m, "post2:treated_dummy", paste0("ATT_", stub, "_fixestSE"))
  list(model = m)
}

# --------------------- weighted SUNAB (save & plot) ----------------------------
sunab_w_and_plot <- function(df_w, y, title, ylab, stub, ref = -1){
  d  <- df_w %>% dplyr::filter(!is.na(.data[[y]]), w_eb > 0)
  es <- fixest::feols(
    as.formula(paste0(y, " ~ sunab(g, year, ref.p = ", ref, ") | place_fips + year")),
    data = d, weights = ~ w_eb, cluster = ~ place_fips
  )
  save_model(es, paste0(stub, "_SUNAB"))
  p <- plot_es(es, title, ylab, paste0(stub, "_ES"))
  # (optional) joint test of pre-treatment leads (−5..−1) = 0
  pre_terms <- rownames(coeftable(es))
  pre_terms <- pre_terms[grepl("^year::-(5|4|3|2|1)$", pre_terms)]
  wald_res  <- if (length(pre_terms)) fixest::wald(es, pre_terms) else NULL
  if (!is.null(wald_res)) save_model(wald_res, paste0(stub, "_SUNAB_pretreat_WALD"))
  list(model = es, plot = p, wald_pre = wald_res)
}

# ============== RUN MAIN WEIGHTED BASELINE (A4R Division→Region + EB) ==========

# We reuse the weight-attached panels you built earlier: pm_w, pm_w2, pc_w, pc_w2

# MAYOR — dem_share
cat("\n\n=== MAIN WEIGHTED BASELINE: MAYOR — dem_share ===\n")
res_mayor_dem_att <- att_w(pm_w,  "dem_share",  stub = "mayor_dem_share_divregEB")
res_mayor_dem_es  <- sunab_w_and_plot(pm_w, "dem_share",
                                      title = "Mayor — Event study (weighted; A4R+EB)",
                                      ylab  = "Effect on dem_share",
                                      stub  = "mayor_dem_share_divregEB")

# MAYOR — turnout_pct
cat("\n\n=== MAIN WEIGHTED BASELINE: MAYOR — turnout_pct ===\n")
res_mayor_turn_att <- att_w(pm_w2, "turnout_pct", stub = "mayor_turnout_divregEB")
res_mayor_turn_es  <- sunab_w_and_plot(pm_w2, "turnout_pct",
                                       title = "Mayor — Event study (weighted; A4R+EB)",
                                       ylab  = "Effect on turnout_pct (pp)",
                                       stub  = "mayor_turnout_divregEB")

# COUNCIL — dem_share
cat("\n\n=== MAIN WEIGHTED BASELINE: COUNCIL — dem_share ===\n")
res_cc_dem_att <- att_w(pc_w,  "dem_share",  stub = "council_dem_share_divregEB")
res_cc_dem_es  <- sunab_w_and_plot(pc_w, "dem_share",
                                   title = "Council — Event study (weighted; A4R+EB)",
                                   ylab  = "Effect on dem_share",
                                   stub  = "council_dem_share_divregEB")

# COUNCIL — votes_per_seat_pct_vap
cat("\n\n=== MAIN WEIGHTED BASELINE: COUNCIL — votes_per_seat_pct_vap ===\n")
res_cc_vps_att <- att_w(pc_w2, "votes_per_seat_pct_vap", stub = "council_vps_divregEB")
res_cc_vps_es  <- sunab_w_and_plot(pc_w2, "votes_per_seat_pct_vap",
                                   title = "Council — Event study (weighted; A4R+EB)",
                                   ylab  = "Effect on votes_per_seat_pct_vap",
                                   stub  = "council_vps_divregEB")

## ============== COMPACT BASELINE ATT TABLE (publication-ready) ================

.att_row <- function(mod, coef_name = "post2:treated_dummy") {
  ct <- fixest::coeftable(mod); stopifnot(coef_name %in% rownames(ct))
  tibble::tibble(
    estimate = unname(ct[coef_name, 1]),
    std_err  = unname(ct[coef_name, 2]),
    t_value  = unname(ct[coef_name, 3]),
    p_fixest = unname(ct[coef_name, 4])   # <-- clustered p (fixest)
  )
}

# Build rows, including donors & ESS from weights actually used in each outcome sample
att_tbl_weighted <- dplyr::bind_rows(
  dplyr::bind_cols(
    panel   = "Mayor",
    outcome = "dem_share",
    cities  = dplyr::n_distinct(pm_w$place_fips[pm_w$w_eb > 0]),
    donors  = donor_ess_info(pm_w, "dem_share")$donors,
    ess     = donor_ess_info(pm_w, "dem_share")$ess,
    .att_row(res_mayor_dem_att$model)
  ),
  dplyr::bind_cols(
    panel   = "Mayor",
    outcome = "turnout_pct",
    cities  = dplyr::n_distinct(pm_w2$place_fips[pm_w2$w_eb > 0]),
    donors  = donor_ess_info(pm_w2, "turnout_pct")$donors,
    ess     = donor_ess_info(pm_w2, "turnout_pct")$ess,
    .att_row(res_mayor_turn_att$model)
  ),
  dplyr::bind_cols(
    panel   = "Council",
    outcome = "dem_share",
    cities  = dplyr::n_distinct(pc_w$place_fips[pc_w$w_eb > 0]),
    donors  = donor_ess_info(pc_w, "dem_share")$donors,
    ess     = donor_ess_info(pc_w, "dem_share")$ess,
    .att_row(res_cc_dem_att$model)
  ),
  dplyr::bind_cols(
    panel   = "Council",
    outcome = "votes_per_seat_pct_vap",
    cities  = dplyr::n_distinct(pc_w2$place_fips[pc_w2$w_eb > 0]),
    donors  = donor_ess_info(pc_w2, "votes_per_seat_pct_vap")$donors,
    ess     = donor_ess_info(pc_w2, "votes_per_seat_pct_vap")$ess,
    .att_row(res_cc_vps_att$model)
  )
)

# --- Formatting helpers (lean, rempsyc-style presentation) ---
stars_from_p <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.01, "***",
                ifelse(p < 0.05, "**",
                       ifelse(p < 0.10, "*", "")
                )))
}
fmt_num <- function(x, digits = 3) number(x, accuracy = 10^-digits, big.mark = ",", trim = TRUE)

# Build the final display table:
# - Stars now use clustered p-values (fixest)
# - Replace "Wild p" column with "Clustered p"
att_tbl_fmt <- att_tbl_weighted %>%
  dplyr::mutate(
    stars   = stars_from_p(p_fixest),
    att_se  = glue("{fmt_num(estimate)}{stars}  ({fmt_num(std_err)})"),
    donors_ess  = glue("{donors} ({fmt_num(ess, 1)})"),
    Clustered.p = fmt_num(p_fixest, 3)
  ) %>%
  dplyr::transmute(
    Panel   = panel,
    Outcome = dplyr::case_when(
      outcome == "dem_share"               ~ "Dem. vote share",
      outcome == "turnout_pct"             ~ "Turnout (%)",
      outcome == "votes_per_seat_pct_vap"  ~ "Votes/seat (% of VAP)",
      TRUE ~ outcome
    ),
    Cities = as.integer(cities),
    `Donor cities (ESS)` = donors_ess,
    `ATT (SE)` = att_se,
    `Clustered p` = Clustered.p
  )

# --- Make a clean, publication-ready table (rempsyc style) ---
tbl_title <- "Average Treatment Effects (ATT) — Division→Region + Entropy Balancing"
tbl_note  <- paste(
  "Models: Two-way FE (city, year); standard errors and p-values clustered by city.",
  "Weights: Treated=1, controls EB-weighted; donors with zero weight excluded.",
  sep = " "
)

# Render & export
tg        <- gridExtra::tableGrob(as.data.frame(att_tbl_fmt), rows = NULL)
title_g   <- grid::textGrob(tbl_title, gp = grid::gpar(fontface = "bold", fontsize = 12))
note_wrapped <- stringr::str_wrap(tbl_note, width = 95)
note_g <- grid::textGrob(label = note_wrapped, x = 0, hjust = 0, gp = grid::gpar(fontsize = 8, col = "#444444"))
note_height <- grid::grobHeight(note_g) + grid::unit(3, "mm")

full_grob <- gridExtra::arrangeGrob(
  gridExtra::arrangeGrob(grid::textGrob(tbl_title, gp = grid::gpar(fontface = "bold", fontsize = 12))),
  gridExtra::tableGrob(as.data.frame(att_tbl_fmt), rows = NULL),
  gridExtra::arrangeGrob(note_g),
  heights = grid::unit.c(
    grid::unit(0.6, "in"),     # title
    grid::unit(1, "null"),     # table
    note_height                # note
  ),
  ncol = 1
)

# Draw & save
grid::grid.newpage(); grid::grid.draw(full_grob)

att_png <- file.path(out_tables, "06_baseline_ATT_weighted_A4R-EB.png")
att_pdf <- file.path(out_tables, "06_baseline_ATT_weighted_A4R-EB.pdf")

ggplot2::ggsave(att_png, ggplotify::as.ggplot(full_grob), width = 7.8, height = 4.6, dpi = 300)

pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
ggplot2::ggsave(
  filename = att_pdf,
  plot     = ggplotify::as.ggplot(full_grob),
  width    = 7.8, height = 4.6, device = pdf_device
)

cat("\nSaved publication-style ATT table to:\n  - ", att_png, "\n  - ", att_pdf, "\n", sep = "")


#### SUNAB Review:
## ==== REVIEW WEIGHTED SUNAB MODELS (A4R Division→Region + EB) ====

suppressPackageStartupMessages({
  for (p in c("here","dplyr","readr","fixest")) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
  library(here); library(dplyr); library(readr); library(fixest)
})

out_root   <- here::here("data_clean","2_output")
out_models <- file.path(out_root, "models")
out_tables <- file.path(out_root, "tables")
dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)

# Helper to load a saved SUNAB model by stub
load_sunab <- function(stub){
  f <- file.path(out_models, paste0(stub, "_SUNAB.rds"))
  if (!file.exists(f)) stop("Missing SUNAB model file: ", f)
  readRDS(f)
}

# Tidy ES (same shape as your plot helper)
es_tidy <- function(es_model) {
  ct <- as.data.frame(fixest::coeftable(es_model))
  ct$term <- rownames(ct)
  keep <- grepl("^year::", ct$term)
  tibble::tibble(
    rel  = as.integer(sub("^year::", "", ct$term[keep])),
    beta = ct$Estimate[keep],
    se   = ct$`Std. Error`[keep],
    t    = ct$`t value`[keep],
    p    = ct$`Pr(>|t|)`[keep]
  ) |>
    dplyr::arrange(rel) |>
    dplyr::mutate(lo = beta - 1.96 * se,
                  hi = beta + 1.96 * se)
}

# Print a concise ES review for console + save a CSV of all ES bins
review_es <- function(stub, nice_name, out_label){
  cat("\n\n=== WEIGHTED SUNAB REVIEW —", nice_name, "===\n")
  m <- load_sunab(stub)
  # 1) Compact model summary header
  cat("Obs:", nobs(m),
      "| Cities:", length(fixest::fixef(m)$place_fips),
      "| Years FE:", length(fixest::fixef(m)$year), "\n")
  if (!is.null(m$collin.var) && length(m$collin.var)) {
    cat("NOTE: Dropped for collinearity:", paste(m$collin.var, collapse=", "), "\n")
  }
  # 2) Joint pre-trend test on leads (-5..-1)
  ct_names <- rownames(coeftable(m))
  pre_terms <- ct_names[grepl("^year::-(5|4|3|2|1)$", ct_names)]
  if (length(pre_terms)) {
    cat("\nJoint test: pre-treatment leads (-5..-1) = 0\n")
    print(fixest::wald(m, pre_terms))
  } else {
    cat("\nNo standard pre-treatment lead bins found in the model.\n")
  }
  # 3) Show a small window around treatment in console
  df <- es_tidy(m)
  window_rel <- c(-5,-4,-3,-2,-1,0,1,2,3,5,7,10)
  show <- df %>% dplyr::filter(rel %in% window_rel)
  if (nrow(show)==0) show <- df %>% dplyr::slice_head(n=12)
  cat("\nSelected event-time bins (beta [lo, hi]):\n")
  print(show %>%
          mutate(ci = sprintf("%.3f [%.3f, %.3f]", beta, lo, hi)) %>%
          select(rel, beta, se, p, ci) %>%
          arrange(rel), n=Inf)
  # 4) Save full ES table to CSV
  csv_path <- file.path(out_tables, paste0(out_label, "_ES_coefs_weighted.csv"))
  readr::write_csv(df, csv_path)
  cat("\nSaved full ES coefficient table to:\n  - ", csv_path, "\n", sep="")
  invisible(m)
}

# ---- Run reviews for all four weighted SUNAB models ----
# Use the same stubs you used when saving:
#   mayor_dem_share_divregEB_SUNAB.rds, etc.
review_es("mayor_dem_share_divregEB",   "Mayor — dem_share",               "mayor_dem_share_divregEB")
review_es("mayor_turnout_divregEB",     "Mayor — turnout_pct",             "mayor_turnout_divregEB")
review_es("council_dem_share_divregEB", "Council — dem_share",             "council_dem_share_divregEB")
review_es("council_vps_divregEB",       "Council — votes_per_seat_pct_vap","council_vps_divregEB")
