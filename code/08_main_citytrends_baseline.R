###############################################################################
# 08_main_citytrends_baseline.R — Baseline with City Trends (A4R Div→Reg + EB)
# - Rebuilds Script 6 weights logic (Division→Region + EB exact→kernel fallback)
# - Estimates ATT and SUNAB with city-specific linear trends
# - Produces publication-style ATT table + wild-cluster bootstrap p-values
###############################################################################

here::i_am("code/08_main_citytrends_baseline.R")


## 0) Packages ------------------------------------------------------------------
pkgs <- c("here","dplyr","tidyr","stringr","purrr","readr","fixest","ebal",
          "ggplot2","ggplotify","gridExtra","glue","scales")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only=TRUE))
theme_set(ggplot2::theme_minimal())

# Wild bootstrap (CRAN, else GH fallback)
if (!requireNamespace("fwildclusterboot", quietly = TRUE)) {
  try(install.packages("fwildclusterboot", repos="https://cloud.r-project.org"))
  if (!requireNamespace("fwildclusterboot", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")
    remotes::install_github("s3alfisc/fwildclusterboot")
  }
}
library(fwildclusterboot)

## 1) Reproducibility -----------------------------------------------------------
set.seed(12345)
if (requireNamespace("dqrng", quietly=TRUE)) dqrng::dqset.seed(12345)

## 2) Load panels (from Script 5 outputs) --------------------------------------
panel_mayor   <- readr::read_rds(here::here("data_clean","1_intermediate","city_panel_analysis_mayor.rds"))
panel_council <- readr::read_rds(here::here("data_clean","1_intermediate","city_panel_analysis_council.rds"))

## 3) Prep helpers (respect stored month-aware post, build g) -------------------
prep_unbalanced <- function(df) {
  df %>%
    mutate(
      post2 = dplyr::case_when(
        "post" %in% names(df)                                 ~ as.integer(!!as.name("post")),
        "rel_year" %in% names(df) & treated_dummy == 1L       ~ as.integer(rel_year >= 0),
        treated_dummy == 1L & !is.na(treat_yr)                ~ as.integer(year >= treat_yr),
        TRUE                                                  ~ 0L
      ),
      g = if_else(treated_dummy == 1L & !is.na(treat_yr), as.integer(treat_yr), 0L),
      year = as.integer(year),
      treated_dummy = as.integer(treated_dummy)
    )
}
pm <- prep_unbalanced(panel_mayor)
pc <- prep_unbalanced(panel_council)

stopifnot(all(c("place_fips","year","treated_dummy","g","post2") %in% names(pm)))
stopifnot(all(c("place_fips","year","treated_dummy","g","post2") %in% names(pc)))

## 4) Census state → Division/Region (same as Script 6) -------------------------
state_div <- tibble::tribble(
  ~state_fips,~division,~region,
  "01",6,3,"02",9,4,"04",8,4,"05",7,3,"06",9,4,"08",8,4,"09",1,1,"10",5,3,"11",5,3,
  "12",5,3,"13",5,3,"15",9,4,"16",8,4,"17",3,2,"18",3,2,"19",4,2,"20",4,2,"21",6,3,
  "22",7,3,"23",1,1,"24",5,3,"25",1,1,"26",3,2,"27",4,2,"28",6,3,"29",4,2,"30",8,4,
  "31",4,2,"32",8,4,"33",1,1,"34",2,1,"35",8,4,"36",2,1,"37",5,3,"38",4,2,"39",3,2,
  "40",7,3,"41",9,4,"42",2,1,"44",1,1,"45",5,3,"46",4,2,"47",6,3,"48",7,3,"49",8,4,
  "50",1,1,"51",5,3,"53",9,4,"54",5,3,"55",3,2,"56",8,4
)

## 5) Build city-level table + donors + EB fallback (clone of Script 6) ---------
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

donor_pool_union <- function(city_tab, vars, min_div = 20, min_reg = 40){
  donors_union <- character(0)
  trt_rows <- city_tab %>% filter(treated_city == 1L)
  
  for (i in seq_len(nrow(trt_rows))){
    this_div <- trt_rows$division[i]
    this_reg <- trt_rows$region[i]
    
    cand_div <- city_tab %>%
      filter(treated_city == 0L, division == this_div) %>%
      drop_na(all_of(vars)) %>% pull(place_fips) %>% unique()
    if (length(cand_div) >= min_div) { donors_union <- union(donors_union, cand_div); next }
    
    cand_reg <- city_tab %>%
      filter(treated_city == 0L, region == this_reg) %>%
      drop_na(all_of(vars)) %>% pull(place_fips) %>% unique()
    if (length(cand_reg) >= min_reg) { donors_union <- union(donors_union, cand_reg); next }
    
    cand_nat <- city_tab %>%
      filter(treated_city == 0L) %>%
      drop_na(all_of(vars)) %>% pull(place_fips) %>% unique()
    donors_union <- union(donors_union, cand_nat)
  }
  donors_union
}

# exact EB; fallback to anisotropic Gaussian kernel (same as 6)
eb_fit_with_fallback <- function(city_tab, vars, donors_union){
  Tset <- city_tab %>% filter(treated_city == 1L) %>% drop_na(all_of(vars))
  Cset <- city_tab %>% filter(treated_city == 0L, place_fips %in% donors_union) %>%
    drop_na(all_of(vars))
  
  target <- colMeans(as.matrix(Tset[, vars, drop=FALSE]), na.rm=TRUE)
  Xc     <- as.matrix(Cset[, vars, drop=FALSE])
  
  eb_try <- try(ebal::ebalance(target.margins = target, X = Xc,
                               base.weights = rep(1, nrow(Xc)),
                               maxit = 10000, constraint.tolerance = 1e-8),
                silent = TRUE)
  
  if (inherits(eb_try, "try-error") || is.null(eb_try$w)) {
    message("EB exact failed; switching to approximate EB (kernel).")
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
  out <- city_tab %>%
    left_join(ctrl_weights, by="place_fips") %>%
    mutate(w_eb = case_when(
      treated_city == 1L            ~ 1,
      is.finite(w_eb) & w_eb > 0    ~ w_eb,
      TRUE                          ~ 0
    ))
  out %>% select(place_fips, treated_city, all_of(vars), w_eb)
}

run_divregion_eb <- function(panel_df, panel_name, outcome){
  cat("\n=== ", panel_name, " | ", outcome, " — Division→Region + EB ===\n", sep="")
  city_tab <- build_city_table(panel_df, y = outcome, pre_years = 1989:1993)
  vars <- c("mean_pre","ln_pop90","poverty_rate_90","unemploy_90","black_city_90")
  donors_union <- donor_pool_union(city_tab, vars, min_div = 20, min_reg = 40)
  wtab <- eb_fit_with_fallback(city_tab, vars, donors_union)
  invisible(list(weights = wtab %>% select(place_fips, w_eb)))
}

attach_weights <- function(panel_df, wtab){
  panel_df %>% left_join(wtab, by="place_fips") %>%
    mutate(w_eb = if_else(treated_dummy==1L, 1, replace_na(w_eb, 0)))
}

## 6) Estimators with CITY TRENDS ------------------------------------------------
est_att_w_trend <- function(df, y, wcol = "w_eb"){
  d <- df %>%
    filter(!is.na(.data[[y]]), .data[[wcol]] > 0) %>%
    mutate(trend = year)
  fixest::feols(
    as.formula(paste0(y, " ~ post2:treated_dummy + i(place_fips, trend) | place_fips + year")),
    data = d, weights = as.formula(paste0("~", wcol)), cluster = ~ place_fips
  )
}

est_es_w_trend <- function(df, y, wcol = "w_eb", ref = -1){
  d <- df %>%
    filter(!is.na(.data[[y]]), .data[[wcol]] > 0) %>%
    mutate(trend = year)
  fixest::feols(
    as.formula(paste0(y, " ~ sunab(g, year, ref.p = ", ref, ") + i(place_fips, trend) | place_fips + year")),
    data = d, weights = as.formula(paste0("~", wcol)), cluster = ~ place_fips
  )
}

## 7) Output dirs ---------------------------------------------------------------
out_root   <- here::here("data_clean","2_output")
out_models <- file.path(out_root, "models")
out_tables <- file.path(out_root, "tables")
out_figs   <- file.path(out_root, "figures")
dir.create(out_models, showWarnings=FALSE, recursive=TRUE)
dir.create(out_tables, showWarnings=FALSE, recursive=TRUE)
dir.create(out_figs,   showWarnings=FALSE, recursive=TRUE)

save_model <- function(obj, name) saveRDS(obj, file = file.path(out_models, paste0(name, ".rds")))
save_df    <- function(df,  name) saveRDS(df,  file = file.path(out_models, paste0(name, "_DATA.rds")))

## 8) Run: build weights (as in 6) → attach → ATT/ES with city trends ----------
# Weights per outcome (like Script 6)
w_may_dem   <- run_divregion_eb(pm, "Mayor",   "dem_share")$weights
w_may_turn  <- run_divregion_eb(pm, "Mayor",   "turnout_pct")$weights
w_cou_dem   <- run_divregion_eb(pc, "Council", "dem_share")$weights
w_cou_vps   <- run_divregion_eb(pc, "Council", "votes_per_seat_pct_vap")$weights

pm_dem  <- attach_weights(pm, w_may_dem)
pm_turn <- attach_weights(pm, w_may_turn)
pc_dem  <- attach_weights(pc, w_cou_dem)
pc_vps  <- attach_weights(pc, w_cou_vps)

# ATT with city trends
m_may_dem   <- est_att_w_trend(pm_dem,  "dem_share")
m_may_turn  <- est_att_w_trend(pm_turn, "turnout_pct")
m_cou_dem   <- est_att_w_trend(pc_dem,  "dem_share")
m_cou_vps   <- est_att_w_trend(pc_vps,  "votes_per_seat_pct_vap")

# ES with city trends (optional plots, same look as 6)
es_tidy <- function(es_model) {
  ct <- as.data.frame(fixest::coeftable(es_model)); ct$term <- rownames(ct)
  keep <- grepl("^year::", ct$term)
  tibble::tibble(
    rel  = as.integer(sub("^year::", "", ct$term[keep])),
    beta = ct$Estimate[keep], se = ct$`Std. Error`[keep]
  ) %>% mutate(lo = beta - 1.96*se, hi = beta + 1.96*se) %>% arrange(rel)
}
plot_es <- function(es_model, title, ylab, file_stub){
  df <- es_tidy(es_model)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=rel,y=beta)) +
    ggplot2::geom_hline(yintercept=0, linetype=2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lo,ymax=hi), alpha=.15) +
    ggplot2::geom_line() + ggplot2::geom_point(size=1.6) +
    ggplot2::labs(title=title, x="Event time ((years from treatment;ref=-1)", y=ylab)
  ggplot2::ggsave(file.path(out_figs, paste0(file_stub, ".png")), plot=p, width=7, height=4.5, dpi=300)
  p
}

es_may_dem  <- est_es_w_trend(pm_dem,  "dem_share")
es_may_turn <- est_es_w_trend(pm_turn, "turnout_pct")
es_cou_dem  <- est_es_w_trend(pc_dem,  "dem_share")
es_cou_vps  <- est_es_w_trend(pc_vps,  "votes_per_seat_pct_vap")

# Save models + data used (so boottest can be re-run later if you want)
save_model(m_may_dem,   "08_CITYTREND_mayor_dem_share_ATT")
save_model(m_may_turn,  "08_CITYTREND_mayor_turnout_pct_ATT")
save_model(m_cou_dem,   "08_CITYTREND_council_dem_share_ATT")
save_model(m_cou_vps,   "08_CITYTREND_council_votes_per_seat_pct_vap_ATT")

save_model(es_may_dem,  "08_CITYTREND_mayor_dem_share_ES")
save_model(es_may_turn, "08_CITYTREND_mayor_turnout_pct_ES")
save_model(es_cou_dem,  "08_CITYTREND_council_dem_share_ES")
save_model(es_cou_vps,  "08_CITYTREND_council_votes_per_seat_pct_vap_ES")

save_df(pm_dem,  "08_CITYTREND_mayor_dem_share")
save_df(pm_turn, "08_CITYTREND_mayor_turnout_pct")
save_df(pc_dem,  "08_CITYTREND_council_dem_share")
save_df(pc_vps,  "08_CITYTREND_council_votes_per_seat_pct_vap")

# ES plots (optional)
plot_es(es_may_dem,  "Mayor — Event study (weighted; city-trend adj)",           "Effect on dem_share",             "08_CITYTREND_mayor_dem_share_ES")
plot_es(es_may_turn, "Mayor — Event study (weighted; city-trend adj)",           "Effect on turnout_pct (pp)",      "08_CITYTREND_mayor_turnout_pct_ES")
plot_es(es_cou_dem,  "Council — Event study (weighted; city-trend adj)",         "Effect on dem_share",             "08_CITYTREND_council_dem_share_ES")
plot_es(es_cou_vps,  "Council — Event study (weighted; city-trend adj)",         "Effect on votes_per_seat_pct_vap","08_CITYTREND_council_vps_ES")

## 9)## ===== Compact ATT table (clustered p) — same style as Script 6 ===============

# deps for saving outputs and pretty tables (no wild bootstrap here)
for (p in c("readr","gridExtra","ggplotify","glue","scales","rempsyc")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(readr); library(gridExtra); library(ggplotify); library(glue); library(scales)

#--- helpers --------------------------------------------------------------------
.att_row <- function(mod, coef_name = "post2:treated_dummy") {
  ct <- fixest::coeftable(mod); stopifnot(coef_name %in% rownames(ct))
  tibble::tibble(
    estimate = unname(ct[coef_name, 1]),
    std_err  = unname(ct[coef_name, 2]),
    t_value  = unname(ct[coef_name, 3]),
    p_fixest = unname(ct[coef_name, 4])   # clustered-by-city p from model
  )
}

donor_ess_info <- function(df_w, y){
  d <- df_w %>% dplyr::filter(!is.na(.data[[y]]), w_eb > 0)
  donors <- d %>% dplyr::filter(treated_dummy == 0L) %>%
    dplyr::distinct(place_fips, w_eb)
  n_donors <- nrow(donors)
  if (n_donors == 0) return(list(donors = 0L, ess = NA_real_))
  w <- donors$w_eb
  ess <- (sum(w)^2) / sum(w^2)
  list(donors = n_donors, ess = ess)
}

stars_from_p <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.01, "***",
                ifelse(p < 0.05, "**",
                       ifelse(p < 0.10, "*", ""))))
}
fmt_num <- function(x, digits = 3)
  number(x, accuracy = 10^-digits, big.mark = ",", trim = TRUE)

#--- build rows from your Script 8 models & panels ------------------------------
# Requires objects already created above: pm_dem, pm_turn, pc_dem, pc_vps
# and models: m_may_dem, m_may_turn, m_cou_dem, m_cou_vps

att_tbl_citytrend <- dplyr::bind_rows(
  dplyr::bind_cols(
    panel   = "Mayor",
    outcome = "dem_share",
    cities  = dplyr::n_distinct(pm_dem$place_fips[pm_dem$w_eb > 0]),
    donors  = donor_ess_info(pm_dem,  "dem_share")$donors,
    ess     = donor_ess_info(pm_dem,  "dem_share")$ess,
    .att_row(m_may_dem)
  ),
  dplyr::bind_cols(
    panel   = "Mayor",
    outcome = "turnout_pct",
    cities  = dplyr::n_distinct(pm_turn$place_fips[pm_turn$w_eb > 0]),
    donors  = donor_ess_info(pm_turn, "turnout_pct")$donors,
    ess     = donor_ess_info(pm_turn, "turnout_pct")$ess,
    .att_row(m_may_turn)
  ),
  dplyr::bind_cols(
    panel   = "Council",
    outcome = "dem_share",
    cities  = dplyr::n_distinct(pc_dem$place_fips[pc_dem$w_eb > 0]),
    donors  = donor_ess_info(pc_dem,  "dem_share")$donors,
    ess     = donor_ess_info(pc_dem,  "dem_share")$ess,
    .att_row(m_cou_dem)
  ),
  dplyr::bind_cols(
    panel   = "Council",
    outcome = "votes_per_seat_pct_vap",
    cities  = dplyr::n_distinct(pc_vps$place_fips[pc_vps$w_eb > 0]),
    donors  = donor_ess_info(pc_vps,  "votes_per_seat_pct_vap")$donors,
    ess     = donor_ess_info(pc_vps,  "votes_per_seat_pct_vap")$ess,
    .att_row(m_cou_vps)
  )
)

#--- format exactly like Script 6 (clustered p + stars, donors with ESS) --------
att_tbl_fmt <- att_tbl_citytrend %>%
  dplyr::mutate(
    stars   = stars_from_p(p_fixest),
    att_se  = glue("{fmt_num(estimate)}{stars}  ({fmt_num(std_err)})"),
    donors_ess  = glue("{donors} ({fmt_num(ess, 1)})"),
    `Clustered p` = fmt_num(p_fixest, 3)
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
    `Clustered p`
  )

#--- draw & save (same look/exports as Script 6) --------------------------------
tbl_title <- "ATT — Division→Region + EB Weights, with City Trends"
tbl_note  <- paste(
  "Models: Two-way FE (city, year) with city-specific linear trends;",
  "SEs & p-values clustered by city.",
  "Weights: Treated=1; controls EB-weighted (A4R donor rules); donors with zero weight excluded."
)

tg        <- gridExtra::tableGrob(as.data.frame(att_tbl_fmt), rows = NULL)
title_g   <- grid::textGrob(tbl_title, gp = grid::gpar(fontface = "bold", fontsize = 12))
note_wrapped <- stringr::str_wrap(tbl_note, width = 95)
note_g    <- grid::textGrob(label = note_wrapped, x = 0, hjust = 0,
                            gp = grid::gpar(fontsize = 8, col = "#444444"))
note_height <- grid::grobHeight(note_g) + grid::unit(3, "mm")

full_grob <- gridExtra::arrangeGrob(
  gridExtra::arrangeGrob(title_g),
  gridExtra::tableGrob(as.data.frame(att_tbl_fmt), rows = NULL),
  gridExtra::arrangeGrob(note_g),
  heights = grid::unit.c(
    grid::unit(0.6, "in"),     # title
    grid::unit(1, "null"),     # table
    note_height                # note
  ),
  ncol = 1
)

grid::grid.newpage(); grid::grid.draw(full_grob)

att_png <- file.path(out_tables, "08_CITYTREND_ATT_weighted_A4R-EB.png")
att_pdf <- file.path(out_tables, "08_CITYTREND_ATT_weighted_A4R-EB.pdf")

ggplot2::ggsave(att_png, ggplotify::as.ggplot(full_grob), width = 7.8, height = 4.6, dpi = 300)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
ggplot2::ggsave(filename = att_pdf,
                plot     = ggplotify::as.ggplot(full_grob),
                width    = 7.8, height = 4.6, device = pdf_device)

cat("\nSaved city-trend ATT table to:\n  - ", att_png, "\n  - ", att_pdf, "\n", sep = "")

#10) Wild-cluster bootstrap p-values (cluster = city) -------------------------
#prepare a named df first 
# Mayor — dem_share
pm_dem_fit <- pm_dem %>% dplyr::filter(!is.na(dem_share), w_eb > 0) %>% dplyr::mutate(trend = year)
m_may_dem  <- fixest::feols(
  dem_share ~ post2:treated_dummy + i(place_fips, trend) | place_fips + year,
  data = pm_dem_fit, weights = ~ w_eb, cluster = ~ place_fips
)

# Mayor- turnout_pct
pm_turn_fit <- pm_turn %>% filter(!is.na(turnout_pct), w_eb > 0) %>% mutate(trend = year)
m_may_turn  <- feols(turnout_pct ~ post2:treated_dummy + i(place_fips, trend) | place_fips + year,
                     data = pm_turn_fit, weights = ~ w_eb, cluster = ~ place_fips)
# City Council - dem_share 
pc_dem_fit  <- pc_dem %>% filter(!is.na(dem_share), w_eb > 0) %>% mutate(trend = year)
m_cou_dem   <- feols(dem_share ~ post2:treated_dummy + i(place_fips, trend) | place_fips + year,
                     data = pc_dem_fit, weights = ~ w_eb, cluster = ~ place_fips)
# City Council - votes_per_seat_pct_vap
pc_vps_fit  <- pc_vps %>% filter(!is.na(votes_per_seat_pct_vap), w_eb > 0) %>% mutate(trend = year)
m_cou_vps   <- feols(votes_per_seat_pct_vap ~ post2:treated_dummy + i(place_fips, trend) | place_fips + year,
                     data = pc_vps_fit, weights = ~ w_eb, cluster = ~ place_fips)

## Make FE-safe columns in each fit dataset
pm_dem_fit  <- pm_dem  %>% dplyr::filter(!is.na(dem_share), w_eb > 0) %>% 
  dplyr::mutate(trend = year, city_fe = factor(place_fips), year_fe = factor(year))
pm_turn_fit <- pm_turn %>% dplyr::filter(!is.na(turnout_pct), w_eb > 0) %>% 
  dplyr::mutate(trend = year, city_fe = factor(place_fips), year_fe = factor(year))
pc_dem_fit  <- pc_dem  %>% dplyr::filter(!is.na(dem_share), w_eb > 0) %>% 
  dplyr::mutate(trend = year, city_fe = factor(place_fips), year_fe = factor(year))
pc_vps_fit  <- pc_vps  %>% dplyr::filter(!is.na(votes_per_seat_pct_vap), w_eb > 0) %>% 
  dplyr::mutate(trend = year, city_fe = factor(place_fips), year_fe = factor(year))

## Mayor — dem_share
m_may_dem <- fixest::feols(
  dem_share ~ post2:treated_dummy + i(city_fe, trend),
  data    = pm_dem_fit,
  fixef   = c("city_fe", "year_fe"),
  weights = ~ w_eb,
  cluster = ~ city_fe
)

## Mayor — turnout_pct
m_may_turn <- fixest::feols(
  turnout_pct ~ post2:treated_dummy + i(city_fe, trend),
  data    = pm_turn_fit,
  fixef   = c("city_fe", "year_fe"),
  weights = ~ w_eb,
  cluster = ~ city_fe
)

## Council — dem_share
m_cou_dem <- fixest::feols(
  dem_share ~ post2:treated_dummy + i(city_fe, trend),
  data    = pc_dem_fit,
  fixef   = c("city_fe", "year_fe"),
  weights = ~ w_eb,
  cluster = ~ city_fe
)

## Council — votes_per_seat_pct_vap
m_cou_vps <- fixest::feols(
  votes_per_seat_pct_vap ~ post2:treated_dummy + i(city_fe, trend),
  data    = pc_vps_fit,
  fixef   = c("city_fe", "year_fe"),
  weights = ~ w_eb,
  cluster = ~ city_fe
)

set.seed(12345); if (requireNamespace("dqrng", quietly=TRUE)) dqrng::dqset.seed(12345)

wb <- function(m) fwildclusterboot::boottest(
  m,
  param       = "post2:treated_dummy",
  clustid     = ~ city_fe,
  B           = 9999,
  type        = "rademacher",
  impose_null = TRUE
)

wb_may_dem  <- wb(m_may_dem)
wb_may_turn <- wb(m_may_turn)
wb_cou_dem  <- wb(m_cou_dem)
wb_cou_vps  <- wb(m_cou_vps)

## ========================= Wild CLUSTER bootstrap — pragmatic block =========================
## Place this after you have:
##   - pm_dem_fit, pm_turn_fit, pc_dem_fit, pc_vps_fit  (FE-safe data frames)
##   - m_may_dem, m_may_turn, m_cou_dem, m_cou_vps      (fixest models with city trends)
## And after you set seeds.

suppressPackageStartupMessages({
  for (p in c("fwildclusterboot","dplyr","purrr","tibble","stringr")) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  }
  library(dplyr); library(purrr); library(tibble); library(stringr)
})

# ------------------------------ Helpers ---------------------------------------

# 1) Quick cluster-level diagnostics for the ATT regressor
wb_sanity <- function(df_fit, yname, cluster = "city_fe",
                      reg_name = "post2:treated_dummy", wcol = "w_eb") {
  
  stopifnot(all(c(cluster, wcol, "post2","treated_dummy") %in% names(df_fit)))
  # the ATT regressor used in the model
  df <- df_fit %>%
    mutate(att_x = as.numeric(post2) * as.numeric(treated_dummy)) %>%
    filter(is.finite(.data[[yname]]), .data[[wcol]] > 0)
  
  # cluster counts + basic variation checks
  cl_tbl <- df %>%
    group_by(.data[[cluster]]) %>%
    summarise(
      n_obs     = n(),
      n_years   = n_distinct(year),
      has_pre   = any(post2 == 0),
      has_post  = any(post2 == 1),
      var_att_x = stats::var(att_x),
      w_sum     = sum(.data[[wcol]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      informative = has_pre & has_post & is.finite(var_att_x) & var_att_x > 0
    )
  
  G          <- nrow(cl_tbl)
  G_inform   <- sum(cl_tbl$informative)
  G_single   <- sum(cl_tbl$n_obs <= 2)                  # very small clusters
  G_no_post  <- sum(!cl_tbl$has_post)
  G_no_pre   <- sum(!cl_tbl$has_pre)
  
  # rules of thumb: wild cluster works ok-ish at ~10–15+ clusters, better >20
  flags <- list(
    ok_G            = G >= 15,
    ok_G_inform     = G_inform >= 10,
    low_singletons  = (G_single / G) <= 0.15
  )
  
  msg <- c(
    sprintf("Clusters (G) = %d; Informative clusters = %d", G, G_inform),
    sprintf("Singleton-ish clusters (<=2 obs) = %d (%.1f%%)", G_single, 100*G_single/max(G,1)),
    sprintf("Clusters with no PRE = %d; no POST = %d", G_no_pre, G_no_post),
    sprintf("Rule-of-thumb checks -> G>=15: %s | Informative>=10: %s | singleton share<=15%%: %s",
            ifelse(flags$ok_G,"OK","WARN"),
            ifelse(flags$ok_G_inform,"OK","WARN"),
            ifelse(flags$low_singletons,"OK","WARN"))
  )
  
  list(
    cluster_table = cl_tbl,
    summary = tibble(
      G = G, G_inform = G_inform, G_single = G_single,
      G_no_pre = G_no_pre, G_no_post = G_no_post,
      ok_G = flags$ok_G, ok_G_inform = flags$ok_G_inform, ok_singletons = flags$low_singletons
    ),
    messages = msg
  )
}

# 2) A safe wrapper around fwildclusterboot::boottest with pragmatic fallbacks
safe_boottest <- function(model, df_fit, cluster_var = "city_fe",
                          param = "post2:treated_dummy",
                          B_main = 9999, type_main = "rademacher",
                          impose_null = TRUE) {
  
  # cluster id vector (avoid formula parsing edge cases)
  cl_vec <- droplevels(df_fit[[cluster_var]])
  
  run_once <- function(B, type, impose_null) {
    fwildclusterboot::boottest(
      model,
      param       = param,
      clustid     = cl_vec,
      B           = B,
      type        = type,        # "rademacher" default; fallback tries "mammen"
      impose_null = impose_null
    )
  }
  
  # Try main
  res <- try(run_once(B_main, type_main, impose_null), silent = TRUE)
  if (!inherits(res, "try-error")) return(list(result = res, variant = "main"))
  
  # If computationally singular or similar, try a lighter/fallback variant:
  # (a) switch weight scheme to "mammen", (b) reduce B a bit to avoid extreme matrices
  res2 <- try(run_once(max(3999, floor(B_main/3)), "mammen", impose_null), silent = TRUE)
  if (!inherits(res2, "try-error")) return(list(result = res2, variant = "fallback:mammen"))
  
  # Last resort: do NOT impose the null (can help with near-singular designs)
  res3 <- try(run_once(max(1999, floor(B_main/5)), "mammen", FALSE), silent = TRUE)
  if (!inherits(res3, "try-error")) return(list(result = res3, variant = "fallback:mammen,no-null"))
  
  # If we get here, it failed
  list(error = TRUE,
       messages = c("boottest failed in all variants.",
                    "Common causes: too few informative clusters, near-perfect collinearity, or all-treated/all-control clusters after weighting."))
}

# 3) A compact printer for sanity + bootstrap outputs
print_wb <- function(label, sanity, wb_out) {
  cat("\n---------------------\n", label, "\n---------------------\n", sep = "")
  cat(paste0(sanity$messages, collapse = "\n"), "\n")
  if (!isTRUE(wb_out$error)) {
    cat("Wildcard run:", wb_out$variant, "\n")
    print(wb_out$result)
  } else {
    cat(paste(wb_out$messages, collapse = " "), "\n")
  }
}

# ------------------------------ Run sanity + WB for all four outcomes ---------

set.seed(12345); if (requireNamespace("dqrng", quietly=TRUE)) dqrng::dqset.seed(12345)

# Label -> (data, model, outcome_col)
WB_CASES <- list(
  "Mayor — dem_share"              = list(df = pm_dem_fit,  m = m_may_dem,  y = "dem_share"),
  "Mayor — turnout_pct"            = list(df = pm_turn_fit, m = m_may_turn, y = "turnout_pct"),
  "Council — dem_share"            = list(df = pc_dem_fit,  m = m_cou_dem,  y = "dem_share"),
  "Council — votes_per_seat_pct_vap" = list(df = pc_vps_fit, m = m_cou_vps,  y = "votes_per_seat_pct_vap")
)

wb_results <- list()
wb_summaries <- list()

for (nm in names(WB_CASES)) {
  df_fit <- WB_CASES[[nm]]$df
  mod    <- WB_CASES[[nm]]$m
  ycol   <- WB_CASES[[nm]]$y
  
  # Step 1: sanity
  sanity <- wb_sanity(df_fit, yname = ycol, cluster = "city_fe", reg_name = "post2:treated_dummy", wcol = "w_eb")
  wb_summaries[[nm]] <- sanity$summary
  
  # Step 2: attempt boottest (+ pragmatic fallbacks)
  wb_out <- safe_boottest(mod, df_fit, cluster_var = "city_fe", param = "post2:treated_dummy",
                          B_main = 9999, type_main = "rademacher", impose_null = TRUE)
  wb_results[[nm]] <- wb_out
  
  # Step 3: print
  print_wb(nm, sanity, wb_out)
}

# Collect a compact sanity table and flag rows that fail rules of thumb
wb_sanity_table <- bind_rows(lapply(names(wb_summaries), function(nm) {
  wb_summaries[[nm]] %>% mutate(model = nm, .before = 1)
}))
readr::write_csv(wb_sanity_table, file.path(out_tables, "08_CITYTREND_WB_sanity_summary.csv"))
cat("\nSaved WB sanity summary to tables/08_CITYTREND_WB_sanity_summary.csv\n")

# Optionally, save each wb object that succeeded
invisible({
  for (nm in names(wb_results)) {
    if (!isTRUE(wb_results[[nm]]$error)) {
      saveRDS(wb_results[[nm]]$result,
              file = file.path(out_models, paste0("08_CITYTREND_", gsub("[^A-Za-z0-9]+","_", nm), "_WB.rds")))
    }
  }
})

## ---------------- Notes / interpretation guidance ----------------
## • If you see WARN in the sanity lines (few clusters, few “informative” clusters,
##   many singletons), treat bootstrap p-values cautiously. Consider reporting
##   clustered-by-city (fixest) p-values as primary and WB as a robustness check.
## • If all WB variants fail with “computationally singular”, typical causes are:
##   (a) the ATT regressor has near-zero variation after absorbing FE + trends,
##   (b) some clusters have only PRE or only POST, or (c) weights zero out most controls.
##   Remedies you can try (domain-appropriate):
##     – Trim to clusters with both PRE and POST (for treated) in a robustness run.
##     – Reduce the trend flexibility (e.g., check if a few city-specific trends were dropped
##       for collinearity; fixest warns about them).
##     – Report WB for outcomes where sanity checks pass and clustered SEs otherwise.



## ================= Timing (E) robustness for turnout =================
## Requires:
##   - pm_turn_fit (has: turnout_pct, city_fe, year_fe, trend, post2, treated_dummy, w_eb, year)
##   - m_may_turn  (your baseline city-trend model)
suppressPackageStartupMessages({
  for (p in c("dplyr","stringr","fixest","tibble")) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  }
  library(dplyr); library(stringr); library(fixest); library(tibble)
})

## ---- Helpers ----

# 1) Create federal-cycle FE and month FE (defensive)
add_timing_fes <- function(df){
  stopifnot(all(c("year","city_fe","post2","treated_dummy") %in% names(df)))
  # Federal election cycle FE: Presidential (%%4==0), Midterm (%%4==2), Odd (%%2==1)
  cycle_fe <- ifelse(df$year %% 2 == 1, "Odd",
                     ifelse(df$year %% 4 == 0, "Presidential", "Midterm"))
  df$cycle_fe <- factor(cycle_fe, levels = c("Presidential","Midterm","Odd"))
  
  # Try to find a month column; if none, make a single-level placeholder
  month_candidates <- c("month_final")
  mcol <- month_candidates[month_candidates %in% names(df)]
  if (length(mcol) >= 1) {
    mm <- as.integer(df[[mcol[1]]])
    if (!all(is.finite(mm))) mm[!is.finite(mm)] <- -1L
    df$month_fe <- factor(mm)
    attr(df, "month_src") <- mcol[1]
  } else {
    df$month_fe <- factor(0L)  # harmless single level if no month info
    attr(df, "month_src") <- NA_character_
  }
  df
}

# 2) Quick counts for sanity
timing_counts <- function(df){
  df %>%
    mutate(tt = if_else(treated_dummy==1 & post2==1, "Treated-POST",
                        if_else(treated_dummy==1 & post2==0, "Treated-PRE",
                                if_else(treated_dummy==0 & post2==1, "Control-POST","Control-PRE")))) %>%
    count(tt, cycle_fe, month_fe, name = "n")
}

# 3) Utility: extract ATT line with clustered p
att_row <- function(m){
  ct <- as.data.frame(fixest::coeftable(m))
  ct$term <- rownames(ct)
  row <- ct[ct$term=="post2:treated_dummy", , drop=FALSE]
  tibble::tibble(
    term = "post2:treated_dummy",
    estimate = row$Estimate,
    se = row$`Std. Error`,
    p_cluster = row$`Pr(>|t|)`
  )
}

## ---- Build timing FEs on the turnout fit data ----
pm_turn_time <- add_timing_fes(pm_turn_fit)

cat("\n[Timing FE check] Source of month variable: ",
    attr(pm_turn_time, "month_src") %||% "none (single-level placeholder)", "\n", sep="")
print(timing_counts(pm_turn_time))

## ========== Variant A: Add cycle + month FEs (replace year FE) ==========
## Keep city FE and city-specific linear trend; drop year FE to avoid collinearity.
m_may_turn_timeFE <- fixest::feols(
  turnout_pct ~ post2:treated_dummy + i(city_fe, trend),
  data    = pm_turn_time,
  fixef   = c("city_fe","cycle_fe","month_fe"),
  weights = ~ w_eb,
  cluster = ~ city_fe
)

cat("\n[Variant A] Turnout with city trends + {cycle, month} FEs (no year FE)\n")
print(summary(m_may_turn_timeFE, cluster = "city_fe"))
print(att_row(m_may_turn_timeFE))

## ========== Variant B: Matched-timing subsample ==========
## Restrict to the cycle(s)/month(s) that treated cities actually use in POST.
treated_post <- pm_turn_time %>% filter(treated_dummy==1, post2==1)
keep_cycles  <- unique(treated_post$cycle_fe)
keep_months  <- unique(treated_post$month_fe)

pm_turn_matched <- pm_turn_time %>%
  filter(cycle_fe %in% keep_cycles, month_fe %in% keep_months)

cat("\n[Matched-timing] Keeping cycles:", paste(levels(keep_cycles), collapse=", "),
    " | months:", paste(levels(keep_months), collapse=", "), "\n")

# If the subset is tiny, warn but still attempt estimation.
cat("[Matched-timing] Rows kept:", nrow(pm_turn_matched),
    " | Cities:", dplyr::n_distinct(pm_turn_matched$city_fe), "\n")

## For comparability with your baseline, we RE-introduce year FE here.
m_may_turn_matched <- fixest::feols(
  turnout_pct ~ post2:treated_dummy + i(city_fe, trend),
  data    = pm_turn_matched,
  fixef   = c("city_fe","year_fe"),
  weights = ~ w_eb,
  cluster = ~ city_fe
)

cat("\n[Variant B] Turnout on matched-timing subsample (baseline FEs: city + year, with city trends)\n")
print(summary(m_may_turn_matched, cluster = "city_fe"))
print(att_row(m_may_turn_matched))

## ---- Optionally save models ----
saveRDS(m_may_turn_timeFE,  file.path(out_models, "08_CITYTREND_mayor_turnout_timingFE.rds"))
saveRDS(m_may_turn_matched, file.path(out_models, "08_CITYTREND_mayor_turnout_matchedTiming.rds"))



