###############################################################################
# 01_treatment_flags.R ── Build 1990-Tract Treatment Flags  ------------------
###############################################################################

  # Builds 1990-tract treatment & control flags exactly the BGK way
  # Inputs :  data_raw/0_admin_raw/blocklevel_admdata.dta
  #           data_raw/0_admin_raw/CBRF.dta
  # Output :  data_clean/1_intermediate/tracts90_flags.rds
  # ---------------------------------------------------------------------------

## 0. Setup ----
source(here::here("code", "00_setup.R"))    
library(haven)                           

## 1. Load Raw Data ----
# These come from BGK FOIA: block-level EZ mapping & “control” (CBRF) files
admin_path <- here("data_clean", "0_admin_raw")
blocks <- read_dta(path(admin_path, "blocklevel_admdata.dta"))
cbrf   <- read_dta(path(admin_path, "CBRF.dta"))


## 2. Inspect for Split Blocks  ----
              # Why? Some 2000 blocks are split across more than one 1990 tract.
h("How many 2000 blocks are split across multiple 1990 tracts?")
cbrf %>% 
  count(blockid00) %>% 
    filter(n > 1) %>% 
    summarise(n_split_blocks = n()) %>% 
    print()


## 3. Crosswalk: Fractional Join, Blocks → 1990 Tracts --
  ## This creates a "split-rows" table (blocks_x) with correct fractional weights.

h("Merging (split-rows) and assigning fractional weights")

blocks_x <- split_join(blocks, cbrf) %>%
  distinct(blockid2000, tract1990, .keep_all = TRUE) %>%      # drops exact (block,tract) duplicate pairs
  group_by(blockid2000) %>% # Quick fix to normalize weights so they sum to 1 per 2000-block:
  mutate(
    weight_block = if_else(is.na(weight_block), 1, weight_block), # Fill in unmatched blocks (so no NAs)
    weight_block = weight_block / sum(weight_block, na.rm = TRUE) # Recompute equal-split weights
  ) %>%
  ungroup()


## 4.  Diagnostics: Check Output ----
# (A) Examples of split (fractional) blocks
h("Show some split examples")

blocks_x %>% 
    filter(weight_block < 1) %>% 
    slice_head(n = 6) %>% 
    select(blockid2000, tract1990, weight_block) %>% 
    print()

# (B) Number of duplicated rows and affected blocks
  # aka # of duplicated rows w/ fractional weight that equals split-block count printed earlier

h("Check how many blocks split and how many rows are duplicated")

dup_rows <- blocks_x %>% filter(weight_block < 1) %>% nrow()
  dup_blocks <- blocks_x %>% filter(weight_block < 1) %>% distinct(blockid2000) %>% nrow()
  cat("Rows with fractional weight:", dup_rows, 
      "\nBlocks affected:", dup_blocks, "\n")

# (C) Weights sum to 1 per 2000 block: Sum-by-Block: 
h("Weights: min/max per block should be 1/1")

      # since each 2000 block now has total weight = 1
        #...we should expect: min_w == 1  &  max_w == 1;max_rows should equal max splits (2 or 3)
blocks_x %>%
  group_by(blockid2000) %>%
  summarise(sum_w = sum(weight_block), n_rows = n()) %>% 
  summarise(min_w = min(sum_w), max_w = max(sum_w),max_rows = max(n_rows)) %>%
  print()


## 5. Label Each Block (safest, conservative labelling) ----
h("Create block labels (ez1, app, future, other)")

  blocks_x <- blocks_x %>%
    mutate(
      ez1    = (assigned_tract == 1),
      app    = (application_tract == 1) & !ez1,
      future = (!ez1) & (program1999 == 1 | program2001 == 1),
      # mutually exclusive label for quick inspection
      label = case_when(
        ez1              ~ "ez1",
        app & future     ~ "future",   # applicant that later wins
        app              ~ "app",      # applicant that never wins
        future           ~ "future",   # non-applicant that later wins (rare)
        TRUE             ~ "other"
      )
    )


## 6.  Aggregate to 1990 Tract Level -------------------------------------------
  # collapse to tract (1990 tracts), create BOTH boolean flags and fractional coverage
  
  h("Aggregate to tract: mutually exclusive booleans + fractional coverage")
  
# (a) boolean “any” flags 
tract_any <- blocks_x %>%
  group_by(tract1990) %>%
  summarise(
    any_ez1    = any(label == "ez1"),
    any_app    = any(label == "app"),
    any_future = any(label == "future"),
    .groups = "drop"
  ) %>%
  mutate(
    # BGK precedence: ez1 takes priority over app/future
    ez1    = any_ez1,
    app    = any_app    & !ez1, # near-miss only if not an EZ1 tract
    future = any_future & !ez1  # future never overrides an EZ1 tract
  ) %>%
  select(tract1990, ez1, app, future)

# (b) Fractional coverage of each tract by label (0..1), using block weights --FOR ROBUSTNESS ONLY.
tract_frac <- blocks_x %>%
  group_by(tract1990) %>%
  summarise(
    w_total     = sum(weight_block, na.rm = TRUE),
    ez1_frac    = sum(weight_block * as.integer(label == "ez1"),    na.rm = TRUE) / w_total,
    app_frac    = sum(weight_block * as.integer(label == "app"),    na.rm = TRUE) / w_total,
    future_frac = sum(weight_block * as.integer(label == "future"), na.rm = TRUE) / w_total,
    .groups = "drop"
  )

# (c) combine; **only** set treat_yr for EZ Round I winners here
tracts90 <- tract_any %>%
  left_join(tract_frac, by = "tract1990") %>%
  mutate(
    pop1990  = 1L,  # placeholder (replaced in Script 02)
    treat_yr = case_when(
      ez1    ~ 1994L,
      TRUE   ~ NA_integer_
    )
  ) %>%
  select(tract1990, ez1, app, future, ez1_frac, app_frac, future_frac, pop1990, treat_yr)

## 7. Diagnostics: Tract Counts & Fraction ranges ----
## quick counts for log ----
h("Counts after aggregation")
tracts90 %>% summarise(
  EZ1      = sum(ez1),
  NearMiss = sum(app),
  Future   = sum(future),
  Others   = sum(!(ez1 | app | future))
) %>% print()

h("Counts after aggregation (booleans) + fractional sanity checks")
tracts90 %>%
  summarise(
    EZ1      = sum(ez1),
    NearMiss = sum(app),
    Future   = sum(future),
    frac_min = min(ez1_frac, app_frac, future_frac, na.rm = TRUE),
    frac_max = max(ez1_frac, app_frac, future_frac, na.rm = TRUE)
  ) %>% print()

## 8. Save Output ----
write_rds(tracts90, here("data_clean", "1_intermediate", "tracts90_flags.rds"))
