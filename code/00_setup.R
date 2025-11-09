###############################################################################
# 00_setup.R ── Initial Project Set-Up --------------------------------------------
###############################################################################

## 0. Anchor the project root --------------------------------------------------
here::i_am("code/00_setup.R")     # <-- THIS ONE LINE is the anchor

## 1. Packages: Install/Load ----
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("fs", quietly = TRUE)) install.packages("fs")
library(here)      # project-root helper- for file paths, works no matter project location
library(fs)        # file path management
library(tidyverse) # tidy data wrangling & for the helper functions below

### Helper Functions ###

## 2. Split-Join Helper: Fractional Block-2-Tract Crosswalk-----
    #fraction-join... so split rows/ duplicates blocks if they map to >1 tract)

split_join <- function(blks, cbrf){
  
  cbrf2 <- cbrf %>% 
    group_by(blockid00) %>% 
    mutate(n_split = n()) %>%
    ungroup()
  
  blks %>% 
    left_join(cbrf2, by = c("blockid2000" = "blockid00"),
              relationship = "many-to-many") %>% 
    
    mutate(weight_block = if_else(is.na(n_split), 1, 1 / n_split)  
    )
}


## 3. Population-Share Helper
    #a pop-weighted share for city-level later on 
pop_share <- function(df, flag){
  with(df, sum(pop1990 * (!!flag)) / sum(pop1990))}
  
# optional helper: header printer
h <- function(txt) cat("\n", strrep("=", 60), "\n", txt, "\n",
                       strrep("=", 60), "\n", sep = "")