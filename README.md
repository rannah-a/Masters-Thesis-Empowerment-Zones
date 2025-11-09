# Empowering Whom? Evaluating the Political Consequences of Federal Empowerment Zones

Replication-based empirical analysis evaluating how U.S. Federal Empowerment Zones (EZs) of the 1990s influenced local electoral outcomes.  
Built from city-level election data (ALEGD), socioeconomic covariates (BGK replication files and NHGIS), and quasi-experimental difference-in-differences methods.

---

##  Abstract
Empowerment Zones were one of the largest U.S. place-based programs, yet their political effects remain unexamined. This project constructs a city–year election panel (1989–2021) merging American Local Government Elections Database (ALEGD) with treatment exposure derived from Busso, Gregory & Kline (2013) and NHGIS 1990 Census data. Using weighted two-way fixed-effects and event-study models, it estimates the political consequences of federal urban investment.

---

##  Methodology Overview
- Design: Difference-in-Differences (with city & year FE + city trends)  
- Weighting: Entropy balancing of treated and donor cities  
- Event Study: Sun & Abraham (2021) staggered-treatment estimator  
- Inference: Clustered SEs (city-level) and wild bootstrap checks  
- Software: R (tidyverse, fixest, ebalance, ggplot2)

---

##  Key Results
- Mayoral Democratic vote share ↑ **0.21 (SE 0.07, p = .003)**
- Mayoral turnout ↓ **6–7 pp (SE 3.65, p ≈ .06)**
- Council outcomes: imprecise, near-zero  
- No pretrend violations detected in event-study checks.

See `/outputs/figures` and `/outputs/tables` for replication figures and tables referenced in the thesis.

---

## Data Availability
Data cannot be redistributed due to licensing restrictions.  
See `/data/README_data.md` for detailed access instructions and citations.

---

## Citation
If you use this work, please cite:

> Ansari, Rannah (2025). Empowering Whom? Evaluating the Political Consequences of Federal Empowerment Zones. Master’s Thesis, University of Barcelona, MSc in Institutions and Political Economy.

---
 
## Author-Rannah Ansari 
MSc in Institutions & Political Economy  
University of Barcelona, 2025  
[LinkedIn](https://linkedin.com/in/rannah-ansari) | [GitHub](https://github.com/rannah-a)
