# Data Access & Availability

This project integrates several restricted and public datasets to construct the city–year analysis panel.

## Primary Sources

1. Busso, Gregory & Kline (2013)  
   “Assessing the Incidence and Efficiency of a Prominent Place-Based Policy,” American Economic Review, 103(2), 897–947.  
   - Data obtained from the AER replication archive (Harvard Dataverse).  
   - Requires approval for access to HUD-linked identifiers.  
   - Used here only for constructing city-level treatment flags (EZ Round I designations).

2. American Local Government Elections Database (ALEGD)
   de Benedictis-Kessner et al. (2023). OSF: [https://osf.io/mv5e6](https://osf.io/mv5e6)  
   - Public access; aggregated to city-year level.  
   - Used to derive turnout, vote share, and council outcomes.

3. National Historical Geographic Information System (NHGIS) 
   IPUMS NHGIS, Version 19.0. [https://doi.org/10.18128/D050.V19.0](https://doi.org/10.18128/D050.V19.0)  
   - Used for 1990 tract populations, place boundaries, and age structure.

---

## Data Restrictions
The replication dataset from Busso, Gregory & Kline (BGK) includes administrative identifiers provided under data-use agreements with HUD.  
These identifiers cannot be publicly shared.  

### Therefore:
- **No raw `.rds` or `.csv` data** are included in this repository.  
- Scripts reference local paths to `.rds` files which must be recreated following the access steps above.  
- Each script is written to run independently once equivalent datasets are placed in `/data_raw/`.

---

## Recreating the Dataset
Steps (once authorized data access is granted):
1. Download BGK replication package and extract tract-level EZ flags.
2. Merge NHGIS 1990 tract shapefiles and demographics.
3. Integrate ALEGD city-level election data (1989–2021).
4. Run `02_city_ez_flags.R` → `04_build_election_panel.R` to rebuild analysis panel.

---

## References
Busso, M., Gregory, J., & Kline, P. (2013). Assessing the Incidence and Efficiency of a Prominent
Place Based Policy. American Economic Review, 103(2), 897–947. https://doi.org/10.1257/AER.103.2.897

de Benedictis-Kessner, J., Lee, D., Velez, Y. R., & Warshaw, C. (2023, May 16). American Local
Government Elections Database. Retrieved from osf.io/mv5e6

Steven Manson, Jonathan Schroeder, David Van Riper, Katherine Knowles, Tracy Kugler, Finn
Roberts, and Steven Ruggles. IPUMS National Historical Geographic Information System:
Version 19.0 [dataset]. Minneapolis, MN: IPUMS. 2024. http://doi.org/10.18128/D050.V19.0 
