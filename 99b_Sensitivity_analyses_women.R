################################################################################
# SENSITIVITY ANALYSIS:
# NFHS4 -> NFHS5 Pre-lockdown -> NFHS5 Post-lockdown
################################################################################

library(haven)
library(survey)
library(tidyverse)

options(survey.lonely.psu = "adjust")

################################################################################
# NFHS5 CHILDREN DATA
################################################################################

child_keep5 <- c(
  "V001","V002","V003","B16",
  "V005","V021","V022","V023",
  "V191A","V190A","V024",
  "B8","B4","B5",
  "HW1",
  "HW2","HW3",
  "HW72","HW57"
)

hh_keep5 <- c(
  "HV001","HV002","HVIDX",
  "SHWEIGHT","HV005","HV021","HV022","HV023",
  "HV006","HV007",
  "HV270","HV270A","HV271","HV271A",
  "HV024","HV025",
  "SH47","SH49"
)

nfhs_child_sub5 <- read_sav(
  "./IAKR7EFL.SAV",
  col_select = all_of(child_keep5)
)

nfhs_hh_child5 <- read_sav(
  "./IAPR7EFL.SAV",
  col_select = all_of(hh_keep5)
)

nfhs_hh_for_child5 <- nfhs_hh_child5 %>%
  rename(
    V001 = HV001,
    V002 = HV002,
    V003 = HVIDX
  )

nfhs_child_analysis <- nfhs_child_sub5 %>%
  left_join(nfhs_hh_for_child5, by = c("V001","V002","V003")) %>%
  rename(
    month_of_interview = HV006,
    year_of_interview  = HV007,
    
    wi_factorscore = HV271,
    wi_factorscore_urbanrural = HV271A,
    wi_combined = HV270,
    wi_combined_urbanrural = HV270A,
    
    child_age_years  = B8,
    child_age_months = HW1,
    child_alive      = B5,
    child_sex        = B4,
    
    child_weight_kgs = HW2,
    child_height_cms = HW3,
    
    child_WFH_SD_WHO   = HW72,
    child_anemia_level = HW57,
    
    child_weight = V005,
    child_psu = V021,
    child_strata = V022,
    child_sample_domain = V023,
    child_wi_factorscore = V191A,
    child_wi_quintile = V190A,
    
    hh_state_weight = SHWEIGHT,
    hh_weight = HV005,
    hh_psu = HV021,
    hh_strata = HV022,
    hh_sample_domain = HV023,
    
    child_region = V024,
    household_region = HV024,
    urban_rural = HV025,
    religion = SH47,
    caste = SH49
  )

################################################################################
# CLEAN NFHS5 VARIABLES
################################################################################

nfhs_child_analysis <- nfhs_child_analysis %>%
  mutate(
    child_weight_kgs = ifelse(child_weight_kgs >= 9990, NA, child_weight_kgs / 10),
    child_height_cms = ifelse(child_height_cms >= 9990, NA, child_height_cms / 10),
    child_WFH_SD_WHO = ifelse(child_WFH_SD_WHO >= 9990, NA, child_WFH_SD_WHO / 100),
    
    hh_state_weight = hh_state_weight / 1e6,
    child_weight = child_weight / 1e6,
    hh_weight = hh_weight / 1e6,
    
    wi_combined_urbanrural = as.numeric(as.character(wi_combined_urbanrural))
  )

nfhs_child_analysis <- nfhs_child_analysis %>%
  filter(child_alive == 1)

################################################################################
# NFHS5 OUTCOMES + PRE/POST VARIABLE
################################################################################

nfhs_child_analysis <- nfhs_child_analysis %>%
  mutate(
    child_wasting = case_when(
      is.na(child_WFH_SD_WHO) ~ NA_real_,
      child_WFH_SD_WHO < -2 ~ 1,
      TRUE ~ 0
    ),
    child_wasting = factor(
      child_wasting,
      levels = c(0, 1),
      labels = c("No wasting", "Wasting")
    ),
    
    child_anemia = case_when(
      is.na(child_anemia_level) ~ NA_real_,
      child_anemia_level %in% c(1, 2, 3) ~ 1,
      TRUE ~ 0
    ),
    child_anemia = factor(
      child_anemia,
      levels = c(0, 1),
      labels = c("Not anemic", "Anemic")
    )
  )

nfhs_child_analysis$date_of_interview <- as.Date(
  paste(
    nfhs_child_analysis$year_of_interview,
    nfhs_child_analysis$month_of_interview,
    "01",
    sep = "-"
  )
)

nfhs_child_analysis$postcovid <- ifelse(
  nfhs_child_analysis$date_of_interview < as.Date("2020-04-01"),
  "Pre-lockdown",
  ifelse(
    nfhs_child_analysis$date_of_interview >= as.Date("2020-06-01"),
    "Post-lockdown",
    NA
  )
)

################################################################################
# 1) LOAD NFHS4 CHILDREN + HOUSEHOLD BIOMARKER DATA
################################################################################

child_keep4 <- c(
  "V001", "V002", "V003", "B16",
  "V005", "V021", "V022", "V023",
  "V191", "V190", "V024",
  "B8", "B4", "B5",
  "HW1",
  "HW2", "HW3",
  "HW72", "HW57"
)

hh_keep4 <- c(
  "HV001", "HV002", "HVIDX",
  "SHV005", "HV005", "HV021", "HV022", "HV023",
  "HV006", "HV007",
  "HV270", "HV271",
  "HV024", "HV025",
  "SH34", "SH36"
)

nfhs_child_sub4 <- read_sav(
  "./99. DHS4/IAKR74FL.SAV",
  col_select = all_of(child_keep4)
)

nfhs_hh_child4 <- read_sav(
  "./99. DHS4/IAPR74FL.SAV",
  col_select = all_of(hh_keep4)
)

################################################################################
# 2) MERGE NFHS4 CHILDREN + HOUSEHOLD DATA
################################################################################

nfhs_hh_for_child4 <- nfhs_hh_child4 %>%
  rename(
    V001 = HV001,
    V002 = HV002,
    V003 = HVIDX
  )

nfhs_child_analysis4 <- nfhs_child_sub4 %>%
  left_join(nfhs_hh_for_child4, by = c("V001", "V002", "V003")) %>%
  rename(
    month_of_interview = HV006,
    year_of_interview  = HV007,
    
    wi_combined = HV270,
    wi_factorscore = HV271,
    
    child_age_years  = B8,
    child_age_months = HW1,
    child_alive      = B5,
    child_sex        = B4,
    
    child_weight_kgs = HW2,
    child_height_cms = HW3,
    child_WFH_SD_WHO = HW72,
    child_anemia_level = HW57,
    
    child_weight = V005,
    child_psu = V021,
    child_strata = V022,
    child_sample_domain = V023,
    child_wi_factorscore = V191,
    child_wi_quintile = V190,
    
    hh_state_weight = SHV005,
    hh_weight = HV005,
    hh_psu = HV021,
    hh_strata = HV022,
    hh_sample_domain = HV023,
    
    child_region = V024,
    household_region = HV024,
    urban_rural = HV025,
    religion = SH34,
    caste = SH36
  )

################################################################################
# 3) CLEAN NFHS4 VARIABLES
################################################################################

nfhs_child_analysis4 <- nfhs_child_analysis4 %>%
  mutate(
    child_weight_kgs = ifelse(child_weight_kgs >= 9990, NA, child_weight_kgs / 10),
    child_height_cms = ifelse(child_height_cms >= 9990, NA, child_height_cms / 10),
    child_WFH_SD_WHO = ifelse(child_WFH_SD_WHO >= 9990, NA, child_WFH_SD_WHO / 100),
    
    wi_factorscore = wi_factorscore / 100000,
    wi_combined_urbanrural = as.numeric(as.character(wi_combined)),
    
    hh_state_weight = hh_state_weight / 1e6,
    child_weight = child_weight / 1e6,
    hh_weight = hh_weight / 1e6
  )

################################################################################
# 4) CREATE NFHS4 OUTCOMES
################################################################################

nfhs_child_analysis4 <- nfhs_child_analysis4 %>%
  mutate(
    child_wasting = case_when(
      is.na(child_WFH_SD_WHO) ~ NA_real_,
      child_WFH_SD_WHO < -2 ~ 1,
      TRUE ~ 0
    ),
    
    child_wasting = factor(
      child_wasting,
      levels = c(0, 1),
      labels = c("No wasting", "Wasting")
    ),
    
    child_anemia = case_when(
      is.na(child_anemia_level) ~ NA_real_,
      child_anemia_level %in% c(1, 2, 3) ~ 1,
      TRUE ~ 0
    ),
    
    child_anemia = factor(
      child_anemia,
      levels = c(0, 1),
      labels = c("Not anemic", "Anemic")
    )
  )

################################################################################
# 5) STATE CODEBOOKS
################################################################################

phase2_states <- c(
  "Uttar Pradesh",
  "Haryana",
  "Punjab",
  "Rajasthan",
  "Uttarakhand",
  "Chhattisgarh",
  "Madhya Pradesh",
  "Jharkhand",
  "Odisha",
  "Arunachal Pradesh",
  "Tamil Nadu",
  "Chandigarh",
  "Delhi",
  "Puducherry"
)

state_codebook_nfhs4 <- tibble::tribble(
  ~region_code4, ~state_name,
  3,  "Arunachal Pradesh",
  6,  "Chandigarh",
  7,  "Chhattisgarh",
  12, "Haryana",
  15, "Jharkhand",
  19, "Madhya Pradesh",
  25, "Delhi",
  26, "Odisha",
  27, "Puducherry",
  28, "Punjab",
  29, "Rajasthan",
  31, "Tamil Nadu",
  33, "Uttar Pradesh",
  34, "Uttarakhand"
)

state_codebook_nfhs5 <- tibble::tribble(
  ~region_code5, ~state_name,
  3,  "Punjab",
  4,  "Chandigarh",
  5,  "Uttarakhand",
  6,  "Haryana",
  7,  "Delhi",
  8,  "Rajasthan",
  9,  "Uttar Pradesh",
  12, "Arunachal Pradesh",
  20, "Jharkhand",
  21, "Odisha",
  22, "Chhattisgarh",
  23, "Madhya Pradesh",
  33, "Tamil Nadu",
  34, "Puducherry"
)

################################################################################
# 6) ADD STATE LABELS TO NFHS4 AND NFHS5
################################################################################

nfhs_child_analysis4 <- nfhs_child_analysis4 %>%
  mutate(child_region_num = as.numeric(child_region)) %>%
  left_join(state_codebook_nfhs4, by = c("child_region_num" = "region_code4"))

nfhs_child_analysis <- nfhs_child_analysis %>%
  mutate(child_region_num = as.numeric(child_region)) %>%
  left_join(state_codebook_nfhs5, by = c("child_region_num" = "region_code5"))

################################################################################
# 7) CREATE NFHS4 PHASE 2 DATASET
################################################################################

nfhs_child_phase2_14_nfhs4 <- nfhs_child_analysis4 %>%
  filter(state_name %in% phase2_states)

################################################################################
# 8) CHECK NFHS5 PHASE 2 STATES BY PRE/POST PERIOD
################################################################################

nfhs_child_analysis %>%
  filter(state_name %in% phase2_states) %>%
  count(state_name, postcovid) %>% print(n=100)

################################################################################
# 9) PREPARE THREE TIMEPOINT DATASETS
################################################################################

## chandigarh has no observations pre-lockdown so exclude it

nfhs4_ready <- nfhs_child_phase2_14_nfhs4 %>%
  filter(state_name != "Chandigarh") %>%
  mutate(
    timepoint = "NFHS4",
    wasting_bin = child_wasting == "Wasting",
    anemia_bin = child_anemia == "Anemic"
  )

nfhs5_pre_ready <- nfhs_child_analysis %>%
  filter(
    state_name %in% phase2_states,
    state_name != "Chandigarh",
    postcovid == "Pre-lockdown"
  ) %>%
  mutate(
    timepoint = "NFHS5 Pre-lockdown",
    wasting_bin = child_wasting == "Wasting",
    anemia_bin = child_anemia == "Anemic"
  )

nfhs5_post_ready <- nfhs_child_analysis %>%
  filter(
    state_name %in% phase2_states,
    state_name != "Chandigarh",
    postcovid == "Post-lockdown"
  ) %>%
  mutate(
    timepoint = "NFHS5 Post-lockdown",
    wasting_bin = child_wasting == "Wasting",
    anemia_bin = child_anemia == "Anemic"
  )

################################################################################
# 10) COMBINE DATASETS
################################################################################

common_vars <- c(
  "timepoint",
  "state_name",
  "child_psu",
  "child_weight",
  "child_strata",
  "hh_state_weight",
  "wi_combined_urbanrural",
  "child_age_months",
  "urban_rural",
  "religion",
  "caste",
  "child_wasting",
  "child_anemia",
  "wasting_bin",
  "anemia_bin"
)

nfhs4_ready <- nfhs4_ready %>%
  mutate(
    religion = as.numeric(haven::zap_labels(religion)),
    caste = as.numeric(haven::zap_labels(caste)),
    urban_rural = as.numeric(haven::zap_labels(urban_rural)),
    wi_combined_urbanrural = as.numeric(haven::zap_labels(wi_combined_urbanrural))
  )

nfhs5_pre_ready <- nfhs5_pre_ready %>%
  mutate(
    religion = as.numeric(haven::zap_labels(religion)),
    caste = as.numeric(haven::zap_labels(caste)),
    urban_rural = as.numeric(haven::zap_labels(urban_rural)),
    wi_combined_urbanrural = as.numeric(haven::zap_labels(wi_combined_urbanrural))
  )

nfhs5_post_ready <- nfhs5_post_ready %>%
  mutate(
    religion = as.numeric(haven::zap_labels(religion)),
    caste = as.numeric(haven::zap_labels(caste)),
    urban_rural = as.numeric(haven::zap_labels(urban_rural)),
    wi_combined_urbanrural = as.numeric(haven::zap_labels(wi_combined_urbanrural))
  )

nfhs_combined_3tp <- bind_rows(
  nfhs4_ready %>% select(all_of(common_vars)),
  nfhs5_pre_ready %>% select(all_of(common_vars)),
  nfhs5_post_ready %>% select(all_of(common_vars))
) %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = c("NFHS4", "NFHS5 Pre-lockdown", "NFHS5 Post-lockdown")
    )
  )

################################################################################
# 11) SURVEY DESIGN
################################################################################

des_3tp <- svydesign(
  ids = ~child_psu,
  strata = ~interaction(child_strata, timepoint),
  weights = ~child_weight,
  data = nfhs_combined_3tp,
  nest = TRUE
)

################################################################################
# 12) SAMPLE SIZE CHECK
################################################################################

nfhs_combined_3tp %>%
  count(timepoint)

################################################################################
# 13) OVERALL WEIGHTED PREVALENCE
################################################################################

wasting_prev_3tp <- svyby(
  ~wasting_bin,
  ~timepoint,
  des_3tp,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

anemia_prev_3tp <- svyby(
  ~anemia_bin,
  ~timepoint,
  des_3tp,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

################################################################################
# 14) WEIGHTED PREVALENCE BY WEALTH QUINTILE
################################################################################

wasting_prev_q_3tp <- svyby(
  ~wasting_bin,
  ~timepoint + wi_combined_urbanrural,
  des_3tp,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

anemia_prev_q_3tp <- svyby(
  ~anemia_bin,
  ~timepoint + wi_combined_urbanrural,
  des_3tp,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

################################################################################
# 15) CLEAN PERCENT TABLES
################################################################################

wasting_prev_3tp_pct <- wasting_prev_3tp %>%
  as_tibble() %>%
  mutate(
    prevalence_pct = wasting_binTRUE * 100,
    se_pct = `se.wasting_binTRUE` * 100,
    ci_l_pct = `ci_l.wasting_binTRUE` * 100,
    ci_u_pct = `ci_u.wasting_binTRUE` * 100
  ) %>%
  select(timepoint, prevalence_pct, se_pct, ci_l_pct, ci_u_pct)

anemia_prev_3tp_pct <- anemia_prev_3tp %>%
  as_tibble() %>%
  mutate(
    prevalence_pct = anemia_binTRUE * 100,
    se_pct = `se.anemia_binTRUE` * 100,
    ci_l_pct = `ci_l.anemia_binTRUE` * 100,
    ci_u_pct = `ci_u.anemia_binTRUE` * 100
  ) %>%
  select(timepoint, prevalence_pct, se_pct, ci_l_pct, ci_u_pct)

wasting_prev_q_3tp_pct <- wasting_prev_q_3tp %>%
  as_tibble() %>%
  mutate(
    prevalence_pct = wasting_binTRUE * 100,
    se_pct = `se.wasting_binTRUE` * 100,
    ci_l_pct = `ci_l.wasting_binTRUE` * 100,
    ci_u_pct = `ci_u.wasting_binTRUE` * 100
  ) %>%
  select(
    timepoint,
    wi_combined_urbanrural,
    prevalence_pct,
    se_pct,
    ci_l_pct,
    ci_u_pct
  )

anemia_prev_q_3tp_pct <- anemia_prev_q_3tp %>%
  as_tibble() %>%
  mutate(
    prevalence_pct = anemia_binTRUE * 100,
    se_pct = `se.anemia_binTRUE` * 100,
    ci_l_pct = `ci_l.anemia_binTRUE` * 100,
    ci_u_pct = `ci_u.anemia_binTRUE` * 100
  ) %>%
  select(
    timepoint,
    wi_combined_urbanrural,
    prevalence_pct,
    se_pct,
    ci_l_pct,
    ci_u_pct
  )

################################################################################
# 16) PRINT PREVALENCE TABLES
################################################################################

wasting_prev_3tp_pct
anemia_prev_3tp_pct

wasting_prev_q_3tp_pct
anemia_prev_q_3tp_pct

################################################################################
# 17A) CREATE NFHS5 PRE VS POST DESIGN FOR ORs
################################################################################

nfhs5_phase2_only <- nfhs_combined_3tp %>%
  filter(timepoint %in% c("NFHS5 Pre-lockdown", "NFHS5 Post-lockdown")) %>%
  mutate(
    timepoint = relevel(
      factor(timepoint),
      ref = "NFHS5 Pre-lockdown"
    ),
    urban_rural = factor(urban_rural),
    religion = factor(religion),
    caste = factor(caste)
  )

des_nfhs5_phase2_only <- svydesign(
  ids = ~child_psu,
  strata = ~interaction(child_strata, timepoint),
  weights = ~child_weight,
  data = nfhs5_phase2_only,
  nest = TRUE
)

or_timepoint <- function(fit) {
  term <- "timepointNFHS5 Post-lockdown"
  est <- coef(fit)[term]
  ci <- confint(fit)[term, ]
  
  c(
    OR = exp(est),
    CI_low = exp(ci[1]),
    CI_high = exp(ci[2])
  )
}

################################################################################
# 17B) ADJUSTED NFHS5 PRE VS POST ODDS RATIOS
################################################################################

fit_wasting_nfhs5 <- svyglm(
  wasting_bin ~ timepoint + child_age_months + caste + religion,
  design = des_nfhs5_phase2_only,
  family = quasibinomial()
)

fit_anemia_nfhs5 <- svyglm(
  anemia_bin ~ timepoint + child_age_months + caste + religion,
  design = des_nfhs5_phase2_only,
  family = quasibinomial()
)

OR_wasting_nfhs5 <- or_timepoint(fit_wasting_nfhs5)
OR_anemia_nfhs5 <- or_timepoint(fit_anemia_nfhs5)

OR_wasting_nfhs5
OR_anemia_nfhs5

################################################################################
# 18) ADJUSTED NFHS5 PRE VS POST ODDS RATIOS BY WEALTH QUINTILE
################################################################################

fit_wasting_q_nfhs5 <- lapply(1:5, function(q) {
  svyglm(
    wasting_bin ~ timepoint + child_age_months + caste + religion,
    design = subset(des_nfhs5_phase2_only,
                    wi_combined_urbanrural == q),
    family = quasibinomial()
  )
})

OR_wasting_q_nfhs5 <- do.call(
  rbind,
  lapply(fit_wasting_q_nfhs5, or_timepoint)
)

rownames(OR_wasting_q_nfhs5) <- paste0("Q", 1:5)

fit_anemia_q_nfhs5 <- lapply(1:5, function(q) {
  svyglm(
    anemia_bin ~ timepoint + child_age_months + caste + religion,
    design = subset(des_nfhs5_phase2_only,
                    wi_combined_urbanrural == q),
    family = quasibinomial()
  )
})

OR_anemia_q_nfhs5 <- do.call(
  rbind,
  lapply(fit_anemia_q_nfhs5, or_timepoint)
)

rownames(OR_anemia_q_nfhs5) <- paste0("Q", 1:5)

OR_wasting_q_nfhs5
OR_anemia_q_nfhs5

################################################################################

################################################################################
# 18b) Adjusted (marginally standardized) prevalence — restricted 13-state sample
################################################################################

# Adjusted prevalence (%) pre and post, via marginal standardization,
# for one design (overall or a quintile subset).
adj_prev_sub_sens <- function(formula, dq, wtvar = "child_weight") {
  fq <- svyglm(formula, design = dq, family = quasibinomial())
  mf <- model.frame(fq)
  w_named <- dq$variables[[wtvar]]
  names(w_named) <- rownames(dq$variables)
  ww <- as.numeric(w_named[rownames(mf)])
  d_pre  <- mf; d_pre$timepoint  <- factor("NFHS5 Pre-lockdown",  levels = levels(mf$timepoint))
  d_post <- mf; d_post$timepoint <- factor("NFHS5 Post-lockdown", levels = levels(mf$timepoint))
  p_pre  <- as.numeric(predict(fq, d_pre,  type = "response"))
  p_post <- as.numeric(predict(fq, d_post, type = "response"))
  c(pre  = 100 * weighted.mean(p_pre,  ww, na.rm = TRUE),
    post = 100 * weighted.mean(p_post, ww, na.rm = TRUE))
}

# --- Overall
adj_was_overall_sens <- adj_prev_sub_sens(
  wasting_bin ~ timepoint + child_age_months + caste + religion, des_nfhs5_phase2_only)
adj_anemia_overall_sens <- adj_prev_sub_sens(
  anemia_bin ~ timepoint + child_age_months + caste + religion, des_nfhs5_phase2_only)

# --- By wealth quintile
adj_was_by_q_sens <- t(sapply(1:5, function(q)
  adj_prev_sub_sens(wasting_bin ~ timepoint + child_age_months + caste + religion,
                    subset(des_nfhs5_phase2_only, wi_combined_urbanrural == q))))
rownames(adj_was_by_q_sens) <- paste0("Q", 1:5)

adj_anemia_by_q_sens <- t(sapply(1:5, function(q)
  adj_prev_sub_sens(anemia_bin ~ timepoint + child_age_months + caste + religion,
                    subset(des_nfhs5_phase2_only, wi_combined_urbanrural == q))))
rownames(adj_anemia_by_q_sens) <- paste0("Q", 1:5)

# --- Print
cat("\n=== RESTRICTED 13-STATE ADJUSTED PREVALENCE (Wasting, children) ===\n")
cat("Overall pre/post:", round(adj_was_overall_sens[["pre"]],1),
    round(adj_was_overall_sens[["post"]],1), "\n")
cat("By quintile (pre, post):\n"); print(round(adj_was_by_q_sens, 1))

cat("\n=== RESTRICTED 13-STATE ADJUSTED PREVALENCE (Anemia, children) ===\n")
cat("Overall pre/post:", round(adj_anemia_overall_sens[["pre"]],1),
    round(adj_anemia_overall_sens[["post"]],1), "\n")
cat("By quintile (pre, post):\n"); print(round(adj_anemia_by_q_sens, 1))

################################################################################