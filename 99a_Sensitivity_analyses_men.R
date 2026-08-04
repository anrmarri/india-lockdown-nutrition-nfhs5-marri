################################################################################
# SENSITIVITY ANALYSIS:
# NFHS4 -> NFHS5 Pre-lockdown -> NFHS5 Post-lockdown
################################################################################

library(haven)
library(survey)
library(tidyverse)

options(survey.lonely.psu = "adjust")

################################################################################
# NFHS5 MEN DATA
################################################################################

men_keep5 <- c(
  "MV001","MV002","MV003",
  "MV005","MV021","MV022","MV023",
  "MV191A","MV190A","MV024"
)

hh_keep5 <- c(
  "HV001","HV002","HVIDX",
  "SHMWEIGHT", "HV005","HV021","HV022","HV023",
  "HV006","HV007",
  "HV270","HV270A","HV271","HV271A",
  "HB1","HB2","HB3","HB35","HB40","HB56","HB57",
  "HV024","HV025","SH47","SH49"
)

nfhs_men_sub5 <- read_sav(
  "./IAMR7EFL.SAV",
  col_select = all_of(men_keep5)
)

nfhs_hh_men5 <- read_sav(
  "./IAPR7EFL.SAV",
  col_select = all_of(hh_keep5)
)

nfhs_hh_for_men5 <- nfhs_hh_men5 %>%
  rename(
    MV001 = HV001,
    MV002 = HV002,
    MV003 = HVIDX
  )

nfhs_men_analysis <- nfhs_men_sub5 %>%
  left_join(nfhs_hh_for_men5, by = c("MV001","MV002","MV003")) %>%
  rename(
    month_of_interview = HV006,
    year_of_interview  = HV007,
    
    wi_combined = HV270,
    wi_combined_urbanrural = HV270A,
    
    man_age = HB1,
    man_weight_kgs = HB2,
    man_height_cms = HB3,
    man_bmi = HB40,
    man_anemia_level = HB57,
    
    man_weight = MV005,
    man_psu = MV021,
    man_strata = MV022,
    
    man_hh_state_weight = SHMWEIGHT,
    
    man_region = MV024,
    urban_rural = HV025,
    religion = SH47,
    caste = SH49
  )

################################################################################
# CLEAN NFHS5 VARIABLES
################################################################################

nfhs_men_analysis <- nfhs_men_analysis %>%
  mutate(
    man_weight_kgs = ifelse(man_weight_kgs >= 9990, NA, man_weight_kgs / 10),
    man_height_cms = ifelse(man_height_cms >= 9990, NA, man_height_cms / 10),
    man_bmi = ifelse(man_bmi >= 9990, NA, man_bmi / 100),
    
    man_hh_state_weight = man_hh_state_weight / 1e6,
    man_weight = man_weight / 1e6,
    
    wi_combined_urbanrural = as.numeric(as.character(wi_combined_urbanrural))
  )

################################################################################
# NFHS5 OUTCOMES + PRE/POST VARIABLE
################################################################################

nfhs_men_analysis <- nfhs_men_analysis %>%
  mutate(
    man_malnourished = case_when(
      is.na(man_bmi) ~ NA_real_,
      man_bmi < 18.5 ~ 1,
      TRUE ~ 0
    ),

    man_malnourished = factor(
      man_malnourished,
      levels = c(0,1),
      labels = c("Not undernourished", "Undernourished")
    ),

    man_anemia = case_when(
      is.na(man_anemia_level) ~ NA_real_,
      man_anemia_level %in% c(1,2,3) ~ 1,
      TRUE ~ 0
    ),

    man_anemia = factor(
      man_anemia,
      levels = c(0,1),
      labels = c("Not anemic", "Anemic")
    )
  )

nfhs_men_analysis$date_of_interview <- as.Date(
  paste(
    nfhs_men_analysis$year_of_interview,
    nfhs_men_analysis$month_of_interview,
    "01",
    sep = "-"
  )
)

nfhs_men_analysis$postcovid <- ifelse(
  nfhs_men_analysis$date_of_interview < as.Date("2020-04-01"),
  "Pre-lockdown",
  ifelse(
    nfhs_men_analysis$date_of_interview >= as.Date("2020-06-01"),
    "Post-lockdown",
    NA
  )
)

################################################################################
# 1) LOAD NFHS4 MEN + HOUSEHOLD BIOMARKER DATA
################################################################################

men_keep4 <- c(
  "MV001", "MV002", "MV003",
  "MV005", "MV021", "MV022", "MV023",
  "MV191", "MV190", "MV024"
)

hh_keep4 <- c(
  "HV001", "HV002", "HVIDX",
  "SHV005", "HV005", "HV021", "HV022", "HV023",
  "HV006", "HV007",
  "HV270", "HV271",
  "HB1", "HB2", "HB3", "HB35", "HB40", "HB56", "HB57",
  "HV024", "HV025", "SH34", "SH36"
)

nfhs_men_sub4 <- read_sav(
  "./99. DHS4/IAMR74FL.SAV",
  col_select = all_of(men_keep4)
)

nfhs_hh_men4 <- read_sav(
  "./99. DHS4/IAPR74FL.SAV",
  col_select = all_of(hh_keep4)
)

################################################################################
# 2) MERGE NFHS4 MEN + HOUSEHOLD DATA
################################################################################

nfhs_hh_for_men4 <- nfhs_hh_men4 %>%
  rename(
    MV001 = HV001,
    MV002 = HV002,
    MV003 = HVIDX
  )

nfhs_men_analysis4 <- nfhs_men_sub4 %>%
  left_join(nfhs_hh_for_men4, by = c("MV001", "MV002", "MV003")) %>%
  rename(
    month_of_interview = HV006,
    year_of_interview  = HV007,
    
    wi_combined = HV270,
    wi_factorscore = HV271,
    
    man_age = HB1,
    man_weight_kgs = HB2,
    man_height_cms = HB3,
    man_smoking = HB35,
    man_bmi = HB40,
    man_hb_gdl_altsmoking = HB56,
    man_anemia_level = HB57,
    
    man_weight = MV005,
    man_psu = MV021,
    man_strata = MV022,
    man_sample_domain = MV023,
    man_wi_factorscore = MV191,
    man_wi_quintile = MV190,
    
    man_hh_state_weight = SHV005,
    hh_weight = HV005,
    hh_psu = HV021,
    hh_strata = HV022,
    hh_sample_domain = HV023,
    
    man_region = MV024,
    household_region = HV024,
    urban_rural = HV025,
    religion = SH34,
    caste = SH36
  )

################################################################################
# 3) CLEAN NFHS4 VARIABLES
################################################################################

nfhs_men_analysis4 <- nfhs_men_analysis4 %>%
  mutate(
    man_weight_kgs = ifelse(man_weight_kgs >= 9990, NA, man_weight_kgs / 10),
    man_height_cms = ifelse(man_height_cms >= 9990, NA, man_height_cms / 10),
    man_bmi = ifelse(man_bmi >= 9990, NA, man_bmi / 100),
    man_hb_gdl_altsmoking = ifelse(
      man_hb_gdl_altsmoking >= 990,
      NA,
      man_hb_gdl_altsmoking / 10
    ),
    
    wi_factorscore = wi_factorscore / 100000,
    wi_combined_urbanrural = as.numeric(as.character(wi_combined)),
    
    man_hh_state_weight = man_hh_state_weight / 1e6,
    man_weight = man_weight / 1e6,
    hh_weight = hh_weight / 1e6
  )

################################################################################
# 4) CREATE NFHS4 OUTCOMES
################################################################################

nfhs_men_analysis4 <- nfhs_men_analysis4 %>%
  mutate(
    anemia_check = case_when(
      is.na(man_hb_gdl_altsmoking) ~ NA_real_,
      man_hb_gdl_altsmoking < 13.0 ~ 1,
      TRUE ~ 0
    )
  )

table(
  anemia_check = nfhs_men_analysis4$anemia_check,
  anemia_level = nfhs_men_analysis4$man_anemia_level,
  useNA = "ifany"
)

nfhs_men_analysis4 <- nfhs_men_analysis4 %>%
  mutate(
    man_malnourished = case_when(
      is.na(man_bmi) ~ NA_real_,
      man_bmi < 18.5 ~ 1,
      TRUE ~ 0
    ),
    man_malnourished = factor(
      man_malnourished,
      levels = c(0, 1),
      labels = c("Not undernourished", "Undernourished")
    ),
    
    man_anemia = case_when(
      is.na(man_anemia_level) ~ NA_real_,
      man_anemia_level %in% c(1, 2, 3) ~ 1,
      TRUE ~ 0
    ),
    man_anemia = factor(
      man_anemia,
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

nfhs_men_analysis4 <- nfhs_men_analysis4 %>%
  mutate(man_region_num = as.numeric(man_region)) %>%
  left_join(state_codebook_nfhs4, by = c("man_region_num" = "region_code4"))

nfhs_men_analysis <- nfhs_men_analysis %>%
  mutate(man_region_num = as.numeric(man_region)) %>%
  left_join(state_codebook_nfhs5, by = c("man_region_num" = "region_code5"))

################################################################################
# 7) CREATE NFHS4 PHASE 2 DATASET
################################################################################

nfhs_men_phase2_14_nfhs4 <- nfhs_men_analysis4 %>%
  filter(state_name %in% phase2_states)

################################################################################
# 8) CHECK NFHS5 PHASE 2 STATES BY PRE/POST PERIOD
################################################################################

nfhs_men_analysis %>%
  filter(state_name %in% phase2_states) %>%
  count(state_name, postcovid) %>% print(n=100)

################################################################################
# 9) PREPARE THREE TIMEPOINT DATASETS
################################################################################

## chandigarh has no observations pre-lockdown so exclude it

nfhs4_ready <- nfhs_men_phase2_14_nfhs4 %>%
  filter(state_name != "Chandigarh") %>%
  mutate(
    timepoint = "NFHS4",
    undernourished_bin = man_malnourished == "Undernourished",
    anemia_bin = man_anemia == "Anemic"
  )

nfhs5_pre_ready <- nfhs_men_analysis %>%
  filter(
    state_name %in% phase2_states,
    state_name != "Chandigarh",
    postcovid == "Pre-lockdown"
  ) %>%
  mutate(
    timepoint = "NFHS5 Pre-lockdown",
    undernourished_bin = man_malnourished == "Undernourished",
    anemia_bin = man_anemia == "Anemic"
  )

nfhs5_post_ready <- nfhs_men_analysis %>%
  filter(
    state_name %in% phase2_states,
    state_name != "Chandigarh",
    postcovid == "Post-lockdown"
  ) %>%
  mutate(
    timepoint = "NFHS5 Post-lockdown",
    undernourished_bin = man_malnourished == "Undernourished",
    anemia_bin = man_anemia == "Anemic"
  )

################################################################################
# 10) COMBINE DATASETS
################################################################################

common_vars <- c(
  "timepoint",
  "state_name",
  "man_psu",
  "man_strata",
  "man_weight",
  "man_hh_state_weight",
  "wi_combined_urbanrural",
  "man_age",
  "urban_rural",
  "religion",
  "caste",
  "man_malnourished",
  "man_anemia",
  "undernourished_bin",
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
  ids = ~man_psu,
  strata = ~interaction(man_strata, timepoint),
  weights = ~man_weight,
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

under_prev_3tp <- svyby(
  ~undernourished_bin,
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

under_prev_q_3tp <- svyby(
  ~undernourished_bin,
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

under_prev_3tp_pct <- under_prev_3tp %>%
  as_tibble() %>%
  mutate(
    prevalence_pct = undernourished_binTRUE * 100,
    se_pct = `se.undernourished_binTRUE` * 100,
    ci_l_pct = `ci_l.undernourished_binTRUE` * 100,
    ci_u_pct = `ci_u.undernourished_binTRUE` * 100
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

under_prev_q_3tp_pct <- under_prev_q_3tp %>%
  as_tibble() %>%
  mutate(
    prevalence_pct = undernourished_binTRUE * 100,
    se_pct = `se.undernourished_binTRUE` * 100,
    ci_l_pct = `ci_l.undernourished_binTRUE` * 100,
    ci_u_pct = `ci_u.undernourished_binTRUE` * 100
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

under_prev_3tp_pct
anemia_prev_3tp_pct

under_prev_q_3tp_pct
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
    caste = factor(caste),
    religion = factor(religion)
  )

des_nfhs5_phase2_only <- svydesign(
  ids = ~man_psu,
  strata = ~interaction(man_strata, timepoint),
  weights = ~man_weight,
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

fit_under_nfhs5 <- svyglm(
  undernourished_bin ~ timepoint + man_age + caste + religion,
  design = des_nfhs5_phase2_only,
  family = quasibinomial()
)

fit_anemia_nfhs5 <- svyglm(
  anemia_bin ~ timepoint + man_age + caste + religion,
  design = des_nfhs5_phase2_only,
  family = quasibinomial()
)

OR_under_nfhs5 <- or_timepoint(fit_under_nfhs5)
OR_anemia_nfhs5 <- or_timepoint(fit_anemia_nfhs5)

OR_under_nfhs5
OR_anemia_nfhs5

################################################################################
# 18) ADJUSTED NFHS5 PRE VS POST ODDS RATIOS BY WEALTH QUINTILE
################################################################################

fit_under_q_nfhs5 <- lapply(1:5, function(q) {
  svyglm(
    undernourished_bin ~ timepoint + man_age + caste + religion,
    design = subset(des_nfhs5_phase2_only,
                    wi_combined_urbanrural == q),
    family = quasibinomial()
  )
})

OR_under_q_nfhs5 <- do.call(
  rbind,
  lapply(fit_under_q_nfhs5, or_timepoint)
)

rownames(OR_under_q_nfhs5) <- paste0("Q", 1:5)

fit_anemia_q_nfhs5 <- lapply(1:5, function(q) {
  svyglm(
    anemia_bin ~ timepoint + man_age + caste + religion,
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

OR_under_q_nfhs5
OR_anemia_q_nfhs5

################################################################################

################################################################################
# 18b) ADJUSTED (marginally standardized) PREVALENCE — restricted 13-state sample
#      Same models as the restricted aORs: timepoint + man_age + caste + religion.
#      Predicts each individual's probability under pre and post (holding their
#      own covariates fixed), then averages over the observed covariate
#      distribution to give composition-standardized adjusted prevalences.
################################################################################

adj_prev_sub_sens <- function(formula, dq, wtvar = "man_weight") {
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

# Overall
adj_under_overall_sens <- adj_prev_sub_sens(
  undernourished_bin ~ timepoint + man_age + caste + religion, des_nfhs5_phase2_only)
adj_anemia_overall_sens <- adj_prev_sub_sens(
  anemia_bin ~ timepoint + man_age + caste + religion, des_nfhs5_phase2_only)

# By wealth quintile
adj_under_by_q_sens <- t(sapply(1:5, function(q)
  adj_prev_sub_sens(undernourished_bin ~ timepoint + man_age + caste + religion,
                    subset(des_nfhs5_phase2_only, wi_combined_urbanrural == q))))
rownames(adj_under_by_q_sens) <- paste0("Q", 1:5)

adj_anemia_by_q_sens <- t(sapply(1:5, function(q)
  adj_prev_sub_sens(anemia_bin ~ timepoint + man_age + caste + religion,
                    subset(des_nfhs5_phase2_only, wi_combined_urbanrural == q))))
rownames(adj_anemia_by_q_sens) <- paste0("Q", 1:5)

# --- Print
cat("\n=== RESTRICTED 13-STATE ADJUSTED PREVALENCE (Underweight, men) ===\n")
cat("Overall pre/post:", round(adj_under_overall_sens[["pre"]],1),
    round(adj_under_overall_sens[["post"]],1), "\n")
cat("By quintile (pre, post):\n"); print(round(adj_under_by_q_sens, 1))

cat("\n=== RESTRICTED 13-STATE ADJUSTED PREVALENCE (Anemia, men) ===\n")
cat("Overall pre/post:", round(adj_anemia_overall_sens[["pre"]],1),
    round(adj_anemia_overall_sens[["post"]],1), "\n")
cat("By quintile (pre, post):\n"); print(round(adj_anemia_by_q_sens, 1))


################################################################################