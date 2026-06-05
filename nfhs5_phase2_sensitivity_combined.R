################################################################################
######################## NFHS-5 Sensitivity Analysis ############################
######################## Phase 2 states: Pre vs Post ############################
################################################################################

# This script combines the men, women, and children sensitivity analyses.
# NFHS-4 sections have been removed. The analysis is restricted to NFHS-5
# Phase 2 states, excluding Chandigarh because it has no pre-lockdown observations.

################################################################################
# Libraries
################################################################################

library(haven)
library(survey)
library(tidyverse)

options(survey.lonely.psu = "adjust")

################################################################################
# Helper functions
################################################################################

print_section <- function(title) {
  cat("\n")
  cat("================================================================================\n")
  cat(title, "\n")
  cat("================================================================================\n")
}

clean_prev_table <- function(prev_obj, estimate_col) {
  prev_obj %>%
    as_tibble() %>%
    mutate(
      prevalence_pct = .data[[estimate_col]] * 100,
      se_pct = se * 100,
      ci_l_pct = ci_l * 100,
      ci_u_pct = ci_u * 100,
      estimate = sprintf("%.1f%% (%.1f, %.1f)", prevalence_pct, ci_l_pct, ci_u_pct)
    )
}

or_postcovid <- function(fit) {
  term <- "postcovidPost-lockdown"
  est <- coef(fit)[term]
  ci <- confint(fit)[term, ]

  tibble(
    aOR = exp(est),
    CI_low = exp(ci[1]),
    CI_high = exp(ci[2]),
    estimate = sprintf("%.2f (%.2f, %.2f)", exp(est), exp(ci[1]), exp(ci[2]))
  )
}

or_by_quintile <- function(design, formula_text) {
  out <- lapply(1:5, function(q) {
    fit <- svyglm(
      as.formula(formula_text),
      design = subset(design, wi_combined_urbanrural == q),
      family = quasibinomial()
    )

    or_postcovid(fit) %>%
      mutate(wealth_quintile = paste0("Q", q), .before = 1)
  })

  bind_rows(out)
}

################################################################################
# Phase 2 state codebook for NFHS-5
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
############################# 1) MEN ANALYSIS ###################################
################################################################################

################################################################################
# 1A) Load NFHS-5 men + household biomarker data
################################################################################

men_keep5 <- c(
  "MV001","MV002","MV003",
  "MV005","MV021","MV022","MV023",
  "MV191A","MV190A","MV024"
)

hh_keep_men5 <- c(
  "HV001","HV002","HVIDX",
  "SHMWEIGHT","HV005","HV021","HV022","HV023",
  "HV006","HV007",
  "HV270","HV270A","HV271","HV271A",
  "HB1","HB2","HB3","HB35","HB40","HB56","HB57",
  "HV024","HV025",
  "SH47","SH49"
)

nfhs_men_sub5 <- read_sav(
  "./IAMR7EFL.SAV",
  col_select = all_of(men_keep5)
)

nfhs_hh_men5 <- read_sav(
  "./IAPR7EFL.SAV",
  col_select = all_of(hh_keep_men5)
)

################################################################################
# 1B) Merge and rename variables
################################################################################

nfhs_hh_for_men5 <- nfhs_hh_men5 %>%
  rename(MV001 = HV001, MV002 = HV002, MV003 = HVIDX)

nfhs_men_analysis <- nfhs_men_sub5 %>%
  left_join(nfhs_hh_for_men5, by = c("MV001", "MV002", "MV003")) %>%
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
# 1C) Clean NFHS-5 men variables
################################################################################

nfhs_men_analysis <- nfhs_men_analysis %>%
  mutate(
    man_weight_kgs = ifelse(man_weight_kgs >= 9990, NA, man_weight_kgs / 10),
    man_height_cms = ifelse(man_height_cms >= 9990, NA, man_height_cms / 10),
    man_bmi = ifelse(man_bmi >= 9990, NA, man_bmi / 100),

    man_hh_state_weight = man_hh_state_weight / 1e6,
    man_weight = man_weight / 1e6,

    wi_combined_urbanrural = as.numeric(haven::zap_labels(wi_combined_urbanrural)),
    man_region_num = as.numeric(haven::zap_labels(man_region)),
    urban_rural = as.numeric(haven::zap_labels(urban_rural)),
    religion = as.numeric(haven::zap_labels(religion)),
    caste = as.numeric(haven::zap_labels(caste))
  )

################################################################################
# 1D) Create outcomes + pre/post variable
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
    ),

    date_of_interview = as.Date(
      paste(year_of_interview, month_of_interview, "01", sep = "-")
    ),

    postcovid = case_when(
      date_of_interview < as.Date("2020-04-01") ~ "Pre-lockdown",
      date_of_interview >= as.Date("2020-06-01") ~ "Post-lockdown",
      TRUE ~ NA_character_
    )
  )

################################################################################
# 1E) Restrict to NFHS-5 Phase 2 states, excluding Chandigarh
################################################################################

men_phase2 <- nfhs_men_analysis %>%
  left_join(state_codebook_nfhs5, by = c("man_region_num" = "region_code5")) %>%
  filter(
    state_name %in% phase2_states,
    state_name != "Chandigarh",
    !is.na(postcovid)
  ) %>%
  mutate(
    postcovid = relevel(factor(postcovid), ref = "Pre-lockdown"),
    caste = factor(caste),
    religion = factor(religion),

    undernourished_bin = case_when(
      is.na(man_malnourished) ~ NA_real_,
      man_malnourished == "Undernourished" ~ 1,
      TRUE ~ 0
    ),
    anemia_bin = case_when(
      is.na(man_anemia) ~ NA_real_,
      man_anemia == "Anemic" ~ 1,
      TRUE ~ 0
    )
  )

################################################################################
# 1F) Men survey design
################################################################################

des_men_phase2 <- svydesign(
  ids = ~man_psu,
  strata = ~interaction(man_strata, postcovid),
  weights = ~man_hh_state_weight,
  data = men_phase2,
  nest = TRUE
)

################################################################################
# 1G) Men prevalence estimates
################################################################################

men_under_prev <- svyby(
  ~undernourished_bin,
  ~postcovid,
  des_men_phase2,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

men_anemia_prev <- svyby(
  ~anemia_bin,
  ~postcovid,
  des_men_phase2,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

men_under_prev_q <- svyby(
  ~undernourished_bin,
  ~postcovid + wi_combined_urbanrural,
  des_men_phase2,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

men_anemia_prev_q <- svyby(
  ~anemia_bin,
  ~postcovid + wi_combined_urbanrural,
  des_men_phase2,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

men_under_prev_pct <- clean_prev_table(men_under_prev, "undernourished_bin") %>%
  select(postcovid, prevalence_pct, se_pct, ci_l_pct, ci_u_pct, estimate)

men_anemia_prev_pct <- clean_prev_table(men_anemia_prev, "anemia_bin") %>%
  select(postcovid, prevalence_pct, se_pct, ci_l_pct, ci_u_pct, estimate)

men_under_prev_q_pct <- clean_prev_table(men_under_prev_q, "undernourished_bin") %>%
  select(postcovid, wi_combined_urbanrural, prevalence_pct, se_pct, ci_l_pct, ci_u_pct, estimate)

men_anemia_prev_q_pct <- clean_prev_table(men_anemia_prev_q, "anemia_bin") %>%
  select(postcovid, wi_combined_urbanrural, prevalence_pct, se_pct, ci_l_pct, ci_u_pct, estimate)

################################################################################
# 1H) Men adjusted ORs
################################################################################

men_under_fit <- svyglm(
  undernourished_bin ~ postcovid + man_age + caste + religion,
  design = des_men_phase2,
  family = quasibinomial()
)

men_anemia_fit <- svyglm(
  anemia_bin ~ postcovid + man_age + caste + religion,
  design = des_men_phase2,
  family = quasibinomial()
)

men_under_aor <- or_postcovid(men_under_fit) %>%
  mutate(outcome = "Undernutrition", group = "Overall dataset", .before = 1)

men_anemia_aor <- or_postcovid(men_anemia_fit) %>%
  mutate(outcome = "Anemia", group = "Overall dataset", .before = 1)

men_under_aor_q <- or_by_quintile(
  des_men_phase2,
  "undernourished_bin ~ postcovid + man_age + caste + religion"
) %>%
  mutate(outcome = "Undernutrition", .before = 1)

men_anemia_aor_q <- or_by_quintile(
  des_men_phase2,
  "anemia_bin ~ postcovid + man_age + caste + religion"
) %>%
  mutate(outcome = "Anemia", .before = 1)

################################################################################
############################ 2) WOMEN ANALYSIS ##################################
################################################################################

################################################################################
# 2A) Load NFHS-5 women + household biomarker data
################################################################################

women_keep5 <- c(
  "V001","V002","V003",
  "V005","V021","V022","V023",
  "V191A","V190A","V024"
)

hh_keep_women5 <- c(
  "HV001","HV002","HVIDX",
  "SHWEIGHT","HV005","HV021","HV022","HV023",
  "HV006","HV007",
  "HV270","HV270A","HV271","HV271A",
  "HA1","HA2","HA3","HA40","HA56","HA57",
  "HV024","HV025",
  "SH47","SH49"
)

nfhs_women_sub5 <- read_sav(
  "./IAIR7EFL.SAV",
  col_select = all_of(women_keep5)
)

nfhs_hh_women5 <- read_sav(
  "./IAPR7EFL.SAV",
  col_select = all_of(hh_keep_women5)
)

################################################################################
# 2B) Merge and rename variables
################################################################################

nfhs_hh_for_women5 <- nfhs_hh_women5 %>%
  rename(V001 = HV001, V002 = HV002, V003 = HVIDX)

nfhs_women_analysis <- nfhs_women_sub5 %>%
  left_join(nfhs_hh_for_women5, by = c("V001", "V002", "V003")) %>%
  rename(
    month_of_interview = HV006,
    year_of_interview  = HV007,

    wi_combined = HV270,
    wi_combined_urbanrural = HV270A,

    woman_age = HA1,
    woman_weight_kgs = HA2,
    woman_height_cms = HA3,
    woman_bmi = HA40,
    woman_anemia_level = HA57,

    woman_weight = V005,
    woman_psu = V021,
    woman_strata = V022,

    hh_state_weight = SHWEIGHT,

    woman_region = V024,
    urban_rural = HV025,
    religion = SH47,
    caste = SH49
  )

################################################################################
# 2C) Clean NFHS-5 women variables
################################################################################

nfhs_women_analysis <- nfhs_women_analysis %>%
  mutate(
    woman_weight_kgs = ifelse(woman_weight_kgs >= 9990, NA, woman_weight_kgs / 10),
    woman_height_cms = ifelse(woman_height_cms >= 9990, NA, woman_height_cms / 10),
    woman_bmi = ifelse(woman_bmi >= 9990, NA, woman_bmi / 100),

    hh_state_weight = hh_state_weight / 1e6,
    woman_weight = woman_weight / 1e6,

    wi_combined_urbanrural = as.numeric(haven::zap_labels(wi_combined_urbanrural)),
    woman_region_num = as.numeric(haven::zap_labels(woman_region)),
    urban_rural = as.numeric(haven::zap_labels(urban_rural)),
    religion = as.numeric(haven::zap_labels(religion)),
    caste = as.numeric(haven::zap_labels(caste))
  )

################################################################################
# 2D) Create outcomes + pre/post variable
################################################################################

nfhs_women_analysis <- nfhs_women_analysis %>%
  mutate(
    woman_malnourished = case_when(
      is.na(woman_bmi) ~ NA_real_,
      woman_bmi < 18.5 ~ 1,
      TRUE ~ 0
    ),
    woman_malnourished = factor(
      woman_malnourished,
      levels = c(0, 1),
      labels = c("Not undernourished", "Undernourished")
    ),

    woman_anemia = case_when(
      is.na(woman_anemia_level) ~ NA_real_,
      woman_anemia_level %in% c(1, 2, 3) ~ 1,
      TRUE ~ 0
    ),
    woman_anemia = factor(
      woman_anemia,
      levels = c(0, 1),
      labels = c("Not anemic", "Anemic")
    ),

    date_of_interview = as.Date(
      paste(year_of_interview, month_of_interview, "01", sep = "-")
    ),

    postcovid = case_when(
      date_of_interview < as.Date("2020-04-01") ~ "Pre-lockdown",
      date_of_interview >= as.Date("2020-06-01") ~ "Post-lockdown",
      TRUE ~ NA_character_
    )
  )

################################################################################
# 2E) Restrict to NFHS-5 Phase 2 states, excluding Chandigarh
################################################################################

women_phase2 <- nfhs_women_analysis %>%
  left_join(state_codebook_nfhs5, by = c("woman_region_num" = "region_code5")) %>%
  filter(
    state_name %in% phase2_states,
    state_name != "Chandigarh",
    !is.na(postcovid)
  ) %>%
  mutate(
    postcovid = relevel(factor(postcovid), ref = "Pre-lockdown"),
    caste = factor(caste),
    religion = factor(religion),

    undernourished_bin = case_when(
      is.na(woman_malnourished) ~ NA_real_,
      woman_malnourished == "Undernourished" ~ 1,
      TRUE ~ 0
    ),
    anemia_bin = case_when(
      is.na(woman_anemia) ~ NA_real_,
      woman_anemia == "Anemic" ~ 1,
      TRUE ~ 0
    )
  )

################################################################################
# 2F) Women survey design
################################################################################

des_women_phase2 <- svydesign(
  ids = ~woman_psu,
  strata = ~interaction(woman_strata, postcovid),
  weights = ~woman_weight,
  data = women_phase2,
  nest = TRUE
)

################################################################################
# 2G) Women prevalence estimates
################################################################################

women_under_prev <- svyby(
  ~undernourished_bin,
  ~postcovid,
  des_women_phase2,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

women_anemia_prev <- svyby(
  ~anemia_bin,
  ~postcovid,
  des_women_phase2,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

women_under_prev_q <- svyby(
  ~undernourished_bin,
  ~postcovid + wi_combined_urbanrural,
  des_women_phase2,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

women_anemia_prev_q <- svyby(
  ~anemia_bin,
  ~postcovid + wi_combined_urbanrural,
  des_women_phase2,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

women_under_prev_pct <- clean_prev_table(women_under_prev, "undernourished_bin") %>%
  select(postcovid, prevalence_pct, se_pct, ci_l_pct, ci_u_pct, estimate)

women_anemia_prev_pct <- clean_prev_table(women_anemia_prev, "anemia_bin") %>%
  select(postcovid, prevalence_pct, se_pct, ci_l_pct, ci_u_pct, estimate)

women_under_prev_q_pct <- clean_prev_table(women_under_prev_q, "undernourished_bin") %>%
  select(postcovid, wi_combined_urbanrural, prevalence_pct, se_pct, ci_l_pct, ci_u_pct, estimate)

women_anemia_prev_q_pct <- clean_prev_table(women_anemia_prev_q, "anemia_bin") %>%
  select(postcovid, wi_combined_urbanrural, prevalence_pct, se_pct, ci_l_pct, ci_u_pct, estimate)

################################################################################
# 2H) Women adjusted ORs
################################################################################

women_under_fit <- svyglm(
  undernourished_bin ~ postcovid + woman_age + caste + religion,
  design = des_women_phase2,
  family = quasibinomial()
)

women_anemia_fit <- svyglm(
  anemia_bin ~ postcovid + woman_age + caste + religion,
  design = des_women_phase2,
  family = quasibinomial()
)

women_under_aor <- or_postcovid(women_under_fit) %>%
  mutate(outcome = "Undernutrition", group = "Overall dataset", .before = 1)

women_anemia_aor <- or_postcovid(women_anemia_fit) %>%
  mutate(outcome = "Anemia", group = "Overall dataset", .before = 1)

women_under_aor_q <- or_by_quintile(
  des_women_phase2,
  "undernourished_bin ~ postcovid + woman_age + caste + religion"
) %>%
  mutate(outcome = "Undernutrition", .before = 1)

women_anemia_aor_q <- or_by_quintile(
  des_women_phase2,
  "anemia_bin ~ postcovid + woman_age + caste + religion"
) %>%
  mutate(outcome = "Anemia", .before = 1)

################################################################################
############################ 3) CHILDREN ANALYSIS ###############################
################################################################################

################################################################################
# 3A) Load NFHS-5 children + household data
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

hh_keep_child5 <- c(
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
  col_select = all_of(hh_keep_child5)
)

################################################################################
# 3B) Merge and rename variables
################################################################################

nfhs_hh_for_child5 <- nfhs_hh_child5 %>%
  rename(V001 = HV001, V002 = HV002, V003 = HVIDX)

nfhs_child_analysis <- nfhs_child_sub5 %>%
  left_join(nfhs_hh_for_child5, by = c("V001", "V002", "V003")) %>%
  rename(
    month_of_interview = HV006,
    year_of_interview  = HV007,

    wi_combined = HV270,
    wi_combined_urbanrural = HV270A,

    child_age_years = B8,
    child_age_months = HW1,
    child_alive = B5,
    child_sex = B4,

    child_weight_kgs = HW2,
    child_height_cms = HW3,

    child_WFH_SD_WHO = HW72,
    child_anemia_level = HW57,

    child_weight = V005,
    child_psu = V021,
    child_strata = V022,

    hh_state_weight = SHWEIGHT,

    child_region = V024,
    urban_rural = HV025,
    religion = SH47,
    caste = SH49
  )

################################################################################
# 3C) Clean NFHS-5 children variables
################################################################################

nfhs_child_analysis <- nfhs_child_analysis %>%
  mutate(
    child_weight_kgs = ifelse(child_weight_kgs >= 9990, NA, child_weight_kgs / 10),
    child_height_cms = ifelse(child_height_cms >= 9990, NA, child_height_cms / 10),
    child_WFH_SD_WHO = ifelse(child_WFH_SD_WHO >= 9990, NA, child_WFH_SD_WHO / 100),

    hh_state_weight = hh_state_weight / 1e6,
    child_weight = child_weight / 1e6,

    wi_combined_urbanrural = as.numeric(haven::zap_labels(wi_combined_urbanrural)),
    child_region_num = as.numeric(haven::zap_labels(child_region)),
    urban_rural = as.numeric(haven::zap_labels(urban_rural)),
    religion = as.numeric(haven::zap_labels(religion)),
    caste = as.numeric(haven::zap_labels(caste))
  ) %>%
  filter(child_alive == 1)

################################################################################
# 3D) Create outcomes + pre/post variable
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
    ),

    date_of_interview = as.Date(
      paste(year_of_interview, month_of_interview, "01", sep = "-")
    ),

    postcovid = case_when(
      date_of_interview < as.Date("2020-04-01") ~ "Pre-lockdown",
      date_of_interview >= as.Date("2020-06-01") ~ "Post-lockdown",
      TRUE ~ NA_character_
    )
  )

################################################################################
# 3E) Restrict to NFHS-5 Phase 2 states, excluding Chandigarh
################################################################################

children_phase2 <- nfhs_child_analysis %>%
  left_join(state_codebook_nfhs5, by = c("child_region_num" = "region_code5")) %>%
  filter(
    state_name %in% phase2_states,
    state_name != "Chandigarh",
    !is.na(postcovid)
  ) %>%
  mutate(
    postcovid = relevel(factor(postcovid), ref = "Pre-lockdown"),
    caste = factor(caste),
    religion = factor(religion),

    wasting_bin = case_when(
      is.na(child_wasting) ~ NA_real_,
      child_wasting == "Wasting" ~ 1,
      TRUE ~ 0
    ),
    anemia_bin = case_when(
      is.na(child_anemia) ~ NA_real_,
      child_anemia == "Anemic" ~ 1,
      TRUE ~ 0
    )
  )

################################################################################
# 3F) Children survey design
################################################################################

des_children_phase2 <- svydesign(
  ids = ~child_psu,
  strata = ~interaction(child_strata, postcovid),
  weights = ~child_weight,
  data = children_phase2,
  nest = TRUE
)

################################################################################
# 3G) Children prevalence estimates
################################################################################

children_wasting_prev <- svyby(
  ~wasting_bin,
  ~postcovid,
  des_children_phase2,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

children_anemia_prev <- svyby(
  ~anemia_bin,
  ~postcovid,
  des_children_phase2,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

children_wasting_prev_q <- svyby(
  ~wasting_bin,
  ~postcovid + wi_combined_urbanrural,
  des_children_phase2,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

children_anemia_prev_q <- svyby(
  ~anemia_bin,
  ~postcovid + wi_combined_urbanrural,
  des_children_phase2,
  svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

children_wasting_prev_pct <- clean_prev_table(children_wasting_prev, "wasting_bin") %>%
  select(postcovid, prevalence_pct, se_pct, ci_l_pct, ci_u_pct, estimate)

children_anemia_prev_pct <- clean_prev_table(children_anemia_prev, "anemia_bin") %>%
  select(postcovid, prevalence_pct, se_pct, ci_l_pct, ci_u_pct, estimate)

children_wasting_prev_q_pct <- clean_prev_table(children_wasting_prev_q, "wasting_bin") %>%
  select(postcovid, wi_combined_urbanrural, prevalence_pct, se_pct, ci_l_pct, ci_u_pct, estimate)

children_anemia_prev_q_pct <- clean_prev_table(children_anemia_prev_q, "anemia_bin") %>%
  select(postcovid, wi_combined_urbanrural, prevalence_pct, se_pct, ci_l_pct, ci_u_pct, estimate)

################################################################################
# 3H) Children adjusted ORs
################################################################################

children_wasting_fit <- svyglm(
  wasting_bin ~ postcovid + child_age_months + caste + religion,
  design = des_children_phase2,
  family = quasibinomial()
)

children_anemia_fit <- svyglm(
  anemia_bin ~ postcovid + child_age_months + caste + religion,
  design = des_children_phase2,
  family = quasibinomial()
)

children_wasting_aor <- or_postcovid(children_wasting_fit) %>%
  mutate(outcome = "Wasting", group = "Overall dataset", .before = 1)

children_anemia_aor <- or_postcovid(children_anemia_fit) %>%
  mutate(outcome = "Anemia", group = "Overall dataset", .before = 1)

children_wasting_aor_q <- or_by_quintile(
  des_children_phase2,
  "wasting_bin ~ postcovid + child_age_months + caste + religion"
) %>%
  mutate(outcome = "Wasting", .before = 1)

children_anemia_aor_q <- or_by_quintile(
  des_children_phase2,
  "anemia_bin ~ postcovid + child_age_months + caste + religion"
) %>%
  mutate(outcome = "Anemia", .before = 1)

################################################################################
# 4) Print outputs
################################################################################

print_section("SAMPLE SIZE CHECKS")

cat("\nMen: unweighted N by pre/post period\n")
print(men_phase2 %>% count(postcovid))

cat("\nWomen: unweighted N by pre/post period\n")
print(women_phase2 %>% count(postcovid))

cat("\nChildren: unweighted N by pre/post period\n")
print(children_phase2 %>% count(postcovid))

print_section("MEN: OVERALL WEIGHTED PREVALENCE")
cat("\nUndernutrition:\n")
print(men_under_prev_pct)
cat("\nAnemia:\n")
print(men_anemia_prev_pct)

print_section("MEN: WEIGHTED PREVALENCE BY WEALTH QUINTILE")
cat("\nUndernutrition:\n")
print(men_under_prev_q_pct)
cat("\nAnemia:\n")
print(men_anemia_prev_q_pct)

print_section("MEN: ADJUSTED ODDS RATIOS")
cat("\nOverall aORs:\n")
print(bind_rows(men_under_aor, men_anemia_aor))
cat("\naORs by wealth quintile:\n")
print(bind_rows(men_under_aor_q, men_anemia_aor_q))

print_section("WOMEN: OVERALL WEIGHTED PREVALENCE")
cat("\nUndernutrition:\n")
print(women_under_prev_pct)
cat("\nAnemia:\n")
print(women_anemia_prev_pct)

print_section("WOMEN: WEIGHTED PREVALENCE BY WEALTH QUINTILE")
cat("\nUndernutrition:\n")
print(women_under_prev_q_pct)
cat("\nAnemia:\n")
print(women_anemia_prev_q_pct)

print_section("WOMEN: ADJUSTED ODDS RATIOS")
cat("\nOverall aORs:\n")
print(bind_rows(women_under_aor, women_anemia_aor))
cat("\naORs by wealth quintile:\n")
print(bind_rows(women_under_aor_q, women_anemia_aor_q))

print_section("CHILDREN: OVERALL WEIGHTED PREVALENCE")
cat("\nWasting:\n")
print(children_wasting_prev_pct)
cat("\nAnemia:\n")
print(children_anemia_prev_pct)

print_section("CHILDREN: WEIGHTED PREVALENCE BY WEALTH QUINTILE")
cat("\nWasting:\n")
print(children_wasting_prev_q_pct)
cat("\nAnemia:\n")
print(children_anemia_prev_q_pct)

print_section("CHILDREN: ADJUSTED ODDS RATIOS")
cat("\nOverall aORs:\n")
print(bind_rows(children_wasting_aor, children_anemia_aor))
cat("\naORs by wealth quintile:\n")
print(bind_rows(children_wasting_aor_q, children_anemia_aor_q))

################################################################################
# End
################################################################################
