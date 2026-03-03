################################################################################
############################ NFHS-5 Women Analysis #############################
########################## Author: anrmarri@bu.edu #############################
################################################################################

# Libraries
library(haven)
library(survey)
library(rineq)
library(viridis)
library(tidyverse)

options(survey.lonely.psu = "adjust")

################################################################################
# 1) Load only required variables
################################################################################

women_keep <- c(
  "V001","V002","V003",          # merge keys
  "V005","V021","V022","V023",   # weights/PSU/strata/domain
  "V191A","V190A","V024"         # wealth + region
)

hh_keep <- c(
  "HV001","HV002","HVIDX",       # merge keys
  "SHWEIGHT","HV005","HV021","HV022","HV023",
  "HV006","HV007",               # month/year
  "HV270","HV270A","HV271","HV271A", # wealth
  "HA1","HA2","HA3","HA40","HA56","HA57", # anthro/biomarkers
  "HV024"                        # region
)

nfhs_women_sub <- read_sav(
  "./IAIR7ESV_women/IAIR7EFL.SAV",
  col_select = all_of(women_keep)
)

nfhs_hh_women <- read_sav(
  "./IAPR7ESV_hh/IAPR7EFL.SAV",
  col_select = all_of(hh_keep)
)

################################################################################
# 2) Merge women + household file and rename variables
################################################################################

nfhs_hh_for_women <- nfhs_hh_women %>%
  rename(V001 = HV001, V002 = HV002, V003 = HVIDX)

nfhs_women_analysis <- nfhs_women_sub %>%
  left_join(nfhs_hh_for_women, by = c("V001","V002","V003")) %>%
  rename(
    month_of_interview = HV006,
    year_of_interview  = HV007,
    
    wi_factorscore = HV271,
    wi_factorscore_urbanrural = HV271A,
    wi_combined = HV270,
    wi_combined_urbanrural = HV270A,
    
    woman_age = HA1,
    woman_weight_kgs = HA2,
    woman_height_cms = HA3,
    woman_bmi = HA40,
    woman_hb_gdl_altsmoking = HA56,
    woman_anemia_level = HA57,
    
    woman_weight = V005,
    woman_psu = V021,
    woman_strata = V022,
    woman_sample_domain = V023,
    woman_wi_factorscore = V191A,
    woman_wi_quintile = V190A,
    
    hh_state_weight = SHWEIGHT,
    hh_weight = HV005,
    hh_psu = HV021,
    hh_strata = HV022,
    hh_sample_domain = HV023,
    
    woman_region = V024,
    household_region = HV024
  )

################################################################################
# 3) Clean special missing codes + rescale to real units
################################################################################

# Special missing codes
nfhs_women_analysis$woman_weight_kgs[nfhs_women_analysis$woman_weight_kgs >= 9990] <- NA
nfhs_women_analysis$woman_height_cms[nfhs_women_analysis$woman_height_cms >= 9990] <- NA
nfhs_women_analysis$woman_bmi[nfhs_women_analysis$woman_bmi >= 9990] <- NA
nfhs_women_analysis$woman_hb_gdl_altsmoking[nfhs_women_analysis$woman_hb_gdl_altsmoking >= 990] <- NA

# Rescale decimals
nfhs_women_analysis$wi_factorscore_urbanrural <- nfhs_women_analysis$wi_factorscore_urbanrural / 100000
nfhs_women_analysis$woman_weight_kgs <- nfhs_women_analysis$woman_weight_kgs / 10
nfhs_women_analysis$woman_height_cms <- nfhs_women_analysis$woman_height_cms / 10
nfhs_women_analysis$woman_bmi <- nfhs_women_analysis$woman_bmi / 100
nfhs_women_analysis$woman_hb_gdl_altsmoking <- nfhs_women_analysis$woman_hb_gdl_altsmoking / 10

# Scale weights
nfhs_women_analysis$hh_state_weight <- nfhs_women_analysis$hh_state_weight / 1e6
nfhs_women_analysis$woman_weight <- nfhs_women_analysis$woman_weight / 1e6
nfhs_women_analysis$hh_weight <- nfhs_women_analysis$hh_weight / 1e6

# Ensure DHS wealth quintile is numeric 1-5
nfhs_women_analysis$wi_combined_urbanrural <- as.numeric(as.character(nfhs_women_analysis$wi_combined_urbanrural))

################################################################################
# 4) Outcomes + exposure
################################################################################

# Undernutrition (BMI < 18.5) using DHS-reported BMI
nfhs_women_analysis$woman_malnourished <- ifelse(
  is.na(nfhs_women_analysis$woman_bmi), NA,
  ifelse(nfhs_women_analysis$woman_bmi < 18.5, 1, 0)
)
nfhs_women_analysis$woman_malnourished <- factor(
  nfhs_women_analysis$woman_malnourished,
  levels = c(0,1),
  labels = c("Not undernourished", "Undernourished")
)

# Anemia (levels 1-3 vs 4)
nfhs_women_analysis$woman_anemia <- ifelse(
  is.na(nfhs_women_analysis$woman_anemia_level), NA,
  ifelse(nfhs_women_analysis$woman_anemia_level %in% c(1,2,3), 1, 0)
)
nfhs_women_analysis$woman_anemia <- factor(
  nfhs_women_analysis$woman_anemia,
  levels = c(0,1),
  labels = c("Not anemic", "Anemic")
)

# Date of interview + pre/post lockdown definition
nfhs_women_analysis$date_of_interview <- as.Date(
  paste(nfhs_women_analysis$year_of_interview, nfhs_women_analysis$month_of_interview, "01", sep = "-")
)

nfhs_women_analysis$postcovid <- ifelse(
  nfhs_women_analysis$date_of_interview < as.Date("2020-04-01"), "Pre-lockdown",
  ifelse(nfhs_women_analysis$date_of_interview >= as.Date("2020-06-01"), "Post-lockdown", NA)
)
nfhs_women_analysis$postcovid <- relevel(factor(nfhs_women_analysis$postcovid), ref = "Pre-lockdown")

################################################################################
# 5) Survey design (women individual weight: V005)
################################################################################

des_women_bio <- svydesign(
  ids = ~woman_psu,
  strata = ~woman_strata,
  weights = ~woman_weight,
  data = nfhs_women_analysis,
  nest = TRUE
)

################################################################################
# 6) Weighted prevalence (overall + by wealth quintile) + chi-square
################################################################################

wprev <- function(des, expr) {
  svymean(as.formula(paste0("~I(", expr, ")")), des, na.rm = TRUE) * 100
}

# Overall prevalence
prev_under_pre  <- wprev(subset(des_women_bio, postcovid == "Pre-lockdown"),  'woman_malnourished=="Undernourished"')
prev_under_post <- wprev(subset(des_women_bio, postcovid == "Post-lockdown"), 'woman_malnourished=="Undernourished"')

prev_an_pre  <- wprev(subset(des_women_bio, postcovid == "Pre-lockdown"),  'woman_anemia=="Anemic"')
prev_an_post <- wprev(subset(des_women_bio, postcovid == "Post-lockdown"), 'woman_anemia=="Anemic"')

# By wealth quintile (weighted % within quintile)
under_by_q_pre <- prop.table(
  svytable(~I(woman_malnourished=="Undernourished") + wi_combined_urbanrural,
           subset(des_women_bio, postcovid=="Pre-lockdown")),
  margin = 2
) * 100

under_by_q_post <- prop.table(
  svytable(~I(woman_malnourished=="Undernourished") + wi_combined_urbanrural,
           subset(des_women_bio, postcovid=="Post-lockdown")),
  margin = 2
) * 100

an_by_q_pre <- prop.table(
  svytable(~I(woman_anemia=="Anemic") + wi_combined_urbanrural,
           subset(des_women_bio, postcovid=="Pre-lockdown")),
  margin = 2
) * 100

an_by_q_post <- prop.table(
  svytable(~I(woman_anemia=="Anemic") + wi_combined_urbanrural,
           subset(des_women_bio, postcovid=="Post-lockdown")),
  margin = 2
) * 100

# Weighted chi-square (overall pre vs post)
des_women_bio <- update(
  des_women_bio,
  undernourished_bin = (woman_malnourished == "Undernourished"),
  anemia_bin         = (woman_anemia == "Anemic")
)

chi_under <- svychisq(~undernourished_bin + postcovid, des_women_bio, na.rm = TRUE)
chi_an    <- svychisq(~anemia_bin + postcovid,         des_women_bio, na.rm = TRUE)

################################################################################
# 7) Weighted odds ratios (overall + by wealth quintile)
# (no forest plots here — you said you'll fill separately)
################################################################################

or_post <- function(fit) {
  term <- "postcovidPost-lockdown"
  est <- coef(fit)[term]
  ci  <- confint(fit)[term, ]
  c(OR = exp(est), CI_low = exp(ci[1]), CI_high = exp(ci[2]))
}

fit_under_all <- svyglm(undernourished_bin ~ postcovid, design = des_women_bio, family = quasibinomial())
fit_an_all    <- svyglm(anemia_bin ~ postcovid,         design = des_women_bio, family = quasibinomial())

OR_under_all <- or_post(fit_under_all)
OR_an_all    <- or_post(fit_an_all)

fit_under_by_q <- lapply(1:5, function(q) {
  svyglm(undernourished_bin ~ postcovid,
         design = subset(des_women_bio, wi_combined_urbanrural == q),
         family = quasibinomial())
})
OR_under_by_q <- do.call(rbind, lapply(fit_under_by_q, or_post))
rownames(OR_under_by_q) <- paste0("Q", 1:5)

fit_an_by_q <- lapply(1:5, function(q) {
  svyglm(anemia_bin ~ postcovid,
         design = subset(des_women_bio, wi_combined_urbanrural == q),
         family = quasibinomial())
})
OR_an_by_q <- do.call(rbind, lapply(fit_an_by_q, or_post))
rownames(OR_an_by_q) <- paste0("Q", 1:5)

################################################################################
# 8) Weighted concentration indices + curve data (same fix as Men template)
################################################################################

make_ci_curve <- function(df, outcome_bin, w, q) {
  
  out <- df %>%
    filter(!is.na(.data[[q]])) %>%
    group_by(.data[[q]]) %>%
    summarise(
      w_pop   = sum(.data[[w]], na.rm = TRUE),
      w_cases = sum(.data[[w]] * (.data[[outcome_bin]] == 1), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(.data[[q]]) %>%
    mutate(
      x = cumsum(w_pop) / sum(w_pop),
      y = cumsum(w_cases) / sum(w_cases)
    )
  
  origin <- out[1, ]
  origin[,] <- 0
  origin[[q]] <- 0
  
  bind_rows(origin, out)
}

df_pre  <- nfhs_women_analysis %>% filter(postcovid == "Pre-lockdown")
df_post <- nfhs_women_analysis %>% filter(postcovid == "Post-lockdown")

df_pre$woman_mal_bin     <- ifelse(df_pre$woman_malnourished == "Undernourished", 1,
                                   ifelse(is.na(df_pre$woman_malnourished), NA, 0))
df_post$woman_mal_bin    <- ifelse(df_post$woman_malnourished == "Undernourished", 1,
                                   ifelse(is.na(df_post$woman_malnourished), NA, 0))
df_pre$woman_anemia_bin  <- ifelse(df_pre$woman_anemia == "Anemic", 1,
                                   ifelse(is.na(df_pre$woman_anemia), NA, 0))
df_post$woman_anemia_bin <- ifelse(df_post$woman_anemia == "Anemic", 1,
                                   ifelse(is.na(df_post$woman_anemia), NA, 0))

curve_under_pre  <- make_ci_curve(df_pre,  "woman_mal_bin",    "woman_weight", "wi_combined_urbanrural")
curve_under_post <- make_ci_curve(df_post, "woman_mal_bin",    "woman_weight", "wi_combined_urbanrural")
curve_an_pre     <- make_ci_curve(df_pre,  "woman_anemia_bin", "woman_weight", "wi_combined_urbanrural")
curve_an_post    <- make_ci_curve(df_post, "woman_anemia_bin", "woman_weight", "wi_combined_urbanrural")

CI_under_pre <- rineq::ci(
  ineqvar  = df_pre$wi_combined_urbanrural,
  outcome  = df_pre$woman_mal_bin,
  weights  = df_pre$woman_weight,
  type     = "CI",
  method   = "direct"
)

CI_under_post <- rineq::ci(
  ineqvar  = df_post$wi_combined_urbanrural,
  outcome  = df_post$woman_mal_bin,
  weights  = df_post$woman_weight,
  type     = "CI",
  method   = "direct"
)

CI_an_pre <- rineq::ci(
  ineqvar  = df_pre$wi_combined_urbanrural,
  outcome  = df_pre$woman_anemia_bin,
  weights  = df_pre$woman_weight,
  type     = "CI",
  method   = "direct"
)

CI_an_post <- rineq::ci(
  ineqvar  = df_post$wi_combined_urbanrural,
  outcome  = df_post$woman_anemia_bin,
  weights  = df_post$woman_weight,
  type     = "CI",
  method   = "direct"
)

################################################################################
# 9) Print results
################################################################################

cat("\n=============================\n")
cat("WEIGHTED PREVALENCE (OVERALL)\n")
cat("=============================\n")
cat("Undernutrition (Pre):\n");  print(prev_under_pre)
cat("Undernutrition (Post):\n"); print(prev_under_post)
cat("Anemia (Pre):\n");          print(prev_an_pre)
cat("Anemia (Post):\n");         print(prev_an_post)

cat("\n====================================\n")
cat("WEIGHTED PREVALENCE BY WEALTH QUINTILE\n")
cat("====================================\n")
cat("Undernutrition (Pre) % by quintile:\n");  print(under_by_q_pre)
cat("Undernutrition (Post) % by quintile:\n"); print(under_by_q_post)
cat("Anemia (Pre) % by quintile:\n");          print(an_by_q_pre)
cat("Anemia (Post) % by quintile:\n");         print(an_by_q_post)

cat("\n=============================\n")
cat("WEIGHTED CHI-SQUARE TESTS\n")
cat("=============================\n")
cat("Undernutrition vs postcovid:\n"); print(chi_under)
cat("Anemia vs postcovid:\n");         print(chi_an)

cat("\n=============================\n")
cat("WEIGHTED ODDS RATIOS\n")
cat("=============================\n")
cat("Overall OR (Undernutrition):\n"); print(OR_under_all)
cat("Overall OR (Anemia):\n");         print(OR_an_all)
cat("\nORs by wealth quintile (Undernutrition):\n"); print(OR_under_by_q)
cat("\nORs by wealth quintile (Anemia):\n");         print(OR_an_by_q)

cat("\n=============================\n")
cat("CONCENTRATION INDICES (WEIGHTED)\n")
cat("=============================\n")
cat("CI Underweight (Pre):\n");  print(CI_under_pre)
cat("CI Underweight (Post):\n"); print(CI_under_post)
cat("CI Anemia (Pre):\n");       print(CI_an_pre)
cat("CI Anemia (Post):\n");      print(CI_an_post)

################################################################################
# 10) Plots and Figures
## Odds ratio: forest plots
# Concentration curves: PRE vs POST overlapping
################################################################################

## Odds ratio: forest plots

# Undernutrition
undernutrition_women <- data.frame(
  Group = c("Q1", "Q2", "Q3", "Q4", "Q5", "Overall dataset"),
  OR = c(1.06, 1.03, 1.00, 0.88, 0.85, 0.95),
  Lower_CI = c(1.02, 0.99, 0.95, 0.84, 0.80, 0.93),
  Upper_CI = c(1.10, 1.07, 1.04, 0.92, 0.90, 0.97)
)

undernutrition_women <- undernutrition_women %>%
  mutate(
    RR_label = sprintf("%.2f (%.2f–%.2f)", OR, Lower_CI, Upper_CI),
    Label_Pos = pmin(Upper_CI + 0.01, 3.3),  # To keep ratios next to respective plots
    Characteristic = factor(Group, levels = rev(Group))  # Preserve order
  )

ggplot(undernutrition_women, aes(y = Group, x = OR)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_text(aes(x = Label_Pos, label = RR_label), hjust = 0, size = 3.5) +
  theme_minimal() +
  labs(title = "Underweight in Women",
       x = "OR and 95% CI",
       y = "Group") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#2E0854") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 12, hjust = 1, color = "black"),
    axis.title.x = element_text(margin = margin(t = 10)),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 70, 10, 10)
  )

# Anemia
anemia_women <- data.frame(
  Group = c("Q1", "Q2", "Q3", "Q4", "Q5", "Overall dataset"),
  OR = c(0.76, 0.73, 0.77, 0.83, 0.94, 0.80),
  Lower_CI = c(0.74, 0.70, 0.74, 0.80, 0.91, 0.79),
  Upper_CI = c(0.79, 0.75, 0.80, 0.86, 0.97, 0.82)
)

#
ggplot(anemia_women, aes(y = Group, x = OR)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  theme_minimal() +
  labs(title = "Anemia in Women",
       x = "OR and 95% CI",
       y = "Group") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#1f9e89") +  
  theme(plot.title = element_text(hjust = 0.5))
#

anemia_women <- anemia_women %>%
  mutate(
    RR_label = sprintf("%.2f (%.2f–%.2f)", OR, Lower_CI, Upper_CI),
    Label_Pos = pmin(Upper_CI + 0.01, 3.3),  # To keep ratios next to respective plots
    Characteristic = factor(Group, levels = rev(Group))  # Preserve order
  )

ggplot(anemia_women, aes(y = Group, x = OR)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_text(aes(x = Label_Pos, label = RR_label), hjust = 0, size = 3.5) +
  theme_minimal() +
  labs(title = "Anemia in Women",
       x = "OR and 95% CI",
       y = "Group") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#1f9e89") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 12, hjust = 1, color = "black"),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 70, 10, 10)
  )

# Underweight curve (Pre + Post on same plot)
plot(curve_under_pre$x, curve_under_pre$y,
     type = "l", lwd = 3, col = "royalblue4",
     main = "Underweight in Women (Pre vs Post Lockdown)",
     xlab = "Cumulative proportion of women (weighted)",
     ylab = "Cumulative proportion underweight (weighted)",
     xlim = c(0,1), ylim = c(0,1)
)
points(curve_under_pre$x, curve_under_pre$y, pch = 16, col = "royalblue4", cex = 1.2)

lines(curve_under_post$x, curve_under_post$y, lwd = 3, col = "firebrick3")
points(curve_under_post$x, curve_under_post$y, pch = 16, col = "firebrick3", cex = 1.2)

abline(0, 1, lwd = 2, lty = 2)

legend("topleft",
       legend = c("Pre-lockdown", "Post-lockdown"),
       col = c("royalblue4", "firebrick3"),
       lwd = 3, pch = 16, bty = "n")

# Anemia curve (Pre + Post on same plot)
plot(curve_an_pre$x, curve_an_pre$y,
     type = "l", lwd = 3, col = "royalblue4",
     main = "Anemia in Women (Pre vs Post Lockdown)",
     xlab = "Cumulative proportion of women (weighted)",
     ylab = "Cumulative proportion anemic (weighted)",
     xlim = c(0,1), ylim = c(0,1)
)
points(curve_an_pre$x, curve_an_pre$y, pch = 16, col = "royalblue4", cex = 1.2)

lines(curve_an_post$x, curve_an_post$y, lwd = 3, col = "firebrick3")
points(curve_an_post$x, curve_an_post$y, pch = 16, col = "firebrick3", cex = 1.2)

abline(0, 1, lwd = 2, lty = 2)

legend("topleft",
       legend = c("Pre-lockdown", "Post-lockdown"),
       col = c("royalblue4", "firebrick3"),
       lwd = 3, pch = 16, bty = "n")

################################################################################
#### END CODE ###
################################################################################