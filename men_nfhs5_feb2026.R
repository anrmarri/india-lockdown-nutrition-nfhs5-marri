################################################################################
############################# NFHS-5 Men Analysis ##############################
########################### Author: anrmarri@bu.edu ############################
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

men_keep <- c(
  "MV001","MV002","MV003",          # merge keys
  "MV005","MV021","MV022","MV023",  # weights/PSU/strata/domain
  "MV191A","MV190A","MV024"         # wealth + region
)

hh_keep <- c(
  "HV001","HV002","HVIDX",          # merge keys
  "SHMWEIGHT", "HV005","HV021","HV022","HV023",
  "HV006","HV007",                  # month/year
  "HV270","HV270A","HV271","HV271A",# wealth
  "HB1","HB2","HB3","HB35","HB40","HB56","HB57", # biomarkers/anthro
  "HV024"                           # region
)

nfhs_men_sub <- read_sav(
  "./IAMR7EDT_men/IAMR7EFL.SAV",
  col_select = all_of(men_keep)
)

nfhs_hh_men <- read_sav(
  "./IAPR7ESV_hh/IAPR7EFL.SAV",
  col_select = all_of(hh_keep)
)

################################################################################
# 2) Merge men + household subsample file and rename variables
################################################################################

nfhs_hh_for_men <- nfhs_hh_men %>%
  rename(MV001 = HV001, MV002 = HV002, MV003 = HVIDX)

nfhs_men_analysis <- nfhs_men_sub %>%
  left_join(nfhs_hh_for_men, by = c("MV001","MV002","MV003")) %>%
  rename(
    month_of_interview = HV006,
    year_of_interview  = HV007,
    wi_factorscore = HV271,
    wi_factorscore_urbanrural = HV271A,
    wi_combined = HV270,
    wi_combined_urbanrural = HV270A,
    
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
    man_wi_factorscore = MV191A,
    man_wi_quintile = MV190A,
    
    man_hh_state_weight = SHMWEIGHT,
    hh_weight = HV005,
    hh_psu = HV021,
    hh_strata = HV022,
    hh_sample_domain = HV023,
    
    man_region = MV024,
    household_region = HV024
  )

################################################################################
# 3) Clean special missing codes + rescale to real units
################################################################################

# Special missing codes
nfhs_men_analysis$man_weight_kgs[nfhs_men_analysis$man_weight_kgs >= 9990] <- NA
nfhs_men_analysis$man_height_cms[nfhs_men_analysis$man_height_cms >= 9990] <- NA
nfhs_men_analysis$man_bmi[nfhs_men_analysis$man_bmi >= 9990] <- NA
nfhs_men_analysis$man_hb_gdl_altsmoking[nfhs_men_analysis$man_hb_gdl_altsmoking >= 990] <- NA

# Rescale decimals
nfhs_men_analysis$wi_factorscore_urbanrural <- nfhs_men_analysis$wi_factorscore_urbanrural / 100000
nfhs_men_analysis$man_weight_kgs <- nfhs_men_analysis$man_weight_kgs / 10
nfhs_men_analysis$man_height_cms <- nfhs_men_analysis$man_height_cms / 10
nfhs_men_analysis$man_bmi <- nfhs_men_analysis$man_bmi / 100
nfhs_men_analysis$man_hb_gdl_altsmoking <- nfhs_men_analysis$man_hb_gdl_altsmoking / 10

# Scale weights
nfhs_men_analysis$man_hh_state_weight <- nfhs_men_analysis$man_hh_state_weight / 1e6
nfhs_men_analysis$man_weight <- nfhs_men_analysis$man_weight / 1e6
nfhs_men_analysis$hh_weight <- nfhs_men_analysis$hh_weight / 1e6

# Ensure DHS wealth quintile is numeric 1-5
nfhs_men_analysis$wi_combined_urbanrural <- as.numeric(as.character(nfhs_men_analysis$wi_combined_urbanrural))


################################################################################
# 4) Outcomes + exposure
################################################################################

# Undernutrition (BMI < 18.5) using DHS-reported BMI
nfhs_men_analysis$man_malnourished <- ifelse(
  is.na(nfhs_men_analysis$man_bmi), NA,
  ifelse(nfhs_men_analysis$man_bmi < 18.5, 1, 0)
)
nfhs_men_analysis$man_malnourished <- factor(
  nfhs_men_analysis$man_malnourished,
  levels = c(0,1),
  labels = c("Not undernourished", "Undernourished")
)

# Anemia (levels 1-3 vs 4)
nfhs_men_analysis$man_anemia <- ifelse(
  is.na(nfhs_men_analysis$man_anemia_level), NA,
  ifelse(nfhs_men_analysis$man_anemia_level %in% c(1,2,3), 1, 0)
)
nfhs_men_analysis$man_anemia <- factor(
  nfhs_men_analysis$man_anemia,
  levels = c(0,1),
  labels = c("Not anemic", "Anemic")
)

# Date of interview + pre/post lockdown definition
nfhs_men_analysis$date_of_interview <- as.Date(
  paste(nfhs_men_analysis$year_of_interview, nfhs_men_analysis$month_of_interview, "01", sep = "-")
)

nfhs_men_analysis$postcovid <- ifelse(
  nfhs_men_analysis$date_of_interview < as.Date("2020-04-01"), "Pre-lockdown",
  ifelse(nfhs_men_analysis$date_of_interview >= as.Date("2020-06-01"), "Post-lockdown", NA)
)
nfhs_men_analysis$postcovid <- relevel(factor(nfhs_men_analysis$postcovid), ref = "Pre-lockdown")

################################################################################
# 5) Survey design (biomarker subsample weight: SHMWEIGHT)
################################################################################

des_men_bio <- svydesign(
  ids = ~man_psu,
  strata = ~man_strata,
  weights = ~man_hh_state_weight,
  data = nfhs_men_analysis,
  nest = TRUE
)


################################################################################
# 6) Weighted prevalence (overall + by wealth quintile) + chi-square
################################################################################

# weighted prevalence (%) for a binary indicator
wprev <- function(des, expr) {
  svymean(as.formula(paste0("~I(", expr, ")")), des, na.rm = TRUE) * 100
}

# Overall prevalence
prev_under_pre  <- wprev(subset(des_men_bio, postcovid == "Pre-lockdown"),  'man_malnourished=="Undernourished"')
prev_under_post <- wprev(subset(des_men_bio, postcovid == "Post-lockdown"), 'man_malnourished=="Undernourished"')

prev_an_pre  <- wprev(subset(des_men_bio, postcovid == "Pre-lockdown"),  'man_anemia=="Anemic"')
prev_an_post <- wprev(subset(des_men_bio, postcovid == "Post-lockdown"), 'man_anemia=="Anemic"')

# By wealth quintile (weighted % within quintile)
under_by_q_pre <- prop.table(
  svytable(~I(man_malnourished=="Undernourished") + wi_combined_urbanrural,
           subset(des_men_bio, postcovid=="Pre-lockdown")),
  margin = 2
) * 100

under_by_q_post <- prop.table(
  svytable(~I(man_malnourished=="Undernourished") + wi_combined_urbanrural,
           subset(des_men_bio, postcovid=="Post-lockdown")),
  margin = 2
) * 100

an_by_q_pre <- prop.table(
  svytable(~I(man_anemia=="Anemic") + wi_combined_urbanrural,
           subset(des_men_bio, postcovid=="Pre-lockdown")),
  margin = 2
) * 100

an_by_q_post <- prop.table(
  svytable(~I(man_anemia=="Anemic") + wi_combined_urbanrural,
           subset(des_men_bio, postcovid=="Post-lockdown")),
  margin = 2
) * 100

# Weighted chi-square (overall pre vs post)
des_men_bio <- update(
  des_men_bio,
  undernourished_bin = (man_malnourished == "Undernourished"),
  anemia_bin         = (man_anemia == "Anemic")
)

chi_under <- svychisq(~undernourished_bin + postcovid, des_men_bio, na.rm = TRUE)
chi_an    <- svychisq(~anemia_bin + postcovid,         des_men_bio, na.rm = TRUE)

################################################################################
# 7) Weighted odds ratios (overall + by wealth quintile)
################################################################################

# OR + 95% CI for postcovid effect from svyglm
or_post <- function(fit) {
  term <- "postcovidPost-lockdown"
  est <- coef(fit)[term]
  ci  <- confint(fit)[term, ]
  c(OR = exp(est), CI_low = exp(ci[1]), CI_high = exp(ci[2]))
}

# Overall ORs
fit_under_all <- svyglm(undernourished_bin ~ postcovid, design = des_men_bio, family = quasibinomial())
fit_an_all    <- svyglm(anemia_bin ~ postcovid,         design = des_men_bio, family = quasibinomial())

OR_under_all <- or_post(fit_under_all)
OR_an_all    <- or_post(fit_an_all)

# By wealth quintile ORs
fit_under_by_q <- lapply(1:5, function(q) {
  svyglm(undernourished_bin ~ postcovid,
         design = subset(des_men_bio, wi_combined_urbanrural == q),
         family = quasibinomial())
})
OR_under_by_q <- do.call(rbind, lapply(fit_under_by_q, or_post))
rownames(OR_under_by_q) <- paste0("Q", 1:5)

fit_an_by_q <- lapply(1:5, function(q) {
  svyglm(anemia_bin ~ postcovid,
         design = subset(des_men_bio, wi_combined_urbanrural == q),
         family = quasibinomial())
})
OR_an_by_q <- do.call(rbind, lapply(fit_an_by_q, or_post))
rownames(OR_an_by_q) <- paste0("Q", 1:5)

################################################################################
# 8) Weighted concentration indices + curve data (FIXED)
################################################################################

# Build weighted concentration curve points (cumulative population share vs cumulative case share)
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
  
  # Add origin row (0,0) without using := on .data
  origin <- out[1, ]
  origin[,] <- 0
  origin[[q]] <- 0
  
  bind_rows(origin, out)
}

df_pre  <- nfhs_men_analysis %>% filter(postcovid == "Pre-lockdown")
df_post <- nfhs_men_analysis %>% filter(postcovid == "Post-lockdown")

df_pre$man_mal_bin     <- ifelse(df_pre$man_malnourished == "Undernourished", 1, ifelse(is.na(df_pre$man_malnourished), NA, 0))
df_post$man_mal_bin    <- ifelse(df_post$man_malnourished == "Undernourished", 1, ifelse(is.na(df_post$man_malnourished), NA, 0))
df_pre$man_anemia_bin  <- ifelse(df_pre$man_anemia == "Anemic", 1, ifelse(is.na(df_pre$man_anemia), NA, 0))
df_post$man_anemia_bin <- ifelse(df_post$man_anemia == "Anemic", 1, ifelse(is.na(df_post$man_anemia), NA, 0))

curve_under_pre  <- make_ci_curve(df_pre,  "man_mal_bin",     "man_hh_state_weight", "wi_combined_urbanrural")
curve_under_post <- make_ci_curve(df_post, "man_mal_bin",     "man_hh_state_weight", "wi_combined_urbanrural")
curve_an_pre     <- make_ci_curve(df_pre,  "man_anemia_bin",  "man_hh_state_weight", "wi_combined_urbanrural")
curve_an_post    <- make_ci_curve(df_post, "man_anemia_bin",  "man_hh_state_weight", "wi_combined_urbanrural")

# Weighted concentration indices (rineq)
CI_under_pre <- rineq::ci(
  ineqvar  = df_pre$wi_combined_urbanrural,
  outcome  = df_pre$man_mal_bin,
  weights  = df_pre$man_hh_state_weight,
  type     = "CI",
  method   = "direct"
)

CI_under_post <- rineq::ci(
  ineqvar  = df_post$wi_combined_urbanrural,
  outcome  = df_post$man_mal_bin,
  weights  = df_post$man_hh_state_weight,
  type     = "CI",
  method   = "direct"
)

CI_an_pre <- rineq::ci(
  ineqvar  = df_pre$wi_combined_urbanrural,
  outcome  = df_pre$man_anemia_bin,
  weights  = df_pre$man_hh_state_weight,
  type     = "CI",
  method   = "direct"
)

CI_an_post <- rineq::ci(
  ineqvar  = df_post$wi_combined_urbanrural,
  outcome  = df_post$man_anemia_bin,
  weights  = df_post$man_hh_state_weight,
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
## Concentration curves: PRE vs POST overlapping
################################################################################

## Odds ratio: forest plots

# Undernutrition
undernutrition_men <- data.frame(
  Group = c("Q1", "Q2", "Q3", "Q4", "Q5", "Overall dataset"),
  OR = c(1.08, 1.30, 1.09, 0.92, 0.94, 1.07),
  Lower_CI = c(0.98, 1.16, 0.97, 0.81, 0.82, 1.01),
  Upper_CI = c(1.20, 1.44, 1.22, 1.04, 1.09, 1.13)
)

undernutrition_men <- undernutrition_men %>%
  mutate(
    RR_label = sprintf("%.2f (%.2f–%.2f)", OR, Lower_CI, Upper_CI),
    Label_Pos = pmin(Upper_CI + 0.01, 3.3),  # To keep ratios next to respective plots
    Characteristic = factor(Group, levels = rev(Group))  # Preserve order
  )

ggplot(undernutrition_men, aes(y = Group, x = OR)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_text(aes(x = Label_Pos, label = RR_label), hjust = 0, size = 3.5) +
  theme_minimal() +
  labs(title = "Underweight in Men",
       x = "OR and 95% CI",
       y = "Group") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#2E0854") +
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

# Anemia
anemia_men <- data.frame(
  Group = c("Q1", "Q2", "Q3", "Q4", "Q5", "Overall dataset"),
  OR = c(0.77, 0.77, 0.76, 0.76, 0.79, 0.77),
  Lower_CI = c(0.70, 0.70, 0.68, 0.68, 0.70, 0.73),
  Upper_CI = c(0.84, 0.86, 0.84, 0.84, 0.88, 0.81)
)

##
ggplot(anemia_men, aes(y = Group, x = OR)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  theme_minimal() +
  labs(title = "Anemia in Men",
       x = "OR and 95% CI",
       y = "Group") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#1f9e89") +  
  theme(plot.title = element_text(hjust = 0.5))
##

anemia_men <- anemia_men %>%
  mutate(
    RR_label = sprintf("%.2f (%.2f–%.2f)", OR, Lower_CI, Upper_CI),
    Label_Pos = pmin(Upper_CI + 0.01, 3.3),  # To keep ratios next to respective plots
    Characteristic = factor(Group, levels = rev(Group))  # Preserve order
  )

ggplot(anemia_men, aes(y = Group, x = OR)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_text(aes(x = Label_Pos, label = RR_label), hjust = 0, size = 3.5) +
  theme_minimal() +
  labs(title = "Anemia in Men",
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

## Concentration curves: PRE vs POST overlapping
# Underweight curve (Pre + Post on same plot)
plot(curve_under_pre$x, curve_under_pre$y,
     type = "l",
     lwd = 3,
     col = "royalblue4",
     main = "Underweight in Men (Pre vs Post Lockdown)",
     xlab = "Cumulative proportion of men (weighted)",
     ylab = "Cumulative proportion underweight (weighted)",
     xlim = c(0,1), ylim = c(0,1)
)

# Add PRE points
points(curve_under_pre$x, curve_under_pre$y,
       pch = 16,          # solid circle
       col = "royalblue4",
       cex = 1.2)

# Add POST line
lines(curve_under_post$x, curve_under_post$y,
      lwd = 3,
      col = "firebrick3")

# Add POST points
points(curve_under_post$x, curve_under_post$y,
       pch = 16,
       col = "firebrick3",
       cex = 1.2)

# Line of equality
abline(0, 1, lwd = 2, lty = 2)

legend("topleft",
       legend = c("Pre-lockdown", "Post-lockdown"),
       col = c("royalblue4", "firebrick3"),
       lwd = 3,
       pch = 16,
       bty = "n")


# Anemia curve (Pre + Post on same plot)
plot(curve_an_pre$x, curve_an_pre$y,
     type = "l",
     lwd = 3,
     col = "royalblue4",
     main = "Anemia in Men (Pre vs Post Lockdown)",
     xlab = "Cumulative proportion of men (weighted)",
     ylab = "Cumulative proportion anemic (weighted)",
     xlim = c(0,1), ylim = c(0,1)
)

# PRE points
points(curve_an_pre$x, curve_an_pre$y,
       pch = 16,
       col = "royalblue4",
       cex = 1.2)

# POST line
lines(curve_an_post$x, curve_an_post$y,
      lwd = 3,
      col = "firebrick3")

# POST points
points(curve_an_post$x, curve_an_post$y,
       pch = 16,
       col = "firebrick3",
       cex = 1.2)

# Line of equality
abline(0, 1, lwd = 2, lty = 2)

legend("topleft",
       legend = c("Pre-lockdown", "Post-lockdown"),
       col = c("royalblue4", "firebrick3"),
       lwd = 3,
       pch = 16,
       bty = "n")

################################################################################
#### END CODE ###
################################################################################