################################################################################
########################### NFHS-5 Children Analysis ###########################
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

child_keep <- c(
  "V001","V002","V003","B16",          # merge keys (+ child line number)
  "V005","V021","V022","V023",         # weights/PSU/strata/domain
  "V191A","V190A","V024",              # wealth + region
  "B8","B4","B5",                      # age years, sex, alive
  "HW1",                               # age months
  "HW2","HW3",                         # weight/height
  "HW72","HW57"                        # WFH z-score WHO, anemia level
)

hh_keep <- c(
  "HV001","HV002","HVIDX",
  "SHWEIGHT","HV005","HV021","HV022","HV023",
  "HV006","HV007",
  "HV270","HV270A","HV271","HV271A",
  "HV024"
)

nfhs_child_sub <- read_sav(
  "./IAKR7ESV_children/IAKR7EFL.SAV",
  col_select = all_of(child_keep)
)

nfhs_hh_children <- read_sav(
  "./IAPR7ESV_hh/IAPR7EFL.SAV",
  col_select = all_of(hh_keep)
)

################################################################################
# 2) Merge children + household file and rename variables
################################################################################

nfhs_hh_for_child <- nfhs_hh_children %>%
  rename(V001 = HV001, V002 = HV002, V003 = HVIDX)

nfhs_child_analysis <- nfhs_child_sub %>%
  left_join(nfhs_hh_for_child, by = c("V001","V002","V003")) %>%
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
    
    child_WFH_SD_WHO    = HW72,
    child_anemia_level  = HW57,
    
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
    household_region = HV024
  )

################################################################################
# 3) Clean special missing codes + rescale to real units
################################################################################

# Special missing codes
nfhs_child_analysis$child_weight_kgs[nfhs_child_analysis$child_weight_kgs >= 9990] <- NA
nfhs_child_analysis$child_height_cms[nfhs_child_analysis$child_height_cms >= 9990] <- NA
nfhs_child_analysis$child_WFH_SD_WHO[nfhs_child_analysis$child_WFH_SD_WHO >= 9990] <- NA

# Rescale decimals
nfhs_child_analysis$wi_factorscore_urbanrural <- nfhs_child_analysis$wi_factorscore_urbanrural / 100000
nfhs_child_analysis$child_weight_kgs <- nfhs_child_analysis$child_weight_kgs / 10
nfhs_child_analysis$child_height_cms <- nfhs_child_analysis$child_height_cms / 10
nfhs_child_analysis$child_WFH_SD_WHO <- nfhs_child_analysis$child_WFH_SD_WHO / 100

# Scale weights
nfhs_child_analysis$hh_state_weight <- nfhs_child_analysis$hh_state_weight / 1e6
nfhs_child_analysis$child_weight    <- nfhs_child_analysis$child_weight / 1e6
nfhs_child_analysis$hh_weight       <- nfhs_child_analysis$hh_weight / 1e6

# Ensure DHS wealth quintile is numeric 1-5
nfhs_child_analysis$wi_combined_urbanrural <- as.numeric(as.character(nfhs_child_analysis$wi_combined_urbanrural))

################################################################################
# 4) Restrict to alive children + outcomes + exposure
################################################################################

# alive only
nfhs_child_analysis <- nfhs_child_analysis %>%
  filter(child_alive == 1)

# Wasting: WFH z-score WHO < -2
nfhs_child_analysis$child_wasting <- ifelse(
  is.na(nfhs_child_analysis$child_WFH_SD_WHO), NA,
  ifelse(nfhs_child_analysis$child_WFH_SD_WHO < -2, 1, 0)
)
nfhs_child_analysis$child_wasting <- factor(
  nfhs_child_analysis$child_wasting,
  levels = c(0,1),
  labels = c("No wasting", "Wasting")
)

# Anemia: levels 1-3 vs 4
nfhs_child_analysis$child_anemia <- ifelse(
  is.na(nfhs_child_analysis$child_anemia_level), NA,
  ifelse(nfhs_child_analysis$child_anemia_level %in% c(1,2,3), 1, 0)
)
nfhs_child_analysis$child_anemia <- factor(
  nfhs_child_analysis$child_anemia,
  levels = c(0,1),
  labels = c("Not anemic", "Anemic")
)

# Date of interview + pre/post lockdown
nfhs_child_analysis$date_of_interview <- as.Date(
  paste(nfhs_child_analysis$year_of_interview, nfhs_child_analysis$month_of_interview, "01", sep = "-")
)

nfhs_child_analysis$postcovid <- ifelse(
  nfhs_child_analysis$date_of_interview < as.Date("2020-04-01"), "Pre-lockdown",
  ifelse(nfhs_child_analysis$date_of_interview >= as.Date("2020-06-01"), "Post-lockdown", NA)
)
nfhs_child_analysis$postcovid <- relevel(factor(nfhs_child_analysis$postcovid), ref = "Pre-lockdown")

################################################################################
# 5) Survey design (children weight: V005)
################################################################################

des_child_bio <- svydesign(
  ids = ~child_psu,
  strata = ~child_strata,
  weights = ~child_weight,
  data = nfhs_child_analysis,
  nest = TRUE
)

################################################################################
# 6) Weighted prevalence (overall + by wealth quintile) + chi-square
################################################################################

wprev <- function(des, expr) {
  svymean(as.formula(paste0("~I(", expr, ")")), des, na.rm = TRUE) * 100
}

prev_was_pre  <- wprev(subset(des_child_bio, postcovid == "Pre-lockdown"),  'child_wasting=="Wasting"')
prev_was_post <- wprev(subset(des_child_bio, postcovid == "Post-lockdown"), 'child_wasting=="Wasting"')

prev_an_pre  <- wprev(subset(des_child_bio, postcovid == "Pre-lockdown"),  'child_anemia=="Anemic"')
prev_an_post <- wprev(subset(des_child_bio, postcovid == "Post-lockdown"), 'child_anemia=="Anemic"')

was_by_q_pre <- prop.table(
  svytable(~I(child_wasting=="Wasting") + wi_combined_urbanrural,
           subset(des_child_bio, postcovid=="Pre-lockdown")),
  margin = 2
) * 100

was_by_q_post <- prop.table(
  svytable(~I(child_wasting=="Wasting") + wi_combined_urbanrural,
           subset(des_child_bio, postcovid=="Post-lockdown")),
  margin = 2
) * 100

an_by_q_pre <- prop.table(
  svytable(~I(child_anemia=="Anemic") + wi_combined_urbanrural,
           subset(des_child_bio, postcovid=="Pre-lockdown")),
  margin = 2
) * 100

an_by_q_post <- prop.table(
  svytable(~I(child_anemia=="Anemic") + wi_combined_urbanrural,
           subset(des_child_bio, postcovid=="Post-lockdown")),
  margin = 2
) * 100

des_child_bio <- update(
  des_child_bio,
  wasting_bin = (child_wasting == "Wasting"),
  anemia_bin  = (child_anemia == "Anemic")
)

chi_was <- svychisq(~wasting_bin + postcovid, des_child_bio, na.rm = TRUE)
chi_an  <- svychisq(~anemia_bin  + postcovid, des_child_bio, na.rm = TRUE)

################################################################################
# 7) Weighted odds ratios (overall + by wealth quintile)
# (no forest plots — you said you'll fill separately)
################################################################################

or_post <- function(fit) {
  term <- "postcovidPost-lockdown"
  est <- coef(fit)[term]
  ci  <- confint(fit)[term, ]
  c(OR = exp(est), CI_low = exp(ci[1]), CI_high = exp(ci[2]))
}

fit_was_all <- svyglm(wasting_bin ~ postcovid, design = des_child_bio, family = quasibinomial())
fit_an_all  <- svyglm(anemia_bin  ~ postcovid, design = des_child_bio, family = quasibinomial())

OR_was_all <- or_post(fit_was_all)
OR_an_all  <- or_post(fit_an_all)

fit_was_by_q <- lapply(1:5, function(q) {
  svyglm(wasting_bin ~ postcovid,
         design = subset(des_child_bio, wi_combined_urbanrural == q),
         family = quasibinomial())
})
OR_was_by_q <- do.call(rbind, lapply(fit_was_by_q, or_post))
rownames(OR_was_by_q) <- paste0("Q", 1:5)

fit_an_by_q <- lapply(1:5, function(q) {
  svyglm(anemia_bin ~ postcovid,
         design = subset(des_child_bio, wi_combined_urbanrural == q),
         family = quasibinomial())
})
OR_an_by_q <- do.call(rbind, lapply(fit_an_by_q, or_post))
rownames(OR_an_by_q) <- paste0("Q", 1:5)

################################################################################
# 8) Weighted concentration indices + curve data
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

df_pre  <- nfhs_child_analysis %>% filter(postcovid == "Pre-lockdown")
df_post <- nfhs_child_analysis %>% filter(postcovid == "Post-lockdown")

df_pre$wasting_bin_num  <- ifelse(df_pre$child_wasting == "Wasting", 1,
                                  ifelse(is.na(df_pre$child_wasting), NA, 0))
df_post$wasting_bin_num <- ifelse(df_post$child_wasting == "Wasting", 1,
                                  ifelse(is.na(df_post$child_wasting), NA, 0))

df_pre$anemia_bin_num  <- ifelse(df_pre$child_anemia == "Anemic", 1,
                                 ifelse(is.na(df_pre$child_anemia), NA, 0))
df_post$anemia_bin_num <- ifelse(df_post$child_anemia == "Anemic", 1,
                                 ifelse(is.na(df_post$child_anemia), NA, 0))

curve_was_pre  <- make_ci_curve(df_pre,  "wasting_bin_num", "child_weight", "wi_combined_urbanrural")
curve_was_post <- make_ci_curve(df_post, "wasting_bin_num", "child_weight", "wi_combined_urbanrural")

curve_an_pre   <- make_ci_curve(df_pre,  "anemia_bin_num",  "child_weight", "wi_combined_urbanrural")
curve_an_post  <- make_ci_curve(df_post, "anemia_bin_num",  "child_weight", "wi_combined_urbanrural")

CI_was_pre <- rineq::ci(
  ineqvar  = df_pre$wi_combined_urbanrural,
  outcome  = df_pre$wasting_bin_num,
  weights  = df_pre$child_weight,
  type     = "CI",
  method   = "direct"
)

CI_was_post <- rineq::ci(
  ineqvar  = df_post$wi_combined_urbanrural,
  outcome  = df_post$wasting_bin_num,
  weights  = df_post$child_weight,
  type     = "CI",
  method   = "direct"
)

CI_an_pre <- rineq::ci(
  ineqvar  = df_pre$wi_combined_urbanrural,
  outcome  = df_pre$anemia_bin_num,
  weights  = df_pre$child_weight,
  type     = "CI",
  method   = "direct"
)

CI_an_post <- rineq::ci(
  ineqvar  = df_post$wi_combined_urbanrural,
  outcome  = df_post$anemia_bin_num,
  weights  = df_post$child_weight,
  type     = "CI",
  method   = "direct"
)

################################################################################
# 9) Print results
################################################################################

cat("\n=============================\n")
cat("WEIGHTED PREVALENCE (OVERALL)\n")
cat("=============================\n")
cat("Wasting (Pre):\n");  print(prev_was_pre)
cat("Wasting (Post):\n"); print(prev_was_post)
cat("Anemia (Pre):\n");   print(prev_an_pre)
cat("Anemia (Post):\n");  print(prev_an_post)

cat("\n====================================\n")
cat("WEIGHTED PREVALENCE BY WEALTH QUINTILE\n")
cat("====================================\n")
cat("Wasting (Pre) % by quintile:\n");  print(was_by_q_pre)
cat("Wasting (Post) % by quintile:\n"); print(was_by_q_post)
cat("Anemia (Pre) % by quintile:\n");   print(an_by_q_pre)
cat("Anemia (Post) % by quintile:\n");  print(an_by_q_post)

cat("\n=============================\n")
cat("WEIGHTED CHI-SQUARE TESTS\n")
cat("=============================\n")
cat("Wasting vs postcovid:\n"); print(chi_was)
cat("Anemia vs postcovid:\n");  print(chi_an)

cat("\n=============================\n")
cat("WEIGHTED ODDS RATIOS\n")
cat("=============================\n")
cat("Overall OR (Wasting):\n"); print(OR_was_all)
cat("Overall OR (Anemia):\n");  print(OR_an_all)
cat("\nORs by wealth quintile (Wasting):\n"); print(OR_was_by_q)
cat("\nORs by wealth quintile (Anemia):\n");  print(OR_an_by_q)

cat("\n=============================\n")
cat("CONCENTRATION INDICES (WEIGHTED)\n")
cat("=============================\n")
cat("CI Wasting (Pre):\n");  print(CI_was_pre)
cat("CI Wasting (Post):\n"); print(CI_was_post)
cat("CI Anemia (Pre):\n");   print(CI_an_pre)
cat("CI Anemia (Post):\n");  print(CI_an_post)

################################################################################
# 10) Plots and Figures
## Odds ratio: forest plots
# Concentration curves: PRE vs POST overlapping
################################################################################

## Odds ratio: forest plots


# Wasting
wasting_children <- data.frame(
  Group = c("Q1", "Q2", "Q3", "Q4", "Q5", "Overall dataset"),
  OR = c(0.85, 0.81, 0.79, 0.75, 0.67, 0.77),
  Lower_CI = c(0.80, 0.76, 0.73, 0.68, 0.58, 0.74),
  Upper_CI = c(0.91, 0.87, 0.86, 0.81, 0.77, 0.80)
)

wasting_children <- wasting_children %>%
  mutate(
    RR_label = sprintf("%.2f (%.2f–%.2f)", OR, Lower_CI, Upper_CI),
    Label_Pos = pmin(Upper_CI + 0.01, 3.3),  # To keep ratios next to respective plots
    Characteristic = factor(Group, levels = rev(Group))  # Preserve order
  )

ggplot(wasting_children, aes(y = Group, x = OR)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_text(aes(x = Label_Pos, label = RR_label), hjust = 0, size = 3.5) +
  theme_minimal() +
  labs(title = "Wasting in Children",
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
anemia_children <- data.frame(
  Group = c("Q1", "Q2", "Q3", "Q4", "Q5", "Overall dataset"),
  OR = c(0.75, 0.77, 0.76, 0.76, 0.91, 0.78),
  Lower_CI = c(0.71, 0.72, 0.71, 0.71, 0.85, 0.75),
  Upper_CI = c(0.80, 0.82, 0.81, 0.82, 0.98, 0.81)
)

#
ggplot(anemia_children, aes(y = Group, x = OR)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  theme_minimal() +
  labs(title = "Anemia in Children",
       x = "OR and 95% CI",
       y = "Group") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#1f9e89") +  
  theme(plot.title = element_text(hjust = 0.5))
#

anemia_children <- anemia_children %>%
  mutate(
    RR_label = sprintf("%.2f (%.2f–%.2f)", OR, Lower_CI, Upper_CI),
    Label_Pos = pmin(Upper_CI + 0.01, 3.3),  # To keep ratios next to respective plots
    Characteristic = factor(Group, levels = rev(Group))  # Preserve order
  )

ggplot(anemia_children, aes(y = Group, x = OR)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2) +
  geom_text(aes(x = Label_Pos, label = RR_label), hjust = 0, size = 3.5) +
  theme_minimal() +
  labs(title = "Anemia in Children",
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

# Wasting curve (Pre + Post)
plot(curve_was_pre$x, curve_was_pre$y,
     type = "l", lwd = 3, col = "royalblue4",
     main = "Wasting in Children Under Five (Pre vs Post Lockdown)",
     xlab = "Cumulative proportion of children (weighted)",
     ylab = "Cumulative proportion wasting (weighted)",
     xlim = c(0,1), ylim = c(0,1)
)
points(curve_was_pre$x, curve_was_pre$y, pch = 16, col = "royalblue4", cex = 1.2)

lines(curve_was_post$x, curve_was_post$y, lwd = 3, col = "firebrick3")
points(curve_was_post$x, curve_was_post$y, pch = 16, col = "firebrick3", cex = 1.2)

abline(0, 1, lwd = 2, lty = 2)

legend("topleft",
       legend = c("Pre-lockdown", "Post-lockdown"),
       col = c("royalblue4", "firebrick3"),
       lwd = 3, pch = 16, bty = "n")

# Anemia curve (Pre + Post)
plot(curve_an_pre$x, curve_an_pre$y,
     type = "l", lwd = 3, col = "royalblue4",
     main = "Anemia in Children Under Five (Pre vs Post Lockdown)",
     xlab = "Cumulative proportion of children (weighted)",
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