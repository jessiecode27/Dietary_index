### data cleaning for covariates ####################################################################
### Major cardiovascular disease (M_CVD: having at least 1 disease = 1; all of these conditions are reported as “No” = 0; 
# all five conditions are missing or coded as “Refused/Don’t know”) defined as having at least one of the following conditions ###

# MCQ160C : Coronary heart disease (CHD)
# MCQ160D : Angina/angina pectoris
# MCQ160B : Congestive heart failure (CHF)
# MCQ160E : Myocardial infarction (heart attack, MI)
# MCQ160F : Stroke

# List NHANES raw variables and new labels
vars_cvd <- c("MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F")
labs_cvd <- c("CHF","CHD","Angina","MI","Stroke")

# Recode: 1 = Yes, 2 = No, 7/9/. = NA
nhanes_merged_0518[labs_cvd] <- lapply(
  nhanes_merged_0518[vars_cvd],
  function(x) ifelse(x == 1, 1L,
                     ifelse(x == 2, 0L, NA_integer_)))

# Define Major CVD
nhanes_merged_0518$M_CVD <- apply(
  nhanes_merged_0518[labs_cvd], 1L,
  function(row) {
    if (all(is.na(row))) NA_integer_                # all rows NA is NA
    else if (any(row == 1L, na.rm = TRUE)) 1L # having at least one of disease
    else 0L                                 
  }
)

### Type 2 diabetes
## Total diabetes = Diagnosed diabetes - self-reported diabetes (DIQ010 == 1) 
#                   OR Fasting plasma glucose (LBXGLU) ≥ 126 mg/dL OR HbA1c (LBXGH) ≥ 6.5%
# AND People with missing values for either fasting glucose or A1C and pregnant women were excluded => will be included in the exclusion criteria 

# Excluded T1D: diabetes diagnosed before the age of 30 years (DID040<30) 
#               AND who was taking only insulin therapy (DIQ050, DID060, DIQ070)
# AND during pregnancy, told you have diabetes (RHQ162==1 only in cycles from 2007) | Gestational diabetes (DIQ175S==28 only in cycles from 2011)

library(dplyr)

nhanes_merged_0518 <- nhanes_merged_0518 %>%
  mutate(
    # Diabetes: DIQ010==1 OR FPG≥126 OR HbA1c≥6.5
    diabetes = ifelse(
      DIQ010 == 1 | LBXGLU >= 126 | LBXGH >= 6.5, 1L,
      ifelse(DIQ010 == 2 & (is.na(LBXGLU) | LBXGLU < 126) & (is.na(LBXGH) | LBXGH < 6.5),
        0L,
        NA_integer_)),
    
    # T1D flag: diabetes==1 & DID040<30 & ONLY insulin (DIQ050==1 AND DID060 in 1:49 AND DIQ070==2)
    t1d_core = (
      ifelse(is.na(diabetes), FALSE, diabetes == 1L) &
        !is.na(DID040) & DID040 < 30 &
        DIQ050 == 1L &                 # now taking insulin
        DIQ070 == 2L &                 # not taking pills
        DID060 %in% 1:49               # used insulin more than 1 month: exclude <1 month (666), refuse/dk (777/999)
    ),
    
    T1D_flag = ifelse(t1d_core, 1L,
                      ifelse(diabetes == 0L, 0L, 0L)),  
    
    # T2D: diabetes==1 AND T1D_flag==0L (NOT T1D) AND RHQ162==2 | DIQ175S==28 (NOT gestational diabetes)
    T2D = ifelse(
      diabetes == 1L & T1D_flag == 0L & (RHQ162 == 2 | DIQ175S !=28), 1L,
      ifelse(diabetes == 0L | RHQ162 == 1 | DIQ175S==28, 0L, NA_integer_)
    )
  ) #Note: If exclude pregnant women from the beginning, code should be refined


### Cancer 
# Recode: 1 = Yes, 2 = No, 7/9/. = NA
nhanes_merged_0518$cancer <- ifelse(
  nhanes_merged_0518$MCQ220 == 1, 1L,
  ifelse(nhanes_merged_0518$MCQ220 == 2, 0L, NA_integer_)
)

### Major chronic disease (combined variable) (MCD) ##################################
# defined as having at least one of the following conditions: CVD, T2D, cancer
nhanes_merged_0518$MCD <- apply(
  nhanes_merged_0518[, c("M_CVD", "T2D", "cancer")], 1,
  function(row) {
    if (all(is.na(row))) NA_integer_
    else if (any(row == 1, na.rm = TRUE)) 1L
    else 0L
  }
)

nhanes_merged_0518 <- nhanes_merged_0518 %>%
  mutate(
    MCD_count = rowSums(across(c(M_CVD, T2D, cancer)), na.rm = TRUE), # Count MCD
    MCD = case_when(
      rowSums(!is.na(across(c(M_CVD, T2D, cancer)))) == 0L ~ NA_integer_, # all NA
      MCD_count >= 1 ~ 1L,  
      TRUE ~ 0L             
    )
  )

### Metabolic syndrome (MetS) #######################################################
# defined as having at least one of the following conditions:
# 1. Abdominal obesity (Waist circumference (cm): BMXWAIST): >40 in (men), >35 in (women) => >101.6 (men), >88.9 (women)
# 2. High blood pressure (mm Hg) 
# BPXSY1:BPXSY4, BPXDI1:BPXDI4 (create average): SBP (Systolic Blood Pressure) ≥130 or DBP (Diastolic Blood Pressure) ≥80
# or BPQ020==1 (Ever told you had high blood pressure)
# 3. Impaired fasting glucose (mg/dL) (LBXGLU) : ≥100
# 4. High triglycerides (mg/dL) (LBXTR): >150 
# 5. Low HDL (HDL cholesterol (mg/dL): LBDHDD) <40 mg/dL (men), <50 mg/dL (women)
# Gender (RIAGENDR): 1 Male, 2 Female

library(dplyr)

nhanes_merged_0518 <- nhanes_merged_0518 %>%
  mutate(
    # --- Abdominal obesity (cm) ---
    wc_high = case_when(
      RIAGENDR == 1L & BMXWAIST > 101.6 ~ 1L,  # men >101.6 cm
      RIAGENDR == 2L & BMXWAIST > 88.9  ~ 1L,  # women >88.9 cm
      RIAGENDR %in% c(1L,2L) & !is.na(BMXWAIST) ~ 0L,               # measured but below threshold
      TRUE ~ NA_integer_                                            # missing waist -> NA
    ),
    
    # --- High blood pressure (use mean of 4 times of measurements) ---
    SBP_mean = {
      x <- rowMeans(as.matrix(across(any_of(c("BPXSY1","BPXSY2","BPXSY3","BPXSY4")))), na.rm = TRUE) #average available SBP readings
      ifelse(is.nan(x), NA_real_, x)   # all-NA -> NA
    },
    DBP_mean = {
      x <- rowMeans(as.matrix(across(any_of(c("BPXDI1","BPXDI2","BPXDI3","BPXDI4")))), na.rm = TRUE)
      ifelse(is.nan(x), NA_real_, x)
    },
    
    bp_high = case_when(
      SBP_mean >= 130 | DBP_mean >= 80  | BPQ020 == 1L ~ 1L,  # or self-reported history of high BP (BPQ020)
      !is.na(SBP_mean) | !is.na(DBP_mean) | BPQ020 %in% c(1L,2L) ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    # --- Impaired fasting glucose ---
    glu_high = case_when(
      LBXGLU >= 100 ~ 1L,
      !is.na(LBXGLU) ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    # --- High triglycerides ---
    tg_high = case_when(
      LBXTR > 150 ~ 1L,
      !is.na(LBXTR) ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    # --- Low HDL ---
    hdl_low = case_when(
      RIAGENDR == 1L & LBDHDD < 40 ~ 1L,  # men
      RIAGENDR == 2L & LBDHDD < 50 ~ 1L,  # women
      RIAGENDR %in% c(1L,2L) & !is.na(LBDHDD) ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    # Final Metabolic Syndrome (≥3 of 5)
    # NA if all five components are NA
    MetS_count = rowSums(across(c(wc_high, bp_high, glu_high, tg_high, hdl_low)), na.rm = TRUE),
    MetS = case_when(
      rowSums(!is.na(across(c(wc_high, bp_high, glu_high, tg_high, hdl_low)))) == 0L ~ NA_integer_,
      MetS_count >= 3 ~ 1L,
      TRUE ~ 0L
    )
  )

View(nhanes_merged_0518)
####### Covariates ############
###############################
library(dplyr)

nhanes_merged_0518 <- nhanes_merged_0518 %>%
  mutate(
    # ---- Age group ----
    RIDAGEYR_CAT = case_when(
      RIDAGEYR >= 20 & RIDAGEYR <= 34 ~ 1L,
      RIDAGEYR >= 35 & RIDAGEYR <= 49 ~ 2L,
      RIDAGEYR >= 50 & RIDAGEYR <= 64 ~ 3L,
      RIDAGEYR >= 65                  ~ 4L,
      TRUE ~ NA_integer_
    ),
    
    # ---- Education (DMDEDUC2) ----
    DMDEDUC2_CAT = case_when(
      DMDEDUC2 %in% c(1,2) ~ 2L,   # <High School
      DMDEDUC2 == 3        ~ 3L,   # High school diploma or GED
      DMDEDUC2 == 4        ~ 4L,   # Some college or Associate’s degree
      DMDEDUC2 == 5        ~ 5L,   # College graduate or above
      DMDEDUC2 %in% c(7,9) ~ NA_integer_,
      TRUE ~ NA_integer_
    ),
    
    # ---- Marital (DMDMARTL) ----
    DMDMARTL_CAT = case_when(
      DMDMARTL == 1           ~ 1L,  # Married
      DMDMARTL %in% 2:6       ~ 0L,  # Not married (includes widowed, divorced, separated, never married, and living with a partner)
      DMDMARTL %in% c(77,99)  ~ NA_integer_,
      TRUE ~ NA_integer_
    ),
    
    # ---- Family Income to Poverty Ratio (INDFMPIR) ----
    INDFMPIR_CAT = case_when(
      INDFMPIR < 1.30   ~ 1L,
      INDFMPIR < 3.50   ~ 2L,
      INDFMPIR >= 3.50  ~ 3L,
      TRUE ~ NA_integer_
    ),
    
    # ---- BMI ----
    BMI_CAT = case_when(
      BMXBMI < 25               ~ 2L,
      BMXBMI >=25 & BMXBMI <30  ~ 3L,
      BMXBMI >=30               ~ 4L,
      TRUE ~ NA_integer_
    ),
    OBESE = case_when(
      BMXBMI >=30 ~ 1L,
      BMXBMI <30  ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    # ---- Abdominal obesity ---- (already have)
    AB_OBESE = case_when(
      RIAGENDR == 1 & BMXWAIST >= 102 ~ 1L,
      RIAGENDR == 2 & BMXWAIST >=  88 ~ 1L,
      !is.na(BMXWAIST)                ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    # ---- Smoking ----
    SMQ_CAT = case_when(
      SMQ020 == 2                          ~ 0L, # Never
      SMQ020 == 1 & SMQ040 == 3            ~ 1L, # Former
      SMQ020 == 1 & SMQ040 %in% c(1,2)     ~ 2L, # Current
      SMQ020 %in% c(7,9) | SMQ040 %in% c(7,9) ~ NA_integer_,
      TRUE ~ NA_integer_
    ),
    
    # ---- Alcohol ----
    ALQ130_CAT = case_when(
      ALQ130 %in% c(777,999) ~ NA_real_,
      !is.na(ALQ130)         ~ as.numeric(ALQ130),
      TRUE                   ~ NA_real_
    ),
    
    # ---- Physical activity ----
    PAQ_MINWEEK = case_when(
      !is.na(PAQ670) & !is.na(PAD675) ~ as.numeric(PAQ670) * as.numeric(PAD675),
      TRUE ~ NA_real_
    ),
    PAQ_CAT = case_when(
      PAQ_MINWEEK == 0                       ~ 0L, # Inactive
      PAQ_MINWEEK > 0 & PAQ_MINWEEK <150     ~ 1L, # Insufficient
      PAQ_MINWEEK >=150 & PAQ_MINWEEK <300   ~ 2L, # Sufficient
      PAQ_MINWEEK >=300                      ~ 3L, # Mod+vigorous
      TRUE ~ NA_integer_
    )
    
    # ---- Supplement use ----
    #DSD010_CAT = case_when(
    #  DSD010 == 1         ~ 1L,
    #  DSD010 == 2         ~ 0L,
    #  DSD010 %in% c(7,9)  ~ NA_integer_,
    #  TRUE ~ NA_integer_
    #)
  )

### label for covariates ###
library(dplyr)

f <- function(x, lv, lab) factor(x, levels = lv, labels = lab)

nhanes_merged_0518 <- nhanes_merged_0518 %>%
  mutate(
    RIAGENDR     = f(RIAGENDR, c(1, 2), c("Male", "Female")),
    RIDAGEYR_CAT = f(RIDAGEYR_CAT, 1:4, c("20-34", "35-49", "50-64", "≥65")),
    RIDRETH1     = f(RIDRETH1, 1:5,
                     c("Mexican American","Other Hispanic","Non-Hispanic White",
                       "Non-Hispanic Black","Other Race")),
    DMDEDUC2_CAT = f(DMDEDUC2_CAT, c(2,3,4,5),
                     c("Less Than High School","High School Graduate/GED",
                       "Some College or Associate's Degree","College Graduate or Above")),
    DMDMARTL_CAT = f(DMDMARTL_CAT, c(1,0), c("Married","Not Married")),
    INDFMPIR_CAT = f(INDFMPIR_CAT, 1:3, c("0–<1.30","1.30–3.49","≥3.50")),
    BMI_CAT      = f(BMI_CAT, c(2,3,4), c("Underweight or Normal Weight","Overweight","Obese"))
  )

################################################################################################
# Check missing data (loading....)
################################################################################################
###### Ex ########################
library(janitor)
table_tabyl <- nhanes_merged_0518 %>%
  tabyl(MCD)
print(table_tabyl)

#nhanes_merged_0518 %>%
#  filter(as.numeric(RIAGENDR) >= 18
#  )

##################################
# visualization of missing data 
#################################
# Check missing data and distribution
library(DataExplorer)

create_report(nhanes_merged_0518,
              output_file = "nhanes_profile_report.html",
              output_dir = getwd())

library(naniar)
library(ggplot2)
library(dplyr)
library(scales)

# all variables
gg_miss_var(nhanes_merged_0518, show_pct = TRUE) + #%
  labs(y = "% missing", x = "Variables") +
  coord_flip()
gg_miss_var(nhanes_merged_0518) #n

# diet + health vars
vars_to_plot <- c("HEI2020_ALL","AHEI_ALL","DASH_ALL","MED_ALL","DII_ALL",
                  "CHF","CHD","Angina","MI","Stroke",
                  "M_CVD","T2D","cancer", "MCD", "MetS")

gg_miss_var(
  nhanes_merged_0518 %>% select(any_of(vars_to_plot)),
  show_pct = TRUE) +
  labs(y = "% missing", x = "Variables") +
  coord_flip()

### fig for missing data (followed order of vars and labeled)  
vars_to_plot <- c("HEI2020_ALL","AHEI_ALL","DASH_ALL","MED_ALL","DII_ALL",
                  "CHF","CHD","Angina","MI","Stroke",
                  "M_CVD","T2D","cancer","MCD","MetS")

# check vars in dataframe
vars_checked <- vars_to_plot[vars_to_plot %in% names(nhanes_merged_0518)]

# calculate n and % missing 
ms <- miss_var_summary(nhanes_merged_0518[vars_checked]) %>%  # -> variable, n_miss, pct_miss (0–100)
  mutate(variable = factor(variable, levels = rev(vars_checked)))   # keep order

max_pct <- max(ms$pct_miss, na.rm = TRUE)

# lollipop fig 
ggplot(ms, aes(y = variable, x = pct_miss)) +
  geom_segment(aes(x = 0, xend = pct_miss, yend = variable), linewidth = 0.7) +
  geom_point(size = 2.3) +
  geom_text(aes(label = sprintf("%.1f%% (n=%s)", pct_miss, scales::comma(n_miss))),
            hjust = -0.05, vjust = 0.5, size = 3.5) +
  scale_x_continuous(
    limits = c(0, max_pct * 1.10),               # leave space for label
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "% missing", y = NULL, title = " ") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# missing plot removed all vars with no missing data
nhanes_merged_0518 %>%
  miss_var_summary() %>%                 # -> variable, n_miss, pct_miss
  filter(pct_miss > 0) %>%               # filtered vars with no missing 
  ggplot(aes(x = pct_miss, y = reorder(variable, pct_miss))) +
  geom_segment(aes(x = 0, xend = pct_miss, yend = variable), linewidth = 0.6) +
  geom_point(size = 2) +
  labs(x = "% missing") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

##################################
# Exclusion criteria: adults >=20, no missing data in ALL dietary scores, not pregnant women
###################################
# filter 
excluded_cols <- c("HEI2020_ALL","AHEI_ALL","DASH_ALL", "MED_ALL", "DII_ALL")
                   #"M_CVD","T2D","cancer")

nhanes_merged_0518_filtered <- nhanes_merged_0518 %>%
  filter(
    RIDAGEYR >= 20,                                        # adults (>=20)
    if_all(all_of(excluded_cols), ~ !is.na(.))  # keep rows with NO missing in ALL excluded_cols
    #if_any(all_of(excluded_cols), ~ !is.na(.)) # (optional) keep rows with AT LEAST ONE non-missing in excluded_cols
    # !(RIAGENDR == 2 & RIDEXPRG == 1)          # exclude women pregnant 
  )

glimpse(nhanes_merged_0518_filtered)

gg_miss_var(nhanes_merged_0518_filtered, show_pct = TRUE) + #%
  labs(y = "% missing", x = "Variables") +
  coord_flip()

# missing plot removed all vars with no missing data
nhanes_merged_0518_filtered %>%
  miss_var_summary() %>%                 # -> variable, n_miss, pct_miss
  filter(pct_miss > 0) %>%               # filtered vars with no missing 
  ggplot(aes(x = pct_miss, y = reorder(variable, pct_miss))) +
  geom_segment(aes(x = 0, xend = pct_miss, yend = variable), linewidth = 0.6) +
  geom_point(size = 2) +
  labs(x = "% missing") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

################################################
# Frequency table (loading...)
################################################
library(dplyr)
library(purrr)
library(gt)

gt_cov_2col <- function(data, vars, digits = 1, title = NULL) {
  stopifnot(is.data.frame(data), is.character(vars))
  
  tab <- map_dfr(vars, function(v) {
    x <- data[[v]]
    if (is.null(x)) stop(sprintf("Variable '%s' not found in data.", v))
    n_total <- length(x)
    
    df <- tibble(level_chr = as.character(x), is_na = is.na(x)) %>%
      mutate(level_chr = ifelse(is_na, "NA", level_chr)) %>%
      count(level_chr, is_na, name = "n") %>%
      arrange(is_na, desc(n)) %>%
      mutate(
        variable = v,
        percent  = 100 * n / n_total
      ) %>%
      transmute(variable, level = level_chr, n, percent)
  })
  
  tab %>%
    gt(groupname_col = "variable") %>%
    cols_label(
      level   = "Level",
      n       = md("**n**"),
      percent = md("**%**")
    ) %>%
    fmt_number(columns = percent, decimals = digits) %>%
    tab_options(table.font.size = px(12), data_row.padding = px(4))
}

# specified vars
covars <- c("MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F",
            "RIAGENDR","RIDAGEYR_CAT","RIDRETH1","DMDEDUC2_CAT","DMDMARTL_CAT",
            "INDFMPIR_CAT","BMI_CAT","OBESE","AB_OBESE","SMQ_CAT","PAQ_CAT",
            "M_CVD","MetS","T2D","MCD")

covars <- c("CHF","CHD","Angina","MI","Stroke",
            "M_CVD","T2D","cancer", "MCD", "MetS")

# table for 
gt_cov_2col(nhanes_merged_0518, covars, digits = 1)

# table for filtered data
health_des <- gt_cov_2col(nhanes_merged_0518_filtered, covars, digits = 1)
health_des
