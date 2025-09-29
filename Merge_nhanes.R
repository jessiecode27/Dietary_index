####################################
# 1) Download & merge health variables
####################################
library(foreign)
library(dplyr)

years   <- c("2005","2007","2009","2011","2013","2015","2017")
suffix  <- c("_D","_E","_F","_G","_H","_I","_J")
cycles  <- c("0506","0708","0910","1112","1314","1516","1718")

# Helper: download XPT and keep SEQN + specified columns
read_keep <- function(url, keep) {
  tf <- tempfile(fileext = ".XPT")
  download.file(url, tf, mode = "wb", quiet = TRUE)
  df <- foreign::read.xport(tf)
  df[, intersect(names(df), c("SEQN", keep)), drop = FALSE]
}

# Groups of variables
demo_cols   <- c("RIAGENDR","RIDAGEYR","SDMVSTRA","SDMVPSU","WTMEC2YR",
                 "DMDEDUC2","RIDRETH1","DMDMARTL","INDFMPIR")
mcq_cols    <- c("MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F","MCQ220")
diq_cols    <- c("DIQ010","DIQ050", "DID060", "DIQ070", "DID040", "DIQ175S")
bpq_cols    <- c("BPQ020")
bpx_cols    <- c("BPXSY1","BPXSY2","BPXSY3","BPXSY4","BPXDI1","BPXDI2","BPXDI3","BPXDI4") # full cycles have four readings
bmx_cols    <- c("BMXBMI","BMXWAIST")
ghb_cols    <- c("LBXGH")
glu_cols    <- c("LBXGLU","LBXIN")     # LBXIN: Insulin (uU/mL) 
trigly_cols <- c("LBDLDL","LBXTR","LBDAPB")
hdl_cols    <- c("LBDHDD")

# covariate components
smq_cols    <- c("SMQ020","SMQ040")      # smoking 
alq_cols    <- c("ALQ130")               # alcohol drinks/day (12 mo)
paq_cols    <- c("PAQ670","PAD675")      # days/wk * minutes/day 
dsd_cols    <- c("DSD010")               # vitamin/mineral supplement use
rhq_cols    <- c("RHQ162")               # gestational diabetes

# Download and merge all health-related files for each cycle
per_cycle <- vector("list", length(cycles))
names(per_cycle) <- cycles

for (i in seq_along(years)) {
  yr  <- years[i]; suf <- suffix[i]
  
  # Short function to build NHANES file URL
  U <- function(stub) sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%s/DataFiles/%s%s.XPT", yr, stub, suf)
  
  # Read each component
  demo_i <- read_keep(U("DEMO"),   demo_cols)
  mcq_i  <- read_keep(U("MCQ"),    mcq_cols)
  diq_i  <- read_keep(U("DIQ"),    diq_cols)
  bpq_i  <- read_keep(U("BPQ"),    bpq_cols)
  bpx_i  <- read_keep(U("BPX"),    bpx_cols)
  bmx_i  <- read_keep(U("BMX"),    bmx_cols)
  ghb_i  <- read_keep(U("GHB"),    ghb_cols)
  glu_i  <- read_keep(U("GLU"),    glu_cols)
  trg_i  <- read_keep(U("TRIGLY"), trigly_cols)
  hdl_i  <- read_keep(U("HDL"),    hdl_cols)
  
  smq_i  <- read_keep(U("SMQ"),    smq_cols)
  alq_i  <- read_keep(U("ALQ"),    alq_cols)
  paq_i  <- read_keep(U("PAQ"),    paq_cols)
  rhq_i  <- read_keep(U("RHQ"),    rhq_cols)
  
  # chưa merge dsd_cols: check lại link path
  
  # Merge into one dataframe for this cycle
  per_cycle[[i]] <- Reduce(
    function(x, y) merge(x, y, by = "SEQN", all.x = TRUE, sort = FALSE),
    list(demo_i, mcq_i, diq_i, bpq_i, bpx_i, bmx_i, ghb_i, glu_i, trg_i, hdl_i, smq_i, alq_i, paq_i, rhq_i)
  )
}

# Bind all cycles and add "cycle" column
demo_health_0518 <- bind_rows(per_cycle, .id = "cycle")
glimpse(demo_health_0518)
View(demo_health_0518)
##############################################
# 2. Dietary index calculation
##############################################
### Read and merge FPED, NUTRIENT, DEMO ###
library(haven)   # read_xpt()/read_sas()

# Set base_path folder 
BASE <- "/Users/jessie/Desktop/Vyy/Projects 2025/Diet score paper - Dr. Linh/DATA"

# Define cycles and their NHANES letter suffixes 
cycles  <- c("0506","0708","0910","1112","1314","1516","1718")
letters <- c("D","E","F","G","H","I","J")

# Create an empty list to store all cycles
NHANES_ALL <- list()

# Loop: read files for each cycle and store as a nested list
for (i in seq_along(cycles)) {
  cyc <- cycles[i]   
  let <- letters[i]  
  
  # Build subfolders for this cycle
  fped_dir <- file.path(BASE, "FPED",     paste0("FPED_", cyc))
  nutr_dir <- file.path(BASE, "NUTRIENT", paste0("NUTRIENT_", cyc))
  demo_dir <- file.path(BASE, "DEMO",     paste0("DEMO_", cyc))
  
  # Read 9 expected tables for the cycle 
  cycle_list <- list(
    DEMO          = read_xpt(file.path(demo_dir, sprintf("DEMO_%s.xpt", let))),             
    
    FPED          = read_sas(file.path(fped_dir, sprintf("fped_dr1tot_%s.sas7bdat", cyc))), # FPED day 1 totals
    NUTRIENT      = read_xpt(file.path(nutr_dir, sprintf("DR1TOT_%s.xpt", let))),           # NUT day 1 totals
    
    FPED2         = read_sas(file.path(fped_dir, sprintf("fped_dr2tot_%s.sas7bdat", cyc))), # FPED day 2 totals
    NUTRIENT2     = read_xpt(file.path(nutr_dir, sprintf("DR2TOT_%s.xpt", let))),           # NUT day 2 totals
    
    NUTRIENT_IND  = read_xpt(file.path(nutr_dir, sprintf("DR1IFF_%s.xpt", let))),           # NUT day 1 individual foods
    FPED_IND      = read_sas(file.path(fped_dir, sprintf("fped_dr1iff_%s.sas7bdat", cyc))), # FPED day 1 individual foods
    
    FPED_IND2     = read_sas(file.path(fped_dir, sprintf("fped_dr2iff_%s.sas7bdat", cyc))), # FPED day 2 individual foods
    NUTRIENT_IND2 = read_xpt(file.path(nutr_dir, sprintf("DR2IFF_%s.xpt", let)))            # NUT day 2 individual foods
  )
  
  # Save this cycle into the outer list, named by its code (e.g., "1718")
  NHANES_ALL[[cyc]] <- cycle_list
}

### Use dietary package ###
library(purrr)
library(dplyr)
library(dietaryindex)

cycles <- c("0506","0708","0910","1112","1314","1516","1718")

# ---------- HEI2020 ----------
HEI2020_0518 <- map_dfr(cycles, function(cyc) {
  dat <- NHANES_ALL[[cyc]]
  out <- HEI2020_NHANES_FPED(
    FPED_PATH      = dat$FPED,
    NUTRIENT_PATH  = dat$NUTRIENT,
    DEMO_PATH      = dat$DEMO,
    FPED_PATH2     = dat$FPED2,
    NUTRIENT_PATH2 = dat$NUTRIENT2
  )
  mutate(out, cycle = cyc)
})

#View(HEI2020_0518)

# ---------- AHEI ----------
AHEI_0518 <- map_dfr(cycles, function(cyc) {
  dat <- NHANES_ALL[[cyc]]
  out <- AHEI_NHANES_FPED(
    FPED_IND_PATH      = dat$FPED_IND,
    NUTRIENT_IND_PATH  = dat$NUTRIENT_IND,
    FPED_IND_PATH2     = dat$FPED_IND2,
    NUTRIENT_IND_PATH2 = dat$NUTRIENT_IND2
  )
  mutate(out, cycle = cyc)
})

#View(AHEI_0518)

# ---------- DASH ----------
DASH_0518 <- map_dfr(cycles, function(cyc) {
  dat <- NHANES_ALL[[cyc]]
  out <- DASH_NHANES_FPED(
    FPED_IND_PATH      = dat$FPED_IND,
    NUTRIENT_IND_PATH  = dat$NUTRIENT_IND,
    FPED_IND_PATH2     = dat$FPED_IND2,
    NUTRIENT_IND_PATH2 = dat$NUTRIENT_IND2
  )
  mutate(out, cycle = cyc)
})

#View(DASH_0518)

# ---------- aMED ----------
AMED_0518 <- map_dfr(cycles, function(cyc) {
  dat <- NHANES_ALL[[cyc]]
  out <- MED_NHANES_FPED(
    FPED_PATH      = dat$FPED,
    NUTRIENT_PATH  = dat$NUTRIENT,
    DEMO_PATH      = dat$DEMO,
    FPED_PATH2     = dat$FPED2,
    NUTRIENT_PATH2 = dat$NUTRIENT2
  )
  mutate(out, cycle = cyc)
})
#View(AMED_0518)

# ---------- DII ----------
DII_0518 <- map_dfr(cycles, function(cyc) {
  dat <- NHANES_ALL[[cyc]]
  out <- DII_NHANES_FPED(
    FPED_PATH      = dat$FPED,
    NUTRIENT_PATH  = dat$NUTRIENT,
    DEMO_PATH      = dat$DEMO,
    FPED_PATH2     = dat$FPED2,
    NUTRIENT_PATH2 = dat$NUTRIENT2
  )
  mutate(out, cycle = cyc)
})
#View(DII_0518)

##############################################
# 3. Mortality data
##############################################
library(readr)
library(dplyr)

# Map cycle 
cyc_years <- c(
  "0506" = "2005_2006",
  "0708" = "2007_2008",
  "0910" = "2009_2010",
  "1112" = "2011_2012",
  "1314" = "2013_2014",
  "1516" = "2015_2016",
  "1718" = "2017_2018"
)

# Function to read mortality
read_mort <- function(dir, cyc) {
  yrs  <- cyc_years[[cyc]]
  path <- file.path(dir, paste0("NHANES_", yrs, "_MORT_2019_PUBLIC.dat"))
  out <- read_fwf(
    file      = path,
    col_types = "iiiiiiii",
    col_positions = fwf_cols(
      SEQN         = c(1, 6),
      eligstat     = c(15, 15),
      mortstat     = c(16, 16),
      ucod_leading = c(17, 19),
      diabetes     = c(20, 20),
      hyperten     = c(21, 21),
      permth_int   = c(43, 45),
      permth_exm   = c(46, 48)
    ),
    na = c("", ".")
  )
  mutate(out, cycle = cyc)
}

# Read and merge mortality data for 2005-2018
MORT_DIR <- "/Users/jessie/Desktop/Vyy/Projects 2025/Diet score paper - Dr. Linh/DATA/Mortality"
cycles   <- names(cyc_years)

mort_0518 <- do.call(
  bind_rows,
  lapply(cycles, function(cyc) read_mort(MORT_DIR, cyc))
)
View(mort_0518)

##############################################
# 4) Merge: HEALTH -> (HEI, AHEI, DASH, aMED, DII) -> MORT
#    Drop duplicate columns from RHS (keep SEQN, cycle)
##############################################
by_keys <- c("SEQN","cycle")

left_join_drop_dups <- function(x, y, by) {
  cols_y <- union(by, setdiff(names(y), names(x)))
  y_sel  <- dplyr::select(y, dplyr::all_of(cols_y))
  dplyr::left_join(x, y_sel, by = by)
}

nhanes_merged_0518 <- Reduce(
  function(x, y) left_join_drop_dups(x, y, by = by_keys),
  list(
    demo_health_0518,
    HEI2020_0518,
    AHEI_0518,
    DASH_0518,
    AMED_0518,
    DII_0518,
    mort_0518
  )
)

glimpse(nhanes_merged_0518)
View(nhanes_merged_0518)

##############################################
# 5) Data Cleaning for Covariates 
##############################################
# Factor RIAGENDR
nhanes_merged_0518$RIAGENDR = factor(nhanes_merged_0518$RIAGENDR, levels = c(1,2), labels = c("Male", "Female"))
# Factor RIDAGEYR_CAT (Age Group): 20-34 (1), 35-49 (2), 50-64 (3), and 65 or older (4)
nhanes_merged_0518$RIDAGEYR_CAT = factor(nhanes_merged_0518$RIDAGEYR_CAT, levels = c(1,2,3,4), labels = c("20-34", "35-49", "50-64", "> 65"))
# Factor RIDRETH1 (Race/Ethicity)
nhanes_merged_0518$RIDRETH1 = factor(nhanes_merged_0518$RIDRETH1, 
                                           levels = c(1, 2, 3, 4, 5), 
                                           labels = c("Mexican American", "Other Hispanic", "Non-Hispanic White", "Non-Hispanic Black", "Other Race")
)
# Factor DMDEDUC2_CAT (Education Level)
nhanes_merged_0518$DMDEDUC2_CAT = factor(
  nhanes_merged_0518$DMDEDUC2_CAT, 
  levels = c(2, 3, 4, 5), 
  labels = c("Less Than High School", "High School Graduate/GED", "Some College or Associate's Degree", "College Graduate or Above")
)
# Factor DMDMARTL_CAT (Marital Status)
nhanes_merged_0518$DMDMARTL_CAT = factor(
  nhanes_merged_0518$DMDMARTL_CAT, 
  levels = c(1, 0), 
  labels = c("Married", "Not Married")
)
# Factor INDFMPIR_CAT (Income to Poverty Ratio)
nhanes_merged_0518$INDFMPIR_CAT = factor(
  nhanes_merged_0518$INDFMPIR_CAT, 
  levels = c(1, 2, 3), 
  labels = c("0-1.3", "1.3-3.5", ">3.5")
)
# Factor BMI_CAT (BMI Category): underweight or normal weight (2), overweight (3), and obese (4)
nhanes_merged_0518$BMI_CAT = factor(
  nhanes_merged_0518$BMI_CAT, 
  levels = c(2, 3, 4), 
  labels = c("Underweight or Normal Weight", "Overweight", "Obese")
)


# update the survey design for combined 14 years NHANES data (2005-2018, combined 14 years)
nhanes_merged_0518_design_14yr = svydesign(
  id = ~SDMVPSU, 
  strata = ~SDMVSTRA, 
  weight = ~WTDR14D, 
  data = nhanes_merged_0518, #set up survey design on the full dataset #can restrict at time of analysis 
  nest = TRUE)

# Check missing data and distribution
library(DataExplorer)

create_report(nhanes_merged_0518,
              output_file = "nhanes_profile_report.html",
              output_dir = getwd())


create_report(PHDI_GHG_FINAL_DATA_0518,
              output_file = "nhanes_james_profile_report.html",
              output_dir = getwd())
