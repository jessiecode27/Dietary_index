# clean memory
rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(gtsummary)
# include COMPLEX SURVEY DESIGN #######################################
library(survey)
options(survey.lonely.psu = "adjust") #accounts for the lonely psu problem from subsetting survey data to small groups

# Load the final merged NHANES dataset
PHDI_GHG_FINAL_DATA_0518 = read_csv("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PHDI_NHANES/NHANES dataset/GHG_PHDI_disease/PHDI_GHG_FINAL_DATA_0518_AHEI_HEI2020_DASH_AMED_DII.csv")

# create the survey design for individual NHANES data (2005-2018)
PHDI_GHG_FINAL_DATA_0518_design = svydesign(
    id = ~SDMVPSU, 
    strata = ~SDMVSTRA, 
    weight = ~WTDR2D, 
    data = PHDI_GHG_FINAL_DATA_0518, #set up survey design on the full dataset #can restrict at time of analysis 
    nest = TRUE)

# update the survey design for combined 14 years NHANES data (2005-2018, combined 14 years)
PHDI_GHG_FINAL_DATA_0518_design_14yr = svydesign(
    id = ~SDMVPSU, 
    strata = ~SDMVSTRA, 
    weight = ~WTDR14D, 
    data = PHDI_GHG_FINAL_DATA_0518, #set up survey design on the full dataset #can restrict at time of analysis 
    nest = TRUE)

# create a function to adjust a variable for total energy intake
adjust_for_energy <- function(var){
  # create a formula for the model
  formula <- as.formula(paste(var, "~ TOTALKCAL_PHDI"))
  # fit the model
  model <- svyglm(formula, PHDI_GHG_FINAL_DATA_0518_design)
  # predict the variable
  predicted <- predict(model, type = "response")
  # return the predicted variable
  return(predicted)
}

# adjust all component scores for total energy intake
PHDI_GHG_FINAL_DATA_0518$PHDI_WGRAIN_kcal_adjusted <- adjust_for_energy("PHDI_WGRAIN")
PHDI_GHG_FINAL_DATA_0518$PHDI_STARCHY_VEG_kcal_adjusted <- adjust_for_energy("PHDI_STARCHY_VEG")
PHDI_GHG_FINAL_DATA_0518$PHDI_VEG_kcal_adjusted <- adjust_for_energy("PHDI_VEG")
PHDI_GHG_FINAL_DATA_0518$PHDI_FRT_kcal_adjusted <- adjust_for_energy("PHDI_FRT")
PHDI_GHG_FINAL_DATA_0518$PHDI_DAIRY_kcal_adjusted <- adjust_for_energy("PHDI_DAIRY")
PHDI_GHG_FINAL_DATA_0518$PHDI_REDPROC_MEAT_kcal_adjusted <- adjust_for_energy("PHDI_REDPROC_MEAT")
PHDI_GHG_FINAL_DATA_0518$PHDI_POULTRY_kcal_adjusted <- adjust_for_energy("PHDI_POULTRY")
PHDI_GHG_FINAL_DATA_0518$PHDI_EGG_kcal_adjusted <- adjust_for_energy("PHDI_EGG")
PHDI_GHG_FINAL_DATA_0518$PHDI_FISH_kcal_adjusted <- adjust_for_energy("PHDI_FISH")
PHDI_GHG_FINAL_DATA_0518$PHDI_NUTS_kcal_adjusted <- adjust_for_energy("PHDI_NUTS")
PHDI_GHG_FINAL_DATA_0518$PHDI_LEGUMES_kcal_adjusted <- adjust_for_energy("PHDI_LEGUMES")
PHDI_GHG_FINAL_DATA_0518$PHDI_SOY_kcal_adjusted <- adjust_for_energy("PHDI_SOY")

# there is no need to adjust the added fat (unsaturated and saturated) and added sugar scores for total energy intake because they are already in the % of total energy format

# calculate the PHDI_ALL scores adjusted by kcal by adding all the component scores
PHDI_GHG_FINAL_DATA_0518$PHDI_ALL_kcal_adjusted = PHDI_GHG_FINAL_DATA_0518$PHDI_WGRAIN_kcal_adjusted + 
    PHDI_GHG_FINAL_DATA_0518$PHDI_STARCHY_VEG_kcal_adjusted + 
    PHDI_GHG_FINAL_DATA_0518$PHDI_VEG_kcal_adjusted + 
    PHDI_GHG_FINAL_DATA_0518$PHDI_FRT_kcal_adjusted + 
    PHDI_GHG_FINAL_DATA_0518$PHDI_DAIRY_kcal_adjusted + 
    PHDI_GHG_FINAL_DATA_0518$PHDI_REDPROC_MEAT_kcal_adjusted + 
    PHDI_GHG_FINAL_DATA_0518$PHDI_POULTRY_kcal_adjusted + 
    PHDI_GHG_FINAL_DATA_0518$PHDI_EGG_kcal_adjusted + 
    PHDI_GHG_FINAL_DATA_0518$PHDI_FISH_kcal_adjusted + 
    PHDI_GHG_FINAL_DATA_0518$PHDI_NUTS_kcal_adjusted + 
    PHDI_GHG_FINAL_DATA_0518$PHDI_LEGUMES_kcal_adjusted + 
    PHDI_GHG_FINAL_DATA_0518$PHDI_SOY_kcal_adjusted + 
    PHDI_GHG_FINAL_DATA_0518$PHDI_ADDED_FAT_UNSAT + 
    PHDI_GHG_FINAL_DATA_0518$PHDI_ADDED_FAT_SAT + 
    PHDI_GHG_FINAL_DATA_0518$PHDI_ADDED_SUGAR

summary(PHDI_GHG_FINAL_DATA_0518$PHDI_ALL)
summary(PHDI_GHG_FINAL_DATA_0518$PHDI_ALL_kcal_adjusted)

# update the survey design for PHDI day 1 with the predicted PHDI_ALL scores
PHDI_GHG_FINAL_DATA_0518_design = svydesign(
    id = ~SDMVPSU, 
    strata = ~SDMVSTRA, 
    weight = ~WTDR2D, 
    data = PHDI_GHG_FINAL_DATA_0518, #set up survey design on the full dataset #can restrict at time of analysis 
    nest = TRUE)

setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PHDI_NHANES/Publication/Results/Main tables")
###################################################################
# Define a function to call svymean and unweighted count
getSummary <- function(varformula, byformula, design){
  # Get mean, stderr, and unweighted sample size
  ## unweighted count
  c <- svyby(varformula, byformula, design, unwtd.count) 
  ## weighted mean by complex survey design
  p <- svyby(varformula, byformula, design, svymean) 
  ## join the two data frames by removing the standard error column in the unweighted count data frame
  outSum <- left_join(select(c, -se), p) 
  ## export results
  outSum
}

# Figure 1 ################################################################
# Calculate weighted mean and se of PHDI_ALL (energy-adjusted) by YEAR status
## By YEAR
PHDI_ALL_survery_summary_kcal_adjusted = getSummary(~PHDI_ALL_kcal_adjusted, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    ## rename the columns
    rename(SE_PHDI = se) %>%
    ## add a column called "group" so that we can use it to plot lines
    mutate(group = "PHDI_ALL_kcal_adjusted")
# Plot raw PHDI_ALL data with error bars
figure_1 = ggplot(PHDI_ALL_survery_summary_kcal_adjusted, aes(x=YEAR, y=PHDI_ALL_kcal_adjusted)) +
  geom_line(aes(group = group), linewidth = 1) +
  geom_point(size = 3, color = "#F8766D") +
  geom_errorbar(aes(ymin = PHDI_ALL_kcal_adjusted - SE_PHDI, ymax = PHDI_ALL_kcal_adjusted + SE_PHDI), width=0.2) +
  theme_bw() +
  labs(x = "NHANES Cycle", y = "Mean PHDI (energy-adjusted)") +
  theme(plot.title = element_text(size=22),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14, angle = 90, hjust = 1),
        axis.text.y = element_text(size=14),
        # legend.text = element_text(size=16),
        # legend.title = element_text(size=18),
        aspect.ratio = 0.5)

# save the figure 1 data
write.csv(PHDI_ALL_survery_summary_kcal_adjusted, "PHDI_ALL_survery_summary_kcal_adjusted.csv")

# Figure S2 ################################################################
# The distribution of PHDI_ALL in the US population
figure_s2 = svyhist(~PHDI_ALL, PHDI_GHG_FINAL_DATA_0518_design_14yr,
        breaks = "Sturges", # or set your own number of breaks
        col = "#619CFF", 
        probability=FALSE,
        xlab = "Planetary Health Diet Index",
        main = "Planetary Health Diet Index Distribution in the US Population",
        xlim = c(0, 140))
    
# print the quantile of PHDI_ALL
print("This is the distribution of PHDI_ALL")
svyquantile(~PHDI_ALL, PHDI_GHG_FINAL_DATA_0518_design_14yr, c(0.1, 0.9))[[1]]

# Figure S3 ################################################################
# The distribution of GHG_TOTAL in the US population
figure_s3 = svyhist(~GHG_TOTAL, PHDI_GHG_FINAL_DATA_0518_design_14yr,
        breaks = "Sturges", # or set your own number of breaks
        col = "#00BA38", 
        border = "#000000",
        probability=FALSE,
        xlab = "Greenhouse Gas Emission (kg CO2) per day",
        main = "Greenhouse Gas Emission Distribution in the US Population",
        xlim = c(0, 40))

# print the quantile of GHG_TOTAL
print("This is the distribution of GHG_TOTAL")
svyquantile(~GHG_TOTAL, PHDI_GHG_FINAL_DATA_0518_design_14yr, c(0.1, 0.9))[[1]]


# Figure S4 ################################################################
# draw the PHDI_ALL trends (raw values) in the US population

# Calculate weighted mean and se of PHDI_ALL by YEAR status
## By YEAR
PHDI_ALL_survery_summary_raw = getSummary(~PHDI_ALL, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    ## rename the columns
    rename(SE_PHDI = se) %>%
    ## add a column called "group" so that we can use it to plot lines
    mutate(group = "PHDI_ALL")
# Plot PHDI_ALL data with error bars 
figure_s4 = ggplot(PHDI_ALL_survery_summary_raw, aes(x=YEAR, y=PHDI_ALL)) +
  geom_line(aes(group=group),linewidth = 1) +
  geom_point(size = 3, color = "#F8766D") +
  geom_errorbar(aes(ymin = PHDI_ALL - SE_PHDI, ymax = PHDI_ALL + SE_PHDI), width=0.2) +
  theme_bw() +
  labs(x = "NHANES Cycle", y = "Mean PHDI") +
  theme(plot.title = element_text(size=22),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14, angle = 90, hjust = 1),
        axis.text.y = element_text(size=14),
        # legend.text = element_text(size=16),
        # legend.title = element_text(size=18),
        aspect.ratio = 0.5) 

# Figure 2 ################################################################
# draw the PHDI dietary component (adjusted by energy) trends in the US population

# draw the PHDI dietary component (adjusted by energy) trends
#' ### Calculate weighted mean and se
#' By YEAR, PHDI_WGRAIN_kcal_adjusted, PHDI_STARCHY_VEG_kcal_adjusted, PHDI_VEG_kcal_adjusted, PHDI_FRT_kcal_adjusted, PHDI_DAIRY_kcal_adjusted, PHDI_REDPROC_MEAT_kcal_adjusted, PHDI_POULTRY_kcal_adjusted, PHDI_EGG_kcal_adjusted, PHDI_FISH_kcal_adjusted, PHDI_NUTS_kcal_adjusted, PHDI_LEGUMES_kcal_adjusted, PHDI_SOY_kcal_adjusted, PHDI_ADDED_FAT_UNSAT_kcal_adjusted, PHDI_ADDED_FAT_SAT_kcal_adjusted, PHDI_ADDED_SUGAR_kcal_adjusted
PHDI_ALL_survery_summary_PHDI_PHDI_WGRAIN_kcal_adjusted = getSummary(~PHDI_WGRAIN_kcal_adjusted, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_WGRAIN_kcal_adjusted, SE_PHDI_component = se) %>%
    ## add a column called "Whole grains"
    mutate(component = "Whole grains")
PHDI_ALL_survery_summary_PHDI_STARTCHY_VEG = getSummary(~PHDI_STARCHY_VEG_kcal_adjusted, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_STARCHY_VEG_kcal_adjusted, SE_PHDI_component = se) %>%
    ## add a column called "Starchy vegetables"
    mutate(component = "Starchy vegetables")
PHDI_ALL_survery_summary_PHDI_PHDI_VEG_kcal_adjusted = getSummary(~PHDI_VEG_kcal_adjusted, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_VEG_kcal_adjusted, SE_PHDI_component = se) %>%
    ## add a column called "Total vegetables"
    mutate(component = "Non-starchy vegetables")
PHDI_ALL_survery_summary_PHDI_PHDI_FRT_kcal_adjusted = getSummary(~PHDI_FRT_kcal_adjusted, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_FRT_kcal_adjusted, SE_PHDI_component = se) %>%
    ## add a column called "Whole fruit"
    mutate(component = "Whole fruit")
PHDI_ALL_survery_summary_PHDI_PHDI_DAIRY_kcal_adjusted = getSummary(~PHDI_DAIRY_kcal_adjusted, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_DAIRY_kcal_adjusted, SE_PHDI_component = se) %>%
    ## add a column called "Dairy"
    mutate(component = "Dairy")
PHDI_ALL_survery_summary_PHDI_PHDI_REDPROC_MEAT_kcal_adjusted = getSummary(~PHDI_REDPROC_MEAT_kcal_adjusted, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_REDPROC_MEAT_kcal_adjusted, SE_PHDI_component = se) %>%
    ## add a column called "Red and processed meat"
    mutate(component = "Red and processed meat")
PHDI_ALL_survery_summary_PHDI_PHDI_POULTRY_kcal_adjusted = getSummary(~PHDI_POULTRY_kcal_adjusted, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_POULTRY_kcal_adjusted, SE_PHDI_component = se) %>%
    ## add a column called "Poultry"
    mutate(component = "Poultry")
PHDI_ALL_survery_summary_PHDI_PHDI_EGG_kcal_adjusted = getSummary(~PHDI_EGG_kcal_adjusted, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_EGG_kcal_adjusted, SE_PHDI_component = se) %>%
    ## add a column called "Eggs"
    mutate(component = "Eggs")
PHDI_ALL_survery_summary_PHDI_PHDI_FISH_kcal_adjusted = getSummary(~PHDI_FISH_kcal_adjusted, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_FISH_kcal_adjusted, SE_PHDI_component = se) %>%
    ## add a column called "Fish"
    mutate(component = "Fish")
PHDI_ALL_survery_summary_PHDI_PHDI_NUTS_kcal_adjusted = getSummary(~PHDI_NUTS_kcal_adjusted, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_NUTS_kcal_adjusted, SE_PHDI_component = se) %>%
    ## add a column called "Nuts and seeds"
    mutate(component = "Nuts and seeds")
PHDI_ALL_survery_summary_PHDI_PHDI_LEGUMES_kcal_adjusted = getSummary(~PHDI_LEGUMES_kcal_adjusted, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_LEGUMES_kcal_adjusted, SE_PHDI_component = se) %>%
    ## add a column called "Nonsoy legumes"
    mutate(component = "Nonsoy legumes")
PHDI_ALL_survery_summary_PHDI_PHDI_SOY_kcal_adjusted = getSummary(~PHDI_SOY_kcal_adjusted, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_SOY_kcal_adjusted, SE_PHDI_component = se) %>%
    ## add a column called "Soy products"
    mutate(component = "Soy products")
PHDI_ALL_survery_summary_PHDI_PHDI_ADDED_FAT_UNSAT_kcal_adjusted = getSummary(~PHDI_ADDED_FAT_UNSAT, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_ADDED_FAT_UNSAT, SE_PHDI_component = se) %>%
    ## add a column called "Added fats (unsaturated)"
    mutate(component = "Added fats (unsaturated)")
PHDI_ALL_survery_summary_PHDI_PHDI_ADDED_FAT_SAT_kcal_adjusted = getSummary(~PHDI_ADDED_FAT_SAT, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_ADDED_FAT_SAT, SE_PHDI_component = se) %>%
    ## add a column called "Added fats (saturated)"
    mutate(component = "Added fats (saturated)")
PHDI_ALL_survery_summary_PHDI_PHDI_ADDED_SUGAR_kcal_adjusted = getSummary(~PHDI_ADDED_SUGAR, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_ADDED_SUGAR, SE_PHDI_component = se) %>%
    ## add a column called "Added sugars"
    mutate(component = "Added sugars")

# combine all data by rows
PHDI_component_survery_summary_kcal_adjusted = rbind(
    PHDI_ALL_survery_summary_PHDI_PHDI_WGRAIN_kcal_adjusted,
    PHDI_ALL_survery_summary_PHDI_STARTCHY_VEG,
    PHDI_ALL_survery_summary_PHDI_PHDI_VEG_kcal_adjusted,
    PHDI_ALL_survery_summary_PHDI_PHDI_FRT_kcal_adjusted,
    PHDI_ALL_survery_summary_PHDI_PHDI_DAIRY_kcal_adjusted,
    PHDI_ALL_survery_summary_PHDI_PHDI_REDPROC_MEAT_kcal_adjusted,
    PHDI_ALL_survery_summary_PHDI_PHDI_POULTRY_kcal_adjusted,
    PHDI_ALL_survery_summary_PHDI_PHDI_EGG_kcal_adjusted,
    PHDI_ALL_survery_summary_PHDI_PHDI_FISH_kcal_adjusted,
    PHDI_ALL_survery_summary_PHDI_PHDI_NUTS_kcal_adjusted,
    PHDI_ALL_survery_summary_PHDI_PHDI_LEGUMES_kcal_adjusted,
    PHDI_ALL_survery_summary_PHDI_PHDI_SOY_kcal_adjusted,
    PHDI_ALL_survery_summary_PHDI_PHDI_ADDED_FAT_UNSAT_kcal_adjusted,
    PHDI_ALL_survery_summary_PHDI_PHDI_ADDED_FAT_SAT_kcal_adjusted,
    PHDI_ALL_survery_summary_PHDI_PHDI_ADDED_SUGAR_kcal_adjusted
    )

PHDI_component_survery_summary_kcal_adjusted$component <- factor(PHDI_component_survery_summary_kcal_adjusted$component)

# Plot weighted mean with error bars for each component
figure_2 = ggplot(PHDI_component_survery_summary_kcal_adjusted, aes(x=YEAR, y=Mean_PHDI_component, color=component)) +
    geom_line(aes(group=component), linewidth = 1) +
    ggrepel::geom_label_repel(
        aes(label = component),
        data = subset(PHDI_component_survery_summary_kcal_adjusted, YEAR == "2005-2006"),
        size = 3,
        box.padding = 0.5,
        point.padding = 0.3
        ) +
    theme_bw() +
    scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0, 10)) + 
    labs(x = "NHANES Cycle", y = "Mean PHDI component (energy-adjusted)") +
    theme(plot.title = element_text(size=22),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14, angle = 90, hjust = 1),
        axis.text.y = element_text(size=14),
        # legend.text = element_text(size=16),
        # legend.title = element_text(size=18),
        # aspect.ratio = 1.5
        ) +
    ## hide the legend
    theme(legend.position = "none")

# save the figure 2 data
write.csv(PHDI_component_survery_summary_kcal_adjusted, "PHDI_component_survery_summary_kcal_adjusted.csv")

# Figure s5 ################################################################
# draw the PHDI dietary component trends (raw values)
#' ### Calculate weighted mean and se of PHDI_ALL by YEAR
#' By YEAR, PHDI_WGRAIN, PHDI_STARCHY_VEG, PHDI_VEG, PHDI_FRT, PHDI_DAIRY, PHDI_REDPROC_MEAT, PHDI_POULTRY, PHDI_EGG, PHDI_FISH, PHDI_NUTS, PHDI_LEGUMES, PHDI_SOY, PHDI_ADDED_FAT_UNSAT, PHDI_ADDED_FAT_SAT, PHDI_ADDED_SUGAR
PHDI_ALL_survery_summary_PHDI_PHDI_WGRAIN = getSummary(~PHDI_WGRAIN, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_WGRAIN, SE_PHDI_component = se) %>%
    ## add a column called "Whole grains"
    mutate(component = "Whole grains")
PHDI_ALL_survery_summary_PHDI_STARTCHY_VEG = getSummary(~PHDI_STARCHY_VEG, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_STARCHY_VEG, SE_PHDI_component = se) %>%
    ## add a column called "Starchy vegetables"
    mutate(component = "Starchy vegetables")
PHDI_ALL_survery_summary_PHDI_PHDI_VEG = getSummary(~PHDI_VEG, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_VEG, SE_PHDI_component = se) %>%
    ## add a column called "Total vegetables"
    mutate(component = "Non-starchy vegetables")
PHDI_ALL_survery_summary_PHDI_PHDI_FRT = getSummary(~PHDI_FRT, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_FRT, SE_PHDI_component = se) %>%
    ## add a column called "Whole fruit"
    mutate(component = "Whole fruit")
PHDI_ALL_survery_summary_PHDI_PHDI_DAIRY = getSummary(~PHDI_DAIRY, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_DAIRY, SE_PHDI_component = se) %>%
    ## add a column called "Dairy"
    mutate(component = "Dairy")
PHDI_ALL_survery_summary_PHDI_PHDI_REDPROC_MEAT = getSummary(~PHDI_REDPROC_MEAT, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_REDPROC_MEAT, SE_PHDI_component = se) %>%
    ## add a column called "Red and processed meat"
    mutate(component = "Red and processed meat")
PHDI_ALL_survery_summary_PHDI_PHDI_POULTRY = getSummary(~PHDI_POULTRY, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_POULTRY, SE_PHDI_component = se) %>%
    ## add a column called "Poultry"
    mutate(component = "Poultry")
PHDI_ALL_survery_summary_PHDI_PHDI_EGG = getSummary(~PHDI_EGG, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_EGG, SE_PHDI_component = se) %>%
    ## add a column called "Eggs"
    mutate(component = "Eggs")
PHDI_ALL_survery_summary_PHDI_PHDI_FISH = getSummary(~PHDI_FISH, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_FISH, SE_PHDI_component = se) %>%
    ## add a column called "Fish"
    mutate(component = "Fish")
PHDI_ALL_survery_summary_PHDI_PHDI_NUTS = getSummary(~PHDI_NUTS, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_NUTS, SE_PHDI_component = se) %>%
    ## add a column called "Nuts and seeds"
    mutate(component = "Nuts and seeds")
PHDI_ALL_survery_summary_PHDI_PHDI_LEGUMES = getSummary(~PHDI_LEGUMES, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_LEGUMES, SE_PHDI_component = se) %>%
    ## add a column called "Nonsoy legumes"
    mutate(component = "Nonsoy legumes")
PHDI_ALL_survery_summary_PHDI_PHDI_SOY = getSummary(~PHDI_SOY, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_SOY, SE_PHDI_component = se) %>%
    ## add a column called "Soy products"
    mutate(component = "Soy products")
PHDI_ALL_survery_summary_PHDI_PHDI_ADDED_FAT_UNSAT = getSummary(~PHDI_ADDED_FAT_UNSAT, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_ADDED_FAT_UNSAT, SE_PHDI_component = se) %>%
    ## add a column called "Added fats (unsaturated)"
    mutate(component = "Added fats (unsaturated)")
PHDI_ALL_survery_summary_PHDI_PHDI_ADDED_FAT_SAT = getSummary(~PHDI_ADDED_FAT_SAT, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_ADDED_FAT_SAT, SE_PHDI_component = se) %>%
    ## add a column called "Added fats (saturated)"
    mutate(component = "Added fats (saturated)")
PHDI_ALL_survery_summary_PHDI_PHDI_ADDED_SUGAR = getSummary(~PHDI_ADDED_SUGAR, ~YEAR, PHDI_GHG_FINAL_DATA_0518_design) %>%
    rename(Mean_PHDI_component = PHDI_ADDED_SUGAR, SE_PHDI_component = se) %>%
    ## add a column called "Added sugars"
    mutate(component = "Added sugars")


# combine all data by rows
PHDI_component_survery_summary_raw = rbind(
    PHDI_ALL_survery_summary_PHDI_PHDI_WGRAIN,
    PHDI_ALL_survery_summary_PHDI_STARTCHY_VEG,
    PHDI_ALL_survery_summary_PHDI_PHDI_VEG,
    PHDI_ALL_survery_summary_PHDI_PHDI_FRT,
    PHDI_ALL_survery_summary_PHDI_PHDI_DAIRY,
    PHDI_ALL_survery_summary_PHDI_PHDI_REDPROC_MEAT,
    PHDI_ALL_survery_summary_PHDI_PHDI_POULTRY,
    PHDI_ALL_survery_summary_PHDI_PHDI_EGG,
    PHDI_ALL_survery_summary_PHDI_PHDI_FISH,
    PHDI_ALL_survery_summary_PHDI_PHDI_NUTS,
    PHDI_ALL_survery_summary_PHDI_PHDI_LEGUMES,
    PHDI_ALL_survery_summary_PHDI_PHDI_SOY,
    PHDI_ALL_survery_summary_PHDI_PHDI_ADDED_FAT_UNSAT,
    PHDI_ALL_survery_summary_PHDI_PHDI_ADDED_FAT_SAT,
    PHDI_ALL_survery_summary_PHDI_PHDI_ADDED_SUGAR
)

PHDI_component_survery_summary_raw$component <- factor(PHDI_component_survery_summary_raw$component)

# Plot weighted mean for each PHDI component
figure_s5 = ggplot(PHDI_component_survery_summary_raw, aes(x=YEAR, y=Mean_PHDI_component, color=component)) +
    geom_line(aes(group=component), linewidth = 1) +
    ggrepel::geom_label_repel(
        aes(label = component),
        data = subset(PHDI_component_survery_summary_raw, YEAR == "2005-2006"),
        size = 3,
        box.padding = 0.5,
        point.padding = 0.3
        ) +
    theme_bw() +
    scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10), limits=c(0, 10)) + 
    labs(x = "NHANES Cycle", y = "Mean PHDI component") +
    theme(plot.title = element_text(size=22),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14, angle = 90, hjust = 1),
        axis.text.y = element_text(size=14),
        # legend.text = element_text(size=16),
        # legend.title = element_text(size=18),
        # aspect.ratio = 1.5
        ) +
    ## hide the legend
    theme(legend.position = "none")

library(jtools)

# update the survey design for combined 14 years NHANES data (2005-2018, combined 14 years)
PHDI_GHG_FINAL_DATA_0518_design_14yr = svydesign(
    id = ~SDMVPSU, 
    strata = ~SDMVSTRA, 
    weight = ~WTDR14D, 
    data = PHDI_GHG_FINAL_DATA_0518, #set up survey design on the full dataset #can restrict at time of analysis 
    nest = TRUE)

# Figure 3: Association between PHDI and GHG using survey design
figure_3 = svyplot(
    GHG_TOTAL~PHDI_ALL, 
    PHDI_GHG_FINAL_DATA_0518_design_14yr, 
    col = "#61faff",
    xlab = "Planetary Health Diet Index",
    ylab = "Greenhouse Gas Emission (kg CO2) per day",
    main = "Association between PHDI and GHG"
)
# add a regression line
abline(svyglm(GHG_TOTAL~PHDI_ALL, design = PHDI_GHG_FINAL_DATA_0518_design_14yr), col = "#00BA38", lwd = 2)
summary(svyglm(GHG_TOTAL~PHDI_ALL, design = PHDI_GHG_FINAL_DATA_0518_design_14yr))

# Given values
estimate = -0.03076
std_error = 0.00252

# Calculate 10 times the beta and its standard error
beta_10 <- 10 * estimate
se_10 <- 10 * std_error

# Calculate the 95% CI
z <- 1.96
ci_lower <- beta_10 - z * se_10
ci_upper <- beta_10 + z * se_10

# Print results
cat("10 times the beta coefficient:", beta_10, "\n")
cat("95% Confidence Interval:", "[", ci_lower, ",", ci_upper, "]\n")

# Figure s6: PHDI association with GHG using survey design (kcal adjusted)
figure_s6 = svyplot(
    GHG_TOTAL~PHDI_ALL_kcal_adjusted, 
    PHDI_GHG_FINAL_DATA_0518_design_14yr, 
    col = "#61faff",
    xlab = "Planetary Health Diet Index",
    ylab = "Greenhouse Gas Emission (kg CO2) per day",
    main = "Association between energy-adjusted PHDI and GHG"
)
# add a regression line
abline(svyglm(GHG_TOTAL~PHDI_ALL_kcal_adjusted, design = PHDI_GHG_FINAL_DATA_0518_design_14yr), col = "#00BA38", lwd = 2)
summary(svyglm(GHG_TOTAL~PHDI_ALL_kcal_adjusted, design = PHDI_GHG_FINAL_DATA_0518_design_14yr))

# Given values
estimate <- -0.05851
std_error <- 0.00713

# Calculate 10 times the beta and its standard error
beta_10 <- 10 * estimate
se_10 <- 10 * std_error

# Calculate the 95% CI
z <- 1.96
ci_lower <- beta_10 - z * se_10
ci_upper <- beta_10 + z * se_10

# Print results
cat("10 times the beta coefficient:", beta_10, "\n")
cat("95% Confidence Interval:", "[", ci_lower, ",", ci_upper, "]\n")


# data cleaning for covariates ###########################################
# convert YEAR to factor
PHDI_GHG_FINAL_DATA_0518$YEAR <- factor(PHDI_GHG_FINAL_DATA_0518$YEAR)

# Create a quintile variable for PHDI_ALL
PHDI_GHG_FINAL_DATA_0518$PHDI_ALL_quintile <- cut(PHDI_GHG_FINAL_DATA_0518$PHDI_ALL, breaks = quantile(PHDI_GHG_FINAL_DATA_0518$PHDI_ALL, probs = seq(0, 1, 0.2), na.rm = TRUE), include.lowest = TRUE, labels = c("1", "2", "3", "4", "5"))
# Factor PHDI_ALL_quintile as Q1, Q2, Q3, Q4, Q5
PHDI_GHG_FINAL_DATA_0518$PHDI_ALL_quintile = factor(PHDI_GHG_FINAL_DATA_0518$PHDI_ALL_quintile, levels = c("1", "2", "3", "4", "5"), labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

# Factor RIAGENDR
PHDI_GHG_FINAL_DATA_0518$RIAGENDR = factor(PHDI_GHG_FINAL_DATA_0518$RIAGENDR, levels = c(1,2), labels = c("Male", "Female"))
# Factor RIDAGEYR_CAT (Age Group): 20-34 (1), 35-49 (2), 50-64 (3), and 65 or older (4)
PHDI_GHG_FINAL_DATA_0518$RIDAGEYR_CAT = factor(PHDI_GHG_FINAL_DATA_0518$RIDAGEYR_CAT, levels = c(1,2,3,4), labels = c("20-34", "35-49", "50-64", "> 65"))
# Factor RIDRETH1 (Race/Ethicity)
PHDI_GHG_FINAL_DATA_0518$RIDRETH1 = factor(PHDI_GHG_FINAL_DATA_0518$RIDRETH1, 
  levels = c(1, 2, 3, 4, 5), 
  labels = c("Mexican American", "Other Hispanic", "Non-Hispanic White", "Non-Hispanic Black", "Other Race")
)
# Factor DMDEDUC2_CAT (Education Level)
PHDI_GHG_FINAL_DATA_0518$DMDEDUC2_CAT = factor(
  PHDI_GHG_FINAL_DATA_0518$DMDEDUC2_CAT, 
  levels = c(2, 3, 4, 5), 
  labels = c("Less Than High School", "High School Graduate/GED", "Some College or Associate's Degree", "College Graduate or Above")
)
# Factor DMDMARTL_CAT (Marital Status)
PHDI_GHG_FINAL_DATA_0518$DMDMARTL_CAT = factor(
  PHDI_GHG_FINAL_DATA_0518$DMDMARTL_CAT, 
  levels = c(1, 0), 
  labels = c("Married", "Not Married")
)
# Factor INDFMPIR_CAT (Income to Poverty Ratio)
PHDI_GHG_FINAL_DATA_0518$INDFMPIR_CAT = factor(
  PHDI_GHG_FINAL_DATA_0518$INDFMPIR_CAT, 
  levels = c(1, 2, 3), 
  labels = c("0-1.3", "1.3-3.5", ">3.5")
)
# Factor BMI_CAT (BMI Category): underweight or normal weight (2), overweight (3), and obese (4)
PHDI_GHG_FINAL_DATA_0518$BMI_CAT = factor(
  PHDI_GHG_FINAL_DATA_0518$BMI_CAT, 
  levels = c(2, 3, 4), 
  labels = c("Underweight or Normal Weight", "Overweight", "Obese")
)

# update the survey design for combined 14 years NHANES data (2005-2018, combined 14 years)
PHDI_GHG_FINAL_DATA_0518_design_14yr = svydesign(
    id = ~SDMVPSU, 
    strata = ~SDMVSTRA, 
    weight = ~WTDR14D, 
    data = PHDI_GHG_FINAL_DATA_0518, #set up survey design on the full dataset #can restrict at time of analysis 
    nest = TRUE)

############################## Table 1: Multivariate logistic regression between PHDI and GHG, disease-related biomarkers or anthropometric measurements ##############################
## Exposure: PHDI_ALL_quintile
## Covariates: TOTALKCAL_PHDI, RIDAGEYR_CAT, RIAGENDR, RIDRETH1, INDFMPIR_CAT, DMDHHSIZ, DMDEDUC2_CAT, DMDMARTL_CAT, SMQ_CAT, ALQ130_CAT, DSD010_CAT, PAQ_CAT, BMI_CAT
## Outcome: BPXSY1 (systolic blood pressure, mmHg)
PHDI_ALL_quintile_BPXSY1 <- svyglm(BPXSY1 ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT + DMDMARTL_CAT + SMQ_CAT + ALQ130_CAT + DSD010_CAT + PAQ_CAT + BMI_CAT,
                        design = PHDI_GHG_FINAL_DATA_0518_design_14yr) #apply survey design

# use tbl_regression to get the regression table
PHDI_ALL_quintile_BPXSY1_tbl = tbl_regression(
    PHDI_ALL_quintile_BPXSY1, 
    exponentiate = FALSE,
    include = PHDI_ALL_quintile,
    label=list(PHDI_ALL_quintile ~ "PHDI quintiles"),
    intercept = TRUE
    ) 

## Exposure: PHDI_ALL_quintile
## Covariates: TOTALKCAL_PHDI, RIDAGEYR_CAT, RIAGENDR, RIDRETH1, INDFMPIR_CAT, DMDHHSIZ, DMDEDUC2_CAT, DMDMARTL_CAT, SMQ_CAT, ALQ130_CAT, DSD010_CAT, PAQ_CAT, BMI_CAT
## Outcome: BPXDI1 (diastolic blood pressure, mm Hg)
PHDI_ALL_quintile_BPXDI1 <- svyglm(BPXDI1 ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT + DMDMARTL_CAT + SMQ_CAT + ALQ130_CAT + DSD010_CAT + PAQ_CAT + BMI_CAT,
                        design = PHDI_GHG_FINAL_DATA_0518_design_14yr) #apply survey design

# use tbl_regression to get the regression table
PHDI_ALL_quintile_BPXDI1_tbl = tbl_regression(
    PHDI_ALL_quintile_BPXDI1, 
    exponentiate = FALSE,
    include = PHDI_ALL_quintile,
    label=list(PHDI_ALL_quintile ~ "PHDI quintiles"),
    intercept = TRUE
    ) 

## Exposure: PHDI_ALL_quintile
## Covariates: TOTALKCAL_PHDI, RIDAGEYR_CAT, RIAGENDR, RIDRETH1, INDFMPIR_CAT, DMDHHSIZ, DMDEDUC2_CAT, DMDMARTL_CAT, SMQ_CAT, ALQ130_CAT, DSD010_CAT, PAQ_CAT, BMI_CAT
## Outcome: HDL cholesterol (mg/dL), LBDHDD
PHDI_ALL_quintile_LBDHDD <- svyglm(LBDHDD ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT + DMDMARTL_CAT + SMQ_CAT + ALQ130_CAT + DSD010_CAT + PAQ_CAT + BMI_CAT,
                        design = PHDI_GHG_FINAL_DATA_0518_design_14yr) #apply survey design

# use tbl_regression to get the regression table
PHDI_ALL_quintile_LBDHDD_tbl = tbl_regression(
    PHDI_ALL_quintile_LBDHDD, 
    exponentiate = FALSE,
    include = PHDI_ALL_quintile,
    label=list(PHDI_ALL_quintile ~ "PHDI quintiles"),
    intercept = TRUE
    )


## Exposure: PHDI_ALL_quintile
## Covariates: TOTALKCAL_PHDI, RIDAGEYR_CAT, RIAGENDR, RIDRETH1, INDFMPIR_CAT, DMDHHSIZ, DMDEDUC2_CAT, DMDMARTL_CAT, SMQ_CAT, ALQ130_CAT, DSD010_CAT, PAQ_CAT, BMI_CAT
## Outcome: Triglyceride (mg/dL), LBXTR
PHDI_ALL_quintile_LBXTR <- svyglm(LBXTR ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT + DMDMARTL_CAT + SMQ_CAT + ALQ130_CAT + DSD010_CAT + PAQ_CAT + BMI_CAT,
                        design = PHDI_GHG_FINAL_DATA_0518_design_14yr) #apply survey design

# use tbl_regression to get the regression table
PHDI_ALL_quintile_LBXTR_tbl = tbl_regression(
    PHDI_ALL_quintile_LBXTR, 
    exponentiate = FALSE,
    include = PHDI_ALL_quintile,
    label=list(PHDI_ALL_quintile ~ "PHDI quintiles"),
    intercept = TRUE
    )

## Exposure: PHDI_ALL_quintile
## Covariates: TOTALKCAL_PHDI, RIDAGEYR_CAT, RIAGENDR, RIDRETH1, INDFMPIR_CAT, DMDHHSIZ, DMDEDUC2_CAT, DMDMARTL_CAT, SMQ_CAT, ALQ130_CAT, DSD010_CAT, PAQ_CAT, BMI_CAT
## Outcome: LDL cholesterol (mg/dL), LBDLDL
PHDI_ALL_quintile_LBDLDL <- svyglm(LBDLDL ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT + DMDMARTL_CAT + SMQ_CAT + ALQ130_CAT + DSD010_CAT + PAQ_CAT + BMI_CAT,
                        design = PHDI_GHG_FINAL_DATA_0518_design_14yr) #apply survey design

# use tbl_regression to get the regression table
PHDI_ALL_quintile_LBDLDL_tbl = tbl_regression(
    PHDI_ALL_quintile_LBDLDL, 
    exponentiate = FALSE,
    include = PHDI_ALL_quintile,
    label=list(PHDI_ALL_quintile ~ "PHDI quintiles"),
    intercept = TRUE
    )

## Exposure: PHDI_ALL_quintile
## Covariates: TOTALKCAL_PHDI, RIDAGEYR_CAT, RIAGENDR, RIDRETH1, INDFMPIR_CAT, DMDHHSIZ, DMDEDUC2_CAT, DMDMARTL_CAT, SMQ_CAT, ALQ130_CAT, DSD010_CAT, PAQ_CAT, BMI_CAT
## Outcome: Glycohemoglobin (%), LBXGH
PHDI_ALL_quintile_LBXGH <- svyglm(LBXGH ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + DMDHHSIZ + DMDEDUC2_CAT + DMDMARTL_CAT + SMQ_CAT + ALQ130_CAT + DSD010_CAT + PAQ_CAT + BMI_CAT,
                        design = PHDI_GHG_FINAL_DATA_0518_design_14yr) #apply survey design

# use tbl_regression to get the regression table
PHDI_ALL_quintile_LBXGH_tbl = tbl_regression(
    PHDI_ALL_quintile_LBXGH, 
    exponentiate = FALSE,
    include = PHDI_ALL_quintile,
    label=list(PHDI_ALL_quintile ~ "PHDI quintiles"),
    intercept = TRUE
    )

## Exposure: PHDI_ALL_quintile
## Covariates: TOTALKCAL_PHDI, RIDAGEYR_CAT, RIAGENDR, RIDRETH1, INDFMPIR_CAT, DMDHHSIZ, DMDEDUC2_CAT, DMDMARTL_CAT, SMQ_CAT, ALQ130_CAT, DSD010_CAT, PAQ_CAT, BMI_CAT
## Outcome: Fasting plasma glucose (mg/dL), LBXGLU
PHDI_ALL_quintile_LBXGLU <- svyglm(LBXGLU ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + DMDHHSIZ + DMDEDUC2_CAT + DMDMARTL_CAT + SMQ_CAT + ALQ130_CAT + DSD010_CAT + PAQ_CAT + BMI_CAT,
                        design = PHDI_GHG_FINAL_DATA_0518_design_14yr) #apply survey design

# use tbl_regression to get the regression table
PHDI_ALL_quintile_LBXGLU_tbl = tbl_regression(
    PHDI_ALL_quintile_LBXGLU, 
    exponentiate = FALSE,
    include = PHDI_ALL_quintile,
    label=list(PHDI_ALL_quintile ~ "PHDI quintiles"),
    intercept = TRUE
    )

## Exposure: PHDI_ALL_quintile
## Covariates: TOTALKCAL_PHDI, RIDAGEYR_CAT, RIAGENDR, RIDRETH1, INDFMPIR_CAT, DMDHHSIZ, DMDEDUC2_CAT, DMDMARTL_CAT, SMQ_CAT, ALQ130_CAT, DSD010_CAT, PAQ_CAT
## Outcome: BMI (kg/m2), BMXBMI
PHDI_ALL_quintile_BMXBMI <- svyglm(BMXBMI ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + DMDHHSIZ + DMDEDUC2_CAT + DMDMARTL_CAT + SMQ_CAT + ALQ130_CAT + DSD010_CAT + PAQ_CAT,
                        design = PHDI_GHG_FINAL_DATA_0518_design_14yr) #apply survey design

# use tbl_regression to get the regression table
PHDI_ALL_quintile_BMXBMI_tbl = tbl_regression(
    PHDI_ALL_quintile_BMXBMI, 
    exponentiate = FALSE,
    include = PHDI_ALL_quintile,
    label = list(PHDI_ALL_quintile ~ "PHDI quintiles"),
    intercept = TRUE
    ) 
## Exposure: PHDI_ALL_quintile
## Covariates: TOTALKCAL_PHDI, RIDAGEYR_CAT, RIAGENDR, RIDRETH1, INDFMPIR_CAT, DMDHHSIZ, DMDEDUC2_CAT, DMDMARTL_CAT, SMQ_CAT, ALQ130_CAT, DSD010_CAT, PAQ_CAT, BMI_CAT
## Outcome: waist circumference (cm), BMXWAIST
PHDI_ALL_quintile_BMXWAIST <- svyglm(BMXWAIST ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + DMDHHSIZ + DMDEDUC2_CAT + DMDMARTL_CAT + SMQ_CAT + ALQ130_CAT + DSD010_CAT + PAQ_CAT,
                        design = PHDI_GHG_FINAL_DATA_0518_design_14yr) #apply survey design

# use tbl_regression to get the regression table
PHDI_ALL_quintile_BMXWAIST_tbl = tbl_regression(
    PHDI_ALL_quintile_BMXWAIST, 
    exponentiate = FALSE,
    include = PHDI_ALL_quintile,
    label = list(PHDI_ALL_quintile ~ "PHDI quintiles"),
    intercept = TRUE
    )


# stack all the tables together using tbl_stack
PHDI_ALL_trend_tbl_list_table_1 <- list(
    PHDI_ALL_quintile_BPXSY1_tbl,
    PHDI_ALL_quintile_BPXDI1_tbl,
    PHDI_ALL_quintile_LBDHDD_tbl,
    PHDI_ALL_quintile_LBXTR_tbl,
    PHDI_ALL_quintile_LBDLDL_tbl,
    PHDI_ALL_quintile_LBXGH_tbl,
    PHDI_ALL_quintile_LBXGLU_tbl,
    PHDI_ALL_quintile_BMXBMI_tbl,
    PHDI_ALL_quintile_BMXWAIST_tbl)

# Final table 3
table_1_biomarker = tbl_stack(
    PHDI_ALL_trend_tbl_list_table_1,
    group_header = c("Systolic blood pressure, mm Hg", "Diastolic blood pressure, mm Hg", "HDL cholesterol, mg/dL", "Triglyceride, mg/dL", "LDL cholesterol, mg/dL", "Glycohemoglobin, %", "Fasting plasma glucose, mg/dL", "BMI, kg/m2", "Waist circumference, cm")
    ) 

# save table 3 as a word file
table_1_biomarker %>%
    as_flex_table() %>%
    flextable::save_as_docx(path = "table_1_disease_biomarker_measurement.docx")


############################## Table 2: Multivariate logistic regression between PHDI and disease ##############################
## Exposure: PHDI_ALL_quintile
## Covariates: TOTALKCAL_PHDI, RIDAGEYR_CAT, RIAGENDR, RIDRETH1, INDFMPIR_CAT, DMDHHSIZ, DMDEDUC2_CAT, DMDMARTL_CAT, SMQ_CAT, ALQ130_CAT, DSD010_CAT, PAQ_CAT
### (BMI category is not adjusted because it is not suitable for assessing obesity)
## Outcome: BMI category, OBESE
PHDI_ALL_quintile_OBESE <- svyglm(OBESE ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT + DMDMARTL_CAT + SMQ_CAT + ALQ130_CAT + DSD010_CAT + PAQ_CAT,
                        design = PHDI_GHG_FINAL_DATA_0518_design_14yr,
                        family = poisson(link = log)
                        ) #apply survey design

# use tbl_regression to get the regression table
PHDI_ALL_quintile_OBESE_tbl = tbl_regression(
    PHDI_ALL_quintile_OBESE, 
    exponentiate = TRUE,
    include = PHDI_ALL_quintile,
    label = list(PHDI_ALL_quintile ~ "PHDI quintiles")
    )

## Exposure: PHDI_ALL_quintile
## Covariates: TOTALKCAL_PHDI, RIDAGEYR_CAT, RIAGENDR, RIDRETH1, INDFMPIR_CAT, DMDHHSIZ, DMDEDUC2_CAT, DMDMARTL_CAT, SMQ_CAT, ALQ130_CAT, DSD010_CAT, PAQ_CAT, BMI_CAT
## Outcome: Abdominal obesity, AB_OBESE
PHDI_ALL_quintile_AB_OBESE <- svyglm(AB_OBESE ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT + DMDMARTL_CAT + SMQ_CAT + ALQ130_CAT + DSD010_CAT + PAQ_CAT,
                        design = PHDI_GHG_FINAL_DATA_0518_design_14yr,
                        family = poisson(link = log)
                        ) #apply survey design

# use tbl_regression to get the regression table
PHDI_ALL_quintile_AB_OBESE_tbl = tbl_regression(
    PHDI_ALL_quintile_AB_OBESE, 
    exponentiate = TRUE,
    include = PHDI_ALL_quintile,
    label = list(PHDI_ALL_quintile ~ "PHDI quintiles")
    )

# stack all the tables together using tbl_stack
PHDI_ALL_trend_tbl_list_table_2_obesity <- list(
    PHDI_ALL_quintile_OBESE_tbl,
    PHDI_ALL_quintile_AB_OBESE_tbl)

# Final table 2
table_2_obesity = tbl_stack(
    tbls = PHDI_ALL_trend_tbl_list_table_2_obesity,
    group_header = c("BMI obesity", "Abdominal obesity")
)

# save table 4 as a word file
table_2_obesity %>%
    as_flex_table() %>%
    flextable::save_as_docx(path = "table_2_obesity.docx")


########################## Table 3 mortality analysis ##############################

################# table 3 full model
# update the survey design for combined 12 years NHANES data (excluding 2017-2018)
PHDI_GHG_FINAL_DATA_0518_design_12yr = svydesign(
    id = ~SDMVPSU, 
    strata = ~SDMVSTRA, 
    weight = ~WTDR12D, 
    data = subset(PHDI_GHG_FINAL_DATA_0518, YEAR != "2017-2018"), #set up survey design on the full dataset #can restrict at time of analysis 
    nest = TRUE)

# Table 3: Hazard ratio
# full model
## covariates: PHDI_ALL_quintile, TOTALKCAL_PHDI, RIDAGEYR_CAT, RIAGENDR, RIDRETH1, INDFMPIR_CAT, DMDHHSIZ, DMDEDUC2_CAT, DMDMARTL_CAT, SMQ_CAT, ALQ130_CAT, DSD010_CAT, PAQ_CAT
table_3 = svycoxph(
    formula = Surv(permth_int, mortstat>0) ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIAGENDR + RIDAGEYR_CAT + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT + DMDMARTL_CAT + SMQ_CAT + ALQ130_CAT + DSD010_CAT + PAQ_CAT,
    design = PHDI_GHG_FINAL_DATA_0518_design_12yr
    )

# use tbl_regression to get the regression table
table_3_tbl = tbl_regression(
    table_3, 
    exponentiate = TRUE,
    include = PHDI_ALL_quintile,
    label = list(PHDI_ALL_quintile ~ "PHDI quintiles")
    ) %>%
    add_n(location = "level") %>%
    add_nevent(location = "level")

# save the table as a word file
table_3_tbl %>%
    as_flex_table() %>%
    flextable::save_as_docx(path = "table_3_mortality_full.docx")

############## table 3 main model
# covariates: PHDI_ALL_quintile (PHDI), TOTALKCAL_PHDI (total energy intake), RIDAGEYR_CAT (age group), RIAGENDR (gender), RIDRETH1 (race/ethinicity), INDFMPIR_CAT (PIR, poverty income ratio)
table_3_main = svycoxph(
    formula = Surv(permth_int, mortstat>0) ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIAGENDR + RIDAGEYR_CAT + RIDRETH1 + INDFMPIR_CAT,
    design = PHDI_GHG_FINAL_DATA_0518_design_12yr
    )

# use tbl_regression to get the regression table
table_3_main_tbl = tbl_regression(
    table_3_main, 
    exponentiate = TRUE,
    include = PHDI_ALL_quintile,
    label = list(PHDI_ALL_quintile ~ "PHDI quintiles")
    ) %>%
    add_n(location = "level") %>%
    add_nevent(location = "level")

# save table 3 as a word file
table_3_main_tbl %>%
    as_flex_table() %>%
    flextable::save_as_docx(path = "table_3_mortality_main.docx")

############## table s7 reduced model adjusted for potential mediator: obesity (OBESE), abdominal obesity (AB_OBESE), HDL cholesterol (LBDHDD), triglyceride (LBXTR), BMI (BMXBMI), waist circumference (BMXWAIST)
# covariates: PHDI_ALL_quintile (PHDI), TOTALKCAL_PHDI (total energy intake), RIDAGEYR_CAT (age group), RIAGENDR (gender), RIDRETH1 (race/ethinicity), INDFMPIR_CAT (PIR, poverty income ratio), OBESE (obesity)
table_s7_reduced_obesity_mediation = svycoxph(
    formula = Surv(permth_int, mortstat>0) ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIAGENDR + RIDAGEYR_CAT + RIDRETH1 + INDFMPIR_CAT + OBESE,
    design = PHDI_GHG_FINAL_DATA_0518_design_12yr
    )

# use tbl_regression to get the regression table
table_s7_reduced_obesity_mediation_tbl = tbl_regression(
    table_s7_reduced_obesity_mediation, 
    exponentiate = TRUE,
    include = PHDI_ALL_quintile,
    label = list(PHDI_ALL_quintile ~ "PHDI quintiles")
    ) %>%
    add_n(location = "level") %>%
    add_nevent(location = "level")

# covariates: PHDI_ALL_quintile (PHDI), TOTALKCAL_PHDI (total energy intake), RIDAGEYR_CAT (age group), RIAGENDR (gender), RIDRETH1 (race/ethinicity), INDFMPIR_CAT (PIR, poverty income ratio), AB_OBESE (abdominal obesity)
table_s7_reduced_abdominal_obesity_mediation = svycoxph(
    formula = Surv(permth_int, mortstat>0) ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIAGENDR + RIDAGEYR_CAT + RIDRETH1 + INDFMPIR_CAT + AB_OBESE,
    design = PHDI_GHG_FINAL_DATA_0518_design_12yr
    )

# use tbl_regression to get the regression table
table_s7_reduced_abdominal_obesity_mediation_tbl = tbl_regression(
    table_s7_reduced_abdominal_obesity_mediation, 
    exponentiate = TRUE,
    include = PHDI_ALL_quintile,
    label = list(PHDI_ALL_quintile ~ "PHDI quintiles")
    ) %>%
    add_n(location = "level") %>%
    add_nevent(location = "level")

# covariates: PHDI_ALL_quintile (PHDI), TOTALKCAL_PHDI (total energy intake), RIDAGEYR_CAT (age group), RIAGENDR (gender), RIDRETH1 (race/ethinicity), INDFMPIR_CAT (PIR, poverty income ratio), LBDHDD (HDL cholesterol)
table_s7_reduced_hdl_cholesterol_mediation = svycoxph(
    formula = Surv(permth_int, mortstat>0) ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIAGENDR + RIDAGEYR_CAT + RIDRETH1 + INDFMPIR_CAT + LBDHDD,
    design = PHDI_GHG_FINAL_DATA_0518_design_12yr
    )

# use tbl_regression to get the regression table
table_s7_reduced_hdl_cholesterol_mediation_tbl = tbl_regression(
    table_s7_reduced_hdl_cholesterol_mediation, 
    exponentiate = TRUE,
    include = PHDI_ALL_quintile,
    label = list(PHDI_ALL_quintile ~ "PHDI quintiles")
    ) %>%
    add_n(location = "level") %>%
    add_nevent(location = "level")

# covariates: PHDI_ALL_quintile (PHDI), TOTALKCAL_PHDI (total energy intake), RIDAGEYR_CAT (age group), RIAGENDR (gender), RIDRETH1 (race/ethinicity), INDFMPIR_CAT (PIR, poverty income ratio), LBXTR (triglyceride)
table_s7_reduced_triglyceride_mediation = svycoxph(
    formula = Surv(permth_int, mortstat>0) ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIAGENDR + RIDAGEYR_CAT + RIDRETH1 + INDFMPIR_CAT + LBXTR,
    design = PHDI_GHG_FINAL_DATA_0518_design_12yr
    )

# use tbl_regression to get the regression table
table_s7_reduced_triglyceride_mediation_tbl = tbl_regression(
    table_s7_reduced_triglyceride_mediation, 
    exponentiate = TRUE,
    include = PHDI_ALL_quintile,
    label = list(PHDI_ALL_quintile ~ "PHDI quintiles")
    ) %>%
    add_n(location = "level") %>%
    add_nevent(location = "level")

# covariates: PHDI_ALL_quintile (PHDI), TOTALKCAL_PHDI (total energy intake), RIDAGEYR_CAT (age group), RIAGENDR (gender), RIDRETH1 (race/ethinicity), INDFMPIR_CAT (PIR, poverty income ratio), BMXBMI (BMI)
table_s7_reduced_BMI_mediation = svycoxph(
    formula = Surv(permth_int, mortstat>0) ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIAGENDR + RIDAGEYR_CAT + RIDRETH1 + INDFMPIR_CAT + BMXBMI,
    design = PHDI_GHG_FINAL_DATA_0518_design_12yr
    )

# use tbl_regression to get the regression table
table_s7_reduced_BMI_mediation_tbl = tbl_regression(
    table_s7_reduced_BMI_mediation, 
    exponentiate = TRUE,
    include = PHDI_ALL_quintile,
    label = list(PHDI_ALL_quintile ~ "PHDI quintiles")
    ) %>%
    add_n(location = "level") %>%
    add_nevent(location = "level")

# covariates: PHDI_ALL_quintile (PHDI), TOTALKCAL_PHDI (total energy intake), RIDAGEYR_CAT (age group), RIAGENDR (gender), RIDRETH1 (race/ethinicity), INDFMPIR_CAT (PIR, poverty income ratio), BMXWAIST (waist circumference)
table_s7_reduced_waist_circumference_mediation = svycoxph(
    formula = Surv(permth_int, mortstat>0) ~ PHDI_ALL_quintile + TOTALKCAL_PHDI + RIAGENDR + RIDAGEYR_CAT + RIDRETH1 + INDFMPIR_CAT + BMXWAIST,
    design = PHDI_GHG_FINAL_DATA_0518_design_12yr
    )

# use tbl_regression to get the regression table
table_s7_reduced_waist_circumference_mediation_tbl = tbl_regression(
    table_s7_reduced_waist_circumference_mediation, 
    exponentiate = TRUE,
    include = PHDI_ALL_quintile,
    label = list(PHDI_ALL_quintile ~ "PHDI quintiles")
    ) %>%
    add_n(location = "level") %>%
    add_nevent(location = "level")

# stack all the tables together using tbl_stack
PHDI_ALL_trend_tbl_list_table_s7 <- list(
    table_s7_reduced_obesity_mediation_tbl,
    table_s7_reduced_abdominal_obesity_mediation_tbl,
    table_s7_reduced_hdl_cholesterol_mediation_tbl,
    table_s7_reduced_triglyceride_mediation_tbl,
    table_s7_reduced_BMI_mediation_tbl,
    table_s7_reduced_waist_circumference_mediation_tbl)

# Final table s7
table_s7 <- tbl_stack(
    tbls = PHDI_ALL_trend_tbl_list_table_s7,
    group_header = c("Obesity", "Abdominal obesity", "HDL cholesterol", "Triglyceride", "BMI", "Waist circumference")
)

setwd("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PHDI_NHANES/Publication/Results/Supplements")
# save table s7 as a word file
table_s7 %>%
    as_flex_table() %>%
    flextable::save_as_docx(path = "table_s7_mortality_mediation.docx")


# PHDI_GHG_FINAL_DATA_0518 <- PHDI_GHG_FINAL_DATA_0518 %>%
#   arrange(YEAR) %>%
#   mutate(PHDI_ALL_quintile_lag1 = lag(PHDI_ALL_quintile, n = 1))

# # install MedSurvey package for mediation analysis of the mortality
# install.packages("lavaan")
# install.packages("/Users/james/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Emory University - Ph.D./PHDI_NHANES/Publication/R script/MedSurvey", repos = NULL, type = "source")
# library(MedSurvey)
# summary(MedData)

# Table s2 ############################################################################################################
# Baseline cahracteristics of the US population by PHDI_ALL quintiles

# update the survey design for combined 14 years NHANES data (2005-2018, combined 14 years)
PHDI_GHG_FINAL_DATA_0518_design_14yr = svydesign(
    id = ~SDMVPSU, 
    strata = ~SDMVSTRA, 
    weight = ~WTDR14D, 
    data = PHDI_GHG_FINAL_DATA_0518, #set up survey design on the full dataset #can restrict at time of analysis 
    nest = TRUE)

# weighted table s2 code for PHDI_ALL quintiles
table_s2 = tbl_svysummary(
  
  data = PHDI_GHG_FINAL_DATA_0518_design_14yr, 
  
  by = PHDI_ALL_quintile, # stratify by PHDI_ALL quintile
  
  include = c("GHG_TOTAL", "TOTALKCAL_PHDI", "RIAGENDR", "RIDAGEYR_CAT", "RIDRETH1", "DMDEDUC2_CAT", "DMDMARTL_CAT", "INDFMPIR_CAT", "BMXBMI", "BMI_CAT", "BMXWAIST", "AB_OBESE"), 
  
  # GHG_TOTAL (Total Greenhouse Gas Emission (kg CO2), TOTALKCAL_PHDI (Total Energy Intake), RIAGENDR (Gender), RIDAGEYR_CAT (Age Group), RIDRETH1 (Race/Ethnicity), DMDEDUC2_CAT (Education Level), DMDMARTL_CAT (Marital Status), INDFMPIR_CAT (Household Poverty to Income Ratio), BMXBMI (BMI), BMI_CAT (BMI Category), BMXWAIST (Waist Circumference), AB_OBESE (abodominal obesity) 
  label = list(GHG_TOTAL ~ "Total Greenhouse Gas Emission (kg CO2) per day", TOTALKCAL_PHDI ~ "Total Energy Intake (kcal)", RIAGENDR ~ "Gender", RIDAGEYR_CAT ~ "Age Group", RIDRETH1 ~ "Race/Ethnicity", DMDEDUC2_CAT ~ "Education Level", DMDMARTL_CAT ~ "Marital Status", INDFMPIR_CAT ~ "Household Poverty to Income Ratio", BMXBMI ~ "Body Mass Index (BMI)", BMI_CAT ~ "BMI Category", BMXWAIST ~ "Waist Circumference", AB_OBESE ~ "Abdominal Obesity"),
  
  type = list(GHG_TOTAL ~ "continuous", TOTALKCAL_PHDI ~ "continuous", RIAGENDR ~ "categorical", RIDAGEYR_CAT ~ "categorical", RIDRETH1 ~ "categorical", DMDEDUC2_CAT ~ "categorical", DMDMARTL_CAT ~ "categorical", INDFMPIR_CAT ~ "categorical", BMXBMI ~ "continuous", BMI_CAT ~ "categorical", BMXWAIST ~ "continuous", AB_OBESE ~ "categorical"),
  
  statistic = list(all_categorical() ~ "{n} ({p}%)", all_continuous() ~ "{mean} ({mean.std.error})"),
  
  # list 1 decimal place for continuous variables
  digits = list(all_continuous() ~ 1, all_categorical() ~ 0),
  
  # tells R to ignore missing values in the summary statistics (e.g., means)
  missing = "no", 
  
  sort = NULL,
  # use column not row percent
  percent = "column")  %>%
  ## add p-values to the output comparing values across groups
  add_p() %>%
  ## add a column with overall summary statistics
  add_overall() %>%
  ## add n
  add_n() %>%
  # adding spanning header
  modify_spanning_header(all_stat_cols() ~ "**PHDI score quintiles in the combined NHANES data**") %>%
  ## Add a custom footnote explaining the summary statistics
  modify_footnote(
    all_stat_cols() ~ "Summary Statistics: n (%) = Weighted Frequency (Weighted Percentage), mean (SD) = Weighted Mean (Weighted Standard Error)"
    ) %>%
  modify_table_styling(
    columns = label,
    rows = label == "BMI Category",
    footnote = "BMI Category: Underweight (<18.5 kg/m^2), Normal weight (18.5-24.9 kg/m^2), Overweight (25.0-29.9 kg/m^2), Obese (≥30.0 kg/m^2)"
  ) %>%
  modify_table_styling(
    columns = label,
    rows = label == "Abdominal Obesity",
    footnote = "Abdominal Obesity Category: Men >= 102 cm (40 in) Waist Circumference, Women >= 88 cm (35 in) Waist Circumference"
  )

# save table s2 as a word file
table_s2 %>%
    as_flex_table() %>%
    flextable::save_as_docx(path = "PHDI_table_s2.docx")


# Table S3 ################################################################
# weighted table for all characteristics over NHANES cycle (weighted n and weighted %)
table_s3 = tbl_svysummary(
  
  data = PHDI_GHG_FINAL_DATA_0518_design, 
  
  by = YEAR, # stratify by YEAR
  
  include = c("PHDI_ALL", "GHG_TOTAL", "TOTALKCAL_PHDI", "RIAGENDR", "RIDAGEYR_CAT", "RIDRETH1", "DMDEDUC2_CAT", "DMDMARTL_CAT", "INDFMPIR_CAT", "BMXBMI", "BMI_CAT", "BMXWAIST", "AB_OBESE"),
  
  # PHDI_ALL (Planetary Health Diet Index Total Score), GHG_TOTAL (Total Greenhouse Gas Emission (kg CO2)), TOTALKCAL_PHDI (Total Energy Intake), RIAGENDR (Gender), RIDAGEYR_CAT (Age Group), RIDRETH1 (Race/Ethnicity), DMDEDUC2_CAT (Education Level), DMDMARTL_CAT (Marital Status), INDFMPIR_CAT (Household Poverty to Income Ratio), BMXBMI (BMI), BMI_CAT (BMI Category), BMXWAIST (Waist Circumference), AB_OBESE (abodominal obesity)
  label = list(PHDI_ALL ~ "Planetary Health Diet Index Total Score", GHG_TOTAL ~ "Total Greenhouse Gas Emission (kg CO2) per day", TOTALKCAL_PHDI ~ "Total Energy Intake (kcal)", RIAGENDR ~ "Gender", RIDAGEYR_CAT ~ "Age Group", RIDRETH1 ~ "Race/Ethnicity", DMDEDUC2_CAT ~ "Education Level", DMDMARTL_CAT ~ "Marital Status", INDFMPIR_CAT ~ "Household Poverty to Income Ratio", BMXBMI ~ "Body Mass Index (BMI)", BMI_CAT ~ "BMI Category", BMXWAIST ~ "Waist Circumference", AB_OBESE ~ "Abdominal Obesity"),
  
  type = list(PHDI_ALL ~ "continuous", GHG_TOTAL ~ "continuous", TOTALKCAL_PHDI ~ "continuous", RIAGENDR ~ "categorical", RIDAGEYR_CAT ~ "categorical", RIDRETH1 ~ "categorical", DMDEDUC2_CAT ~ "categorical", DMDMARTL_CAT ~ "categorical", INDFMPIR_CAT ~ "categorical", BMI_CAT ~ "categorical", BMXWAIST ~ "continuous", AB_OBESE ~ "categorical"),
  
  statistic = list(all_categorical() ~ "{n} ({p}%)", all_continuous() ~ "{mean} ({mean.std.error})"),
  
  # list 1 decimal place for continuous variables
  digits = list(all_continuous() ~ 1, all_categorical() ~ 0),
  
  missing = "no", #tells R to ignore missing values in the summary statistics (e.g., means)
  
  sort = NULL,
  # use column not row percent
  percent = "column")  %>%
  ## add p-values to the output comparing values across groups
  add_p() %>%
  ## add a column with overall summary statistics
  add_overall() %>%
  ## add n
  add_n() %>%
  # adding spanning header
  modify_spanning_header(all_stat_cols() ~ "**NHANES cycle**") %>%
  # Add a custom footnote explaining the summary statistics
  modify_footnote(
    all_stat_cols() ~ "Summary Statistics: n (%) = Weighted Frequency (Weighted Percentage), mean (SD) = Weighted Mean (Weighted Standard Error)"
    ) %>%
  modify_table_styling(
    columns = label,
    rows = label == "BMI Category",
    footnote = "BMI Category: Underweight (<18.5 kg/m^2), Normal weight (18.5-24.9 kg/m^2), Overweight (25.0-29.9 kg/m^2), Obese (≥30.0 kg/m^2)"
  ) %>%
  modify_table_styling(
    columns = label,
    rows = label == "Abdominal Obesity",
    footnote = "Abdominal Obesity Category: Men >= 102 cm (40 in) Waist Circumference, Women >= 88 cm (35 in) Waist Circumference"
  )

# save table S2 as a word file
table_s3 %>%
    as_flex_table() %>%
    flextable::save_as_docx(path = "PHDI_table_s3.docx")

######################################## Table s4: Covariate-adjusted PHDI trend, NHANES 2005-2018 ########################################
# Add a new variable, YEAR_mid
## if YEAR = 2005-2006, YEAR_mid = 2005.5
## if YEAR = 2007-2008, YEAR_mid = 2007.5
## if YEAR = 2009-2010, YEAR_mid = 2009.5
## if YEAR = 2011-2012, YEAR_mid = 2011.5
## if YEAR = 2013-2014, YEAR_mid = 2013.5
## if YEAR = 2015-2016, YEAR_mid = 2015.5
## if YEAR = 2017-2018, YEAR_mid = 2017.5
PHDI_GHG_FINAL_DATA_0518 = PHDI_GHG_FINAL_DATA_0518 %>%
    mutate(YEAR_mid = case_when(
        YEAR == "2005-2006" ~ 2005.5,
        YEAR == "2007-2008" ~ 2007.5,
        YEAR == "2009-2010" ~ 2009.5,
        YEAR == "2011-2012" ~ 2011.5,
        YEAR == "2013-2014" ~ 2013.5,
        YEAR == "2015-2016" ~ 2015.5,
        YEAR == "2017-2018" ~ 2017.5
    ))


# update the survey design for including the new variable YEAR_mid
PHDI_GHG_FINAL_DATA_0518_design = svydesign(
    id = ~SDMVPSU, 
    strata = ~SDMVSTRA, 
    weight = ~WTDR2D, 
    data = PHDI_GHG_FINAL_DATA_0518, #set up survey design on the full dataset #can restrict at time of analysis 
    nest = TRUE)

# Covariate-adjusted PHDI score estimated by multivariate linear regression analysis
# PHDI_ALL with YEAR_mid being exposure and all other covariates being confounders

## PHDI_ALL ~ YEAR_mid + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT
PHDI_ALL_trend_total_0 <- svyglm(PHDI_ALL ~ YEAR_mid + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT,
                        PHDI_GHG_FINAL_DATA_0518_design) #apply survey design 
# use tbl_regression to get the regression table
PHDI_ALL_trend_total_0_tbl = tbl_regression(
    PHDI_ALL_trend_total_0, 
    exponentiate = FALSE,
    include=YEAR_mid,
    label=list(YEAR_mid ~ "PHDI trend all")
    )

# PHDI_ALL with YEAR_mid being exposure and all other covariates, except RIDAGEYR_CAT, being confounders
## PHDI_ALL ~ YEAR_mid + TOTALKCAL_PHDI + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT
### apply survey design and subset RIDAGEYR_CAT == '20-34'
PHDI_ALL_trend_total_1 <- svyglm(PHDI_ALL ~ YEAR_mid + TOTALKCAL_PHDI + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT,
                        subset(PHDI_GHG_FINAL_DATA_0518_design, RIDAGEYR_CAT == '20-34')
                        ) 
# use tbl_regression to get the regression table
PHDI_ALL_trend_total_1_tbl = tbl_regression(
    PHDI_ALL_trend_total_1, 
    exponentiate = FALSE,
    include=YEAR_mid,
    label=list(YEAR_mid ~ "PHDI trend age group 20-34")
    )

### apply survey design and subset RIDAGEYR_CAT == '35-49'
PHDI_ALL_trend_total_2 <- svyglm(PHDI_ALL ~ YEAR_mid + TOTALKCAL_PHDI + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT,
                        subset(PHDI_GHG_FINAL_DATA_0518_design, RIDAGEYR_CAT == '35-49')
                        )
# use tbl_regression to get the regression table
PHDI_ALL_trend_total_2_tbl = tbl_regression(
    PHDI_ALL_trend_total_2, 
    exponentiate = FALSE,
    include=YEAR_mid,
    label=list(YEAR_mid ~ "PHDI trend age group 35-49")
    )

### apply survey design and subset RIDAGEYR_CAT == `50-64`
PHDI_ALL_trend_total_3 <- svyglm(PHDI_ALL ~ YEAR_mid + TOTALKCAL_PHDI + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT,
                        subset(PHDI_GHG_FINAL_DATA_0518_design, RIDAGEYR_CAT == '50-64')
                        )
# use tbl_regression to get the regression table
PHDI_ALL_trend_total_3_tbl = tbl_regression(
    PHDI_ALL_trend_total_3, 
    exponentiate = FALSE,
    include=YEAR_mid,
    label=list(YEAR_mid ~ "PHDI trend age group 50-64")
    )

### apply survey design and subset RIDAGEYR_CAT == '>65'
PHDI_ALL_trend_total_4 <- svyglm(PHDI_ALL ~ YEAR_mid + TOTALKCAL_PHDI + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT,
                        subset(PHDI_GHG_FINAL_DATA_0518_design, RIDAGEYR_CAT == '> 65')
                        )
# use tbl_regression to get the regression table
PHDI_ALL_trend_total_4_tbl = tbl_regression(
    PHDI_ALL_trend_total_4, 
    exponentiate = FALSE,
    include=YEAR_mid,
    label=list(YEAR_mid ~ "PHDI trend age group > 65")
    )

# PHDI_ALL with YEAR_mid being exposure and all other covariates, except RIAGENDR, being confounders
## PHDI_ALL ~ YEAR_mid + RIDAGEYR_CAT + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT
### apply survey design and subset RIAGENDR == 'Male'
PHDI_ALL_trend_total_5 <- svyglm(PHDI_ALL ~ YEAR_mid + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT,
                        subset(PHDI_GHG_FINAL_DATA_0518_design, RIAGENDR == 'Male')
                        )
# Use tbl_regression to get the regression table for males
PHDI_ALL_trend_total_5_tbl <- tbl_regression(
    PHDI_ALL_trend_total_5, 
    exponentiate = FALSE,
    include = YEAR_mid,
    label = list(YEAR_mid ~ "PHDI trend Male")
)

### apply survey design and subset RIAGENDR == 'Female'
PHDI_ALL_trend_total_6 <- svyglm(PHDI_ALL ~ YEAR_mid + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT,
                        subset(PHDI_GHG_FINAL_DATA_0518_design, RIAGENDR == 'Female')
                        )
# Use tbl_regression to get the regression table for females
PHDI_ALL_trend_total_6_tbl <- tbl_regression(
    PHDI_ALL_trend_total_6, 
    exponentiate = FALSE,
    include = YEAR_mid,
    label = list(YEAR_mid ~ "PHDI trend Female")
)

# Stratification based on RIDRETH1
# Assuming you have the levels for RIDRETH1 as "Mexican American", "Other Hispanic", "Non-Hispanic White", "Non-Hispanic Black", "Other Race"
levels_RIDRETH1 <- levels(PHDI_GHG_FINAL_DATA_0518$RIDRETH1)
PHDI_ALL_trend_RIDRETH1_tbl_list <- lapply(levels_RIDRETH1, function(eth_group) {
  model <- svyglm(PHDI_ALL ~ YEAR_mid + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT,
                  subset(PHDI_GHG_FINAL_DATA_0518_design, RIDRETH1 == eth_group))
  tbl_regression(model, exponentiate = FALSE, include = YEAR_mid, label = list(YEAR_mid ~ paste("PHDI trend", eth_group)))
})

# Stratification based on INDFMPIR_CAT
levels_INDFMPIR_CAT <- levels(PHDI_GHG_FINAL_DATA_0518$INDFMPIR_CAT)
PHDI_ALL_trend_INDFMPIR_CAT_tbl_list <- lapply(levels_INDFMPIR_CAT, function(pir_group) {
  model <- svyglm(PHDI_ALL ~ YEAR_mid + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + DMDHHSIZ + DMDEDUC2_CAT,
                  subset(PHDI_GHG_FINAL_DATA_0518_design, INDFMPIR_CAT == pir_group))
  tbl_regression(model, exponentiate = FALSE, include = YEAR_mid, label = list(YEAR_mid ~ paste("PHDI trend PIR", pir_group)))
})

# Stratification based on DMDEDUC2_CAT
levels_DMDEDUC2_CAT <- levels(PHDI_GHG_FINAL_DATA_0518$DMDEDUC2_CAT)
PHDI_ALL_trend_DMDEDUC2_CAT_tbl_list <- lapply(levels_DMDEDUC2_CAT, function(edu_group) {
  model <- svyglm(PHDI_ALL ~ YEAR_mid + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ,
                  subset(PHDI_GHG_FINAL_DATA_0518_design, DMDEDUC2_CAT == edu_group))
  tbl_regression(model, exponentiate = FALSE, include = YEAR_mid, label = list(YEAR_mid ~ paste("PHDI trend", edu_group)))
})

# Stratification based on BMI_CAT
levels_BMI_CAT <- levels(PHDI_GHG_FINAL_DATA_0518$BMI_CAT)
PHDI_ALL_trend_BMI_CAT_tbl_list <- lapply(levels_BMI_CAT, function(bmi_group) {
  model <- svyglm(PHDI_ALL ~ YEAR_mid + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT,
                  subset(PHDI_GHG_FINAL_DATA_0518_design, BMI_CAT == bmi_group))
  tbl_regression(model, exponentiate = FALSE, include = YEAR_mid, label = list(YEAR_mid ~ paste("PHDI trend", bmi_group)))
})

# Stratification based on AB_OBESE
levels_AB_OBESE <- levels(PHDI_GHG_FINAL_DATA_0518$AB_OBESE)
PHDI_ALL_trend_AB_OBESE_tbl_list <- lapply(levels_AB_OBESE, function(ab_group) {
  model <- svyglm(PHDI_ALL ~ YEAR_mid + TOTALKCAL_PHDI + RIDAGEYR_CAT + RIAGENDR + RIDRETH1 + INDFMPIR_CAT + DMDHHSIZ + DMDEDUC2_CAT,
                  subset(PHDI_GHG_FINAL_DATA_0518_design, AB_OBESE == ab_group))
  tbl_regression(model, exponentiate = FALSE, include = YEAR_mid, label = list(YEAR_mid ~ paste("PHDI trend abdominal obesity", ab_group)))
})

# stack all the tables together using tbl_stack
# Merge first 6 tables in a list
PHDI_ALL_trend_tbl_list <- list(
  PHDI_ALL_trend_total_0_tbl,
  PHDI_ALL_trend_total_1_tbl,
  PHDI_ALL_trend_total_2_tbl,
  PHDI_ALL_trend_total_3_tbl,
  PHDI_ALL_trend_total_4_tbl,
  PHDI_ALL_trend_total_5_tbl,
  PHDI_ALL_trend_total_6_tbl
)

# Final linear trend table
table_s3_linear_trend = tbl_stack(c(PHDI_ALL_trend_tbl_list, PHDI_ALL_trend_RIDRETH1_tbl_list, PHDI_ALL_trend_INDFMPIR_CAT_tbl_list, PHDI_ALL_trend_DMDEDUC2_CAT_tbl_list, PHDI_ALL_trend_BMI_CAT_tbl_list, PHDI_ALL_trend_AB_OBESE_tbl_list))

# save table s3 as a word file
table_s3_linear_trend %>%
    as_flex_table() %>%
    flextable::save_as_docx(path = "PHDI_table_s4_covariate_adjusted_PHDI_linear_trend.docx")

# Supplementary table 5: Correlation between PHDI and other dietary indexes
# Compute Pearson correlation
correlation_PHDI_AHEI_HEI2020_DASH_AMED_DII <- svycor(~ PHDI_ALL + AHEI_ALL + HEI2020_ALL + DASH_ALL + MED_ALL + DII_ALL, design = PHDI_GHG_FINAL_DATA_0518_design_14yr)
print(correlation_PHDI_AHEI_HEI2020_DASH_AMED_DII)

