# Figure 3: Association between PHDI quintiles and GHG emission and mortality
# This script creates a dual-axis plot showing HR for mortality and GHG emissions

library(survey)
library(ggplot2)
library(dplyr)

# Assuming 'design' is already loaded with survey design object
# and nhanes_merged_0518_filtered contains the data

# Step 1: Calculate GHG means and SE by PHDI quintiles
ghg_by_quintile <- svyby(
  formula = ~GHG_TOTAL, 
  by = ~PHDI_cat, 
  design = design, 
  FUN = svymean, 
  na.rm = TRUE
)

# Add quintile labels and calculate 95% CI
ghg_data <- data.frame(
  quintile = c("Q1", "Q2", "Q3", "Q4", "Q5"),
  ghg_mean = ghg_by_quintile$GHG_TOTAL,
  ghg_se = ghg_by_quintile$se,
  ghg_lower = ghg_by_quintile$GHG_TOTAL - 1.96 * ghg_by_quintile$se,
  ghg_upper = ghg_by_quintile$GHG_TOTAL + 1.96 * ghg_by_quintile$se
)

# Step 2: Calculate HR for mortality by PHDI quintiles
# Fit Cox proportional hazards model with Q1 as reference
formula_cox <- as.formula(paste(
  "Surv(RIDAGEYR, mortstat > 0) ~ PHDI_cat + RIAGENDR + RIDRETH1 +",
  "DMDEDUC2_CAT + DMDMARTL_CAT + INDFMPIR_CAT + SMQ_CAT +",
  "TOTALKCAL_PHDI_mean + ALQ130_CAT + DMDHHSIZ + DSD010"
))

model_cox <- svycoxph(formula_cox, design = design)

# Extract HR and CI for each quintile
coef_summary <- summary(model_cox)
phdi_rows <- grep("^PHDI_cat", rownames(coef_summary$coefficients))

# Create HR data frame (Q1 is reference with HR=1)
hr_data <- data.frame(
  quintile = c("Q1", "Q2", "Q3", "Q4", "Q5"),
  hr = c(1, exp(coef_summary$coefficients[phdi_rows, "coef"])),
  hr_lower = c(1, coef_summary$conf.int[phdi_rows, "lower .95"]),
  hr_upper = c(1, coef_summary$conf.int[phdi_rows, "upper .95"])
)

# Step 3: Merge data for plotting
plot_data <- merge(ghg_data, hr_data, by = "quintile")
plot_data$quintile_num <- 1:5

# Step 4: Create dual-axis plot
# Scale factor to align the two y-axes
ghg_range <- range(c(plot_data$ghg_lower, plot_data$ghg_upper))
hr_range <- range(c(plot_data$hr_lower, plot_data$hr_upper))

# Calculate scaling factor
scale_factor <- diff(ghg_range) / diff(hr_range)
offset <- ghg_range[1] - hr_range[1] * scale_factor

# Create plot
p <- ggplot(plot_data, aes(x = quintile_num)) +
  # GHG line and error bars (blue)
  geom_line(aes(y = ghg_mean), color = "#00BFFF", linewidth = 1.2) +
  geom_point(aes(y = ghg_mean), color = "#00BFFF", size = 3) +
  geom_errorbar(aes(ymin = ghg_lower, ymax = ghg_upper), 
                color = "#00BFFF", width = 0.2, linewidth = 0.8) +
  
  # HR line and error bars (red) - scaled to GHG axis
  geom_line(aes(y = hr * scale_factor + offset), color = "#E74C3C", linewidth = 1.2) +
  geom_point(aes(y = hr * scale_factor + offset), color = "#E74C3C", size = 3) +
  geom_errorbar(aes(ymin = hr_lower * scale_factor + offset, 
                    ymax = hr_upper * scale_factor + offset), 
                color = "#E74C3C", width = 0.2, linewidth = 0.8) +
  
  # Primary y-axis (HR)
  scale_y_continuous(
    name = "HR for all-cause mortality",
    limits = c(ghg_range[1], ghg_range[2]),
    sec.axis = sec_axis(
      trans = ~(. - offset) / scale_factor,
      name = "Total Greenhouse Gas Emissions (kg CO2) per person per day",
      breaks = seq(0.5, 2, by = 0.25)
    )
  ) +
  
  # X-axis
  scale_x_continuous(
    breaks = 1:5,
    labels = paste0("Q", 1:5),
    name = "Quintiles of PHDI"
  ) +
  
  # Theme
  theme_classic() +
  theme(
    axis.title.y.left = element_text(color = "#E74C3C", size = 11, face = "bold"),
    axis.text.y.left = element_text(color = "#E74C3C", size = 10),
    axis.title.y.right = element_text(color = "#00BFFF", size = 11, face = "bold"),
    axis.text.y.right = element_text(color = "#00BFFF", size = 10),
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
  )

# Add legend manually
p <- p + 
  annotate("text", x = 2.5, y = ghg_range[2] * 0.95, 
           label = "HR", color = "#E74C3C", size = 4, fontface = "bold") +
  annotate("text", x = 3.5, y = ghg_range[2] * 0.95, 
           label = "Greenhouse Gas Emissions (kg CO2) per person per day", 
           color = "#00BFFF", size = 4, fontface = "bold")

print(p)

# Save the plot
ggsave("Figure3_PHDI_GHG_Mortality.png", 
       plot = p, 
       width = 10, 
       height = 6, 
       dpi = 400)

# Print summary statistics
cat("\n=== GHG by PHDI Quintile ===\n")
print(ghg_data)

cat("\n=== HR for Mortality by PHDI Quintile ===\n")
print(hr_data)
