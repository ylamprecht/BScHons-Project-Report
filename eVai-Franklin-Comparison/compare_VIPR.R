# This R script creates a Boxplot showing the distribution of VIPR pathogenicity
# scores across each of eVai's/ Franklin's ACMG-AMP variant classifications

# Load required libraries
library(ggplot2)
library(dplyr)

tool <- "Franklin"  # Set the tool to either "eVai" or "Franklin"

# Conditional processing based on the tool
if (tool == "eVai") {
  # Read the data for eVai
  data <- read.csv("eVai_VIPR_cohort5.csv")
  
  # Recode the Classification column for eVai to standardise labels
  data$Classification <- case_when(
    data$Classification == "Uncertain significance" ~ "VUS",
    TRUE ~ data$Classification  # Keep all other values unchanged
  )
  
} else if (tool == "Franklin") {
  # Read the data for Franklin
  data <- read.csv("Franklin_VIPR_cohort5.csv")
  
  # Recode the Genoox_Classification column for Franklin to standarise labels
  data$Classification <- case_when(
    data$Genoox_Classification == "PATHOGENIC" ~ "Pathogenic",
    data$Genoox_Classification == "LIKELY_PATHOGENIC" ~ "Likely pathogenic",
    data$Genoox_Classification %in% c("UNCERTAIN_SIGNIFICANCE", "POSSIBLY_BENIGN", 
                                      "POSSIBLY_PATHOGENIC_LOW", "POSSIBLY_PATHOGENIC_MODERATE") ~ "VUS",
    data$Genoox_Classification == "LIKELY_BENIGN" ~ "Likely benign",
    data$Genoox_Classification == "BENIGN" ~ "Benign",
    TRUE ~ NA_character_  
  )
}

# Filter to keep only rows with non-missing VIPR_Pathogenicity scores
data_filtered <- data %>%
  filter(!is.na(VIPR_Pathogenicity))

# Ensure the classification is in the specified order (not alphabetical)
data_filtered$Classification <- factor(data_filtered$Classification, 
                                       levels = c("Pathogenic", "Likely pathogenic", "VUS", "Likely benign", "Benign"))

# Generate the boxplot to visualise VIPR pathogenicity scores by classification
ggplot(data_filtered, aes(x = Classification, y = VIPR_Pathogenicity)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)), vjust = -0.5, color = "red") +
  labs(x = "Franklin Classification", y = "VIPR Pathogenicity Score", 
       title = "Distribution of VIPR Pathogenicity Scores\n by Franklin Classification") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.margin = margin(t = 30, r = 20, b = 20, l = 20)
  )
