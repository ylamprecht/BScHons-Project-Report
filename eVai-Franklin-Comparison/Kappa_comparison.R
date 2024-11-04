# R script to calculate Cohen’s weighted Kappa as an index of interrater agreement
# between eVai and Franklin on ordinal data (variant classifications)

# Load the required packages
library(irr)
library(ggplot2)

# List of files containing variant classifications
files <- c("CVD1_common_variants.csv", "CVD2_common_variants.csv", 
           "CVD3_common_variants.csv", "CVD4_common_variants.csv", 
           "CVD12_common_variants.csv", "CVD20_common_variants.csv", 
           "CVD21_common_variants.csv", "CVD22_common_variants.csv", 
           "CVD23_common_variants.csv", "CVD26_common_variants.csv", 
           "CVD27_common_variants.csv", "CVD29_common_variants.csv", 
           "CVD30_common_variants.csv", "CVD38_common_variants.csv", 
           "CVD39_common_variants.csv", "CVD40_common_variants.csv", 
           "CVD42_common_variants.csv", "CVD43_common_variants.csv", 
           "CVD45_common_variants.csv", "CVD46_common_variants.csv", "CVD_common_variants.csv")

sample_names <- c("CVD1", "CVD2", "CVD3", "CVD4", "CVD12", "CVD20", 
                  "CVD21", "CVD22", "CVD23", "CVD26", "CVD27", "CVD29", 
                  "CVD30", "CVD38", "CVD39", "CVD40", "CVD42", "CVD43", 
                  "CVD45", "CVD46", "overall")

# Initialise a vector to store Kappa values for each sample
kappa_values <- numeric(length(files))

#  Loop through each file, calculate Kappa, and store the result
for (i in seq_along(files)) {
  
  df <- read.csv(files[i])
  
  # Map Franklin (Genoox_Classification) values to numeric codes
  df$Franklin_Code <- with(df, ifelse(Genoox_Classification == "BENIGN", 1,
                                      ifelse(Genoox_Classification == "LIKELY_BENIGN", 2,
                                             ifelse(Genoox_Classification %in% c("UNCERTAIN_SIGNIFICANCE", "POSSIBLY_BENIGN", 
                                                                                 "POSSIBLY_PATHOGENIC_LOW", "POSSIBLY_PATHOGENIC_MODERATE"), 3,
                                                    ifelse(Genoox_Classification == "LIKELY_PATHOGENIC", 4,
                                                           ifelse(Genoox_Classification == "PATHOGENIC", 5, NA))))))
  
  # Map eVai_Classification values to numeric codes
  df$eVai_Code <- with(df, ifelse(eVai_Classification == "Benign", 1,
                                  ifelse(eVai_Classification == "Likely benign", 2,
                                         ifelse(eVai_Classification == "Uncertain significance", 3,
                                                ifelse(eVai_Classification == "Likely pathogenic", 4,
                                                       ifelse(eVai_Classification == "Pathogenic", 5, NA))))))
  
  # Create the new n*2 dataframe
  df_classifications <- df[, c("Franklin_Code", "eVai_Code")]
  
  # Calculate the weighted Kappa statistic using squared (Fleiss-Cohen) weights
  kappa_result <- kappa2(df_classifications, weight = "squared")
  
  # Store the Kappa value in the vector
  kappa_values[i] <- kappa_result$value
  
  # Print out the Kappa result for each sample
  cat(sample_names[i], ": Kappa =", kappa_result$value, "\n")
}

# Create a data frame for plotting
plot_df <- data.frame(Sample = sample_names, Kappa = kappa_values)

# Create a scatter plot of Kappa values across samples
ggplot(plot_df, aes(x = Sample, y = Kappa)) +
  geom_point(size = 3) +
  labs(title = "Weighted Kappa Values Across Samples",
       x = "Sample", y = "Kappa Value") +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1))

