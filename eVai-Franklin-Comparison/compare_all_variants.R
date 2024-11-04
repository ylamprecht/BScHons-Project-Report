# This R script generates Venn Diagrams to:
# 1) Compare the genes found by eVai and Franklin
# 2) Compare the variants found by eVai and Franklin

# Load required libraries
library(dplyr)
library(VennDiagram)
library(grid)

# Read the data
franklin_data <- read.csv("clean_raw_Franklin_outputs/combined_franklin_compare.csv", stringsAsFactors = FALSE)
evai_data <- read.csv("clean_eVai_outputs/combined_eVai_compared.csv", stringsAsFactors = FALSE)

# Standardise the Franklin data for direct comparison to eVai data
# Replace 'n.' with 'c.' in the Nucleotide field
franklin_data$Nucleotide <- gsub("^n\\.", "c.", franklin_data$Nucleotide)

# Create unique identifiers for Franklin data
franklin_data <- franklin_data %>%
  mutate(Variant_ID = paste(Gene, Nucleotide, sep = "_"))

# Function to update HGVS_Coding in eVai data where Nucleotide is a substring
update_evai_hgvs <- function(franklin_df, evai_df) {
  updated_rows <- evai_df
  # Loop through each row in Franklin data
  for (i in seq_len(nrow(franklin_df))) {
    gene <- franklin_df$Gene[i]
    nucleotide <- franklin_df$Nucleotide[i]
    updated_rows <- updated_rows %>%
      mutate(
        HGVS_Coding = if_else(
          Gene == gene & grepl(nucleotide, HGVS_Coding, fixed = TRUE),
          nucleotide,  # Update HGVS_Coding if condition is met
          HGVS_Coding  # Keep original if condition is not met
        )
      )
  }
  return(updated_rows)
}

# Update eVai data using the function defined above
updated_evai_data <- update_evai_hgvs(franklin_data, evai_data)

# Create unique identifiers for updated eVai data
updated_evai_data <- updated_evai_data %>%
  mutate(Variant_ID = paste(Gene, HGVS_Coding, sep = "_"))

# Find common variants based on Variant_ID
common_variants <- franklin_data %>%
  inner_join(updated_evai_data, by = "Variant_ID") %>%
  rename(eVai_Classification = Classification)

# Find unique variants in Franklin
unique_franklin <- franklin_data %>%
  anti_join(common_variants %>% select(Variant_ID), by = "Variant_ID")

# Find unique variants in eVai
unique_evai <- updated_evai_data %>%
  anti_join(common_variants %>% select(Variant_ID), by = "Variant_ID")

# Create a Venn diagram for the variants found
venn_data <- list(
  Franklin = franklin_data$Variant_ID,
  eVai = updated_evai_data$Variant_ID
)

venn_plot <- venn.diagram(
  x = venn_data,
  category.names = c("Franklin", "eVai"),
  filename = NULL,
  fill = c("#FF9999", "#99CCFF"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  main = "All Variants",
  main.cex = 2,
  main.fontfamily = "Arial",
  show.percent = TRUE,
  print.mode = c("percent", "raw"),
  fontfamily = "Arial",
  cat.fontfamily = "Arial", 
  cat.dist = c(0.01, 0.01)
)

grid.draw(venn_plot)

# Create a Venn diagram for the genes found
venn_data_genes <- list(
  Franklin = unique(franklin_data$Gene),
  eVai = unique(updated_evai_data$Gene)
)

venn_plot_genes <- venn.diagram(
  x = venn_data_genes,
  category.names = c("Franklin", "eVai"),
  filename = NULL,
  fill = c("#FF9999", "#99CCFF"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  main = "All Genes",
  main.cex = 2,
  main.fontfamily = "Arial",
  show.percent = TRUE,
  print.mode = c("percent", "raw"),
  fontfamily = "Arial",
  cat.fontfamily = "Arial",
  cat.dist = c(0.01, 0.01)
)

grid.draw(venn_plot_genes)

# Write results to new CSV files
write.csv(common_variants, "comparison_results/CVD_common_variants.csv", row.names = FALSE)
write.csv(unique_franklin, "comparison_results/CVD_unique_franklin_variants.csv", row.names = FALSE)
write.csv(unique_evai, "comparison_results/CVD_unique_evai_variants.csv", row.names = FALSE)
