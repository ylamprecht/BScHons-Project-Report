# This R script generates Venn Diagrams to compare the:
# 1) Pathogenic/ Likely Pathogenic variants, 
# 2) Variants of Uncertain Significance (VUS) and 
# 3) Benign/ Likely Benign variants found by eVai and Franklin

# Load required libraries
library(dplyr)
library(VennDiagram)

# Read the data
common_variants <- read.csv("CVD_common_variants.csv", stringsAsFactors = FALSE)

# Filter Pathogenic/Likely Pathogenic variants
pathogenic_variants <- common_variants %>%
  filter(
    Genoox_Classification %in% c("PATHOGENIC", "LIKELY_PATHOGENIC") |
      eVai_Classification %in% c("Pathogenic", "Likely pathogenic")
  )

# Filter VUS variants
vus_variants <- common_variants %>%
  filter(
    Genoox_Classification %in% c("UNCERTAIN_SIGNIFICANCE", "POSSIBLY_BENIGN", "POSSIBLY_PATHOGENIC_LOW", "POSSIBLY_PATHOGENIC_MODERATE") |
      eVai_Classification %in% c("Uncertain significance")
  )

# Filter Benign/Likely Benign variants
benign_variants <- common_variants %>%
  filter(
    Genoox_Classification %in% c("LIKELY_BENIGN", "BENIGN") |
      eVai_Classification %in% c("Likely benign", "Benign")
  )

# Pathogenic/Likely Pathogenic Venn data
venn_data_pathogenic <- list(
  Franklin = pathogenic_variants %>% filter(Genoox_Classification %in% c("PATHOGENIC", "LIKELY_PATHOGENIC")) %>% pull(Variant_ID),
  eVai = pathogenic_variants %>% filter(eVai_Classification %in% c("Pathogenic", "Likely pathogenic")) %>% pull(Variant_ID)
)

# VUS Venn data
venn_data_vus <- list(
  Franklin = vus_variants %>% filter(Genoox_Classification %in% c("UNCERTAIN_SIGNIFICANCE", "POSSIBLY_BENIGN", "POSSIBLY_PATHOGENIC_LOW", "POSSIBLY_PATHOGENIC_MODERATE")) %>% pull(Variant_ID),
  eVai = vus_variants %>% filter(eVai_Classification %in% c("Uncertain significance")) %>% pull(Variant_ID)
)

# Benign/Likely Benign Venn data
venn_data_benign <- list(
  Franklin = benign_variants %>% filter(Genoox_Classification %in% c("LIKELY_BENIGN", "BENIGN")) %>% pull(Variant_ID),
  eVai = benign_variants %>% filter(eVai_Classification %in% c("Likely benign", "Benign")) %>% pull(Variant_ID)
)

# Function to create and display Venn diagrams
create_venn_diagram <- function(venn_data, title) {
  venn_plot <- venn.diagram(
    x = venn_data,
    category.names = c("Franklin", "eVai"),
    filename = NULL,
    fill = c("#FF9999", "#99CCFF"),
    alpha = 0.5,
    cex = 2,
    cat.cex = 2,
    show.percent = FALSE,
    print.mode = c("percent", "raw"),
    fontfamily = "arial",
    cat.fontfamily = "arial",
    cat.dist = c(0.01, 0.01)
  )
  grid.draw(venn_plot)
}

# Generate and display Venn diagrams for each variant classification
create_venn_diagram(venn_data_pathogenic, "Pathogenic/ Likely Pathogenic Variants")
create_venn_diagram(venn_data_vus, "VUS Variants")
create_venn_diagram(venn_data_benign, "Benign/ Likely Benign Variants")
