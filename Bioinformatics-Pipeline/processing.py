"""
This script is the third step in the bioinformatics pipeline: Cohort Comparison.
It merges and refines CSV files, counts variant occurrences, and combines cohort data
to identify candidate host genetic variants associated with severe COVID-19 phenotypes.

The script consists of three main functions:
1. clean_franklin: Cleans and merges input CSV files for each sample by removing
   unnecessary columns and duplicate rows.
2. count_variants: Counts and categorises variants by classification for each cohort,
   generating summary files by gene and variant.
3. combine_cohorts: Combines variant count data across multiple cohorts to provide
   a consolidated view of variant occurrences across cohorts.

Output:
- Generates eight output files, one for each classification, summarising variant 
  distribution across cohorts.

Author: Yentl Lamprecht (24952818@sun.ac.za)
Date: 2024/11/04
"""

import csv
import os
from collections import defaultdict


def clean_franklin(input_file_path_default, input_file_path_UTR, output_file_path):
    """
    Cleans and merges two Franklin input CSV files for a sample by removing unnecessary columns
    and duplicate rows. Writes the cleaned data to an output CSV file.

    Parameters:
    - input_file_path_default (str): Path to the first input CSV file (default variants).
    - input_file_path_UTR (str): Path to the second input CSV file (UTR variants).
    - output_file_path (str): Path to save the cleaned and merged output CSV file.

    Returns:
    - None: Writes the cleaned data to the specified output file.
    """
    # Columns to retain in the output file for easier analysis
    columns_to_keep = ["Gene", "Nucleotide", "Genoox_Classification", "Zygosity", "Inheritance_Model"]

    # Read and merge rows from both input files
    merged_rows = []
    with open(input_file_path_default, 'r', encoding='utf-8') as infile1, open(input_file_path_UTR, 'r', encoding='utf-8') as infile2:
        reader1 = csv.DictReader(infile1)
        reader2 = csv.DictReader(infile2)
        
        for row in reader1:
            merged_rows.append(row)
        for row in reader2:
            merged_rows.append(row)

    # Write the cleaned and filtered data to a new file
    with open(output_file_path, 'w', newline='', encoding='utf-8') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=columns_to_keep)
        writer.writeheader()

        seen_entries = set()
        duplicate_entries = []

        # Process and clean each row
        for row in merged_rows:
            unique_id = (row["Gene"], row["Nucleotide"])
            if unique_id in seen_entries:
                duplicate_entries.append(row)
            else:
                seen_entries.add(unique_id)
                cleaned_row = {key: row[key].replace('""', '').strip('"') for key in columns_to_keep}
                writer.writerow(cleaned_row)

        if duplicate_entries:
            print(f"Duplicate entries found in {output_file_path}:")
            for entry in duplicate_entries:
                print(entry)
        else:
            print(f"No duplicates found in {output_file_path}.")
    print("Data cleaning complete.")

  
def count_variants(path_to_csv_files):
    """
    Counts and categorises variants by classification for each cohort.
    For each classification, it groups data by Gene_Nucleotide and Gene, creating two summary
    CSV files for further analysis.

    Parameters:
    - path_to_csv_files (str): Path to the directory containing CSV files for the cohort samples.

    Returns:
    - None: Generates CSV files summarising variant counts by Gene_Nucleotide and Gene.
    """

    # Extract cohort name from the path (e.g., 'Cohort_1' from 'Cohort_1/')
    cohort_name = os.path.basename(path_to_csv_files)

    classifications = ['PATHOGENIC', 'LIKELY_PATHOGENIC', 'UNCERTAIN_SIGNIFICANCE', 
                      'POSSIBLY_PATHOGENIC_LOW', 'POSSIBLY_PATHOGENIC_MODERATE', 
                      'POSSIBLY_BENIGN', 'LIKELY_BENIGN', 'BENIGN']
    
    gene_nucleotide_data = {classification: defaultdict(list) for classification in classifications}
    gene_data = {classification: defaultdict(set) for classification in classifications}

    # Process each .csv file in the specified directory
    for filename in os.listdir(path_to_csv_files):
        if filename.endswith('.csv'):
            sample_name = filename.replace('.csv', '')
            file_path = os.path.join(path_to_csv_files, filename)
            with open(file_path, mode='r') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    gene = row['Gene']
                    nucleotide = row['Nucleotide']
                    genoox_classification = row['Genoox_Classification']
                    gene_nucleotide = f"{gene}_{nucleotide}"

                    # Group data by classification
                    if genoox_classification in gene_nucleotide_data:
                        gene_nucleotide_data[genoox_classification][gene_nucleotide].append(sample_name)
                        gene_data[genoox_classification][gene].add(sample_name)

    # Output results for each classification
    for classification in classifications:
        gene_nucleotide_output = [
            {
                'Gene_Nucleotide': gene_nucleotide,
                'Sample_Count': len(samples),
                'Samples': ', '.join(samples)
            }
            for gene_nucleotide, samples in gene_nucleotide_data[classification].items()
        ]
        gene_nucleotide_output.sort(key=lambda x: x['Sample_Count'], reverse=True)
        # Write output files to the variant_classifications directory
        output_dir = 'variant_classifications'
        os.makedirs(output_dir, exist_ok=True)
        cohort_dir = os.path.join(output_dir, cohort_name)
        os.makedirs(cohort_dir, exist_ok=True)

        with open(os.path.join(cohort_dir, f'{cohort_name}_{classification}_gene_nucleotide.csv'), mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=['Gene_Nucleotide', 'Sample_Count', 'Samples'])
            writer.writeheader()
            writer.writerows(gene_nucleotide_output)

    print("Variant counts complete.")


def combine_cohorts(classifications, cohorts):
    """
    Combines variant count data across multiple cohorts for a set of classifications.
    Each cohort's data is merged based on the Gene_Nucleotide identifier to provide an overview
    of the variant counts across all cohorts.

    Parameters:
    - classifications (list): List of variant classifications to include in the merged data.
    - cohorts (list): List of cohort names to process and merge.

    Returns:
    - None: Outputs a CSV file for each classification showing the combined counts across cohorts.
    """

    # Helper function to read data from CSV for each cohort
    def read_csv_to_dict(filepath, cohort_name):
        data = {}
        with open(filepath, newline='', mode='r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                gene_nucleotide = row['Gene_Nucleotide']
                sample_count = int(row['Sample_Count'])
                data[gene_nucleotide] = {cohort_name: sample_count}
        return data

    # Merge data across cohorts for a given classification
    def merge_cohort_data(classification):
        merged_data = {}
        for cohort in cohorts:
            file_path = os.path.join("variant_classifications", cohort, f"{cohort}_{classification}_gene_nucleotide.csv")
            # Check if the file exists before reading
            if os.path.exists(file_path):
                cohort_data = read_csv_to_dict(file_path, cohort)
                for gene_nucleotide, counts in cohort_data.items():
                    if gene_nucleotide not in merged_data:
                        merged_data[gene_nucleotide] = {cohort: 0 for cohort in cohorts}
                    merged_data[gene_nucleotide].update(counts)
            else:
                print(f"Warning: The file {file_path} does not exist. Skipping...")

        output_dir = "variant_classifications/all_cohorts"
        os.makedirs(output_dir, exist_ok=True)

        # Output file for merged data
        output_file = os.path.join(output_dir, f"combined_cohorts_{classification}_gene_nucleotide.csv")
        with open(output_file, mode='w', newline='') as csvfile:
            fieldnames = ['Gene_Nucleotide'] + cohorts + ['Total']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for gene_nucleotide, counts in merged_data.items():
                total = counts['Cohort_1'] + counts['Cohort_2'] + counts['Cohort_3'] + counts['Cohort_4']
                writer.writerow({
                    'Gene_Nucleotide': gene_nucleotide,
                    'Cohort_1': counts['Cohort_1'],
                    'Cohort_2': counts['Cohort_2'],
                    'Cohort_3': counts['Cohort_3'],
                    'Cohort_4': counts['Cohort_4'],
                    'Total': total,
                    'Cohort_5': counts['Cohort_5']
                })

    for classification in classifications:
        merge_cohort_data(classification)
    print("Merging complete.")


# Main function to run the entire pipeline
def run_pipeline():
    # Define the list of cohorts to process
    cohort_folders = ['Cohort_1_raw', 'Cohort_2_raw', 'Cohort_3_raw', 'Cohort_4_raw', 'Cohort_5_raw']
    
    # Define corresponding output folders
    cohort_output_folders = ['Cohort_1', 'Cohort_2', 'Cohort_3', 'Cohort_4', 'Cohort_5']

    # Loop through each cohort folder
    for cohort_folder, cohort_output_folder in zip(cohort_folders, cohort_output_folders):
        # List all files in the current cohort folder
        files = os.listdir(cohort_folder)
        os.makedirs(cohort_output_folder, exist_ok=True)

        # Check if the cohort folder is empty
        if not files:
            print(f"Skipping {cohort_folder}: folder is empty.")
            continue
        
        # Process each sample (two input files per sample)
        samples = set()
        for file in files:
            if 'single_snp_variants' in file:
                # Extract sample name
                sample_name = file.split('_single_snp_variants')[0]
                samples.add(sample_name)
        
        # Check if there are any valid samples to process
        if not samples:
            print(f"No valid samples found in {cohort_folder}.")
            continue
        
        # For each sample, clean and merge the input files
        for sample in samples:
            file_default = os.path.join(cohort_folder, f"{sample}_single_snp_variants.csv")
            file_UTR = os.path.join(cohort_folder, f"{sample}_single_snp_variants (1).csv")
            output_file = os.path.join(cohort_output_folder, f"{sample}.csv")

            # Clean the files for the current sample
            clean_franklin(file_default, file_UTR, output_file)
        
        # Count variants for the current cohort
        count_variants(cohort_output_folder)

    # Combine data from different cohorts
    classifications = [
        "BENIGN", "LIKELY_BENIGN", "LIKELY_PATHOGENIC", "PATHOGENIC",
        "POSSIBLY_BENIGN", "POSSIBLY_PATHOGENIC_LOW", "POSSIBLY_PATHOGENIC_MODERATE", "UNCERTAIN_SIGNIFICANCE"
    ]
    
    # Define cohort names for combining data
    cohorts = ['Cohort_1', 'Cohort_2', 'Cohort_3', 'Cohort_4', 'Cohort_5']
    
    # Combine cohorts
    combine_cohorts(classifications, cohorts)
    print("Pipeline complete.")

# Run the pipeline
if __name__ == "__main__":
    run_pipeline()  
