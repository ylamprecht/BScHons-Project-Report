# This script merges and refines Franklin output CSV files, retaining only the
# Gene, Nucleotide, Classification, Zygosity, and Inheritance Model columns 
# to streamline subsequent analysis.

import csv

# Define file paths
input_file_path_default = 'Cohort_3_raw/CVD930_single_snp_variants.csv'
input_file_path_UTR = 'Cohort_3_raw/CVD930_single_snp_variants (1).csv'
output_file_path = 'Cohort_3/CVD930.csv'

# Define the columns to keep
columns_to_keep = ["Gene", "Nucleotide", "Genoox_Classification", "Zygosity", "Inheritance_Model"]

# Merge the two input files
merged_rows = []
with open(input_file_path_default, 'r', encoding='utf-8') as infile1, open(input_file_path_UTR, 'r', encoding='utf-8') as infile2:
    reader1 = csv.DictReader(infile1)
    reader2 = csv.DictReader(infile2)
    
    # Append rows from the first file
    for row in reader1:
        merged_rows.append(row)
    
    # Append rows from the second file
    for row in reader2:
        merged_rows.append(row)

# Open a new file to save the cleaned and filtered data
with open(output_file_path, 'w', newline='', encoding='utf-8') as outfile:
    writer = csv.DictWriter(outfile, fieldnames=columns_to_keep)

    # Write the header
    writer.writeheader()

    # Track seen entries to detect duplicates
    seen_entries = set()
    duplicate_entries = []

    # Process each row in the merged data
    for row in merged_rows:
        # Create a unique identifier based on Gene and Nucleotide
        unique_id = (row["Gene"], row["Nucleotide"])

        # Check if the entry is a duplicate
        if unique_id in seen_entries:
            duplicate_entries.append(row)
        else:
            seen_entries.add(unique_id)
            # Clean up the data by removing extra quotes and keeping only the columns of interest
            cleaned_row = {key: row[key].replace('""', '').strip('"') for key in columns_to_keep}
            writer.writerow(cleaned_row)

    # Print out any duplicates
    if duplicate_entries:
        print("Duplicate entries found:")
        for entry in duplicate_entries:
            print(entry)
    else:
        print("No duplicates found.")

print("Data merging, cleaning, and duplicate check complete.")
