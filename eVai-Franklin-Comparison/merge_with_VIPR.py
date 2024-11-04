# This script merges eVai or Franklin outputs with VIPR Pathogenicity scores.
# It allows for comparison between the classifications provided by   
# eVai or Franklin and the VIPR Pathogenicity scores.

import os
import csv

def read_csv_to_dict(file_path, key_fields):
    """
    Reads a CSV file and returns a dictionary where each entry's key is created
    from specified fields, allowing easy look-up for merging data.
    
    Parameters:
    - file_path (str): Path to the CSV file to read.
    - key_fields (list): List of fields to use as the dictionary keys.
    
    Returns:
    - dict: Dictionary of rows keyed by values in `key_fields`.
    """
    data_dict = {}
    with open(file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # Rename 'Nucleotide' field to 'HGVS_Coding' for consistency
            if 'Nucleotide' in row:
                row['HGVS_Coding'] = row.pop('Nucleotide')
            # Create a tuple key based on specified key_fields
            key = tuple(row[field] for field in key_fields)
            data_dict[key] = row  # Add row to dictionary with the tuple key
    return data_dict

def merge_csv_files(input_folder, key_fields):
    """
    Merges all sample CSV files, removing duplicates based on
    the Gene and Nucleotide key.
    
    Parameters:
    - input_folder (str): Directory path containing CSV files to merge.
    - key_fields (list): Fields used to identify unique entries.
    
    Returns:
    - dict: Merged data from all CSV files, with unique entries based on key_fields.
    """
    merged_data = {}  # Dictionary to store unique merged data
    all_files = [f for f in os.listdir(input_folder) if f.endswith('.csv')]

    for file in all_files:
        file_path = os.path.join(input_folder, file)
        with open(file_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                # Create a tuple key based on key_fields to identify unique entries
                key = tuple(row[field] for field in key_fields)
                # Only add row if the key is not already in merged_data (remove duplicates)
                if key not in merged_data:
                    merged_data[key] = row
    
    return merged_data

def merge_with_prioritised(merged_data, prioritised_file, key_fields, output_file):
    """
    Merges the previously merged CSV data with the VIPR output file to include the
    VIPR_Pathogenicity score where applicable.
    
    Parameters:
    - merged_data (dict): Dictionary of data from merged CSV files.
    - prioritised_file (str): Path to the prioritised CSV file with VIPR_Pathogenicity.
    - key_fields (list): Fields to match entries between files.
    - output_file (str): Path to the output CSV file.
    """
    # Load prioritised data from CSV, using key_fields to create unique keys
    prioritised_data = read_csv_to_dict(prioritised_file, key_fields)
    
    # Prepare the list of fieldnames for the output file, adding 'VIPR_Pathogenicity' column
    fieldnames = list(next(iter(merged_data.values())).keys()) + ['VIPR_Pathogenicity']
    
    # Write the final merged data to the output CSV file
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()  # Write the header row for the CSV
        
        for key, row in merged_data.items():
            # Check if key is present in prioritised_data to add VIPR_Pathogenicity
            if key in prioritised_data:
                row['VIPR_Pathogenicity'] = prioritised_data[key]['VIPR_Pathogenicity']
            else:
                row['VIPR_Pathogenicity'] = ''  # If no match, leave VIPR_Pathogenicity blank
            writer.writerow(row)  # Write the row to the CSV

# Example usage:
input_folder = 'Franklin_output/Cohort_5_clean_eVai_outputs'  # Folder with CSV files to merge
output_file = 'eVai_VIPR_cohort5.csv'  # Final output file with merged data and VIPR_Pathogenicity
prioritised_file = 'extracted_prioritised_file.csv'  # File with VIPR_Pathogenicity information
key_fields = ['Gene', 'HGVS_Coding']  # Fields used to uniquely identify and merge rows

# Step 1: Merge all CSV files in the input folder, removing duplicates based on key_fields
merged_data = merge_csv_files(input_folder, key_fields)

# Step 2: Merge the merged data with the prioritised file to add VIPR_Pathogenicity information
merge_with_prioritised(merged_data, prioritised_file, key_fields, output_file)
