# This script refines eVai output TXT files, retaining only the
# Gene, HGVS_Coding, Classification, Zygosity, Condition Inheritance  
# and Pathogenicity Socre columns to streamline subsequent analysis.

import csv

input_file_path = 'eVai_outputs/CVD46_eVai.txt'
output_file_path = 'clean_eVai_outputs/CVD46_eVai_compared.csv'

columns_to_keep = ["Gene", "HGVS_Coding", "Classification", "Sample_Zygosity", "Condition_Inheritance", "Pathogenicity_Score"]

def clean_data(input_file_path, output_file_path, columns_to_keep):
    with open(input_file_path, 'r') as infile, open(output_file_path, 'w', newline='') as outfile:
        headers_found = False
        headers = []
        seen_entries = set()
        
        writer = csv.writer(outfile)
        
        for line in infile:
            line = line.strip()
            if line.startswith("##"):
                continue  # Skip metadata lines
            
            # Process the header line
            if not headers_found:
                if line.startswith("Gene"):
                    headers = line.split(",")
                    writer.writerow(columns_to_keep)  # Write the header to the output file
                    headers_found = True
                continue
            
            # Process data lines
            # Split the line into gene and the rest of the data
            if line:
                # Extract gene name and the rest of the line
                gene_name, rest_of_line = line.split(',', 1)
                
                # Handle rest of the line: Split by '","' and clean entries
                line_data = rest_of_line.split('","')
                line_data = [gene_name.strip('"')] + [entry.strip('"') for entry in line_data]
                
                # Create a dictionary from the line data
                data_dict = dict(zip(headers, line_data))
                
                # Check for duplicates
                gene = data_dict.get("Gene", "")
                hgvs_coding = data_dict.get("HGVS_Coding", "")
                entry_key = (gene, hgvs_coding)
                
                if entry_key not in seen_entries:
                    seen_entries.add(entry_key)
                    
                    # Clean and select columns to keep
                    cleaned_row = [data_dict.get(key, "") for key in columns_to_keep]
                    
                    # Write the cleaned row to the output file
                    writer.writerow(cleaned_row)

# Call the function
clean_data(input_file_path, output_file_path, columns_to_keep)
