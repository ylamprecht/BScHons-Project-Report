# This script processes the output .txt file from VIPR, to create a .csv file 
# retaining only the genes, nucleotide variants, and VIPR Pathogenicity scores.

import csv

def extract_c_dot(sequence, start_str='c.', end_chars=[':', ';']):
    """
    Extract all substrings from the input sequence that start with 'c.'
    and end with the first occurrence of specified end characters.

    Parameters:
    sequence (str): The input string from which to extract substrings.
    start_str (str): The starting string to look for (default is 'c.').
    end_chars (list): A list of characters that denote the end of the substring.

    Returns:
    list: A list of extracted substrings.
    """
    extracted_list = []  # List to store extracted substrings
    start_idx = sequence.find(start_str)  # Find the first occurrence of 'c.'
    
    while start_idx != -1:
        # Find the earliest occurrence of an end character
        end_idx = len(sequence)
        for char in end_chars:
            idx = sequence.find(char, start_idx)  # Search for end character
            if idx != -1:
                end_idx = min(end_idx, idx)  # Update end index to the earliest found
        
        # Append the extracted part from 'c.' to ':' or ';' or end of the line
        extracted_list.append(sequence[start_idx:end_idx])
        start_idx = sequence.find(start_str, end_idx)  # Look for more 'c.' in the sequence

    return extracted_list

def process_tsv(input_file, output_file):
    """
    Process the input TSV file to extract Gene, Nucleotide and VIPR Pathogenicity score
    and save it to a CSV file.

    Parameters:
    input_file (str): The path to the input .txt file in TSV format.
    output_file (str): The path to the output .csv file.
    """
    with open(input_file, 'r') as tsvfile, open(output_file, 'w', newline='') as csvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        fieldnames = ['Gene', 'Nucleotide', 'VIPR_Pathogenicity']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        # Process each row in the input file
        for row in reader:
            gene = row['Gene.refGene']  # Extract gene name
            vipr_pathogenicity = row['.pred_P_LP']  # Extract VIPR Pathogenicity score
            gene_detail = row['GeneDetail.refGene'] 
            aa_change = row['AAChange.refGene']
            
            # Determine which field to extract c_dot values from
            if gene_detail == ".":
                c_dots = extract_c_dot(aa_change, start_str='c.', end_chars=[':'])
            elif aa_change == ".":
                c_dots = extract_c_dot(gene_detail, start_str='c.', end_chars=[';', '\n'])
            else:
                c_dots = []

            # Write a new row for each extracted c_dot as "Nucleotide"
            for c_dot in c_dots:
                writer.writerow({
                    'Gene': gene,
                    'Nucleotide': c_dot,
                    'VIPR_Pathogenicity': vipr_pathogenicity
                })

# Define input and output file paths
input_file = 'prioritised_file.txt'
output_file = 'extracted_prioritised_file_2.csv'

# Run the processing function
process_tsv(input_file, output_file)
