from Bio import SeqIO
import pandas as pd

# Define the paths to the input GBFF file, CSV file, and output GBFF file
gbff_input_path = r'input\merged_gbff_v2.gbff'
tsv_file_path = r'input\DeepECv2_result.tsv'
gbff_output_path = r'output\merged_gbff_v4.gbff'
i = 0
j = 0

# Read data from CSV
df_EC = pd.read_csv(tsv_file_path, sep='\t')

# Read the GBFF file
seq_record = SeqIO.read(gbff_input_path, 'genbank')

# Iterate through features
for feature in seq_record.features:
    if feature.type == 'CDS':
        locus_tag = feature.qualifiers.get('locus_tag', [])[0]  # This variable is to get the locus tag
        # Check if the gene position is in the CSV data
        if locus_tag in df_EC['locus_tag'].values:
            deepEC_number = df_EC.loc[df_EC['locus_tag'] == locus_tag, 'EC'].values[0]
            # Check if 'EC_number' key exists
            if 'EC_number' in feature.qualifiers:
                existing_EC_number = feature.qualifiers['EC_number'][0]
                # Check if 'EC_number' is different from 'deepEC_number'
                if existing_EC_number != deepEC_number:
                    print(f'1{locus_tag}')
                    i += 1
                    feature.qualifiers['DeepEC_number'] = [deepEC_number]
            else:
                print(f'2{locus_tag}')
                j += 1
                # If 'EC_number' doesn't exist, directly assign 'DeepEC_number'
                feature.qualifiers['DeepEC_number'] = [deepEC_number]

# Write the modified SeqRecord to a new GBFF file
with open(gbff_output_path, 'w') as output_gbff:
    SeqIO.write(seq_record, output_gbff, 'genbank')

print(f'Modified GBFF file "{gbff_output_path}" has been created.')
print(f'Iterations: {i}, Breaks: {j}')
