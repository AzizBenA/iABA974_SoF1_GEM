from Bio import SeqIO
import csv
import pandas as pd

# Define the paths to the input GBFF file, CSV file, and output GBFF file
gbff_input_path = 'merged_gbff_EC_Protein.gbff'
csv_file_path = 'Genome_annotation_prokka.csv'
gbff_output_path = 'merged_gbff_EC_Protein_prokka.gbff'
i = 0
# Read data from CSV
df = pd.read_csv(csv_file_path)

# Read the GBFF file
seq_record = SeqIO.read(gbff_input_path, 'genbank')

# Iterate through features
for feature in seq_record.features:
    # Check if the feature is a CDS
    if feature.type == 'CDS' and 'hypothetical protein' in feature.qualifiers.get('product', []):
        # Get the gene position from the location
        gene_position = int(feature.location.start) + 1  # Convert to integer
        # Check if the gene position is in the CSV data
        if gene_position in df['cds_position'].values:
            protein_function = df.loc[df['cds_position'] == gene_position, 'protein_function'].values[0]
            if not pd.isna(protein_function) and protein_function != 'hypothetical protein':
                print(f'gene_position: {gene_position}')
                print(f'gene_product: {protein_function}')
                i += 1
                feature.qualifiers['product'] = [protein_function]

# Write the modified SeqRecord to a new GBFF file
with open(gbff_output_path, 'w') as output_gbff:
    SeqIO.write(seq_record, output_gbff, 'genbank')

print(f'Modified GBFF file "{gbff_output_path}" has been created.')
print(i)
