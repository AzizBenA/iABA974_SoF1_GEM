from Bio import SeqIO
import csv
import pandas as pd

# Path to the GBFF file
gbff_input_path = r'input\merged_gbff_v2.gbff'

# Path to the output files
output_csv_path = r'output\SoF1_merged_csv.csv'
output_excel_path = r'output\SoF1_excel_merged.xlsx'

input_handle = open(gbff_input_path, "r")

# fixing the columns of the dataframe

columns = ['sequence_id', 'type', 'start', 'stop', 'strand', 'locus_tag', 'gene', 'product', 'db_xref', 'EC',
           'protein_translation']


i=1
# Create an empty list to store rows
rows = []
for seq_record in SeqIO.parse(input_handle, "genbank"):
    for seq_feature in seq_record.features:
        if seq_feature.type != "gene" and seq_feature.type != "source":
            # Create a new row dictionary
            row = {}

            row['sequence_id'] = f'SoF1_gene_{i}'
            row['type'] = seq_feature.type
            row['start'] = int(seq_feature.location.start) + 1
            row['stop'] = int(seq_feature.location.end) + 1
            row['strand'] = seq_feature.location.strand
            locus_tag = seq_feature.qualifiers.get('locus_tag', [''])[0]
            row['locus_tag'] = locus_tag if locus_tag else ''
            row['gene'] = seq_feature.qualifiers.get('gene', [''])[0]
            row['product'] = seq_feature.qualifiers['product'][0] if 'product' in seq_feature.qualifiers else "Unknown"
            row['db_xref'] = ', '.join(seq_feature.qualifiers.get('db_xref', ['']))
            row['EC'] = seq_feature.qualifiers.get('EC_number', [''])[0]
            row['protein_translation'] = translation = seq_feature.qualifiers['translation'][
                0] if 'translation' in seq_feature.qualifiers else ""

            # Append the row to the list
            rows.append(row)

            i += 1

# Convert the list of rows to a DataFrame
df = pd.DataFrame(rows, columns=columns)
df.to_csv(output_csv_path, index=False)

# Export to Excel
df.to_excel(output_excel_path, index=False)