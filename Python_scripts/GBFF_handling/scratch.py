from Bio import SeqIO
import csv
import pandas as pd

csv_file_path = r'input\SoF1_merged_csv.csv'
tsv_file_path = r'input\DeepECv2_result.tsv'

# Path to the output files
output_csv_path = r'output\SoF1_merged_csv.csv'
output_excel_path = r'output\SoF1_excel_merged.xlsx'

df_EC = pd.read_csv(tsv_file_path, sep='\t')
i = 0

df_csv = pd.read_csv(csv_file_path)

# Add a new 'deepEC' column with NaN values
df_csv['deepEC'] = float('nan')

for value in df_EC['locus_tag']:
    deepEC_number = df_EC.loc[df_EC['locus_tag'] == value, 'EC'].values[0]
    EC_number = df_csv.loc[df_csv['locus_tag'] == value, 'EC'].values[0]
    deepEC_number_purified = deepEC_number.replace('EC=', '')

    if deepEC_number_purified not in str(EC_number):
        print(deepEC_number, EC_number)
        i += 1
        # Set the 'deepEC' value for the specific row
        df_csv.loc[df_csv['locus_tag'] == value, 'deepEC'] = deepEC_number

# Insert the 'deepEC' column at index 10
df_csv.insert(10, 'deepEC', df_csv.pop('deepEC'))

print(i)

df_csv.to_csv(output_csv_path, index=False)

# Export to Excel
df_csv.to_excel(output_excel_path, index=False)