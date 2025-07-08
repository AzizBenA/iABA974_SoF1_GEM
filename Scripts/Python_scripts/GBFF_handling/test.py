import pandas as pd

def find_EC(file_path):
    # Read the tab-delimited file into a DataFrame
    df_EGGNOG = pd.read_csv(file_path, sep='\t')
    df_EGGNOG.index = range(1, len(df_EGGNOG) + 1)

    tsv_data = {'start': [], 'EC': []}  # Initialize an empty dictionary with lists

    for index, row in df_EGGNOG.iterrows():
        # Iterate over columns
        for column_name, value in row.items():
            if 'em_EC=' in str(value):
                tsv_data['start'].append(row['start'])
                tsv_data['EC'].append(value.replace('em_EC=', 'EC='))

    df_EC = pd.DataFrame(tsv_data)
    df_EC.to_csv('EGGNOG_EC.csv', index=False)
    return print('Done')

find_EC('TSV_EGGNOG.txt')
