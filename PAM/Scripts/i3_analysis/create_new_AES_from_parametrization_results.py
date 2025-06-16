import pandas as pd
import os
import numpy as np

from PAModelpy import CatalyticEvent


RESULT_PARAMETRIZATION_FILE = os.path.join('Results', '2_parametrization', 'pam_parametrizer_diagnostics_iML1515_randomized.xlsx')
NEW_AES_SUFFIX = 'iMl1515_241209'
SECTOR_PARAM_FILE = os.path.join('Results','2_parametrization','proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx')


def search_index_in_parameter_file(df:pd.DataFrame, protein:str, reaction:str, direction:str):
    all_protein_rows = df[df.enzyme_id == protein]
    all_rxn_and_protein_rows = all_protein_rows[all_protein_rows.rxn_id == reaction]
    the_row = all_rxn_and_protein_rows[all_rxn_and_protein_rows.direction == direction]
    return the_row.index

def create_new_aes_parameter_file(old_param_file:str = SECTOR_PARAM_FILE,
                                  result_file_path:str = RESULT_PARAMETRIZATION_FILE,
                                  new_aes_suffix: str = NEW_AES_SUFFIX):

    parameter_files = pd.read_excel(old_param_file, sheet_name=None)
    aes_parameter_file =parameter_files['ActiveEnzymes']

    parametrization_results = pd.read_excel(
        result_file_path, sheet_name='Best_Individuals').drop_duplicates(
        ['rxn_id', 'direction', 'enzyme_id'], keep='last')

    # extract reaction id from catalytic event id
    parametrization_results['rxn_id'] = [CatalyticEvent._extract_reaction_id_from_catalytic_reaction_id(id) for id in
                                         parametrization_results['rxn_id']]

    for index, row in parametrization_results.iterrows():
        aes_index = search_index_in_parameter_file(aes_parameter_file, row.enzyme_id, row.rxn_id, row.direction)
        aes_parameter_file.loc[aes_index, 'kcat_values'] = row['kcat[s-1]']

    write_mode = 'w'
    kwargs = {}
    result_file = os.path.join('Results', '3_analysis', 'parameter_files',
                               f'proteinAllocationModel_EnzymaticData_{new_aes_suffix}.xlsx')
    if os.path.isfile(result_file):
        write_mode = 'a'
        kwargs = {'if_sheet_exists': 'replace'}

    with pd.ExcelWriter(result_file,
            mode=write_mode, engine='openpyxl', **kwargs) as writer:
        aes_parameter_file.to_excel(writer, sheet_name='ActiveEnzymes', index=False)
        for sheet, df in parameter_files.items():
            if sheet != 'ActiveEnzymes':
                df.to_excel(writer, sheet_name = sheet, index=False)

def gaussian(x, mean, amplitude, standard_deviation):
    return amplitude * np.exp( - (x - mean)**2 / (2*standard_deviation ** 2))

if __name__ == '__main__':
    # other_files = [os.path.join('Results', '3_analysis', 'parameter_files',
    #                            'proteinAllocationModel_EnzymaticData_iML1515_241009.xlsx')]
    #
    for file_nmbr in range(1,11):
        suffix = f'iML1515_{file_nmbr}'
        result_file = os.path.join('Results', '2_parametrization', 'diagnostics', f'pam_parametrizer_diagnostics_{file_nmbr}.xlsx')
        create_new_aes_parameter_file(result_file_path= result_file,
                                      new_aes_suffix= suffix)
    # result_file = os.path.join('Results', '2_parametrization', 'diagnostics',
    #                            f'pam_parametrizer_diagnostics_mciML1515.xlsx')
    # create_new_aes_parameter_file(result_file_path= result_file,
    #                                   new_aes_suffix= 'mciML1515')

