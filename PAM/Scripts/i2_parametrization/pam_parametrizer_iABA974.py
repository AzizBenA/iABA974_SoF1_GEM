import os
import sys
import pandas as pd
import numpy as np
import warnings
from typing import List
warnings.filterwarnings("ignore")



from PAModelpy.configuration import Config

from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults
from Modules.PAM_parametrizer import PAMParametrizer
from Scripts.pam_generation import setup_sof1_pam
from PAModelpy.utils.pam_generation import increase_kcats_in_parameter_file
from Modules.utils.pamparametrizer_setup import set_up_sector_config



MAX_SUBSTRATE_UPTAKE_RATE = 0.5
MIN_SUBSTRATE_UPTAKE_RATE = -3

def set_up_validation_data() -> list[ValidationData]:

    c_uptake_id = 'EX_co2_e'
    validation_fluxes = pd.read_excel(r"C:\Users\adam-1y9pxi6q95ggpmx\Documents\Python Scripts\PAM_Parametrization\Data\SoF1_phenotypes\sof1_phenotypes.xlsx",
        sheet_name ='fluxes'
    )[['Growth', 'EX_co2_e', 'EX_h2_e', 'EX_o2_e']]
    validation_fluxes[f'{c_uptake_id}_ub'] = validation_fluxes[c_uptake_id]
    validation_data = ValidationData(valid_data = validation_fluxes,
                                     id = 'EX_co2_e',
                                     validation_range = [
                                         validation_fluxes[c_uptake_id].min()-1, 0
                                     ])

    validation_data._reactions_to_plot = ['Growth', 'EX_co2_e', 'EX_h2_e', 'EX_o2_e']
    validation_data._reactions_to_validate = ['Growth', 'EX_co2_e', 'EX_h2_e', 'EX_o2_e']
    return [validation_data]

def set_up_hyperparameter(processes: int,
                          gene_flow_events:int,
                          filename_extension:str,
                          num_kcats_to_mutate:int = 20,
                          threshold_iteration:int = 10):
    hyperparams = HyperParameters
    hyperparams.threshold_iteration = threshold_iteration
    hyperparams.number_of_kcats_to_mutate = num_kcats_to_mutate
    hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_iABA974'
    hyperparams.filename_extension = filename_extension
    hyperparams.genetic_algorithm_filename_base += filename_extension

    hyperparams.genetic_algorithm_hyperparams['time_limit'] = 60000
    hyperparams.genetic_algorithm_hyperparams['processes'] = processes
    hyperparams.genetic_algorithm_hyperparams['number_gene_flow_events'] = gene_flow_events
    hyperparams.genetic_algorithm_hyperparams['number_generations'] = 5
    hyperparams.genetic_algorithm_filename_base = 'genetic_algorithm_run_iABA974_'
    hyperparams.genetic_algorithm_hyperparams['print_progress'] = True
    return hyperparams


def set_up_pamparametrizer(min_substrate_uptake_rate:float, max_substrate_uptake_rate: float,
                           pam_info_file: str = os.path.join(
                                             'Results','1_preprocessing',
                                             'proteinAllocationModel_iABA974_EnzymaticData_250523.xlsx'),
                           processes: int =4,
                           gene_flow_events: int = 4,
                           filename_extension:str = 'iABA974',
                           num_kcats_to_mutate: int =20,
                           threshold_iteration:int =10,
                           kcat_increase_factor: int = 1):
    pam_info_file_path_out = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'Results', '2_parametrization',
                                 'proteinAllocationModel_iABA974_EnzymaticData_multi.xlsx'))
    increase_kcats_in_parameter_file(kcat_increase_factor,
                                         pam_info_file_path_ori=pam_info_file,
                                         pam_info_file_path_out=pam_info_file_path_out)

    pputida_pam = setup_sof1_pam(pam_info_file)
    #close off all exchanges not necessary for medium, as pputida doesn't exchange anything
    for exchange in pputida_pam.exchanges:
        if exchange.id not in pputida_pam.medium:
            pputida_pam.change_reaction_bounds(exchange.id, 0,0)

    pputida_pam.GLUCOSE_EXCHANGE_RXNID = 'EX_glc__D_e'


    validation_data = set_up_validation_data()
    hyperparameters = set_up_hyperparameter(processes,
                                            gene_flow_events,
                                            filename_extension,
                                            num_kcats_to_mutate,
                                            threshold_iteration)

    sector_configs = set_up_sector_config(
        pam_info_file = pam_info_file_path_out,
        sectors_not_related_to_growth = ['TranslationalProteinSector','UnusedEnzymeSector']
    )

    return PAMParametrizer(pamodel=pputida_pam,
                           validation_data=validation_data,
                           hyperparameters=hyperparameters,
                           sector_configs = sector_configs,
                           substrate_uptake_id='EX_co2_e',
                           max_substrate_uptake_rate=max_substrate_uptake_rate,
                           min_substrate_uptake_rate=min_substrate_uptake_rate
                           )

def run_parametrizations(n_iterations:int=5,
                         pam_info_file: str = os.path.join(
                                             'Results', '1_preprocessing',
                                             'proteinAllocationModel_iABA974_EnzymaticData_250523.xlsx')
                         ) -> None:
    for i in range(1, n_iterations+1):
        print('Working on iteration number', i, 'out of ',n_iterations)
        print('------------------------------------------------------------------------------------------------')
        pam_parametrizer = set_up_pamparametrizer(MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE,
                                                  pam_info_file = pam_info_file,
                                                  filename_extension = f'iABA974_{i}',
                                                  kcat_increase_factor= 1)
        #
        pam_parametrizer.run(remove_subruns=True, binned='False')

if __name__ == "__main__":
    # pam_parametrizer = set_up_pamparametrizer(MIN_SUBSTRATE_UPTAKE_RATE, MAX_SUBSTRATE_UPTAKE_RATE,
    #                      c_sources = ['Glycerol', 'Glucose', 'Succinate', 'Fructose','m-Xylene','Toluene','Benzoate', 'Octanoate'])
    # #
    # pam_parametrizer.run(remove_subruns=True, binned = 'False')
    #pam_info_file = os.path.join('Results', '1_preprocessing',
                                 #'proteinAllocationModel_iABA974_EnzymaticData_250523.xlsx')
    pam_info_file = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'Results', '1_preprocessing',
                                                 'proteinAllocationModel_iABA974_EnzymaticData_250523.xlsx'))
    if len(sys.argv)>1:
        pam_info_file = sys.argv[1]
    # set_up_validation_data(pam_info_file=pam_info_file)
    run_parametrizations(pam_info_file=pam_info_file)
#for running:
#python -m Scripts.i2_parametrization.pam_parametrizer_iABA974