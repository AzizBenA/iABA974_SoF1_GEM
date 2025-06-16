import math
import os
import pandas as pd
from typing import Callable
from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer

from PAModelpy.utils import set_up_pam

MIN_SUBSTRATE_UPTAKE_RATE = -11
MAX_SUBSTRATE_UPTAKE_RATE = -0.1
MAX_MU_ECOLI = 0.7 #estimate
MAX_MU_ECOLI_ALE = math.log(2)*1.63#from lenski experiment: ln(2)*1.63(doublings per hour)

UE_0 = 0.111925# intercept for the unused protein to mu relation (unused enzymes at zero growth) in g/gDW

max_mu_for_UEmu_determination =[MAX_MU_ECOLI, #max growth
                                MAX_MU_ECOLI+(MAX_MU_ECOLI_ALE-MAX_MU_ECOLI)/2, #intermediate timepoint
                                MAX_MU_ECOLI_ALE, #max growth after ALE
                                MAX_MU_ECOLI_ALE+(MAX_MU_ECOLI_ALE-MAX_MU_ECOLI)/2] #to check if there is still some unused protein left after ALE
pam_info_file = os.path.join(
        'Results', '1_preprocessing', 'proteinAllocationModel_iML1515_EnzymaticData_250423.xlsx')

def scan_unused_enzyme_sector_max_growth_rates(max_mu_list:list,
                                               intercept_unused_enzymes:float,
                                               kcat_increase_factor:int,
                                               result_df_file_path: str,
                                               num_replicate_simulations:int = 3,
                                               setup_pamparametrizer_function: Callable = set_up_pamparametrizer,
                                               pam_parametrizer_kwargs: dict = {'c_sources': ['Glucose'],
                                                                                'threshold_iteration': 5,
                                                                                'processes':2,
                                                                                'gene_flow_events': 2},
                                               ):
    result_df = pd.DataFrame(columns=['maxmu', 'rsquared', 'iteration'])
    for i, max_mu in enumerate(max_mu_list):
        slope = -intercept_unused_enzymes / max_mu  # dy/dx (y is protein concentration, x is growth rate)
        print('\n\n----------------------------------------------------------------------')
        print(f'Starting optimization using {max_mu} h-1 as growth rate without unused enzymes')
        print(
            f'The unused protein sector represents the following relation:\n\t\t P_UE[g/gDW] = {UE_0}[g/gDW] - {abs(slope)}[g/gDW/h] * mu[h-1]')
        for iteration in range(num_replicate_simulations):
            print(f'\n===========\nStarting iteration number {iteration}')
            pam_parametrizer = setup_pamparametrizer_function(MIN_SUBSTRATE_UPTAKE_RATE,
                                                              MAX_SUBSTRATE_UPTAKE_RATE,
                                                              pam_info_file = pam_info_file,
                                                              kcat_increase_factor=kcat_increase_factor,
                                                              filename_extension=f'iML1515_UE_{i}_{iteration}', # change the filename to save the results in a recognizable name
                                                              **pam_parametrizer_kwargs)
            pam_parametrizer.sector_configs['UnusedEnzymeSector']['slope'] = slope
            #reset sector configuration
            pam_parametrizer.calculate_sector_parameters_for_multiple_csources(reset=True)

            pam_parametrizer.run(remove_subruns=True, binned='False')

            result_df.loc[len(result_df)] = [max_mu, pam_parametrizer.final_error, iteration]

            reset_pam_parametrizer(pam_parametrizer)

    result_df.to_excel(result_df_file_path)
    print("Ended the scan of unused enzyme sector parameters")


def reset_pam_parametrizer(parametrizer):
    # need to reset best individual and computational performance df
    parametrizer.parametrization_results.best_individuals = pd.DataFrame(
        columns=['run_id', 'enzyme_id', 'direction', 'rxn_id', 'kcat[s-1]', 'ga_error'])
    parametrizer.parametrization_results.computational_time = pd.DataFrame(columns=['run_id', 'time_s', 'time_h'])

    # reset final errors for correct saving
    parametrizer.parametrization_results.final_errors = pd.DataFrame(columns=['run_id', 'r_squared'])
    ecoli_pam = set_up_pam(pam_info_file = pam_info_file)
    parametrizer._pamodel = ecoli_pam

if __name__ == "__main__":
    #example usage iML1515
    scan_unused_enzyme_sector_max_growth_rates(max_mu_for_UEmu_determination,
                                               intercept_unused_enzymes= UE_0,
                                               kcat_increase_factor=3,
                                               result_df_file_path= os.path.join('Results',
                                                                                 '1_preprocessing',
                                                                                 'iML1515_UES_intercept_screen.xlsx'))