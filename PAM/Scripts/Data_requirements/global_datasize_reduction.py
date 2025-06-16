import pandas as pd
import numpy as np
import os
from datetime import date

from Scripts.i2_parametrization.pam_parametrizer_performance_analysis import get_statistics_from_df, save_pam_parametrizer_results_to_df

from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer
RESULT_FILE_PATH = os.path.join('Results', f'pam_parametrizer_statistics_datasize_reduction_{date.today()}.xlsx')


def run_parametrization_workflow(iteration, iterations,
                                 processes, frac_data,
                                 gene_flow_events, num_kcats_to_mutate,
                                 best_individual_df, computational_performance_df,
                                 param_max_iteration: int=5,
                                 min_substrate_uptake = -11, max_substrate_uptake = -0.1,
                                 kcat_increase_factor =1):
    print('\n\n###################################################################################')
    print('starting with iteration number ', iteration, ' out of ', iterations, ' iterations with ',
          frac_data*1e2, '% of the available datapoints\n')
    print('------------------------------------------------------------------------------------------------')
    percentage_data_str = str(int(round(frac_data,2)*100))
    parametrizer = set_up_pamparametrizer(min_substrate_uptake, max_substrate_uptake, processes=processes,
                                          threshold_iteration = param_max_iteration,
                                          gene_flow_events=gene_flow_events,
                                          filename_extension= f'datareduc_{percentage_data_str}_{iteration}',
                                          num_kcats_to_mutate = num_kcats_to_mutate,
                                          kcat_increase_factor = kcat_increase_factor)

    nmrb_rows_to_sample = int(frac_data*len(parametrizer.validation_data.get_by_id('EX_glc__D_e').valid_data))

    parametrizer.validation_data.get_by_id('EX_glc__D_e').valid_data = parametrizer.validation_data.get_by_id('EX_glc__D_e').valid_data.sample(nmrb_rows_to_sample)

    #need to reset best individual, computational performance, and final errors df
    parametrizer.parametrization_results.best_individuals = pd.DataFrame(columns=['run_id', 'enzyme_id', 'direction', 'rxn_id', 'kcat[s-1]', 'ga_error'])
    parametrizer.parametrization_results.computational_time = pd.DataFrame(columns=['run_id', 'time_s', 'time_h'])
    parametrizer.parametrization_results.final_errors = pd.DataFrame(columns=['run_id', 'r_squared'])

    parametrizer.run(binned='False')

    results_file_path = os.path.join('Results','2_parametrization','diagnostics',
                                     f'pam_parametrizer_diagnostics_datareduc_{percentage_data_str}_{iteration}.xlsx')
    best_individual_df, computational_performance_df = save_pam_parametrizer_results_to_df(iteration,str(nmrb_rows_to_sample),
                                                   best_individual_df,computational_performance_df,
                                                   results_file_path)
    return best_individual_df, computational_performance_df

def analyse_parametrizer_performance():
    min_substrate = -11
    max_substrate = -0.1
    iterations = 4
    processes = 4
    kcat_increase_factor = 3
    gene_flow_events = processes

    # 0. Initialize result dataframes
    best_individual_df = pd.DataFrame(columns=['iteration', 'binned',
                                               'run_id', 'enzyme_id', 'rxn_id', 'kcat[s-1]',
                                               'ga_error', 'r_squared'])
    computational_performance_df = pd.DataFrame(columns=['iteration', 'binned', 'run_id', 'time_s', 'time_h'])

    # 1. Run different configuations of the toy model parametrization for n iterations and save the results
    for frac_data in np.arange(0.1,1,0.1):
        for iteration in range(iterations):
            best_individual_df, computational_performance_df = run_parametrization_workflow(
                iteration+1, iterations, processes,frac_data, gene_flow_events, 10,
                best_individual_df, computational_performance_df,
                min_substrate_uptake=min_substrate, max_substrate_uptake=max_substrate,
                kcat_increase_factor=kcat_increase_factor
            )

    # 2. Calculate mean error and stdev for each method and for a single iteration
    error_per_iteration_config = get_statistics_from_df(best_individual_df,
                                                        group_by=['iteration', 'binned'],
                                                        columns=['ga_error', 'r_squared'])
    error_per_config = get_statistics_from_df(best_individual_df,
                                              group_by=['binned'],
                                              columns=['ga_error', 'r_squared'])
    computational_performance_per_config = get_statistics_from_df(computational_performance_df,
                                                                  group_by=['binned'],
                                                                  columns=['time_s'])
    # 3. save the results to excel
    with pd.ExcelWriter(RESULT_FILE_PATH) as writer:
        # Write each DataFrame to a specific sheet
        best_individual_df.to_excel(writer, sheet_name='Best_Individuals', index=False)
        computational_performance_df.to_excel(writer, sheet_name='Computational_Time',
                                              index=False)
        error_per_config.to_excel(writer, sheet_name='stats_error_per_config', index=False)
        error_per_iteration_config.to_excel(writer, sheet_name='stats_error_per_iteration', index=False)
        computational_performance_per_config.to_excel(writer, sheet_name='stats_computation_time', index=False)



if __name__ == '__main__':
    analyse_parametrizer_performance()