import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import matplotlib.ticker as ticker
from typing import Callable

from Modules.utils.pamparametrizer_visualization import plot_simulation, plot_valid_data
from Modules.utils.pamparametrizer_analysis import set_up_pam_parametrizer_and_get_substrate_uptake_rates
from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer as set_up_pamparametrizer_ecoli
from Scripts.i2_parametrization.pam_parametrizer_iJN1463 import set_up_pamparametrizer as set_up_pamparametrizer_putida
# from Scripts.i2_parametrization.mcpam import set_up_pamparametrizer as set_up_pamparametrizer_mcpam



FONTSIZE = 16

def recreate_progress_plot(best_individual_df:pd.DataFrame,
                           fig_file_path:str,
                           return_error_df = False,
                           set_up_parametrizer: Callable = None,
                           pamparam_kwargs: dict = {'max_substrate_uptake_rate': -0.1,
                                                    'min_substrate_uptake_rate': -11,
                                                    'kcat_increase_factor': 3}
                           ) -> None:
    FIGWIDTH = 12
    FIGHEIGHT = 12
    FONTSIZE = 20

    if set_up_parametrizer is None:
        set_up_parametrizer = set_up_pamparametrizer_ecoli

    error_df = pd.DataFrame(columns=['run_id', 'error'])
    parametrizer, substrate_rates = set_up_pam_parametrizer_and_get_substrate_uptake_rates(
        set_up_parametrizer=set_up_parametrizer,
        parametrizer_kwargs = pamparam_kwargs
    )


    fig, axs = plot_valid_data(parametrizer, fontsize=FONTSIZE)
    print('Run reference simulations')
    # fluxes = run_simulations(pamodel, substrate_rates)
    fluxes, substrates = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                   substrate_rates = substrate_rates,
                                                                   sensitivity = False)

    fig, axs = plot_simulation(fig, axs, fluxes,
                               substrates.keys(),
                               parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                               iteration=0, color='black')

    groups = best_individual_df.groupby('run_id')
    for j, group in groups:
        print('\nRun number ', j + 1)
        for i, row in group.iterrows():
            kcat_dict = {row['rxn_id']: {row['direction']: row['kcat[s-1]']}}
            parametrizer.pamodel_no_sensitivity.change_kcat_value(enzyme_id=row['enzyme_id'], kcats=kcat_dict)
        # fluxes = run_simulations(pamodel, substrate_rates)
        fluxes, substrates = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                       substrate_rates=substrate_rates,
                                                                       sensitivity=False)
        fig, axs = plot_simulation(fig, axs, fluxes, substrates.keys(),
                                   parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                                   iteration=j + 1, max_iteration=len(groups))
        for flux, rate in zip(fluxes, substrates.keys()):
            parametrizer.parametrization_results.add_fluxes_from_fluxdict(flux_dict=flux,
                                                                          bin_id='final',
                                                                          substrate_reaction_id= parametrizer.substrate_uptake_id,
                                                                          substrate_uptake_rate=rate,
                                                                          fluxes_abs=False)
        error = parametrizer.calculate_final_error()
        error_df.loc[len(error_df)] = [j, error]


    lines, labels = fig.axes[1].get_legend_handles_labels()

    fig.set_figwidth(FIGWIDTH)
    fig.set_figheight(FIGHEIGHT)
    fig.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.1, 0.99), frameon=False,
               fontsize=FONTSIZE - 5)
    fig.tight_layout()
    plt.savefig(fig_file_path)
    plt.close(fig)
    if return_error_df: return error_df

def create_empty_plot(plot_reference_simulations:bool = True):
    parametrizer =  set_up_pamparametrizer_ecoli(-12,-0.1,
                                                 pam_info_file = os.path.join(
                                                     'Results','1_preprocessing','proteinAllocationModel_iML1515_EnzymaticData_250423.xlsx'
                                                 ))
    parametrizer._init_results_objects()
    substrate_rates = parametrizer._init_validation_df([parametrizer.min_substrate_uptake_rate,
                                                        parametrizer.max_substrate_uptake_rate])['EX_glc__D_e']

    fig, axs = plot_valid_data(parametrizer, fontsize=FONTSIZE)

    if plot_reference_simulations:
        substrate_rates = sorted(substrate_rates)
        fluxes, substrates = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                  substrate_rates=substrate_rates,
                                                                  sensitivity=False)
        fig, axs = plot_simulation(fig, axs, fluxes,
                                   substrates.keys(),
                                   parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                                   iteration=0, color='black')

    fig.set_size_inches(18.5, 10.5)
    return fig, axs

def main_ecoli():
    result_file = os.path.join('Results', '2_parametrization',
                               'diagnostics', 'pam_parametrizer_diagnostics_2.xlsx')

    best_indiv_df = pd.read_excel(result_file, sheet_name='Best_Individuals')
    #
    fig_file_path = os.path.join('Results', '3_analysis', 'pam_parametrizer_progess_cleaned_iML1515_2.png')
    recreate_progress_plot(best_indiv_df, fig_file_path)

def main_putida():
    result_file = os.path.join('Results', '2_parametrization', 'diagnostics', 'pam_parametrizer_diagnostics_iJN1463_1.xlsx')

    best_indiv_df = pd.read_excel(result_file, sheet_name='Best_Individuals')
    #
    fig_file_path = os.path.join('Results', '3_analysis', 'pam_parametrizer_progess_cleaned_iJN1463_1.png')
    recreate_progress_plot(best_indiv_df,
                           fig_file_path,
                           set_up_parametrizer=set_up_pamparametrizer_putida)

def main_mcecoli():
    result_file = os.path.join('Results', '2_parametrization', 'diagnostics',
                               'pam_parametrizer_diagnostics_mciML1515_1.xlsx')

    best_indiv_df = pd.read_excel(result_file, sheet_name='Best_Individuals')
    #
    fig_file_path = os.path.join('Results', '3_analysis', 'pam_parametrizer_progess_cleaned_mcpam.png')
    recreate_progress_plot(best_indiv_df,
                           fig_file_path,
                           set_up_parametrizer=set_up_pamparametrizer_mcpam)

if __name__ == '__main__':
    create_empty_plot()

