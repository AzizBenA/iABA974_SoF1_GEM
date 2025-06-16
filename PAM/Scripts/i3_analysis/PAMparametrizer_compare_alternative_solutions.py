import pandas as pd
import os
import matplotlib.pyplot as plt

from typing import Tuple, Callable

from PAModelpy.utils.pam_generation import set_up_pam
from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer
# from Scripts.i2_parametrization.mcpam import set_up_pamparametrizer as set_up_pamparametrizer_mcpam

from Modules.utils.pam_generation import create_pamodel_from_diagnostics_file
from Modules.utils.pamparametrizer_analysis import set_up_pam_parametrizer_and_get_substrate_uptake_rates
from Modules.utils.pamparametrizer_visualization import plot_valid_data, plot_simulation



def recreate_progress_plot(best_indiv_files:list[str],
                           labels:list[str],
                           fig_file_path:str,
                           legend:bool = True
                           ):
    FIGWIDTH = 12
    FIGHEIGHT = 12
    FONTSIZE = 20

    j=0

    parametrizer, substrate_rates = set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates()

    substrate_rates = sorted(substrate_rates)
    fig, axs = plot_valid_data(parametrizer, fontsize=FONTSIZE)
    print('Run reference simulations')
    # fluxes = run_simulations(pamodel, substrate_rates)
    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                   substrate_rates = substrate_rates,
                                                                   sensitivity = False)
    fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                               parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                               iteration=0, color='black')

    for file, label in zip(best_indiv_files, labels):
        j +=1
        print('\nAlternative ', label, ' from file ', file)
        parametrizer.pamodel = create_pamodel_from_diagnostics_file(file,
                                                                    parametrizer._pamodel.copy(copy_with_pickle = True))
        fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                       substrate_rates=substrate_rates,
                                                                       sensitivity=False)
        fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                                   parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                                   iteration=j + 1, max_iteration=len(best_indiv_files), label = label)


    lines, labels = fig.axes[1].get_legend_handles_labels()

    fig.set_figwidth(FIGWIDTH)
    fig.set_figheight(FIGHEIGHT)
    if legend:
        fig.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.1, 0.99), frameon=False,
                   fontsize=FONTSIZE - 5)
    fig.tight_layout()
    plt.savefig(fig_file_path)
    plt.close(fig)

def visualize_progress_plot(alternative_param_files:list[str],
                           labels:list[str],
                           fig_file_path:str,
                            set_up_parametrizer: Callable = set_up_pamparametrizer
                           ):
    FIGWIDTH = 12
    FIGHEIGHT = 12
    FONTSIZE = 20

    j=0

    parametrizer, substrate_rates = set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates(set_up_parametrizer)
    fig, axs = plot_valid_data(parametrizer, fontsize=FONTSIZE)
    pam = parametrizer._pamodel.copy(copy_with_pickle = True)
    print('Run reference simulations')
    # fluxes = run_simulations(pamodel, substrate_rates)
    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                   substrate_rates = substrate_rates,
                                                                   sensitivity = False)
    fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                               parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                               iteration=0, color='black')

    for file, label in zip(alternative_param_files, labels):
        j +=1
        print('\nAlternative ', label, ' from file ', file)

        parametrizer.pamodel = pam.copy(copy_with_pickle=True)

        fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                                       substrate_rates=substrate_rates,
                                                                       sensitivity=False)
        fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                                   parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                                   iteration=j + 1, max_iteration=len(alternative_param_files), label = label)


    lines, labels = fig.axes[1].get_legend_handles_labels()

    fig.set_figwidth(FIGWIDTH)
    fig.set_figheight(FIGHEIGHT)
    fig.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.1, 0.99), frameon=False,
               fontsize=FONTSIZE - 5)
    fig.tight_layout()
    plt.savefig(fig_file_path)
    plt.close(fig)

def set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates() -> Tuple:
    kwargs = {'min_substrate_uptake_rate':-12,
              'max_substrate_uptake_rate': -0.1,
              'kcat_increase_factor': 3}
    return set_up_pam_parametrizer_and_get_substrate_uptake_rates(set_up_pamparametrizer,
                                                           kwargs)

    return parametrizer, substrate_rates

def main_ecoli():
    NUM_MODELS = 8
    PAM_KCAT_FILES = [os.path.join('Results', '2_parametrization', 'diagnostics',
                                   f'pam_parametrizer_diagnostics_{file_nmbr}.xlsx') for file_nmbr in
                      range(1, NUM_MODELS+1)]
    labels = [f'Alternative {i}' for i in range(1, NUM_MODELS+1)]
    fig_file_path = os.path.join('Results', '3_analysis', 'pam-parametrizer_alternatives_cleaned_figure.png')

    alternative_param_files = PAM_KCAT_FILES
    recreate_progress_plot(alternative_param_files, labels, fig_file_path, legend=True)

# def main_mcecoli():
#     PAM_DATA_FILE_PATH = os.path.join('Results', '2_parametrization','proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx')
#     # PAM_KCAT_FILES = [os.path.join('Results', '3_analysis', 'parameter_files',
#     #                                f'proteinAllocationModel_EnzymaticData_iML1515_{file_nmbr}.xlsx') for file_nmbr in
#     #                   range(1, 11)]
#     PAM_KCAT_FILES = [os.path.join('Results', '2_parametrization', 'diagnostics',
#                                    f'pam_parametrizer_diagnostics_mcML1515_1.xlsx')]
#     labels = [f'Alternative 1' ]
#     fig_file_path = os.path.join('Results', '3_analysis', 'pam-parametrizer_alternatives_cleaned_figure_mcpam.png')
#
#     alternative_param_files = PAM_KCAT_FILES
#     recreate_progress_plot(alternative_param_files,
#                            labels, fig_file_path,
#                            legend=True,
#                            setup_parametrizer = set_up_pamparametrizer_mcpam)

if __name__ == '__main__':
    main_ecoli()
