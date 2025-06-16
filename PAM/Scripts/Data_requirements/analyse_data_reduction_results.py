from typing import Iterable, Tuple, Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import seaborn as sns

from Scripts.i3_analysis.PAMparametrizer_compare_alternative_solutions import set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates

from Modules.utils.error_calculation import nanaverage, calculate_smape_for_reaction
from Modules.utils.pam_generation import create_pamodel_from_diagnostics_file

NUM_DATAPOINTS = 53

def run_simulations_and_calculate_error(parametrizer, substrate_rates) -> Tuple[dict[str, float], dict[str, float]]:
    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                     substrate_rates=substrate_rates,
                                                     sensitivity=False)
    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                     substrate_rates=substrate_rates,
                                                     sensitivity=False)

    for flux, rate in zip(fluxes, substrate_rates):
        parametrizer.parametrization_results.add_fluxes_from_fluxdict(flux_dict=flux,
                                                                      bin_id='final',
                                                                      substrate_reaction_id=parametrizer.substrate_uptake_id,
                                                                      substrate_uptake_rate=rate,
                                                                      fluxes_abs=False)
    smape = get_smape_for_parametrization_experiment(parametrizer)
    error = get_rsquared_for_parametrization_experiment(parametrizer)
    return error, smape

def get_smape_for_parametrization_experiment(parametrizer) -> dict[str, float]:
    error = {}
    for valid_data in parametrizer.validation_data:
        substrate_uptake_id = valid_data.id
        reactions_to_validate = valid_data._reactions_to_validate
        validation_df = parametrizer._get_validation_data_to_validate(valid_data)

        flux_df = parametrizer.parametrization_results.flux_results.get_by_id(substrate_uptake_id).fluxes_df

        if len(flux_df) == 0:  # means model is infeasible
            return -1

        for rxn in reactions_to_validate: #+ self.validation_data._get_biomass_reactions():
            # only select the rows which are filled with data
            validation_data = validation_df.dropna(axis=0, subset=rxn)
            # if there are no reference data points, continue to the next reaction
            if len(validation_data) == 0:
                error[rxn] = np.nan
                continue

            error[rxn] = calculate_smape_for_reaction(rxn, validation_df, substrate_uptake_id,
                                                               flux_df)

    return error

def get_rsquared_for_parametrization_experiment(parametrizer) -> dict[str, float]:
    error = {}
    for valid_data in parametrizer.validation_data:
        substrate_uptake_id = valid_data.id
        reactions_to_validate = valid_data._reactions_to_validate
        validation_df = parametrizer._get_validation_data_to_validate(valid_data)
        error_list = parametrizer._calculate_error_for_reactions(substrate_uptake_id,
                                                            validation_df,
                                                            reactions_to_validate)
        error = {rxn: error for rxn, error in zip(reactions_to_validate, error_list)}

    return error


def get_error_for_parametrization_experiment(parametrizer,
                                             best_indiv_file_path:str,
                                             substrate_rates: Iterable
                                             ) -> Tuple[dict[str,float], dict[str,float]]:
    pamodel = parametrizer.pamodel.copy(copy_with_pickle=True)
    parametrizer.pamodel = create_pamodel_from_diagnostics_file(best_indiv_file_path, pamodel)
    return run_simulations_and_calculate_error(parametrizer, substrate_rates)

def compute_final_error_on_full_dataset_for_all_experiments(base_file_path: str,
                                                            datasizes: Iterable,
                                                            num_replicates:float) -> pd.DataFrame:
    parametrizer, substrate_rates = set_up_ecoli_pam_parametrizer_and_get_substrate_uptake_rates()
    #ensure all the errors are calculated based pn the same dataset
    parametrizer.validation_data.get_by_id('EX_glc__D_e').sampled_valid_data = parametrizer.validation_data.get_by_id('EX_glc__D_e').valid_data
    smape_columns = [f'smape_{rxn}' for rxn in ['mean']+ parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_validate]
    rsquared_columns = [f'rsquared_{rxn}' for rxn in ['mean']+ parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_validate]
    final_errors = pd.DataFrame(columns = ['perc_data', 'sample']+rsquared_columns+smape_columns)

    for datasize in datasizes:
        print('\n------------------------------------------------------')
        print(f"Analyzing parametrization with {datasize}% of the total amount of data to train")
        for sample in range(1,num_replicates+1):
            print(f'\tReplicate {sample}')
            file_path = f'{base_file_path}{datasize}_{sample}.xlsx'
            if os.path.exists(file_path):
                error, smape = get_error_for_parametrization_experiment(parametrizer, file_path, substrate_rates)

                final_errors = final_errors.append({'perc_data':datasize,
                                     'sample':sample,
                                     'rsquared_mean': nanaverage(list(error.values())),
                                     'smape_mean': nanaverage(list(smape.values())),
                                     **{f'rsquared_{rxn}':r for rxn, r in error.items()},
                                     **{f'smape_{rxn}': sm for rxn, sm in smape.items()}},
                                    ignore_index=True)
                print(final_errors)
    return final_errors

def plot_progression_of_errors(final_errors: pd.DataFrame,
                               metrics: Literal['smape', 'r_squared'] ='rsquared',
                               fig_file_path: str = os.path.join('Results', 'data_reduction_results', 'error_per_reaction_num_datapoints.png'),
                               fontsize:int = 16,
                               ax = None,
                               legend=True
                               ) -> None:
    metrics_mapper = {'rsquared': r'$R^{2}$',
                      'smape': 'Symmetric Mean Absolute Percentage Error'}
    rxn_mapper = {'BIOMASS_Ec_iML1515_core_75p37M': 'Growth rate',
                  'EX_co2_e': 'CO2 excretion',
                  'EX_o2_e': 'O2 uptake',
                  'EX_ac_e': 'Acetate excretion',
                  'mean': 'mean'}
    # Extract all columns that start with "rsquared_" (assuming reaction columns follow this pattern)
    reaction_columns = [col for col in final_errors.columns if col.startswith(f"{metrics}_")]

    model_colors = sns.color_palette("colorblind6", n_colors=len(reaction_columns) - 1)
    cmap = dict(zip(reaction_columns[1:], model_colors))
    cmap[f'{metrics}_mean'] = 'black'

    # Melt the DataFrame to make it long-form for easier plotting
    long_data = final_errors.melt(
        id_vars=["perc_data"],
        value_vars=reaction_columns,
        var_name="reaction",
        value_name=metrics
    )

    # Group by `perc_data` and `Reaction` to calculate the mean, min, and max
    summary_stats = long_data.groupby(["perc_data", "reaction"])[metrics].agg(['mean', 'min', 'max']).reset_index()

    # Plot
    if ax is None:
        plt.figure(figsize=(12, 6))
        plt.tight_layout()

    # Plot each reaction with error bars
    for reaction in summary_stats['reaction'].unique():
        reaction_data = summary_stats[summary_stats["reaction"] == reaction]
        plt.errorbar(
            [perc / 100 * NUM_DATAPOINTS for perc in reaction_data["perc_data"]],
            reaction_data["mean"],
            yerr=[reaction_data["mean"] - reaction_data["min"], reaction_data["max"] - reaction_data["mean"]],
            label=rxn_mapper[reaction.replace(f'{metrics}_', '')],
            capsize=3,
            color=cmap[reaction],
            linewidth=2
        )

    # Customize the plot
    # plt.title("R_squared Values by Reaction", fontsize=16)
    plt.xlabel("Number of datapoints", fontsize=fontsize)
    plt.ylabel(metrics_mapper[metrics], fontsize=fontsize)
    if legend:
        plt.legend(
            loc='upper center',
            bbox_to_anchor=(0.5, -0.2),
            ncol=3,
            fontsize=fontsize,
            frameon=False  # removes the border box
        )

    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    # plt.ylim([40,100])
    #
    plt.savefig(fig_file_path)

if __name__ == '__main__':
    # diagnostic_file_path_base = os.path.join('Results', 'data_reduction_results', 'diagnostics', 'pam_parametrizer_diagnostics_datareduc_')
    # final_errors = compute_final_error_on_full_dataset_for_all_experiments(diagnostic_file_path_base,
    #                                                                        datasizes= np.arange(10,100,10),
    #                                                                        num_replicates=8)
    # final_errors.to_excel(os.path.join('Results', 'data_reduction_results', 'r_squared_for_analysis.xlsx'), index=False)
    final_errors = pd.read_excel(os.path.join('Results', 'data_reduction_results', 'r_squared_for_analysis.xlsx'))
    plot_progression_of_errors(final_errors, metrics = 'smape')

