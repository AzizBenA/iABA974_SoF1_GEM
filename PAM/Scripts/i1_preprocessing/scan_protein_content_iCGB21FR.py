import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Modules.utils.sector_config_functions import perform_linear_regression

from Scripts.i2_parametrization.pam_parametrizer_iCGB21FR import set_up_pamparametrizer, MAX_SUBSTRATE_UPTAKE_RATE
from Modules.utils.pamparametrizer_analysis import set_up_pam_parametrizer_and_get_substrate_uptake_rates
from Modules.utils.pamparametrizer_visualization import plot_valid_data, plot_simulation

FIGWIDTH = 12
FIGHEIGHT = 12
FONTSIZE = 20
F_ACTIVE_PROTEINS = 0.5#g_p/g_CDW
MAX_GROWTH_ALE = 0.94#1/h; EVO5 strain from Graf et al 2019, growth rate obtained from Matamouros et al (2023)
UNUSED_PROTEIN_INTERCEPT = 0.37 #g_unusedprotein/g_protein; 37% (Bruggeman et al (2020)) of all the proteins which can be measured in E.coli

def main():
    DATA_FILE_PATH = os.path.join('Data', 'Cglutanicum_phenotypes',
                                  'cglutanicum_proteomics.xlsx') # measurements to setup translational_protein_sector
    df_proteomics = pd.read_excel(DATA_FILE_PATH, sheet_name='ribosome_fraction')

    total_protein_range = np.arange(0.51,0.7,0.01) #g_p/g_CDW based on experimental quatitative proteomics measurements
    parametrizer_kwargs = {
        'max_substrate_uptake_rate':MAX_SUBSTRATE_UPTAKE_RATE,
        'pam_info_file': os.path.join(
            'Results', '1_preprocessing','proteinAllocationModel_iCGB21FR_EnzymaticData_250221.xlsx'
        )
    }
    parametrizer, substrate_rates = set_up_pam_parametrizer_and_get_substrate_uptake_rates(set_up_pamparametrizer,
                                                                                           parametrizer_kwargs)
    fig, axs = plot_valid_data(parametrizer, fontsize=FONTSIZE)
    print('Run reference simulations')
    # fluxes = run_simulations(pamodel, substrate_rates)
    fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                     substrate_rates=substrate_rates,
                                                     sensitivity=False)

    fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                               parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                               iteration=0, color='black', label = 'original')

    for i, total_protein in enumerate(total_protein_range):
        metabolic_protein_fraction = total_protein*F_ACTIVE_PROTEINS

        #translational_sector
        slope, intercept = perform_linear_regression(x=df_proteomics['growth_rate'][df_proteomics['growth_rate'] < 0.4],
                                                     y=df_proteomics['ribosomal_protein_fraction'][
                                                           df_proteomics['growth_rate'] < 0.4] * total_protein)
        trans_sector = parametrizer.pamodel_no_sensitivity.sectors.get_by_id('TranslationalProteinSector')

        parametrizer.pamodel_no_sensitivity.change_sector_parameters(trans_sector, slope,
                                         intercept , lin_rxn_id='Growth')
        #unused enzymes sector
        trans_sector = parametrizer.pamodel_no_sensitivity.sectors.get_by_id('UnusedEnzymeSector')
        parametrizer.pamodel_no_sensitivity.change_sector_parameters(
            trans_sector,
            slope=UNUSED_PROTEIN_INTERCEPT * metabolic_protein_fraction / MAX_GROWTH_ALE,
            intercept=UNUSED_PROTEIN_INTERCEPT * metabolic_protein_fraction,
            lin_rxn_id='Growth'
        )

        #total protein
        parametrizer.pamodel_no_sensitivity.change_total_protein_constraint(metabolic_protein_fraction)

        #run simulations
        fluxes, _ = parametrizer.run_simulations_to_plot(substrate_uptake_id='EX_glc__D_e',
                                                         substrate_rates=substrate_rates,
                                                         sensitivity=False)

        #visualize
        fig, axs = plot_simulation(fig, axs, fluxes, [abs(rate) for rate in substrate_rates],
                                   parametrizer.validation_data.get_by_id('EX_glc__D_e')._reactions_to_plot,
                                   iteration=i, max_iteration=len(total_protein_range), label = round(total_protein,2))


    lines, labels = fig.axes[1].get_legend_handles_labels()

    fig.set_figwidth(FIGWIDTH)
    fig.set_figheight(FIGHEIGHT)
    fig.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.1, 0.99), frameon=False,
                   fontsize=FONTSIZE - 5)
    fig.tight_layout()
    plt.savefig(os.path.join('Results', '1_preprocessing', 'protein_content_scan_cglutanicum.png'))
    plt.close(fig)


if __name__ == '__main__':
    main()

