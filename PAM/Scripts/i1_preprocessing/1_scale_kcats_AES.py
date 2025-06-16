import pandas as pd
import numpy as np
import os
from typing import Union, Callable

from PAModelpy import PAModel
import matplotlib.pyplot as plt

from Modules.utils.pam_generation import setup_pputida_pam, setup_yeast_pam, setup_cglutanicum_pam
from PAModelpy.utils import set_up_pam
from Scripts.i2_parametrization.pam_parametrizer_iML1515 import set_up_pamparametrizer
from Scripts.i2_parametrization.pam_parametrizer_iJN1463 import set_up_pamparametrizer as set_up_pamparam_pputida
from Scripts.i2_parametrization.pam_parametrizer_iCGB21FR import set_up_pamparametrizer as set_up_pamparam_cglutanicum
# from Scripts.i2_parametrization.pam_parametrizer_yeast9 import set_up_pamparametrizer as set_up_pamparam_yeast9




VALID_DATA_PATH = os.path.join('Data', 'Ecoli_phenotypes', 'Ecoli_phenotypes_py_rev.xls')
valid_data_df = pd.read_excel(VALID_DATA_PATH, sheet_name='Yields').sort_values('BIOMASS_Ec_iML1515_core_75p37M',
                                                                                    ascending = False)
max_mu = valid_data_df.at[0,'BIOMASS_Ec_iML1515_core_75p37M']

def change_kcat_with_factor(pam:PAModel,
                            multiplication_factor: int) -> None:
    for enzyme in pam.enzymes:
        kcats = enzyme.rxn2kcat.copy()
        for rxn, kcat_dict in kcats.items():
            # if the enzyme is part of a complex it has only zero kcats
            if all([val == 0 for val in kcat_dict.values()]):
                continue
            for dir, kcat in kcat_dict.items():
                kcats[rxn][dir] = kcat * multiplication_factor
        pam.change_kcat_value(enzyme.id, kcats)

def scan_kcat_factors(max_factor:int,
                      scan_figure_file_path: str,
                      min_factor:int=1,
                      stepsize:Union[float, int] = 1,
                      setup_pam_function: Callable = set_up_pam,
                      pam_info_file: str = os.path.join(
                          'Results', '1_preprocessing', 'proteinAllocationModel_iML1515_EnzymaticData_250225.xlsx'
                      ),
                      setup_pamparametrizer_function: Callable = set_up_pamparametrizer,
                      pam_parametrizer_kwargs: dict = {'c_sources': ['Glucose']},
                      substrate_reaction_id = 'EX_glc__D_e',
                      substrate_uptake_rates: np.arange = np.arange(-11,1,1)) -> None:

    factors_to_scan = np.arange(min_factor,max_factor+stepsize,stepsize)

    pam_parametrizer_kwargs['threshold_iteration'] = len(factors_to_scan)
    pam_parametrizer = setup_pamparametrizer_function(min_substrate_uptake_rate = substrate_uptake_rates[0],
                                                      max_substrate_uptake_rate =substrate_uptake_rates[-1],
                                                      **pam_parametrizer_kwargs)
    fig, axs = pam_parametrizer.plot_valid_data()

    for mult_factor in factors_to_scan:
        print("\n----------------------------------------------------------------------")
        print(f"Multiplying the kcat values with {mult_factor}")
        pam_parametrizer.iteration += 1

        pam = setup_pam_function(pam_info_file,
                                 sensitivity=False)
        pam.change_reaction_bounds(substrate_reaction_id, -1e3, 0)
        change_kcat_with_factor(pam, mult_factor)

        pam_parametrizer.pamodel= pam
        fig = pam_parametrizer.plot_simulation(fig,axs, cbar_label='Scaling factor')
        pam.change_reaction_bounds(substrate_reaction_id, -1e3, 0)
        pam.optimize()
        print(f'The model with these parameters is *{pam.solver.status}* and the maximal possible growth rate is {pam.objective.value} h-1') #TODO why is it not giving output?

    plt.savefig(scan_figure_file_path)
    print("\n-------------------------------------------------------------------\nDone scanning the kcat multiplication factors")

def scan_kcat_factors_iML1515():
    scan_kcat_factors(10, os.path.join('Results', '1_preprocessing','figures', 'multifactor_scan_iML1515.png'))


def scan_kcat_factors_pputida():
    scan_kcat_factors(max_factor=10,
                      min_factor=1,
                      stepsize=1,
                      scan_figure_file_path= os.path.join('Results','1_preprocessing','figures','multifactor_scan_iJN1463.png'),
                      setup_pam_function=setup_pputida_pam,
                      pam_info_file=os.path.join(
                          'Results', '1_preprocessing', 'proteinAllocationModel_iJN1463_EnzymaticData_250225.xlsx'
                      ),
                      setup_pamparametrizer_function=set_up_pamparam_pputida,
                      substrate_reaction_id= 'EX_glc__D_e',
                      substrate_uptake_rates= np.arange(-15,1,1))

def scan_kcat_factors_cglutanicum():
    scan_kcat_factors(max_factor=10,
                      min_factor=1,
                      stepsize=1,
                      scan_figure_file_path= os.path.join('Results','1_preprocessing','figures','multifactor_scan_iCGB21FR.png'),
                      setup_pam_function=setup_cglutanicum_pam,
                      pam_info_file=os.path.join(
                          'Results', '1_preprocessing','proteinAllocationModel_iCGB21FR_EnzymaticData_250227.xlsx'
                      ),
                      setup_pamparametrizer_function=set_up_pamparam_cglutanicum,
                      substrate_reaction_id= 'EX_glc__D_e',
                      substrate_uptake_rates= np.arange(-20,1,1))

# def scan_kcat_factors_yeast9():
#     scan_kcat_factors(max_factor=10,
#                       min_factor=1,
#                       stepsize=1,
#                       scan_figure_file_path= os.path.join('Results','1_preprocessing','figures','multifactor_scan_yeast9.png'),
#                       setup_pam_function=setup_yeast_pam,
#                       pam_info_file=os.path.join(
#                           'Results', '1_preprocessing','proteinAllocationModel_yeast9_EnzymaticData_TurnUp.xlsx'
#                       ),
#                       setup_pamparametrizer_function=set_up_pamparam_yeast9,
#                       substrate_reaction_id= 'r_1714',
#                       substrate_uptake_rates= np.arange(-20,1,1))

if __name__ == '__main__':
    # scan_kcat_factors_iML1515()
    scan_kcat_factors_pputida()
    # scan_kcat_factors_cglutanicum()
    # scan_kcat_factors_yeast9()