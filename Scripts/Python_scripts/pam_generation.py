import cobra
import pandas as pd
import os
from typing import Union, Optional
import numpy as np

# load PAModelpy modules
from PAModelpy.PAModel import PAModel
from PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from PAModelpy.configuration import Config
from PAModelpy.utils.pam_generation import set_up_pam, parse_reaction2protein, _order_enzyme_complex_id


'Function library for making Protein Allocation Models as described in the publication'

import os
os.chdir("C:/Users/adam-1y9pxi6q95ggpmx/Documents/PhD/SoF1-GEM/")

def setup_sof1_pam(pam_info_file:str= os.path.join(
                                             'Data','PAM_data',
                                             'proteinAllocationModel_EnzymaticData_iABA974_101025.xlsx'),
                     model:str = os.path.join('Data','Models','250924_iABA974_BLG.sbml'),
                     total_protein: Union[bool, float] = None, active_enzymes: bool = True,
                    translational_enzymes: bool = True, unused_enzymes: bool = True, sensitivity = True,
                    gam:Optional[float] = None, ngam:Optional[float] = None):
    config = Config()
    config.reset()
    config.BIOMASS_REACTION = 'Growth'

    if total_protein is None:
        total_protein = 0.67836*0.56 #measured mass fraction *percentage used in metabolism based on E.coli measurements
    #make sure default enzyme ids for locus tags without enzyme name are parsed correctly
    config.ENZYME_ID_REGEX += r'|Enzyme_*|Enzyme_AAFOLC_[0-9]+'
    pam = set_up_pam(pam_info_file, model, config,
                     total_protein, active_enzymes, translational_enzymes,
                     unused_enzymes, sensitivity = sensitivity)
    #open up exchanges to properly study effect of limiting reaction
    for sub_rxn in ['EX_co2_e', 'EX_o2_e', 'EX_h2_e']:
        pam.change_reaction_bounds(sub_rxn, -1e2, 1e2)

    #change ngam and gam
    if ngam is not None:
        pam.change_reaction_bounds('ATPM', ngam)
        print(fr'Changing non-growth associated maintenance to {ngam} mmol_ATP/gCDW/h')
    if gam is not None:
        m_to_coeff = {pam.metabolites.atp_c:-gam,
                  **{pam.metabolites.get_by_id(mid): gam for mid in ['adp_c', 'pi_c', 'h_c']
                     }
                  }
        pam.reactions.Growth.add_metabolites(m_to_coeff)

        print(fr'Changing growth associated maintenance to {gam} mmol_ATP/gCDW')
    return pam



def parse_coefficients(pamodel):
    Ccsc = list()

    for csc in ['flux_ub', 'flux_lb', 'enzyme_max', 'enzyme_min', 'proteome', 'sector']:
        Ccsc += pamodel.capacity_sensitivity_coefficients[
            pamodel.capacity_sensitivity_coefficients['constraint'] == csc].coefficient.to_list()

    Cesc = pamodel.enzyme_sensitivity_coefficients.coefficient.to_list()

    return Ccsc, Cesc

def parse_esc(pamodel):
    return pamodel.enzyme_sensitivity_coefficients.coefficient.to_list()

# if __name__ == '__main__':
#     import numpy as np
#
#     pam = setup_ecolicore_pam()
#     for i in np.arange(0, 12, 1):
#         pam.change_reaction_bounds('EX_glc__D_e', lower_bound=-i, upper_bound=0)
#         pam.optimize()
#         print(pam.objective.value)