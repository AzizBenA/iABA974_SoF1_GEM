import cobra
import pandas as pd
import os
from typing import Union
import numpy as np

# load PAModelpy modules
from PAModelpy.PAModel import PAModel
from PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from PAModelpy.configuration import Config
from PAModelpy.utils.pam_generation import set_up_pam, parse_reaction2protein, _order_enzyme_complex_id


from .toy_ec_pam import build_toy_gem, build_active_enzyme_sector, build_translational_protein_sector, build_unused_protein_sector

'Function library for making Protein Allocation Models as described in the publication'



def setup_sof1_pam(pam_info_file:str= os.path.join(
                                             'SoF1','Results',
                                             'proteinAllocationModel_iABA974_EnzymaticData_250523.xlsx'),
                     model:str = r"C:\Users\adam-1y9pxi6q95ggpmx\Documents\Python Scripts\PAM_Parametrization\Models\250521_iABA974_kegg_bjp.sbml",
                     total_protein: Union[bool, float] = 0.3, active_enzymes: bool = True,
                    translational_enzymes: bool = True, unused_enzymes: bool = True, sensitivity = True):
    config = Config()
    config.reset()
    config.BIOMASS_REACTION = 'Growth'
    #make sure default enzyme ids for locus tags without enzyme name are parsed correctly
    config.ENZYME_ID_REGEX += r'|Enzyme_*|Enzyme_AAFOLC_[0-9]+'
    pam = set_up_pam(pam_info_file, model, config,
                     total_protein, active_enzymes, translational_enzymes,
                     unused_enzymes, sensitivity = sensitivity)
    return pam

def setup_toy_pam(sensitivity:bool =True, kcat_fwd:list =[1, 0.5, 1, 0.5 ,0.45, 1.5]) -> PAModel:
    config = Config()
    #setting the configuration for the toy model
    config.BIOMASS_REACTION = 'R7'
    config.GLUCOSE_EXCHANGE_RXNID = 'R1'
    config.CO2_EXHANGE_RXNID = 'R8'
    config.ACETATE_EXCRETION_RXNID = 'R9'
    config.ENZYME_ID_REGEX = r'E([0-9]|[1-9][0-9])'

    Etot = 0.6*1e-3
    model = build_toy_gem()
    active_enzyme = build_active_enzyme_sector(config, kcat_fwd)
    unused_enzyme = build_unused_protein_sector(config)
    unused_enzyme.id = 'UnusedEnzymeSector'
    translation_enzyme = build_translational_protein_sector(config)
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
                      translational_sector=translation_enzyme,
                      unused_sector=unused_enzyme, p_tot=Etot, sensitivity=sensitivity)
    pamodel.objective = 'R7'
    config.reset()
    return pamodel


def setup_ecolicore_pam(total_protein:bool = True, active_enzymes: bool = True,
                        translational_enzymes:bool = True, unused_enzymes:bool = True, sensitivity =True,
                        pam_data_file_path:str = ''):
    # Setting the relative paths
    DATA_DIR = os.path.join(os.getcwd(), 'Data')
    MODEL_DIR = os.path.join(os.getcwd(), 'Models')
    if len(pam_data_file_path)<1:
        PAM_DATA_FILE_PATH = os.path.join(DATA_DIR, 'proteinAllocationModel_iML1515_EnzymaticData_230503.xls')
    else:
        PAM_DATA_FILE_PATH = pam_data_file_path

    # some other constants
    BIOMASS_REACTION = 'BIOMASS_Ecoli_core_w_GAM'
    TOTAL_PROTEIN_CONCENTRATION = 0.16995  # [g_prot/g_cdw]

    config = Config()
    config.reset()

    # load the genome-scale information
    model = cobra.io.load_json_model(os.path.join(MODEL_DIR, 'e_coli_core.json'))

    #load example data for the E.coli iML1515 model
    if active_enzymes:
        # load active enzyme sector information
        enzyme_db = pd.read_excel(PAM_DATA_FILE_PATH, sheet_name='ActiveEnzymes')
        enzyme_db = enzyme_db.set_index('rxnID')
        # correct reaction IDS
        for idx in enzyme_db.index.to_list():
            # transprt reactions<

            if 'pp' in idx:
                idx_new = idx.replace('pp', '')
                if idx_new not in enzyme_db.index:
                    enzyme_db.rename(index={idx: idx_new}, inplace=True)
            if 'ex' in idx:
                idx_new = idx.replace('ex', '')
                if idx_new not in enzyme_db.index:
                    enzyme_db.rename(index={idx: idx_new}, inplace=True)

                    # replace NaN values with unique identifiers
        # replace NaN enzyme ids with a dummy enzyme identifier
        # select the NaN values
        nan_values = enzyme_db['EC_nmbr'].isnull()
        # make a list with unique ids
        nan_ids = [f'E{i}' for i in range(nan_values.sum())]
        # replace nan values by unique id
        enzyme_db.loc[nan_values, 'EC_nmbr'] = nan_ids

        # create enzyme objects for each gene-associated reaction

        # parse the enzyme information (kcat values, identifiers and molmasses)
        kcats_dict = {rxn:kcat for rxn, kcat in enzyme_db.loc[:, 'kcat'].to_dict().items() if rxn in model.reactions}
        ec_dict = {rxn:kcat for rxn, kcat in enzyme_db.loc[:, 'EC_nmbr'].to_dict().items() if rxn in model.reactions}
        molmass_dict = {rxn:kcat for rxn, kcat in enzyme_db.loc[:, 'molMass'].to_dict().items() if rxn in model.reactions}

        kcats = {}
        # save fwd and bckw kcats separately in the form of: {rxnID: {'f': kcat_f, 'b': kcat_b}}
        for rxn, kcat in kcats_dict.items():
            # reversible reaction
            if rxn[-2:] == '_f' or rxn[-2:] == '_b':
                direction = rxn[-1]
                # check if the reaction already exists in the kcat dictionary
                try:
                    kcats[rxn[:-2]][direction] = kcat
                except:
                    kcats[rxn[:-2]] = {direction: kcat}
            # irreversible reaction
            else:
                kcats[rxn] = {'f': kcat}

        rxn2ec = {}
        # parse the enzyme identifiers for the reactions
        for rxn, ec in ec_dict.items():
            if rxn[-2:] == '_f' or rxn[-2:] == '_b':
                rxn = rxn[:-2]
            for enz in str(ec).split(','):
                rxn2ec[rxn] = enz.strip()

        molmass = {}
        # parse the enzyme molmasses for the reactions
        for rxn, mw in molmass_dict.items():
            if rxn[-2:] == '_f' or rxn[-2:] == '_b':
                rxn = rxn[:-2]
            molmass[rxn] = mw

        rxn2protein = {}
        for rxn, ec in rxn2ec.items():
            ec_dict = {**kcats[rxn], **{'molmass': molmass[rxn]}}
            # add enzyme to enzymes related to reaction if these are already stored
            if rxn in rxn2protein.keys():
                rxn2protein[rxn] = {**rxn2protein[rxn], **{ec: ec_dict}}
            # if not create new reaction entry
            else:
                rxn2protein[rxn] = {ec: ec_dict}

        active_enzyme_sector = ActiveEnzymeSector(rxn2protein=rxn2protein, configuration=config)

    else:
        active_enzyme_sector = None

    if translational_enzymes:
        # translational protein sector parameter (substrate dependent)
        id_list_tps = ['EX_glc__D_e']
        tps_0 = [0.04992]  # g/gDW
        tps_mu = [-0.002944]  # g h/gDW -> transformed to match glucose uptake variable
        molmass_tps = [405903.94]  # g/mol

        # translational protein sector
        translation_enzyme_sector = TransEnzymeSector(
            id_list=id_list_tps,
            tps_0=tps_0,
            tps_mu=tps_mu,
            mol_mass=molmass_tps,
            configuration=config
        )
    else:
        translation_enzyme_sector = None

    if unused_enzymes:
        id_list_ups = [BIOMASS_REACTION]
        ups_0 = [0.0407]  # g/gDW
        ups_mu = [-0.0214]  # g h/gDW -> negative relation with growth rate
        molmass_ups = [405903.94]  # g/mol

        unused_enzyme_sector = UnusedEnzymeSector(
            id_list=id_list_ups,
            ups_0=ups_0,
            ups_mu=ups_mu,
            mol_mass=molmass_ups,
            configuration=config,
        )
        unused_enzyme_sector.id = 'UnusedEnzymeSector'
    else:
        unused_enzyme_sector = None

    if total_protein: total_protein = TOTAL_PROTEIN_CONCENTRATION

    pa_model = PAModel(id_or_model=model, p_tot=total_protein, sensitivity=sensitivity, configuration=config,
                       active_sector=active_enzyme_sector, translational_sector=translation_enzyme_sector, unused_sector=unused_enzyme_sector)
    pa_model.BIOMASS_REACTION = BIOMASS_REACTION
    return pa_model

def setup_ecoli_pam(total_protein: Union[bool, float] = True, active_enzymes: bool = True,
                     translational_enzymes: bool = True, unused_enzymes:bool = True, sensitivity = True,
                     pam_data_file_path: str = ''):
    config = Config()
    config.reset()
    # Setting the relative paths
    DATA_DIR = os.path.join(os.getcwd(), 'Data')
    MODEL_DIR = os.path.join(os.getcwd(), 'Models')
    if len(pam_data_file_path) < 1:
        PAM_DATA_FILE_PATH = os.path.join(DATA_DIR, 'proteinAllocationModel_iML1515_EnzymaticData_230503.xls')
    else:
        PAM_DATA_FILE_PATH = pam_data_file_path

    # some other constants
    BIOMASS_REACTION = 'BIOMASS_Ec_iML1515_core_75p37M'
    TOTAL_PROTEIN_CONCENTRATION = 0.258  # [g_prot/g_cdw]

    #setup the gem ecoli iML1515 model
    model = cobra.io.read_sbml_model(os.path.join(MODEL_DIR, 'iML1515.xml'))

    #check if a different total protein concentration is given
    if isinstance(total_protein, float):
        TOTAL_PROTEIN_CONCENTRATION = total_protein

    # load example data for the E.coli iML1515 model
    if active_enzymes:
        # load active enzyme sector information
        active_enzyme_info_old = pd.read_excel(PAM_DATA_FILE_PATH, sheet_name='ActiveEnzymes')

        # replace NaN values with unique identifiers
        # select the NaN values
        nan_values = active_enzyme_info_old['EC_nmbr'].isnull()
        # make a list with unique ids
        nan_ids = [f'E{i}' for i in range(nan_values.sum())]
        # replace nan values by unique id
        active_enzyme_info_old.loc[nan_values, 'EC_nmbr'] = nan_ids

        # parse the enzyme information (kcat values, identifiers and molmasses)
        kcats_dict = active_enzyme_info_old.set_index(keys='rxnID').loc[:, 'kcat'].to_dict()
        ec_dict = active_enzyme_info_old.set_index(keys='rxnID').loc[:, 'EC_nmbr'].to_dict()
        molmass_dict = mol_mass = active_enzyme_info_old.set_index(keys='rxnID').loc[:, 'molMass'].to_dict()

        kcats = {}
        # save fwd and bckw kcats separately in the form of: {rxn_id: {'f': kcat_f, 'b': kcat_b}}
        for rxn, kcat in kcats_dict.items():
            # reversible reaction
            if rxn[-2:] == '_f' or rxn[-2:] == '_b':
                direction = rxn[-1]
                # check if the reaction already exists in the kcat dictionary
                try:
                    kcats[rxn[:-2]][direction] = kcat
                except:
                    kcats[rxn[:-2]] = {direction: kcat}
            # irreversible reaction
            else:
                kcats[rxn] = {'f': kcat}

        rxn2ec = {}
        # parse the enzyme identifiers for the reactions
        for rxn, ec in ec_dict.items():
            if rxn[-2:] == '_f' or rxn[-2:] == '_b':
                rxn = rxn[:-2]
            for enz in str(ec).split(','):
                rxn2ec[rxn] = enz.strip()

        molmass = {}
        # parse the enzyme molmasses for the reactions
        for rxn, mw in molmass_dict.items():
            if rxn[-2:] == '_f' or rxn[-2:] == '_b':
                rxn = rxn[:-2]
            molmass[rxn] = mw

        rxn2protein = {}
        for rxn, ec in rxn2ec.items():
            ec_dict = {**kcats[rxn], **{'molmass': molmass[rxn]}}
            # add enzyme to enzymes related to reaction if these are already stored
            if rxn in rxn2protein.keys():
                rxn2protein[rxn] = {**rxn2protein[rxn], **{ec: ec_dict}}
            # if not create new reaction entry
            else:
                rxn2protein[rxn] = {ec: ec_dict}

        active_enzyme_sector = ActiveEnzymeSector(rxn2protein=rxn2protein)

    else:
        active_enzyme_sector = None

    if translational_enzymes:
        translational_info = pd.read_excel(PAM_DATA_FILE_PATH, sheet_name='Translational')
        translation_enzyme_sector = TransEnzymeSector(
            id_list=[translational_info[translational_info.Parameter == 'id_list'].loc[0, 'Value']],
            tps_0=[translational_info[translational_info.Parameter == 'tps_0'].loc[1, 'Value']],
            tps_mu=[translational_info[translational_info.Parameter == 'tps_mu'].loc[2, 'Value']],
            mol_mass=[translational_info[translational_info.Parameter == 'mol_mass'].loc[3, 'Value']])
    else:
        translation_enzyme_sector = None

    if unused_enzymes:
        unused_protein_info = pd.read_excel(PAM_DATA_FILE_PATH, sheet_name='ExcessEnzymes')
        ups_0 = unused_protein_info[unused_protein_info.Parameter == 'ups_0'].loc[2, 'Value']

        try:
            ups_mu = [unused_protein_info[unused_protein_info.Parameter == 'ups_mu'].loc[1, 'Value']]
        except:
            smax = unused_protein_info[unused_protein_info.Parameter == 's_max_uptake'].loc[1, 'Value']
            ups_mu =[ups_0 / smax]

        unused_protein_sector = UnusedEnzymeSector(
            id_list=[unused_protein_info[unused_protein_info.Parameter == 'id_list'].loc[0, 'Value']],
            ups_mu=ups_mu,
            ups_0=[ups_0],
            mol_mass=[unused_protein_info[unused_protein_info.Parameter == 'mol_mass'].loc[3, 'Value']])
        unused_protein_sector.id = 'UnusedProteinSector'
    else:
        unused_protein_sector = None

    if total_protein: total_protein = TOTAL_PROTEIN_CONCENTRATION

    pamodel = PAModel(id_or_model=model, p_tot=total_protein,
                       active_sector=active_enzyme_sector, translational_sector=translation_enzyme_sector,
                       unused_sector=unused_protein_sector, sensitivity=sensitivity)

    pamodel.BIOMASS_REACTION = BIOMASS_REACTION
    return pamodel

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