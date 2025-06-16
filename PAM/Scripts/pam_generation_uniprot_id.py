import cobra
import pandas as pd
import os
from typing import Union
import json

# load PAMpy modules
from PAModelpy.PAModel import PAModel
from PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector, CustomSector
from PAModelpy.configuration import Config
from PAModelpy.utils.pam_generation import set_up_pam, parse_reaction2protein
from cobra.io.sbml import read_sbml_model

'Function library for making Protein Allocation Models as described in the publication'

def setup_ecolicore_pam(total_protein:bool = True,
                         active_enzymes: bool = True,
                         translational_enzymes:bool = True,
                         unused_enzymes:bool = True,
                         sensitivity:bool =True):
    # Setting the relative paths
    PAM_DATA_FILE_PATH = os.path.join('Results',
                                      '1_preprocessing',
                                      'proteinAllocationModel_iML1515_EnzymaticData_250225.xlsx')

    config = Config()
    config.reset()

    # some other constants
    BIOMASS_REACTION = 'BIOMASS_Ecoli_core_w_GAM'
    TOTAL_PROTEIN_CONCENTRATION = 0.16995  # [g_prot/g_cdw]
    config.BIOMASS_REACTION = BIOMASS_REACTION

    # load the genome-scale information
    model = cobra.io.load_json_model(os.path.join('Models', 'e_coli_core.json'))

    #load example data for the E.coli iML1515 model
    if active_enzymes:
        # load active enzyme sector information
        enzyme_db = parse_enzyme_db(pd.read_excel(PAM_DATA_FILE_PATH, sheet_name='ActiveEnzymes'))
        # create enzyme objects for each gene-associated reaction
        rxn2protein, protein2gene = parse_reaction2protein(enzyme_db, model)

        active_enzyme_sector = ActiveEnzymeSector(rxn2protein=rxn2protein, protein2gene=protein2gene,
                                                  configuration=config)

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
        )
    else:
        unused_enzyme_sector = None

    if total_protein: total_protein = TOTAL_PROTEIN_CONCENTRATION

    pa_model = PAModel(id_or_model=model, p_tot=total_protein, sensitivity=sensitivity,
                       active_sector=active_enzyme_sector, translational_sector=translation_enzyme_sector,
                       unused_sector=unused_enzyme_sector, configuration=config)
    return pa_model

def set_up_ecoli_pam(pam_info_file:str= os.path.join('Results',
                                      '1_preprocessing',
                                      'proteinAllocationModel_iML1515_EnzymaticData_250225.xlsx'),
                     model:str = 'iML1515.xml', config:Config = None,
                     total_protein: Union[bool, float] = True, active_enzymes: bool = True,
                     translational_enzymes: bool = True, unused_enzymes: bool = True, sensitivity = True,
                     enzyme_db:pd.DataFrame = None):


    if config is None:
        config = Config()
        config.reset()

    # some other constants
    TOTAL_PROTEIN_CONCENTRATION = 0.258  # [g_prot/g_cdw]

    #setup the gem ecoli iML1515 model
    model = cobra.io.read_sbml_model(os.path.join('Models', model))

    #check if a different total protein concentration is given
    if isinstance(total_protein, float):
        TOTAL_PROTEIN_CONCENTRATION = total_protein

    # load example data for the E.coli iML1515 model
    if active_enzymes:
        # load active enzyme sector information
        if enzyme_db is None:
            enzyme_db = parse_enzyme_db(pd.read_excel(pam_info_file, sheet_name='ActiveEnzymes'))
        # create enzyme objects for each gene-associated reaction
        rxn2protein, protein2gene = parse_reaction2protein(enzyme_db, model)

        active_enzyme_sector = ActiveEnzymeSector(rxn2protein=rxn2protein, protein2gene=protein2gene,
                                                  configuration=config)

    else:
        active_enzyme_sector = None

    if translational_enzymes:
        translational_info = pd.read_excel(pam_info_file, sheet_name='Translational')
        translation_enzyme_sector = TransEnzymeSector(
            id_list=[translational_info[translational_info.Parameter == 'id_list'].loc[0, 'Value']],
            tps_0=[translational_info[translational_info.Parameter == 'tps_0'].loc[1, 'Value']],
            tps_mu=[translational_info[translational_info.Parameter == 'tps_mu'].loc[2, 'Value']],
            mol_mass=[translational_info[translational_info.Parameter == 'mol_mass'].loc[3, 'Value']],
            configuration = config)
    else:
        translation_enzyme_sector = None

    if unused_enzymes:
        unused_protein_info = pd.read_excel(pam_info_file, sheet_name='UnusedEnzyme').set_index('Parameter')

        ups_0 = unused_protein_info.at['ups_0', 'Value']
        ups_mu = unused_protein_info.at['ups_mu', 'Value']

        unused_protein_sector = UnusedEnzymeSector(
            id_list=[unused_protein_info.at['id_list', 'Value']],
            ups_mu=[ups_mu],
            ups_0=[ups_0],
            mol_mass=[unused_protein_info.at['mol_mass', 'Value']],
            configuration = config)
    else:
        unused_protein_sector = None

    if total_protein: total_protein = TOTAL_PROTEIN_CONCENTRATION

    pamodel = PAModel(id_or_model=model, p_tot=total_protein,
                       active_sector=active_enzyme_sector, translational_sector=translation_enzyme_sector,
                       unused_sector=unused_protein_sector, sensitivity=sensitivity, configuration = config
                      )
    return pamodel

def parse_enzyme_db(enzyme_db: pd.DataFrame) -> None:
    try:
        enzyme_db = enzyme_db.drop(['kegg_id', 'Reactants', 'Products', 'EC', 'Length'],axis=1).pivot_table(
                          index=['rxn_id', 'GPR', 'gene', 'uniprot_id', 'molMass'],
                          columns='direction',
                          values='kcat_values').reset_index().rename({'f': 'kcat_f', 'b':'kcat_b'}, axis=1)
    except:
        enzyme_db = enzyme_db.pivot_table(
            index=['rxn_id', 'GPR', 'gene', 'uniprot_id', 'molMass'],
            columns='direction',
            values='kcat_values').reset_index().rename({'f': 'kcat_f', 'b': 'kcat_b'}, axis=1)
    return enzyme_db

def get_rxn2kcat_protein2gene_dict(param_file):
    # create enzyme objects for each gene-associated reaction
    ecoli_pam = set_up_ecoli_pam(param_file, sensitivity=False)
    enzyme_db = pd.read_excel(param_file, sheet_name='ActiveEnzymes').iloc[:,1:]
    rxn2protein, protein2gene = parse_reaction2protein(enzyme_db, ecoli_pam)

    ae_sector = ecoli_pam.sectors.ActiveEnzymeSector
    new_rxn2prot = ae_sector.rxn2protein.copy()
    for rxn, enz_dict in ae_sector.rxn2protein.items():
        if rxn[:2]!='CE':
            for enzyme_id, enzyme_dict in enz_dict.items():
                protein_reaction = enzyme_dict['protein_reaction_association']
                if ae_sector._enzyme_is_enzyme_complex(protein_reaction, enzyme_id):
                    for pr in protein_reaction:
                        if len(pr) > 1:
                            enzyme_complex_id = '_'.join(pr)
                            new_rxn2prot[rxn] = {**new_rxn2prot[rxn],
                                                 **{enzyme_complex_id: enzyme_dict}}
    return new_rxn2prot, protein2gene

def pickle_rxn2kcat_protein2gene_dicts(pam_info_file: str = os.path.join('Data', 'proteinAllocationModel_iML1515_EnzymaticData_240730.xlsx')):
    # setup the gem ecoli iML1515 model
    model = cobra.io.read_sbml_model(os.path.join('Models', 'iML1515.xml'))
    enzyme_db = pd.read_excel(pam_info_file, sheet_name='ActiveEnzymes')

    # create enzyme objects for each gene-associated reaction
    rxn2kcat, protein2gene = get_rxn2kcat_protein2gene_dict(enzyme_db, model)

    # with open('Models/rxn2protein_iML1515_uniprot.json', 'w') as file:
    #     json.dump(rxn2protein, file, indent=4, sort_keys=True)
    with open('Models/rxn2kcat_iML1515_uniprot.json', 'w') as file:
        json.dump(rxn2kcat, file, indent=4, sort_keys=True)
    with open('Models/protein2gene_iML1515_uniprot.json', 'w') as file:
        json.dump(protein2gene, file, indent=4, sort_keys=True)


def increase_kcats_in_parameter_file(kcat_increase_factor: int,
                                     pam_info_file_path_ori: str,
                                     pam_info_file_path_out: str = os.path.join(
    'Data', 'proteinAllocationModel_iML1515_EnzymaticData_240730_multi.xlsx')):

    pam_info_file_ori = pd.read_excel(pam_info_file_path_ori, sheet_name='ActiveEnzymes')
    transprot = pd.read_excel(pam_info_file_path_ori, sheet_name='Translational')
    unusedenz = pd.read_excel(pam_info_file_path_ori, sheet_name='UnusedEnzyme')

    pam_info_file_new = pam_info_file_ori.rename({'kcat_values': 'kcat_values_ori'})
    pam_info_file_new['kcat_values'] = pam_info_file_ori['kcat_values']*kcat_increase_factor

    write_mode = 'w'
    kwargs = {}
    if os.path.isfile(pam_info_file_path_out):
        write_mode = 'a'
        kwargs = {'if_sheet_exists': 'replace'}

    with pd.ExcelWriter(pam_info_file_path_out, engine = 'openpyxl',
                        mode = write_mode, **kwargs) as writer:
        pam_info_file_new.to_excel(writer, sheet_name = 'ActiveEnzymes', index=False)
        transprot.to_excel(writer, sheet_name='Translational', index=False)
        unusedenz.to_excel(writer, sheet_name='UnusedEnzyme', index=False)




if __name__ == '__main__':
    # pam = setup_yeast_pam(os.path.join(
    #     'Results', '1_preprocessing',  'yeast9','proteinAllocationModel_yeast9_EnzymaticData_TurnUp_multi.xlsx'),
    #     sensitivity=True)
    # # pam.change_total_protein_constraint(1)
    # pam.change_reaction_bounds('r_1714', -1e3,0)
    # pam.optimize()
    # print(pam.summary())
    # print(pam.enzyme_sensitivity_coefficients.sort_values('coefficient', ascending=False))
    # print(pam.capacity_sensitivity_coefficients.sort_values('coefficient', ascending=False))

    pam_info_file: str = os.path.join('Results', '1_preprocessing', 'proteinAllocationModel_iML1515_EnzymaticData_241209.xlsx')
    model = cobra.io.read_sbml_model(os.path.join('Models', 'iML1515.xml'))

    pam = set_up_pam(pam_info_file, model)
    print(pam.catalytic_events.CE_ICYSDS.enzymes)
    # enzyme_db = pd.read_excel(pam_info_file, sheet_name='ActiveEnzymes')

    # create enzyme objects for each gene-associated reaction
    # rxn2kcat, protein2gene = get_rxn2kcat_protein2gene_dict(enzyme_db, model)


