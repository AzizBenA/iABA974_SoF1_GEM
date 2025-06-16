import pandas as pd
from cobra import Configuration
from cobra import Model, Reaction, Metabolite
import os
from scipy.stats import linregress

import numpy as np

#importing the tools from the PAModelpy package
from PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from PAModelpy.PAModel import PAModel
from PAModelpy.configuration import Config

Config.BIOMASS_REACTION = 'R7'
Config.GLUCOSE_EXCHANGE_RXNID = 'R1'
Config.CO2_EXHANGE_RXNID = 'R8'
Config.ACETATE_EXCRETION_RXNID = 'R9'

#need to have gurobipy installed

BASE_DIR = os.path.split(os.getcwd())[0]
DATA_DIR = os.path.join(BASE_DIR, 'Data')
RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')


#global variables:
global metabolites, n, m, Etot
#global variables:
global metabolites, n, m, Etot
metabolites = ['Substrate', 'ATP', 'CO2', 'Precursor', 'Biomass', 'Byproduct', 'Intermediate']
n = 9
m = 7
Etot = 0.6*1e-3 #will be adjusted in the model with 1e3

#functions:
def build_toy_gem():
    '''
    Rebuild the toymodel as in the MATLAB script.
    sub int byp atp co2 pre bio
R1 = [ 1,  0,  0,  0,  0,  0,  0];
R2 = [-1,  1,  0,  0,  1,  0,  0];
R3 = [ 0, -1,  1,  1,  0,  0,  0];
R3r= -R3;
R4 = [ 0, -1,  0,  2,  1,  0,  0];
R5 = [ 0, -1,  0,  0,  0,  1,  0];
R6 = [ 0,  0,  0, -1,  0, -1,  1];
R7 = [ 0,  0,  0,  0,  0,  0, -1];
R8 = [ 0,  0,  0,  0, -1,  0,  0];
R9 = [ 0,  0, -1,  0,  0,  0,  0];
S  = [R1;R2;R3;R3r;R4;R5;R6;R7;R8;R9]';

    :return: Cobrapy model instance as model
    '''
    #set up model basics
    model = Model('toy_model')
    cobra_config = Configuration()
    cobra_config.solver = 'gurobi'
    for i in range(1, n + 1):
        rxn = Reaction('R' + str(i))
        lower_bound = 0
        upper_bound = 1e6
        #force flux through the system
        if i == 1:
            lower_bound = 1
        #reversible reactions 3, 5 and 9
        if i ==3 or i==5 or 1==9:
            lower_bound =  -1e6
        #constrain nutrient (substrate or byproduct) uptake rate
        if i != 1 or i != 9:
            upper_bound = 100
        else:
            upper_bound = 10

        rxn.lower_bound = lower_bound
        rxn.upper_bound = upper_bound
        model.add_reactions([rxn])


    # add metabolites to the reactions:
    # R1:
    r1 = model.reactions.get_by_id('R1')
    r1.add_metabolites({Metabolite('Substrate'): 1})
    # R2:
    r2 = model.reactions.get_by_id('R2')
    r2.add_metabolites({Metabolite('Substrate'): -1, Metabolite('Intermediate'): 1, Metabolite('CO2'): 1})
    # R3:
    r3 = model.reactions.get_by_id('R3')
    r3.add_metabolites({Metabolite('Intermediate'): -1, Metabolite('Byproduct'):1, Metabolite('ATP'):1})
    # R4:
    r4 = model.reactions.get_by_id('R4')
    r4.add_metabolites({Metabolite('Intermediate'): -1, Metabolite('ATP'): 2, Metabolite('CO2'):1})
    # R5:
    r5 = model.reactions.get_by_id('R5')
    r5.add_metabolites({Metabolite('Intermediate'): -1, Metabolite('Precursor'): 1})
    # R6:
    r6 = model.reactions.get_by_id('R6')
    r6.add_metabolites({Metabolite('ATP'): -1, Metabolite('Precursor'): -1, Metabolite('Biomass'): 1})
    # Exchange reactions
    # R7:
    r7 = model.reactions.get_by_id('R7')
    r7.add_metabolites({Metabolite('Biomass'): -1})
    # R8:
    r8 = model.reactions.get_by_id('R8')
    r8.add_metabolites({Metabolite('CO2'): -1})
    # R9:
    r9 = model.reactions.get_by_id('R9')
    r9.add_metabolites({Metabolite('Byproduct'): -1})

    return model

def build_active_enzyme_sector(Config, kcat_fwd: list  =[1, 0.5, 1, 0.5 ,0.45, 1.5] ):
    kcat_rev = [kcat for kcat in kcat_fwd]
    rxn2kcat = {}
    for i in range(n-3): # all reactions have an enzyme, except excretion reactions
        rxn_id = f'R{i+1}'
        # 1e-6 to correct for the unit transformation in the model (meant to make the calculations preciser for different unit dimensions)
        #dummy molmass like in MATLAB script
        rxn2kcat = {**rxn2kcat, **{rxn_id: {f'E{i+1}':{'f': kcat_fwd[i]/(3600*1e-6), 'b': kcat_rev[i]/(3600*1e-6), 'molmass': 1e6}}}}

    return ActiveEnzymeSector(rxn2protein = rxn2kcat, configuration=Config)

def build_unused_protein_sector(Config):
    return UnusedEnzymeSector(id_list = ['R1'], ups_mu=[-0.01*1e-3], ups_0=[0.1*1e-3], mol_mass= [1], configuration=Config)

def build_translational_protein_sector(Config):
    return TransEnzymeSector(id_list = ['R7'], tps_mu=[0.01*1e-3], tps_0=[0.01*1e-3], mol_mass= [1], configuration=Config)
def run_simulations(pamodel, substrate_rates):
    result_df = pd.DataFrame(columns= ['R1_ub','R1', 'R7', 'R8', 'R9'])

    for substrate in substrate_rates:
        pamodel.change_reaction_bounds(rxn_id='R1',
                                       lower_bound=0, upper_bound=substrate)
        print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
        pamodel.optimize()
        if pamodel.solver.status == 'optimal' and pamodel.objective.value>0:
            results_row = []
            for rxn in ['R1', 'R7', 'R8', 'R9']:
                results_row += [pamodel.reactions.get_by_id(rxn).flux]

            result_df.loc[len(result_df)] = [substrate] + results_row
    return result_df

def evaluate_toy_model_fitness(toy_model: PAModel, substrate_rates = [0.001, 0.091],
                               reference_data_file_path:str = 'Scripts/i2_parametrization/Data/toy_model_simulations_ga.csv') -> float:
    """
    Evaluate the fitness of the toymodel compared to the reference dataset generated using kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5]
    :return: float: error average difference of validation and result for the total of substrate uptake range and available reactiosn
    """
    validation_results = pd.read_csv(reference_data_file_path)
    simulation_results = run_simulations(toy_model, substrate_rates)

    error = []
    for rxn in validation_results.columns[2:]:
        line = linregress(x=[abs(substrate_rates[0]), abs(substrate_rates[-1])],
                              y=[simulation_results[rxn].iloc[0], simulation_results[rxn].iloc[-1]])
        ref_data_rxn = validation_results.assign(
            simulation=lambda x: line.intercept + line.slope * x['R1_ub'])

        # simulation mean
        data_average = ref_data_rxn[rxn].mean()
        # error: squared difference
        ref_data_rxn = ref_data_rxn.assign(error=lambda x: (x[rxn] - x['simulation']) ** 2)

        # calculate R^2:
        residual_ss = sum(ref_data_rxn.error)
        total_ss = sum([(data - data_average) ** 2 for data in ref_data_rxn[rxn]])
        # calculating r_squared is only feasible of the numerator and the denomenator are both nonzero
        if (residual_ss == 0) | (total_ss == 0):
            r_squared = 0
        else:
            r_squared = 1 - residual_ss / total_ss

        error += [r_squared]
    return sum(error)/len(error)



if __name__ == "__main__":
    model = build_toy_gem()
    #[1, 0.5, 1, 0.5, 0.45, 1.5] the initial dataset
    # kcat_fwd = [0.000714801759902067,0.0031403256767839,1, 0.0204382115642068,0.00142103969403549, 1.5]
    kcat_fwd = [1, 0.5, 5, 0.1, 0.25, 1.5]  # the 'final' dataset
    active_enzyme = build_active_enzyme_sector(Config, kcat_fwd=kcat_fwd)
    unused_enzyme = build_unused_protein_sector(Config)
    translation_enzyme = build_translational_protein_sector(Config)
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
                      translational_sector = translation_enzyme,
                      unused_sector = unused_enzyme, p_tot=Etot, configuration=Config)


    #optimize biomass formation
    pamodel.objective={pamodel.reactions.get_by_id('R7') :1}
    print('final error: ',evaluate_toy_model_fitness(pamodel, substrate_rates=list(np.arange(1e-3, 1e-1, 1e-2))),'\n\n')

    substrate_rates = np.arange(1e-3, 1e-1, 1e-2)

    simulation_results = run_simulations(pamodel, substrate_rates)
    # simulation_results.to_csv(RESULT_DF_FILE)
    print(simulation_results.to_markdown())
    #
    #
    # pamodel.change_kcat_value('E4', {'R4':{'f':10, 'b':10}})
    # simulation_results = run_simulations(pamodel, substrate_rates)
    # # simulation_results.to_csv(RESULT_DF_FILE)
    # print(simulation_results.to_markdown())
    #
    # print(evaluate_toy_model_fitness(pamodel))
