
import cobra
import pandas as pd
from cobra.io import load_model
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
from cobra import Model, Reaction, Metabolite, Gene
import numpy as np
from pathlib import Path
import requests
import time


def load_model(subfolder, filename):
    """Load an SBML model from a given subfolder and filename."""
    base_dir = Path.home() / "Documents" / "PhD" / "10-19 Research" / "11 Data" / "11.09_Models"
    model_path = base_dir / subfolder / filename
    return read_sbml_model(str(model_path))

def metabolite_summary_df(model, summary, file_path):
    """
    Create a dataframe and an excel file for the metabolite.summary function of cobrapy

    Parameters:
    - model: COBRApy model object
    - summary: COBRApy summary function output
    - metabolite: metabolite object
    """

    # Extracting the data for file naming
    year = time.strftime("%Y")[2:4]
    timestr = time.strftime("%m%d")

    # Extract the producing and consuming fluxes
    producing_flux_df = summary.producing_flux
    consuming_flux_df = summary.consuming_flux

    # Add a 'type' column to distinguish between producing and consuming
    producing_flux_df['type'] = 'producing'
    consuming_flux_df['type'] = 'consuming'

    # Combine both DataFrames
    df_summary = pd.concat([producing_flux_df, consuming_flux_df])
    df_summary_sorted = df_summary.sort_values(by='flux', key=abs, ascending=False)

    df_summary_sorted['reaction_name'] = df_summary_sorted.index.map(lambda x: model.reactions.get_by_id(x).name)

    # Add a new column for the reaction formulas
    df_summary_sorted['reaction_formula'] = df_summary_sorted.index.map(lambda x: model.reactions.get_by_id(x).reaction)
    df_summary_sorted['reaction_genes'] = df_summary_sorted.index.map(lambda x: [gene.id for gene in model.reactions.get_by_id(x).genes])

    return df_summary_sorted.to_excel(file_path, index=False)

def search_reaction_info(model, search_word):
    """
    Searches for reactions in the model whose ID contains the search word.
    Prints detailed information including name, bounds, genes, and annotations.

    Parameters:
    - model: COBRApy model object
    - search_word: string to match part of the reaction ID
    """
    found = False
    for reaction in model.reactions:
        if search_word in reaction.id:
            found = True
            print(f"\n🔎 Reaction ID: {reaction.id}")
            print(f" - Name        : {reaction.name}")
            print(f" - Equation    : {reaction.reaction}")
            print(f" - Lower Bound : {reaction.lower_bound}")
            print(f" - Upper Bound : {reaction.upper_bound}")
            print(f" - Annotations : {reaction.annotation}")
            print(f" - Genes       : {[gene.id for gene in reaction.genes]}")
    
    if not found:
        print(f"No reactions found with ID containing '{search_word}'")
   

def search_metabolite_info(model, search_word):
    """
    Searches for metabolites in the model whose ID contains the search_word.
    Prints detailed information including related reactions and associated genes.
    
    Parameters:
    - model: COBRApy model object
    - search_word: string to search for in metabolite IDs
    """
    for metabolite in model.metabolites:
        if search_word in metabolite.id:
            # Print metabolite basic info
            print(f"\n🧪 Metabolite: {metabolite.id}")
            print(f" - Name       : {metabolite.name}")
            print(f" - Formula    : {metabolite.formula}")
            print(f" - Charge     : {metabolite.charge}")
            print(f" - Compartment: {metabolite.compartment}")
            print(f" - Annotations: {metabolite.annotation}")
            
            # Get related reactions
            related_reactions = list(metabolite.reactions)

            # Get genes related to the reactions
            related_genes = set()
            for reaction in related_reactions:
                related_genes.update(reaction.genes)

            # Print related reaction and gene IDs
            print(f" - Related Reactions: {', '.join(r.id for r in related_reactions)}")
            print(f" - Associated Genes : {', '.join(g.id for g in related_genes)}")


def check_reaction_flux(objective_rxn, model, reaction_id):
    """
    Optimizes the model and checks the flux through a specified reaction.

    Parameters:
    - model: COBRApy model object
    - reaction_id: string, ID of the reaction to check

    Returns:
    - flux_value: float or None
    """
    model.objective = objective_rxn
    solution = model.optimize()

    if solution.status == 'optimal':
        # Check if the reaction ID exists in the solution
        if reaction_id not in solution.fluxes:
            print(f"Reaction '{reaction_id}' not found in the model.")
            return None

        # Retrieve flux value
        flux_value = solution.fluxes[reaction_id]
        if flux_value != 0:
            print(f"The reaction '{reaction_id}' has a flux of {flux_value:.3f} mmol/gDW/hr.")
        else:
            print(f"There is no flux through the reaction '{reaction_id}'.")
        
        return flux_value
    else:
        print("Optimization was not successful. Check model constraints or objective function.")
        return None
    

def get_active_reactions_for_metabolite(model, solution, metabolite_id, flux_threshold=0.5):
    """
    Returns a list of reactions involving a specific metabolite that carry significant flux.

    Parameters:
    - model: COBRApy model object
    - solution: COBRApy solution object (after optimization)
    - metabolite_id: string, ID of the metabolite (e.g. 'h2_c')
    - flux_threshold: float, minimum absolute flux to include (default 0.5)

    Returns:
    - List of [reaction name, reaction ID, equation, flux] for each active reaction
    """
    try:
        metabolite = model.metabolites.get_by_id(metabolite_id)
    except KeyError:
        print(f" Metabolite '{metabolite_id}' not found in model.")
        return []

    active_reactions = [
        [rxn.name, rxn.id, rxn.reaction, solution.fluxes[rxn.id]]
        for rxn in metabolite.reactions
        if abs(solution.fluxes[rxn.id]) > flux_threshold
    ]

    return active_reactions


def print_reaction_details(model, reaction_id):
    """
    Prints detailed information about a reaction in the model.

    Parameters:
    - model: COBRApy model object
    - reaction_id: string, ID of the reaction (e.g., 'NOR_syn_1')
    """
    try:
        reaction = model.reactions.get_by_id(reaction_id)
    except KeyError:
        print(f"Reaction '{reaction_id}' not found in the model.")
        return

    print(f"\nReaction ID   : {reaction.id}")
    print(f"Name          : {reaction.name}")
    print(f"Equation      : {reaction.reaction}")
    print(f"Lower Bound   : {reaction.lower_bound}")
    print(f"Upper Bound   : {reaction.upper_bound}")

    print("\nMetabolites and Stoichiometry:")
    net_charge = 0
    for metabolite, coefficient in reaction.metabolites.items():
        print(f" - {metabolite.id:<15}: {coefficient:>6}  "
              f"Name: {metabolite.name:<50} "
              f"Charge: {metabolite.charge:>3}  "
              f"Formula: {metabolite.formula}")
        net_charge += coefficient * (metabolite.charge if metabolite.charge else 0)

    print(f"\nSum of Charge Balance: {net_charge}")

    print("\nAssociated Genes:")
    if reaction.genes:
        for gene in reaction.genes:
            print(f" - {gene.id}")
    else:
        print(" - None")


def print_gene_reaction_info(model, gene_id):
    """
    Prints information about a gene and all reactions associated with it.

    Parameters:
    - model: COBRApy model object
    - gene_id: string, ID of the gene (e.g., 'AAFOLC_02510')
    """
    try:
        gene = model.genes.get_by_id(gene_id)
    except KeyError:
        print(f"Gene '{gene_id}' not found in the model.")
        return

    print(f"\nGene ID   : {gene.id}")
    print(f"Name      : {gene.name}\n")

    print("Reactions Associated with the Gene:")
    for reaction in gene.reactions:
        print(f"\n- Reaction ID   : {reaction.id}")
        print(f"  Name          : {reaction.name}")
        print(f"  Equation      : {reaction.reaction}")
        print(f"  Bounds        : [{reaction.lower_bound}, {reaction.upper_bound}]")
        
        print("  Metabolites:")
        for met, coeff in reaction.metabolites.items():
            print(f"    {met.id:<15}: {coeff:>6}")


def check_model_mass_balance(model):
    """
    Checks for mass or charge imbalances in the reactions of a COBRA model.

    Parameters:
    - model: COBRApy model object

    Prints:
    - The number of imbalanced reactions
    - Each reaction's ID and its imbalance details
    """
    imbalanced_reactions = check_mass_balance(model)

    print(f"\nNumber of imbalanced reactions: {len(imbalanced_reactions)}")

    for reaction, imbalance in imbalanced_reactions.items():
        print(f"Reaction ID: {reaction.id}, Imbalance: {imbalance}")

def update_biomass_from_dataframe(model, biomass_df, reaction_id='Growth', id_column='Biomass_equation', coeff_column='normalized_mmol/g'):
    """
    Updates metabolite coefficients in a biomass reaction based on values in a given DataFrame.

    Parameters:
    - model: COBRApy model object
    - biomass_df: pandas DataFrame with metabolite IDs and updated coefficients
    - reaction_id: string, ID of the biomass reaction in the model (default: 'Growth')
    - id_column: string, column name in the DataFrame containing metabolite IDs
    - coeff_column: string, column name in the DataFrame containing new coefficient values

    Returns:
    - int: number of metabolites successfully updated
    """
    reaction = model.reactions.get_by_id(reaction_id)
    update_count = 0

    for element in biomass_df[id_column]:
        found = False
        for metabolite, coefficient in reaction.metabolites.items():
            if element == metabolite.id:
                found = True
                # Get the updated coefficient from the DataFrame
                new_coeff = biomass_df.loc[biomass_df[id_column] == element, coeff_column].values[0]

                print(f'\n✔ Updating {element}')
                print(f' - Original coefficient in model : {coefficient}')
                print(f' - New coefficient from DataFrame: {new_coeff}')

                # Update metabolite coefficient (replace, not add)
                reaction.add_metabolites({metabolite: new_coeff}, combine=False)
                updated_coefficient = reaction.metabolites[metabolite]

                print(f' - Updated coefficient           : {updated_coefficient}')

                # Test optimization
                solution = model.slim_optimize()
                print(f' - Growth rate after update      : {solution}')

                update_count += 1
                break

        if not found:
            print(f'⚠ Metabolite {element} not found in biomass reaction.')

    print(f"\n✅ Total metabolites updated: {update_count}")
    return update_count


def find_missing_biomass_metabolites(model, biomass_df, reaction_id='Growth',
                                     id_column='Biomass_equation',
                                     coeff_column='normalized_mmol/g',
                                     name_column='Unnamed: 1',
                                     output_path=None):
    """
    Finds metabolites from a biomass DataFrame that are not present in the model's biomass reaction.

    Parameters:
    - model: COBRApy model object
    - biomass_df: pandas DataFrame containing biomass components
    - reaction_id: string, ID of the biomass reaction (default: 'Growth')
    - id_column: string, column name for metabolite IDs
    - coeff_column: string, column name for coefficient values
    - name_column: string, column name for metabolite names
    - output_path: optional string path to save the result as Excel

    Returns:
    - DataFrame of missing metabolites (metabolite_id, metabolite_name, coefficient)
    """
    reaction = model.reactions.get_by_id(reaction_id)
    missing_data = []

    for element in biomass_df[id_column]:
        if element not in [m.id for m in reaction.metabolites]:
            try:
                normalized_val = biomass_df.loc[biomass_df[id_column] == element, coeff_column].values[0]
                name_val = biomass_df.loc[biomass_df[id_column] == element, name_column].values[0]
                missing_data.append({
                    'metabolite_id': element,
                    'metabolite_name': name_val,
                    'coefficient': normalized_val
                })
            except IndexError:
                print(f"⚠ Could not find values for metabolite '{element}' in the DataFrame.")

    df_missing = pd.DataFrame(missing_data)

    print(f"\n🔎 Total missing metabolites: {len(df_missing)}")

    if output_path:
        df_missing.to_excel(output_path, index=False)
        print(f"📁 Missing metabolite list exported to: {output_path}")

    return df_missing

def find_unlisted_biomass_metabolites(model, biomass_df,
                                      reaction_id='Growth',
                                      id_column='Biomass_equation',
                                      output_path=None):
    """
    Finds metabolites that are in the model's biomass reaction but not in the given biomass DataFrame.

    Parameters:
    - model: COBRApy model object
    - biomass_df: pandas DataFrame containing biomass components
    - reaction_id: string, ID of the biomass reaction (default: 'Growth')
    - id_column: string, column name for metabolite IDs in the DataFrame
    - output_path: optional string path to save the result as Excel

    Returns:
    - DataFrame of unlisted metabolites (metabolite_id, metabolite_name, coefficient)
    """
    reaction = model.reactions.get_by_id(reaction_id)
    model_metabolite_ids = {met.id for met in reaction.metabolites}
    listed_metabolite_ids = set(biomass_df[id_column])

    missing_in_df = model_metabolite_ids - listed_metabolite_ids
    unlisted_data = []

    for met in reaction.metabolites:
        if met.id in missing_in_df:
            unlisted_data.append({
                'metabolite_id': met.id,
                'metabolite_name': met.name,
                'coefficient': reaction.metabolites[met]
            })

    df_unlisted = pd.DataFrame(unlisted_data)

    print(f"\n🔎 Metabolites in model biomass but not in DataFrame: {len(df_unlisted)}")

    if output_path:
        df_unlisted.to_excel(output_path, index=False)
        print(f"📁 Result saved to: {output_path}")

    return df_unlisted


def compare_biomass_metabolites_between_models(model_a, model_b,
                                               reaction_id_a='Growth',
                                               reaction_id_b='BIOMASS_KT2440_WT3'):
    """
    Compares biomass reactions from two models and returns metabolites present in model_b but missing in model_a.

    Parameters:
    - model_a: COBRApy model object (reference model)
    - model_b: COBRApy model object (comparison model)
    - reaction_id_a: biomass reaction ID in model_a (default: 'Growth')
    - reaction_id_b: biomass reaction ID in model_b (default: 'BIOMASS_KT2440_WT3')

    Returns:
    - DataFrame containing the list of metabolites not found in model_a's biomass reaction
    """
    reaction_a = model_a.reactions.get_by_id(reaction_id_a)
    reaction_b = model_b.reactions.get_by_id(reaction_id_b)

    metabolites_not_found = []

    for met_b in reaction_b.metabolites:
        if met_b.id not in {met.id for met in reaction_a.metabolites}:
            metabolites_not_found.append({
                'metabolite_id': met_b.id,
                'metabolite_name': met_b.name
            })

    df_missing = pd.DataFrame(metabolites_not_found)
    print(f"Number of metabolites not found in '{reaction_id_a}': {len(df_missing)}")
    return df_missing

def update_GAM_in_biomass(model, reaction_id='Growth', gam_dict_1=None, gam_dict_2=None, verbose=True):
    """
    Updates the GAM-related metabolite coefficients in the biomass reaction and optimizes the model.

    Parameters:
    - model: COBRApy model object
    - reaction_id: string, ID of the biomass reaction (default: 'Growth')
    - gam_dict_1: dict, metabolites and coefficients (usually negative side)
    - gam_dict_2: dict, metabolites and coefficients (usually positive side)
    - verbose: bool, whether to print details during update

    Returns:
    - float: optimized objective value (growth rate)
    """
    if gam_dict_1 is None:
        gam_dict_1 = {'atp_c': -49.0, 'h2o_c': -49.0}
    if gam_dict_2 is None:
        gam_dict_2 = {'adp_c': 49.0, 'pi_c': 49.0, 'h_c': 49.0}

    reaction = model.reactions.get_by_id(reaction_id)

    for gam_dict in [gam_dict_1, gam_dict_2]:
        for met_id, new_coeff in gam_dict.items():
            metabolite = model.metabolites.get_by_id(met_id)
            if metabolite in reaction.metabolites:
                old_coeff = reaction.metabolites[metabolite]
                reaction.add_metabolites({metabolite: new_coeff - old_coeff})
                if verbose:
                    print(f'Updated {met_id}: {old_coeff} → {new_coeff}')
            else:
                if verbose:
                    print(f"⚠️ Metabolite {met_id} not found in biomass reaction.")

    # Re-optimize the model
    solution = model.slim_optimize()
    if verbose:
        print(f"\nOptimal growth rate after GAM update: {solution}")

    return solution

def update_metabolite_charges_from_excel(model, excel_path, id_column='metabolite changed', charge_column='New charge'):
    """
    Updates metabolite charges in a COBRA model based on values from an Excel file.

    Parameters:
    - model: COBRApy model object
    - excel_path: str or Path, path to the Excel file
    - id_column: str, column name for metabolite IDs in the Excel file
    - charge_column: str, column name for new charge values in the Excel file

    Returns:
    - int: number of successfully updated metabolites
    """
    df = pd.read_excel(excel_path)
    updated_count = 0

    for met_id in df[id_column]:
        if met_id in model.metabolites:
            metabolite = model.metabolites.get_by_id(met_id)
            new_charge = int(df.loc[df[id_column] == met_id, charge_column].values[0])
            old_charge = metabolite.charge
            metabolite.charge = new_charge
            print(f"✔ Modified metabolite: {met_id} | Old charge: {old_charge} → New charge: {new_charge}")
            updated_count += 1
        else:
            print(f"⚠ Metabolite '{met_id}' not found in the model.")

    print(f"\n🔁 Total metabolites updated: {updated_count}")
    return updated_count


def get_bigg_models_for_imbalanced_reactions(model, limit=None, verbose=True):
    """
    Retrieves BiGG model IDs for reactions that are mass/charge imbalanced in the given model.

    Parameters:
    - model: COBRApy model object
    - limit: optional int, max number of reactions to check (for testing or speed)
    - verbose: bool, whether to print progress

    Returns:
    - dict: mapping of reaction ID → model ID from BiGG database
    """
    imbalanced = check_mass_balance(model)
    reaction_model_mapping = {}

    for i, (reaction, balance_dict) in enumerate(imbalanced.items()):
        if limit and i >= limit:
            break
        reaction_id = reaction.id
        url = f'http://bigg.ucsd.edu/api/v2/universal/reactions/{reaction_id}'

        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            models = data.get('models_containing_reaction', [])
            if models:
                reaction_model_mapping[reaction_id] = models[0]['bigg_id']
                if verbose:
                    print(f"{reaction_id} → {models[0]['bigg_id']}")
            else:
                if verbose:
                    print(f"{reaction_id} found, but no models listed.")
        else:
            print(f"Failed to retrieve {reaction_id}. Status code: {response.status_code}")

    print(f"\n✅ Retrieved model info for {len(reaction_model_mapping)} imbalanced reactions.")
    return reaction_model_mapping

def update_metabolite_charges_from_bigg(model, reaction_model_mapping):
    """
    Updates zero-charge metabolites in reactions using charge values from the BiGG database.

    Parameters:
    - model: COBRApy model object
    - reaction_model_mapping: dict mapping reaction IDs to BiGG model names

    Returns:
    - int: number of metabolites whose charges were updated
    """
    updated_count = 0

    for reaction_id, bigg_model_name in reaction_model_mapping.items():
        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError:
            print(f"Reaction '{reaction_id}' not found in model.")
            continue

        for metabolite, coefficient in reaction.metabolites.items():
            if metabolite.charge == 0:
                url = f'http://bigg.ucsd.edu/api/v2/models/{bigg_model_name}/metabolites/{metabolite.id}'
                response = requests.get(url)

                if response.status_code == 200:
                    data = response.json()
                    bigg_charge = data.get('charge')
                    if bigg_charge is not None and bigg_charge != metabolite.charge:
                        print(f"Updating metabolite '{metabolite.id}'")
                        print(f" - Current charge: {metabolite.charge}")
                        print(f" - BiGG charge   : {bigg_charge}")
                        metabolite.charge = bigg_charge
                        updated_count += 1
                else:
                    print(f"Failed to retrieve data for '{metabolite.id}'. Status code: {response.status_code}")

    print(f"\nTotal metabolites updated: {updated_count}")
    return updated_count


def update_metabolite_info(non_curated_model, curated_model):
    """
    Update the charge and formula of metabolites in the non-curated model based on the curated model.

    Parameters:
    - non_curated_model: cobra.Model
      The non-curated model to update.
    - curated_model: cobra.Model
      The curated model with correct metabolite information.
      
    Returns:
    - updated_model: cobra.Model
      The non-curated model with updated metabolite information.
    """
    # Create a copy of the non-curated model
    copy_of_model = non_curated_model.copy()
    # Create a copy of the curated model
    copy_of_model1 = curated_model.copy()

    # Iterate through metabolites in the non-curated model
    for metabolite in copy_of_model.metabolites:
        # Check if the metabolite exists in the curated model
        if metabolite.id in copy_of_model1.metabolites:
            # Find the corresponding metabolite in the curated model
            curated_metabolite = copy_of_model1.metabolites.get_by_id(metabolite.id)
            # Update the charge and formula of the non-curated metabolite
            metabolite.charge = curated_metabolite.charge
            metabolite.formula = curated_metabolite.formula
        else:
            print(f"Metabolite '{metabolite.id}' not found in the curated model.")
    
    return copy_of_model

def add_missing_reactions(non_curated_model, curated_model):
    """
    Add missing reactions from the curated model to the non-curated model.

    Parameters:
    - non_curated_model: cobra.Model
      The non-curated model to update.
    - curated_model: cobra.Model
      The curated model with the complete set of reactions.
      
    Returns:
    - updated_model: cobra.Model
      The non-curated model with added reactions from the curated model.
    """
    # Create a copy of the non-curated model
    updated_model = non_curated_model.copy()

    # Iterate through reactions in the curated model
    for reaction in curated_model.reactions:
        # Check if the reaction exists in the non-curated model
        if reaction.id not in updated_model.reactions:
            # Add the reaction to the non-curated model
            updated_model.add_reactions([reaction])
            print(f"Reaction '{reaction.id}' was added to the non-curated model.")
        else:
            pass

    # Optimize the updated model
    solution = updated_model.optimize()
    print(f'Biomass reaction flux after adding reactions: {solution.objective_value}')

    # Verify the addition of reactions
    for reaction in curated_model.reactions:
        if reaction.id not in updated_model.reactions:
            print(f"Reaction '{reaction.id}' was not added correctly.")
    
    return updated_model

def add_missing_metabolites(non_curated_model, curated_model):
    """
    Add missing metabolites from the curated model to the non-curated model.

    Parameters:
    - non_curated_model: cobra.Model
      The non-curated model to update.
    - curated_model: cobra.Model
      The curated model with the complete set of metabolites.
      
    Returns:
    - updated_model: cobra.Model
      The non-curated model with added metabolites from the curated model.
    """
    # Create a copy of the non-curated model
    updated_model = non_curated_model.copy()

    # Iterate through metabolites in the curated model
    for metabolite in curated_model.metabolites:
        # Check if the metabolite exists in the non-curated model
        if metabolite.id not in updated_model.metabolites:
            # Create a new metabolite object
            new_metabolite = Metabolite(
                id=metabolite.id,             # Metabolite ID
                name=metabolite.name,         # Metabolite Name
                formula=metabolite.formula,   # Chemical Formula
                compartment=metabolite.compartment, # Compartment (e.g., cytoplasm)
            )

            # Add additional properties
            new_metabolite.charge = metabolite.charge
            new_metabolite.annotation = metabolite.annotation

            # Add the new metabolite to the non-curated model
            updated_model.add_metabolites([new_metabolite])
            print(f"Metabolite '{metabolite.id}' was added to the non-curated model.")
        else:
            pass

    return updated_model



def export_reactions_to_excel(model, file_path):
    """
    Export reaction details from the model to an Excel file.

    Parameters:
    - model: cobra.Model
      The metabolic model from which reactions will be extracted.
    - file_path: str
      The file path where the reaction data will be saved as an Excel file.
      
    Returns:
    - None
    """
    # Create an empty DataFrame to store reaction data
    df_reactions = pd.DataFrame(columns=[
        'reaction_id', 'reaction_name', 'reaction_equation', 
        'reaction_lower_bound', 'reaction_upper_bound', 
        'reaction_gene', 'reaction_annotation'
    ])

    # Iterate over all reactions in the model and populate the DataFrame
    for reaction in model.reactions:
        # Create a DataFrame for the current reaction
        reaction_data = pd.DataFrame([{
            "reaction_id": reaction.id,
            "reaction_name": reaction.name,
            "reaction_equation": reaction.reaction,
            "reaction_lower_bound": reaction.lower_bound,
            "reaction_upper_bound": reaction.upper_bound,
            "reaction_gene": reaction.gene_reaction_rule,
            "reaction_annotation": reaction.annotation
        }])
        
        # Concatenate the current reaction's DataFrame to the main DataFrame
        df_reactions = pd.concat([df_reactions, reaction_data], ignore_index=True)

    # Export the DataFrame to an Excel file
    df_reactions.to_excel(file_path, index=False)
    print(f"Reaction data has been exported to {file_path}")

    
def export_metabolites_to_excel(model, file_path):
    """
    Export metabolite details from the model to an Excel file.

    Parameters:
    - model: cobra.Model
      The metabolic model from which metabolites will be extracted.
    - file_path: str
      The file path where the metabolite data will be saved as an Excel file.
      
    Returns:
    - None
    """
    # Create an empty DataFrame to store metabolite data
    df_metabolites = pd.DataFrame(columns=[
        'metabolite_id', 'metabolite_name', 'metabolite_formula', 
        'metabolite_charge', 'related_reactions'
    ])

    # Iterate over all metabolites in the model and populate the DataFrame
    for metabolite in model.metabolites:
        # Create a DataFrame for the current metabolite
        metabolite_data = pd.DataFrame([{
            "metabolite_id": metabolite.id,
            "metabolite_name": metabolite.name,
            "metabolite_formula": metabolite.formula,
            "metabolite_charge": metabolite.charge,
            "related_reactions": list(metabolite.reactions),
        }])
        
        # Concatenate the current metabolite's DataFrame to the main DataFrame
        df_metabolites = pd.concat([df_metabolites, metabolite_data], ignore_index=True)

    # Export the DataFrame to an Excel file
    df_metabolites.to_excel(file_path, index=False)
    print(f"Metabolite data has been exported to {file_path}")

def export_metabolites_from_reaction(model, reaction_id, file_path):
    """
    Export metabolite details from a specific reaction in the model to an Excel file.

    Parameters:
    - model: cobra.Model
      The metabolic model containing the reaction.
    - reaction_id: str
      The ID of the reaction to extract metabolites from.
    - file_path: str
      The file path where the metabolite data will be saved as an Excel file.
      
    Returns:
    - None
    """
    # Create an empty DataFrame to store metabolite data
    df_reaction = pd.DataFrame(columns=['met_id', 'met_name', 'met_coeff', 'met_charge'])

    # Get the reaction from the model
    try:
        reaction = model.reactions.get_by_id(reaction_id)
        
        # Iterate through metabolites and coefficients in the reaction
        for metabolite, coefficient in reaction.metabolites.items():
            print(f"{metabolite.id:<15}: {coefficient:>25}  Name: {metabolite.name:<60}  Charge: {metabolite.charge:>3}  Formula: {metabolite.formula}")
            
            # Create a DataFrame for the current metabolite
            reaction_data = pd.DataFrame([{
                "met_id": metabolite.id,
                "met_name": metabolite.name,
                "met_coeff": coefficient,
                "met_charge": metabolite.charge,
            }])

            # Concatenate the current metabolite's DataFrame to the main DataFrame
            df_reaction = pd.concat([df_reaction, reaction_data], ignore_index=True)

        # Export the DataFrame to an Excel file
        df_reaction.to_excel(file_path, index=False)
        print(f"Metabolite data from reaction '{reaction_id}' has been exported to {file_path}")
    
    except KeyError:
        print(f"Reaction '{reaction_id}' not found in the model.")