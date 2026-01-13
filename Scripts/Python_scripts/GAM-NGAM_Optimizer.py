import os
import numpy as np
import pandas as pd
from cobra.io import load_model
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
from cobra import Model, Reaction, Metabolite, Gene
from sklearn.metrics import r2_score 
from datetime import datetime



def compute_r2_score(excel_data):
    """
    Computes the R² (coefficient of determination) between experimental and simulated growth rates.
    """
    if 'dilRate' not in excel_data.columns or 'growth_rate' not in excel_data.columns:
        raise ValueError("DataFrame must contain 'dilRate' and 'growth_rate' columns.")
    
    # Remove rows with NaN in either column
    df_clean = excel_data.dropna(subset=['dilRate', 'growth_rate'])

    if df_clean.empty:
        raise ValueError("No valid data left after dropping NaNs.")

    r2_value = r2_score(df_clean['dilRate'], df_clean['growth_rate'])


    return r2_value

def compute_mean_difference(excel_data):
    """
    Computes the mean absolute difference between experimental and simulated growth rates.

    Parameters:
    - excel_data: Pandas DataFrame containing 'dilRate' (experimental) and 'growth_rate' (simulated).

    Returns:
    - mean_difference: The average absolute difference between experimental and simulated values.
    """
    
    # Ensure columns exist
    if 'dilRate' not in excel_data.columns or 'growth_rate' not in excel_data.columns:
        raise ValueError("DataFrame must contain 'dilRate' and 'growth_rate' columns.")
    
    # Compute absolute difference for each row
    differences = abs(excel_data['dilRate'] - excel_data['growth_rate'])
    # Compute the mean difference
    mean_difference = differences.mean()

    return mean_difference

def compute_growth_rates(model, excel_data):
    """
    Computes growth rates by adjusting H2, O2, and CO2 fluxes in the model and optimizing.

    Parameters:
    - model: COBRApy model object
    - excel_data: Pandas DataFrame containing columns 'H2 Flux theoretic' and 'CO2 flux theoretic'

    Returns:
    - Updated DataFrame with calculated growth rates and model-predicted fluxes.
    """
    model.objective = 'Growth'


    # Add columns for growth rate and modeled fluxes
    excel_data['growth_rate'] = None
    excel_data['H2 Flux model'] = None
    excel_data['O2 Flux model'] = None
    excel_data['CO2 Flux model'] = None

    # Get the exchange reactions
    h2_reaction = model.reactions.get_by_id('EX_h2_e')
    o2_reaction = model.reactions.get_by_id('EX_o2_e')
    co2_reaction = model.reactions.get_by_id('EX_co2_e')

    # Iterate over each row in the DataFrame
    for index, row in excel_data.iterrows():
        
        # Set bounds for H2 flux
        h2_reaction.bounds = (-row['H2 Flux theoretic'], 0.0)

        # Set bounds for O2 flux
        o2_reaction.bounds = (-100.0, 0.0)

        # Set bounds for CO2 flux
        co2_reaction.bounds = (-100.0, 0.0)

        # Optimize the model
        solution = model.slim_optimize()
        #print(solution)

        # Store results in DataFrame
        excel_data.at[index, 'growth_rate'] = solution
        excel_data.at[index, 'H2 Flux model'] = h2_reaction.flux
        excel_data.at[index, 'O2 Flux model'] = o2_reaction.flux
        excel_data.at[index, 'CO2 Flux model'] = co2_reaction.flux
    return excel_data

def calculate_difference_for_gam(model, gam_value):
    """
    Updates the Growth-Associated Maintenance (GAM) coefficients in the biomass reaction.
    
    Parameters:
    - model: COBRApy model object
    - gam_value: integer, the new GAM coefficient value

    Returns:
    - solution: The optimized solution after modifying the biomass reaction
    """

    # Define the GAM metabolite changes
    GAM_metabolites_1 = {'atp_c': -gam_value, 'h2o_c': -gam_value}
    GAM_metabolites_2 = {'adp_c': gam_value, 'pi_c': gam_value, 'h_c': gam_value}

    # Get the biomass reaction
    reaction = model.reactions.get_by_id('Growth')

    # Function to update metabolite coefficients in the reaction
    def update_coefficients(metabolites):
        for element, new_coefficient in metabolites.items():
            metabolite = model.metabolites.get_by_id(element)
            if metabolite in reaction.metabolites:
                old_coefficient = reaction.metabolites[metabolite]
                reaction.add_metabolites({metabolite: new_coefficient - old_coefficient})

    # Update coefficients for both metabolite groups
    update_coefficients(GAM_metabolites_1)
    update_coefficients(GAM_metabolites_2)


    compute_growth_rates(model,excel_data)

    mean_diff = compute_mean_difference(excel_data)
    r2_value  = compute_r2_score(excel_data)

    
    return mean_diff,r2_value





if __name__ == "__main__":
    model_name = input("Insert the name of the model : ").strip()
    type_data = input("Do you want to include all the date points ? Yes or No: ").strip()
    removal_data = input("Do you want to remove the outlier data points ? Yes or No: ").strip()

    
    model_file_path = os.path.join("Models", model_name)
    
    model = read_sbml_model(model_file_path)

    
    if removal_data.lower() == 'yes' and type_data.lower() == "no" :
        exp_data_path =os.path.join("Data","240823_growth_data.xlsx" ) # the average point per sample
        excel_data = pd.read_excel(exp_data_path)
        excel_data = excel_data.drop([4,5])
    elif removal_data.lower() == 'no' and type_data.lower() == "no" :
        exp_data_path =os.path.join("Data","240823_growth_data.xlsx" ) # the average point per sample
        excel_data = pd.read_excel(exp_data_path)
    else:
        exp_data_path =os.path.join("Data","240716_exp_gas_consumption.xlsx" ) # all data points are included
        excel_data = pd.read_excel(exp_data_path)



    # Define the H2 uptake range and step size
    h2_range = (8.7, 14.7)  # Start and end values
    step_size = 0.1  # Step size for iteration

    ngam_values = []  # Store ATPM flux values

    # list to store result
    results = []
    # initilize the ATPM bounds to avoid infeasible solution in the first itteration
    model.reactions.get_by_id('ATPM').bounds = -1000.0 , 1000.0
    # Iterate over the defined range
    for uptake in np.arange(h2_range[0], h2_range[1] + step_size, step_size):
        # Get the reaction
        reaction = model.reactions.get_by_id('EX_h2_e')
        reaction.bounds = (-uptake, 1000.0)  # Set H2 uptake flux
        
        model.objective = 'ATPM'  # Set ATP maintenance as the objective
        solution = model.optimize()  # Optimize model
        
        #ngam_values.append(solution.objective_value)  # Store ATPM values

        reaction = model.reactions.get_by_id('ATPM') 
        reaction.bounds = solution.objective_value, 1000.0 #change ATPM bounds

        print('uptake :',uptake)
        print('flux through atpm ',solution.objective_value)
        for gam_value in np.arange(10,150,10):
            mean_difference, r2_value = calculate_difference_for_gam(model, gam_value)  # Get both values
            print(f"GAM Value: {gam_value}, Mean Difference: {mean_difference}")
            # Store the result in a dictionary
            results.append({
                "H2_uptake": uptake,
                "ATPM_flux": solution.objective_value,
                "GAM_value": gam_value,
                "Mean_Difference": mean_difference,
                "R2_Score": r2_value  # Store the R² score
            })
    date_str = datetime.now().strftime("%y%m%d")
    # Convert results into a Pandas DataFrame
    results_df = pd.DataFrame(results)
    # Save to an Excel or CSV file for further analysis
    results_df.to_excel(os.path.join("Results","NGAM_GAM_optimizer",f"{date_str}_gam_analysis_results.xlsx"), index=False)
    


    # Find the row with the minimum Mean Difference
    result_low_mean = results_df.loc[results_df['Mean_Difference'].idxmin()]
    result_r2_square = results_df.loc[results_df['R2_Score'].idxmax()]

    # Print the best GAM value with lowest mean difference
    print("Best GAM Value with Lowest Mean Difference:")
    print(result_low_mean)
    print("Best GAM Value with highest r2 score:")
    print(result_r2_square)

    # to run the script : python Scripts\Python_scripts\GAM-NGAM_Optimizer.py