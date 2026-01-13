import os
import numpy as np
import pandas as pd
from datetime import datetime
from cobra.io import read_sbml_model
from sklearn.metrics import r2_score


# this version is the  update by CHATGPT it also include r2 score for the substrate uptakes and uses the while Model when changing the GAM .





# ===============================
# DATA PROCESSING HELPERS
# ===============================

def load_experimental_data(include_all, remove_outliers):
    if include_all:
        path = os.path.join("Data", "240716_exp_gas_consumption.xlsx")
        df = pd.read_excel(path)
    else:
        path = os.path.join("Data", "240823_growth_data.xlsx")
        df = pd.read_excel(path)
        if remove_outliers:
            df = df.drop([4, 5])
    return df


# ===============================
# MODEL HELPERS
# ===============================

def set_gam(model, gam):
    """Apply GAM without cumulative drift."""
    biomass = model.reactions.get_by_id('Growth')
    
    # Stoichiometry for GAM:
    gam_rxn = {
        'atp_c': -gam,
        'h2o_c': -gam,
        'adp_c': gam,
        'pi_c': gam,
        'h_c': gam,
    }
    
    # Reset and apply
    for met_id, coeff in gam_rxn.items():
        met = model.metabolites.get_by_id(met_id)
        old = biomass.metabolites.get(met, 0)
        biomass.add_metabolites({met: coeff - old})


def compute_simulation(model, df):
    """Simulate all experimental points once for a given GAM."""
    h2 = model.reactions.get_by_id('EX_h2_e')
    o2 = model.reactions.get_by_id('EX_o2_e')
    co2 = model.reactions.get_by_id('EX_co2_e')

    growth_rates = []
    h2_model = []
    o2_model = []
    co2_model = []

    for h2_flux in df['H2 Flux theoretic']:
        with model:
            h2.bounds = (-h2_flux, 0)
            o2.bounds = (-100, 0)
            co2.bounds = (-100, 0)
            sol = model.optimize()
            growth_rates.append(sol.objective_value)
            h2_model.append(h2.flux)
            o2_model.append(o2.flux)
            co2_model.append(co2.flux)

    return pd.DataFrame({
        'growth_rate': growth_rates,
        'H2 Flux model': h2_model,
        'O2 Flux model': o2_model,
        'CO2 Flux model': co2_model
    })


# ===============================
# METRICS
# ===============================

def eval_metrics(exp, sim):
    exp_clean = exp.dropna(subset=['dilRate'])

    mean_diff = np.mean(np.abs(exp_clean['dilRate'] - sim['growth_rate']))
    r2_growth = r2_score(exp_clean['dilRate'], sim['growth_rate'])
    r2_o2 = r2_score(exp['O2 flux theoretic'], abs(sim['O2 Flux model']))
    r2_co2 = r2_score(exp['CO2 flux theoretic'], abs(sim['CO2 Flux model']))

    return mean_diff, r2_growth, r2_o2, r2_co2


# ===============================
# MAIN OPTIMIZER
# ===============================

def optimize_gam_ngam(model, df, h2_range, gam_range, step):
    results = []

    atpm = model.reactions.get_by_id('ATPM')
    h2_ex = model.reactions.get_by_id('EX_h2_e')

    for uptake in np.arange(h2_range[0], h2_range[1] + 1e-9, step):
        with model:
            # NGAM determination
            h2_ex.bounds = (-uptake, 1000)
            model.objective = atpm
            sol = model.optimize()
            ngam = sol.objective_value

        # Fix NGAM permanently for GAM sweeps
        atpm.bounds = (ngam, 1000)

        for gam in gam_range:
            with model:
                set_gam(model, gam)
                model.objective = 'Growth'
                sim = compute_simulation(model, df)
            
            mean_diff, r2_g, r2_o2, r2_co2 = eval_metrics(df, sim)

            results.append({
                'H2_uptake': uptake,
                'NGAM_flux': ngam,
                'GAM_value': gam,
                'Mean_Difference': mean_diff,
                'R2_Growth': r2_g,
                'R2_O2': r2_o2,
                'R2_CO2': r2_co2
            })

            print(f"[OK] H2={uptake}, NGAM={ngam:.1f}, GAM={gam}, mean={mean_diff:.3f}")

    return pd.DataFrame(results)


# ===============================
# CLI WRAPPER
# ===============================

def main():
    model_name = input("Model filename: ").strip()
    use_all = input("Include all datapoints? (yes/no): ").lower() == 'yes'
    drop_outliers = input("Remove outliers? (yes/no): ").lower() == 'yes'

    model = read_sbml_model(os.path.join("Models", model_name))
    df = load_experimental_data(use_all, drop_outliers)
    model.reactions.get_by_id('ATPM').bounds = -1000.0 , 1000.0

    results = optimize_gam_ngam(
        model=model,
        df=df,
        h2_range=(8.7, 14.7),
        gam_range=np.arange(10, 150, 10),
        step=0.1
    )

    date = datetime.now().strftime("%y%m%d")
    out_path = os.path.join("Results", "NGAM_GAM_optimizer", f"{date}_gam_results.xlsx")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    results.to_excel(out_path, index=False)

    print("\n✔ Results saved:", out_path)
    print("\nBEST MODELS:")
    print("Lowest mean diff:\n", results.loc[results['Mean_Difference'].idxmin()])
    print("Best R2 Growth:\n", results.loc[results['R2_Growth'].idxmax()])
    print("Best R2 O2:\n", results.loc[results['R2_O2'].idxmax()])
    print("Best R2 CO2:\n", results.loc[results['R2_CO2'].idxmax()])


if __name__ == "__main__":
    main()
