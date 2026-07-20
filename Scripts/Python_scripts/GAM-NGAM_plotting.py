import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mplt
from typing import Iterable, List, Union, Optional, Dict

from cobra.io import read_sbml_model
from cobra import Model
from PAModelpy import PAModel


def get_model_fluxes(substrate_rates: Iterable,
                     model: Union[PAModel, Model],
                     h2_ex_id: str = "EX_h2_e") -> List[pd.Series]:
    """
    For each H2 uptake rate (positive number in input), constrain EX_h2_e to [-rate, 0]
    and optimize the model objective (whatever is currently set).
    Returns a list of flux Series (sol.fluxes).
    """
    fluxes = []
    for sub in substrate_rates:
        model.reactions.get_by_id(h2_ex_id).lower_bound = -float(sub)
        model.reactions.get_by_id(h2_ex_id).upper_bound = 0.0
        sol = model.optimize()
        fluxes.append(sol.fluxes if sol.status == "optimal" else pd.Series(dtype=float))
    return fluxes


def change_atp_maintenance(model: Union[Model, PAModel],
                           gam: Optional[float] = None,
                           ngam: Optional[float] = None,
                           biomass_rxn_id: str = 'Growth'
                           ) -> Union[PAModel, Model]:
    """
    - NGAM: set ATPM lower bound
    - GAM: modify biomass reaction ATP/ADP/Pi/H coefficients (adds to existing coeffs)
    NOTE: Calling this multiple times on the same model accumulates coefficients.
          Use a fresh copy each time (pickle copy) like you do below.
    """
    if ngam is not None:
        if isinstance(model, PAModel):
            model.change_reaction_bounds('ATPM', ngam)
        else:
            model.reactions.ATPM.lower_bound = float(ngam)
        print(f"Changing NGAM to {ngam} mmol_ATP/gCDW/h")

    if gam is not None:
        gam = float(gam)
        m_to_coeff = {
            model.metabolites.atp_c: -gam,
            **{model.metabolites.get_by_id(mid): gam for mid in ['adp_c', 'pi_c', 'h_c']}
        }
        model.reactions.get_by_id(biomass_rxn_id).add_metabolites(m_to_coeff, combine=True)
        print(f"Changing GAM to {gam} mmol_ATP/gCDW")

    return model


if __name__ == '__main__':

    ENERGY_FIG_FILE = ('Results/NGAM_GAM_optimizer/effect_of_energy_requirements_H2_2.png')

    # ---------- Load one GEM ----------
    
    model = read_sbml_model(os.path.join("Models", "250808_iABA974.sbml"))
    model_id = "iABA974"
    model.reactions.get_by_id('ATPM').bounds = -1000.0 , 1000.0 # initilize the ATPM bounds to avoid infeasible solution in the first itteration

    m_to_coeff = {
            model.metabolites.atp_c: 0.0,
            **{model.metabolites.get_by_id(mid): 0.0 for mid in ['adp_c', 'pi_c', 'h_c']}
        }
    model.reactions.get_by_id("Growth").add_metabolites(m_to_coeff, combine=True)
    print(f"Changing GAM to {0.0} mmol_ATP/gCDW")
    # ---------- Experimental data ----------
    chemostat = pd.read_excel('Data/240716_exp_gas_grouped.xlsx')

    # Choose the experimental columns you want to compare against
    # IMPORTANT: Replace these with your real column names in the Excel file.
    exp_cols = {
        "EX_o2_e":  "O2 flux theoretic",        # e.g. mmol/gCDW/h (put exact column name)
        "EX_co2_e": "CO2 flux theoretic",       # e.g. mmol/gCDW/h
        "Growth":   "dilRate",        # if dilution rate is your growth proxy
    }

    # x-axis = experimental H2 uptake (positive)
    x_exp = np.abs(chemostat["H2 Flux theoretic"].to_numpy(dtype=float))

    # ---------- H2 scan range for model lines ----------
    # You can either:
    # A) scan the same x values as experimental (recommended for clean overlay)
    substrate_rates = x_exp
    # or B) scan a dense range:
    # substrate_rates = np.arange(0, 20.0, 0.5)

    # ---------- GAM/NGAM sets ----------
    gam_ngam_sets = [
        (168.71, 11.55), # exp data without outliers
        (75.72, 15.27), # exp data with outliers
        (40.0, 21.2), # output of the optimizer
    ]

    colors = [mplt.colormaps['Set2'](i / max(1, (len(gam_ngam_sets)-1))) for i in range(len(gam_ngam_sets))]

    # ---------- Plot ----------
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
    axs = axs.flatten()

    # set objective once (growth)
    model.objective = "Growth"

    # baseline model (no changes)
    # p0 = pickle.dumps(model)
    # fluxes_base = get_model_fluxes(substrate_rates=substrate_rates, model=pickle.loads(p0), h2_ex_id="EX_h2_e")


    # build flux lists per condition using fresh copies (avoid accumulating GAM edits)
    flux_lists = []
    labels = []


    for (gam, ngam) in gam_ngam_sets:
        
        p = pickle.dumps(model)
        m_i = change_atp_maintenance(model=pickle.loads(p), ngam=ngam, gam=gam, biomass_rxn_id="Growth")
        m_i.objective = "Growth"
        flux_i = get_model_fluxes(substrate_rates=substrate_rates, model=m_i, h2_ex_id="EX_h2_e")
        flux_lists.append(flux_i)
        labels.append(f"{model_id}  GAM={gam:.3g}, NGAM={ngam:.3g}")
        print(f"for the GAM {gam} and NGAM {ngam} the analysis was done")

    # # include baseline if you want it
    # flux_lists = [fluxes_base] + flux_lists
    # labels = [f"{model_id} Draft"] + labels
    # colors = ["k"] + colors


    # what to plot in each panel (rxn_id, exp_column)
    panel_specs = [
        ("EX_o2_e",  exp_cols["EX_o2_e"]),
        ("EX_co2_e", exp_cols["EX_co2_e"]),
        ("Growth",   exp_cols["Growth"]),
    ]

    for ax, (rxn, exp_col) in zip(axs, panel_specs):
        # experimental y
        y_exp = np.abs(chemostat[exp_col].to_numpy(dtype=float))

        ax.scatter(x_exp, y_exp, color="black", s=35, zorder=3, label="Experiment")

        # simulated lines
        for color, fluxlist, lab in zip(colors, flux_lists, labels):
            # y_sim from model fluxes; Growth is special: use biomass reaction id "Growth"
            if rxn == "Growth":
                y_sim = np.array([abs(f.get("Growth", np.nan)) for f in fluxlist], dtype=float)
            else:
                y_sim = np.array([abs(f.get(rxn, np.nan)) for f in fluxlist], dtype=float)

            ax.plot(substrate_rates, y_sim, color=color, linewidth=2, label=lab)

        ax.set_xlabel("H2 uptake rate [mmol/gCDW/h]")
        ax.set_ylabel(rxn)
        ax.grid(True, alpha=0.2)

    # one legend
    handles, leg_labels = axs[0].get_legend_handles_labels()
    fig.legend(handles=handles, labels=leg_labels, loc="lower center",
               bbox_to_anchor=(0.5, 0.0), ncols=2, frameon=False)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.18)
    plt.savefig(ENERGY_FIG_FILE, dpi=300)
    plt.show()