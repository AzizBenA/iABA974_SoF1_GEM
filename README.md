# iABA974: genome-scale metabolic model of *Xanthobacter* sp. SoF1

iABA974 is a manually curated genome-scale metabolic reconstruction of the hydrogen-oxidizing bacterium *Xanthobacter* sp. SoF1. It provides a framework for studying growth on CO2 and H2, energy metabolism, carbon-source utilization, and the metabolic demands of recombinant protein production. A protein allocation model (PAM) extends the reconstruction with enzyme capacity and proteome allocation constraints.

## Reconstruction

The reconstruction followed five main stages:

1. **Genome annotation.** Bakta, eggNOG, RAST, and Prokka annotations were combined with DeepEC Transformer predictions to improve functional assignments and enzyme coverage.
2. **Draft reconstruction.** CarveMe generated an initial network containing 2,258 reactions and 1,534 metabolites.
3. **Manual curation.** Gene associations, reaction stoichiometry, and mass and charge balance were reviewed. Curation addressed the Calvin cycle, hydrogen oxidation and the electron transport chain, and pathways supplying biomass precursors. Unsupported exchange reactions were removed, and the biomass equation incorporated measured SoF1 composition and reference-model information.
4. **Calibration and validation.** ATP maintenance requirements were calibrated against chemostat measurements, and carbon-source growth predictions were evaluated using Biolog phenotype data. MEMOTE reports document model quality at successive reconstruction stages.
5. **Protein allocation.** Enzyme turnover numbers, protein masses, and active, translational, and unused protein sectors were integrated into a PAM to investigate metabolic capacity and growth limitations.

The repository supports flux balance and variability analysis, gene deletion studies, gas-uptake phenotype phase planes, and recombinant protein production analyses. The manuscript identifies respiratory efficiency and absolute gas-uptake predictions as areas for further refinement.

## Models

Models are provided in SBML format with gene-protein-reaction associations and cytosolic, periplasmic, and extracellular compartments.

| File in `Models/` | Description | Reactions | Metabolites | Genes |
| --- | --- | ---: | ---: | ---: |
| `260302_iABA974.sbml` | Most recent dated base-model snapshot | 2,389 | 1,614 | 974 |
| `250924_iABA974_BLG.sbml` | Beta-lactoglobulin production variant | 2,392 | 1,616 | 974 |

Counts above were checked directly against the SBML files. The totals in the manuscript's Table 2 match the BLG variant; archived snapshots can differ. Earlier models and reaction/metabolite spreadsheets retain the reconstruction's development record.

## Repository structure

| Path | Contents |
| --- | --- |
| [`Models/`](Models/) | Draft and curated SBML models, production variants, and model information spreadsheets. |
| [`data/PAM_data/`](data/PAM_data/) | Enzyme parameter workbooks for protein allocation modeling. |
| [`data/Escher_data/`](data/Escher_data/) | Escher pathway maps, a JSON model, and flux overlays. |
| [`Scripts/`](Scripts/) | Jupyter notebooks for model curation, FVA, strain design, PAM parametrization, and visualization; shared utilities in `jupyter_utils.py` and `model_utils.py`. |
| [`Scripts/Python_scripts/`](Scripts/Python_scripts/) | Scripts for maintenance-energy optimization, phase-plane analysis, PAM construction, and annotation processing. |
| [`Results/`](Results/) | Saved simulation tables and figures, including FVA, gene knockouts, maintenance fitting, and phenotype phase planes. |
| [`Memote/`](Memote/) | Archived HTML model-quality reports. |

The `data/` folder retains only PAM and Escher resources. Experimental measurements and annotation datasets have been removed from this folder. Historical notebooks that load those inputs require them to be supplied separately. Some scripts also contain machine-specific paths or use `Data/` instead of `data/`; adjust these before running, particularly on case-sensitive systems. PAM workflows additionally require PAModelpy and, for parameter fitting, PAMparametrizer.

## Quick start

Install COBRApy in your Python environment:

```bash
python -m pip install cobra
```

From the repository root, load the base model and solve it using its stored objective and constraints:

```python
from cobra.io import read_sbml_model

model = read_sbml_model("Models/260302_iABA974.sbml")
solution = model.optimize()
print("Solver status:", solution.status)
print("Objective value:", solution.objective_value)
```

For a specific growth experiment, configure the medium, gas uptake bounds, and ATP maintenance parameters before interpreting the prediction. The stored model is a starting point; this example does not reproduce the manuscript's calibrated experiments.
