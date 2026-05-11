#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import jupyter_utils as ju


def main():
    # Paths
    base_dir = Path.home() / "Documents" / "PhD" / "SOF1-GEM"

    # Load model
    model = ju.load_model(base_dir, "Models", "250808_iABA974.sbml")

    # Parameters
    pi_value = 0.014
    step = 0.00001
    growth_threshold = 1e-6

    # Conversion constants
    M_PI = 94.97  # g/mol
    DS = 2.424    # %
    X = DS * 10   # gDW/L

    # Data storage
    pi_value_data = []

    # Loop
    while pi_value >= 0:
        model_copy = model.copy()
        rxn = model_copy.reactions.get_by_id('EX_pi_e')
        rxn.bounds = (-pi_value, 0.0)

        solution = model_copy.optimize()

        if solution.status != "optimal" or solution.objective_value < growth_threshold:
            print(f"Stopped at PI uptake: {pi_value:.6f}, "
                  f"status: {solution.status}, growth: {solution.objective_value}")
            break

        # Convert to g/L/h
        pi_rate = pi_value * M_PI * 1e-3 * X

        print(f"PI uptake: {pi_value:.6f}, "
              f"Growth rate: {solution.objective_value:.6f}, "
              f"PI rate: {pi_rate:.6f} g/L/h")

        pi_value_data.append((pi_value, solution.objective_value, pi_rate))

        pi_value -= step

    # Save results
    df = pd.DataFrame(pi_value_data, columns=['PI Flux', 'Growth Rate', 'PI_rate_gLh'])

    output_path = base_dir / "Results" / "260511_pi_flux_growth_rate.csv"
    df.to_csv(output_path, index=False)

    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()