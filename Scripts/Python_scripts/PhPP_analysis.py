import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# ----------------------------
# 1) Load production envelope (this comes from the build in fuction of cobrapy)
# ----------------------------
envelope_path = "Results/PhPP/production_envelope_basemodel_h2_o2.csv"  
df = pd.read_csv(envelope_path)

# Pick the 2 axes (these are the two varied reactions in your CSV)
x_col = "EX_h2_e"
y_col = "EX_o2_e"

# IMPORTANT: exchange reactions are often negative for uptake.
# If you want "uptake rates" as positive numbers on the plot, flip the sign here.
USE_POSITIVE_UPTAKE = True
if USE_POSITIVE_UPTAKE:
    df["x_plot"] = -df[x_col].astype(float)
    df["y_plot"] = -df[y_col].astype(float)
else:
    df["x_plot"] = df[x_col].astype(float)
    df["y_plot"] = df[y_col].astype(float)

# Objective columns produced by cobrapy.production_envelope
# flux_maximum is typically the maximum objective at that (x,y) grid point
obj_col = "flux_maximum"

# Keep only feasible points (some grids can have NaN)
df_feas = df.dropna(subset=["x_plot", "y_plot", obj_col]).copy()

# ----------------------------
# 2) Build "line of optimality"
# ----------------------------
# A common definition in a 2D envelope grid:
# For each x value, pick the y giving the maximum objective (flux_maximum).
# (You can also do it per y, depending on how you define optimality.)
df_feas["x_round"] = df_feas["x_plot"].round(10)
df_feas["obj_round"] = df_feas[obj_col].round(5)
opt_points = (
    df_feas
      # compute the max growth at each x
      .assign(max_obj=df_feas.groupby("x_round")["obj_round"].transform("max"))
      # keep only rows that achieve that max (ties allowed)
      .query(f"obj_round == max_obj")
      # among ties, pick the lowest O2 uptake (smallest y_plot)
      .sort_values(["x_round", "y_plot"], ascending=[True, True])
      .groupby("x_round", as_index=False)
      .first()
      .sort_values("x_round")
)
print(opt_points)
x_opt = opt_points["x_plot"].to_numpy()
y_opt = opt_points["y_plot"].to_numpy()

# Fit a straight line to the optimality locus: y = m_opt*x + b_opt
opt_model = LinearRegression()
opt_model.fit(x_opt.reshape(-1, 1), y_opt)
m_opt = float(opt_model.coef_[0])
b_opt = float(opt_model.intercept_)


# ----------------------------
# 3) Load experimental points
# ----------------------------

# Put your file path here:
exp_path = 'Data/240716_exp_gas_consumption.xlsx'
exp = pd.read_excel(exp_path)

x_col_exp = "H2 Flux theoretic"
y_col_exp = "O2 flux theoretic"
if x_col_exp in exp.columns and y_col_exp in exp.columns:
    exp_x_raw = exp[x_col_exp].astype(float).to_numpy()
    exp_y_raw = exp[y_col_exp].astype(float).to_numpy()
else:
    raise ValueError(
        f"Experimental CSV must contain either ({x_col},{y_col}) or (x,y) columns."
    )



# ----------------------------
# 4) Linear regression on experimental points + R²
# ----------------------------
exp_model = LinearRegression()
exp_model.fit(exp_x_raw.reshape(-1, 1), exp_y_raw)
m_exp = float(exp_model.coef_[0])
b_exp = float(exp_model.intercept_)

y_exp_pred = exp_model.predict(exp_x_raw.reshape(-1, 1))
r2_exp_fit = r2_score(exp_y_raw, y_exp_pred)

# ----------------------------
# 5) R² of experimental points vs optimality line
# ----------------------------
y_opt_pred_at_exp = opt_model.predict(exp_y_raw.reshape(-1, 1))


# ----------------------------
# 6) Plot: envelope + optimality + exp points + exp regression
# ----------------------------
fig, ax = plt.subplots(figsize=(7.5, 6))

# Envelope "cloud" (grid points)
ax.scatter(
    df_feas["x_plot"], df_feas["y_plot"],
    s=10, alpha=0.25, label="Production envelope grid"
)

# Optimality line (from fitted optimal points)
x_line = np.linspace(df_feas["x_plot"].min(), df_feas["x_plot"].max(), 200)
y_line_opt = m_opt * x_line + b_opt
ax.plot(x_line, y_line_opt, linewidth=2, label=f"Line of optimality: y={m_opt:.3f}x+{b_opt:.3f}")

# Optional: show the actual optimal points used to fit the line
ax.scatter(x_opt, y_opt, s=20, alpha=0.9, label="Optimality locus points")

# Experimental points
ax.scatter(exp_x_raw, exp_y_raw, s=45, marker="o", label="Experimental points")

# Experimental regression line
y_line_exp = m_exp * x_line + b_exp
ax.plot(
    x_line, y_line_exp, linestyle="--", linewidth=2,
    label=f"Exp regression: y={m_exp:.3f}x+{b_exp:.3f} (R²={r2_exp_fit:.3f})"
)



ax.set_xlabel(f"{x_col} ({'uptake (+)' if USE_POSITIVE_UPTAKE else 'raw'})")
ax.set_ylabel(f"{y_col} ({'uptake (+)' if USE_POSITIVE_UPTAKE else 'raw'})")
ax.set_title("2D Production Envelope with Optimality Line + Experimental Fit")
ax.legend()
ax.grid(True, alpha=0.25)

plt.savefig(f"Results/PhPP/{x_col}_{y_col}_PhPP_2d_GEM.png")
plt.tight_layout()
plt.show()



