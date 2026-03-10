import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# ----------------------------
# Global plotting style
# ----------------------------
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.linewidth": 1.0,
})

# ----------------------------
# Load model data
# ----------------------------
df = pd.read_csv("Results/PhPP/production_envelope_basemodel_h2_o2.csv")

obj_col = "flux_maximum"
x_col = "EX_h2_e"
y_col = "EX_o2_e"

# ----------------------------
# Load experimental data
# ----------------------------
exp_path = "Data/240716_exp_gas_consumption.xlsx"
exp = pd.read_excel(exp_path)

exp = exp[["dilRate", "H2 Flux theoretic", "O2 flux theoretic", "CO2 flux theoretic"]].copy()

# group every 3 rows = one condition with 3 replicates
exp_grouped = exp.groupby(exp.index // 3).mean(numeric_only=True)
exp_std = exp.groupby(exp.index // 3).std(numeric_only=True)

exp_grouped["dilRate_std"] = exp_std["dilRate"]
exp_grouped["h2_std"] = exp_std["H2 Flux theoretic"]
exp_grouped["o2_std"] = exp_std["O2 flux theoretic"]
exp_grouped["co2_std"] = exp_std["CO2 flux theoretic"]

# Experimental x/y for plotting
x_col_exp = "O2 flux theoretic"
y_col_exp = "H2 Flux theoretic"

exp_x = exp_grouped[x_col_exp].astype(float).to_numpy()
exp_y = exp_grouped[y_col_exp].astype(float).to_numpy()
exp_xerr = exp_grouped["o2_std"].astype(float).to_numpy()
exp_yerr = exp_grouped["h2_std"].astype(float).to_numpy()

# ----------------------------
# Experimental regression
# ----------------------------
exp_model = LinearRegression()
exp_model.fit(exp_x.reshape(-1, 1), exp_y)

m_exp = float(exp_model.coef_[0])
b_exp = float(exp_model.intercept_)
y_exp_pred = exp_model.predict(exp_x.reshape(-1, 1))
r2_exp_fit = r2_score(exp_y, y_exp_pred)

# ----------------------------
# Convert model uptake to positive values for plotting
# ----------------------------
df["x_plot"] = -df[x_col].astype(float)   # H2 uptake
df["y_plot"] = -df[y_col].astype(float)   # O2 uptake

df_feas = df.dropna(subset=["x_plot", "y_plot", obj_col]).copy()

# ----------------------------
# Build line of optimality
# For each H2 uptake, keep max growth;
# if tied, keep lowest O2 uptake
# ----------------------------
df_feas["x_round"] = df_feas["x_plot"].round(10)
df_feas["obj_round"] = df_feas[obj_col].round(5)

opt_points = (
    df_feas
    .assign(max_obj=df_feas.groupby("x_round")["obj_round"].transform("max"))
    .query("obj_round == max_obj")
    .sort_values(["x_round", "y_plot"], ascending=[True, True])
    .groupby("x_round", as_index=False)
    .first()
    .sort_values("x_round")
)

x_opt = opt_points["x_plot"].to_numpy()   # H2 uptake
y_opt = opt_points["y_plot"].to_numpy()   # O2 uptake

opt_model = LinearRegression()
opt_model.fit(x_opt.reshape(-1, 1), y_opt)

m_opt = float(opt_model.coef_[0])
b_opt = float(opt_model.intercept_)

# ----------------------------
# Build heatmap grid
# columns = O2, index = H2
# ----------------------------
pivot = df.pivot(
    index="EX_h2_e",
    columns="EX_o2_e",
    values=obj_col
)

X = -pivot.columns.values   # O2 uptake positive
Y = -pivot.index.values     # H2 uptake positive
Z = pivot.values

# ----------------------------
# Plot
# ----------------------------
fig, ax = plt.subplots(figsize=(6.3, 5.4))

mesh = ax.pcolormesh(
    X, Y, Z,
    shading="auto",
    cmap="viridis",
    edgecolors=(0, 0, 0, 0.35),   # softer than pure black
    linewidth=0.25,
    antialiased=True
)

# Line of optimality
# plot x = O2, y = H2
h2_line = np.linspace(x_opt.min(), x_opt.max(), 300)
o2_line_opt = m_opt * h2_line + b_opt

ax.plot(
    o2_line_opt, h2_line,
    color="white",
    linewidth=3.2,
    alpha=0.95,
    zorder=4
)
ax.plot(
    o2_line_opt, h2_line,
    color="crimson",
    linewidth=2.0,
    zorder=5,
    label="Line of optimality"
)


# Experimental error bars
ax.errorbar(
    exp_x, exp_y,
    xerr=exp_xerr,
    yerr=exp_yerr,
    fmt="none",
    ecolor="black",
    elinewidth=0.9,
    capsize=2,
    alpha=0.9,
    zorder=6
)

# Experimental points
ax.scatter(
    exp_x, exp_y,
    color="black",
    s=34,
    zorder=7,
    label="Experimental data"
)

# Experimental regression
o2_line_exp = np.linspace(exp_x.min(), exp_x.max(), 300)
h2_line_exp = m_exp * o2_line_exp + b_exp

ax.plot(
    o2_line_exp, h2_line_exp,
    linestyle="--",
    color="black",
    linewidth=1.8,
    zorder=6,
    label="Experimental fit"
)

# Colorbar
cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
cbar.set_label("Growth rate (h$^{-1}$)")

# Labels
ax.set_xlabel("O$_2$ uptake (mmol gDW$^{-1}$ h$^{-1}$)")
ax.set_ylabel("H$_2$ uptake (mmol gDW$^{-1}$ h$^{-1}$)")

# Clean frame
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Compact annotation box
annotation = (
    f"Optimality line:\n"
    f"O$_2$ = {m_opt:.3f}·H$_2$ {b_opt:+.3f}\n\n"
    f"Experimental fit:\n"
    f"H$_2$ = {m_exp:.3f}·O$_2$ {b_exp:+.3f}\n"
    f"R$^2$ = {r2_exp_fit:.3f}"
)

ax.text(
    0.03, 0.97, annotation,
    transform=ax.transAxes,
    va="bottom", ha="right",
    fontsize=6,
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", alpha=0.85)
)

# Legend
ax.legend(frameon=False, loc="lower right")

plt.tight_layout()

# Save for paper
# plt.savefig("Results/PhPP/PhPP_heatmap_optimality_exp.pdf", bbox_inches="tight")
# plt.savefig("Results/PhPP/PhPP_heatmap_optimality_exp.png", bbox_inches="tight")

plt.show()