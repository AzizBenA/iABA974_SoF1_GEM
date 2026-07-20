import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# ----------------------------
# Global plotting style
# ----------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.default": "regular",   # keeps subscripts/superscripts non-italic, matches body text
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "font.size": 10,
    "axes.labelsize": 12,
    "axes.labelweight": "medium",
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
})
# ----------------------------
# Load model data
# ----------------------------
df_phpp = pd.read_csv("Results/PhPP/production_envelope_basemodel_h2_o2.csv")

growth_col = "flux_maximum"
y_col = "EX_h2_e"
x_col = "EX_o2_e"

# ----------------------------
# Load experimental data
# ----------------------------
exp_path = "Data/240716_exp_gas_consumption.xlsx"
exp_raw = pd.read_excel(exp_path)

exp = exp_raw[["dilRate", "H2 Flux theoretic", "O2 flux theoretic", "CO2 flux theoretic"]].copy()


# ----------------------------
# Convert model uptake to positive values for plotting
# ----------------------------
df_phpp["x_plot"] = -df_phpp["EX_h2_e"].astype(float)   # H2 uptake
df_phpp["y_plot"] = -df_phpp["EX_o2_e"].astype(float)   # O2 uptake

df_feas = df_phpp.dropna(subset=["x_plot", "y_plot", growth_col]).copy()

def group_experimental_data(exp_df):
    # group every 3 rows = one condition with 3 replicates
    #the idea is to average the replicates and get the standard deviation for error bars
    exp_grouped = exp_df.groupby(exp_df.index // 3).mean(numeric_only=True)
    exp_std = exp_df.groupby(exp_df.index // 3).std(numeric_only=True)

    exp_grouped["dilRate_std"] = exp_std["dilRate"]
    exp_grouped["h2_std"] = exp_std["H2 Flux theoretic"]
    exp_grouped["o2_std"] = exp_std["O2 flux theoretic"]
    exp_grouped["co2_std"] = exp_std["CO2 flux theoretic"]

    return exp_grouped

def calculate_standard_error_for_plotting(exp_grouped):
    # Calculate standard error for each column
    # Experimental x/y for plotting
    x_col_exp = "O2 flux theoretic"
    y_col_exp = "H2 Flux theoretic"

    exp_x = exp_grouped[x_col_exp].astype(float).to_numpy()
    exp_y = exp_grouped[y_col_exp].astype(float).to_numpy()
    exp_xerr = exp_grouped["o2_std"].astype(float).to_numpy()
    exp_yerr = exp_grouped["h2_std"].astype(float).to_numpy()

    return exp_x, exp_y, exp_xerr, exp_yerr


def get_exp_linear_regression(exp_x, exp_y):
    # Experimental regression
    exp_model = LinearRegression()
    exp_model.fit(exp_x.reshape(-1, 1), exp_y)

    m_exp = float(exp_model.coef_[0])
    b_exp = float(exp_model.intercept_)
    y_exp_pred = exp_model.predict(exp_x.reshape(-1, 1))
    r2_exp_fit = r2_score(exp_y, y_exp_pred)

    return m_exp, b_exp, r2_exp_fit

def calculate_optimality_line(df_feas, growth_col):
    # ----------------------------
    # Build line of optimality
    # For each H2 uptake, keep max growth;
    # if tied, keep lowest O2 uptake
    # ----------------------------

    df_feas["x_round"] = df_feas["x_plot"].round(10)
    df_feas["obj_round"] = df_feas[growth_col].round(5)

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
    opt_model.fit(y_opt.reshape(-1, 1), x_opt)

    m_opt = float(opt_model.coef_[0])
    b_opt = float(opt_model.intercept_)

    return x_opt, y_opt, m_opt, b_opt

def build_heatmap_grid(df_phpp, growth_col, x_col, y_col):
    # ----------------------------
    # Build heatmap grid
    # columns = O2, index = H2
    # ----------------------------
    pivot = df_phpp.pivot(
        index=y_col,
        columns=x_col,
        values=growth_col
    )

    X_centers = -pivot.columns.values   # O2 uptake, positive
    Y_centers = -pivot.index.values     # H2 uptake, positive
    Z = pivot.values

    # pivot.columns/index are ascending in the raw (negative) flux values,
    # so negating them leaves X_centers/Y_centers descending. Sort ascending
    # and reorder Z to match before computing bin edges.
    x_sort = np.argsort(X_centers)
    y_sort = np.argsort(Y_centers)
    X_centers = X_centers[x_sort]
    Y_centers = Y_centers[y_sort]
    Z = Z[np.ix_(y_sort, x_sort)]

    def centers_to_edges(centers):
        #Convert an array of bin centers to bin edges (assumes near-uniform spacing).
        centers = np.asarray(centers, dtype=float)
        d = np.diff(centers)
        edges = np.empty(len(centers) + 1)
        edges[1:-1] = centers[:-1] + d / 2
        edges[0] = centers[0] - d[0] / 2
        edges[-1] = centers[-1] + d[-1] / 2
        return edges

    X_edges = centers_to_edges(X_centers)   # O2 uptake bin edges
    Y_edges = centers_to_edges(Y_centers)   # H2 uptake bin edges

    return X_edges, Y_edges, Z




# ----------------------------
# Main execution
# ----------------------------



grouped_exp = group_experimental_data(exp)
exp_x, exp_y, exp_xerr, exp_yerr = calculate_standard_error_for_plotting(grouped_exp)
m_exp, b_exp, r2_exp_fit = get_exp_linear_regression(exp_x, exp_y)
x_opt, y_opt, m_opt, b_opt = calculate_optimality_line(df_feas, growth_col)
X_edges, Y_edges, Z = build_heatmap_grid(df_phpp, growth_col, x_col, y_col)




# ----------------------------
# Plot
# ----------------------------
fig, ax = plt.subplots(figsize=(6.3, 5.4))
mesh = ax.pcolormesh(
    X_edges, Y_edges, Z,
    shading="flat",              # X/Y are now true cell edges -> mesh fills to the boundary
    cmap="viridis",
    edgecolors=(0, 0, 0, 0.35),  # softer than pure black
    linewidth=0.25,
    antialiased=True
)
# ----------------------------
# Full-plot grid, aligned to the same bin edges as the heatmap
# so it's continuous across both filled and blank (NaN) regions
# ----------------------------
ax.set_xticks(X_edges, minor=True)
ax.set_yticks(Y_edges, minor=True)
ax.grid(which="minor", color=(0, 0, 0, 0.35), linewidth=0.25)
ax.tick_params(which="minor", length=0)   # no visible tick marks, just gridlines
# ----------------------------
# Trim y-limits to the first/last row that actually has data
# (skip all-NaN rows so the axis starts right at the first colored cell)
# ----------------------------
row_has_data = ~np.isnan(Z).all(axis=1)   # True for rows with at least one valid cell
first_row = np.argmax(row_has_data)        # index of first data-containing row
ax.set_ylim(Y_edges[first_row],20)


# Line of optimality
# plot x = O2, y = H2
o2_line = np.linspace(y_opt.min(), y_opt.max(), 300)   # O2 range (renamed correctly)
h2_line_opt = m_opt * o2_line + b_opt                    # predicted H2 = f(O2)

ax.plot(
    o2_line, h2_line_opt,
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

# Experimental regression line (kept available, currently unused)
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
cbar = fig.colorbar(mesh, ax=ax, pad=0.02, fraction=0.046)
cbar.set_label(r"Growth rate ($\mathrm{h^{-1}}$)", labelpad=8)
cbar.ax.tick_params(labelsize=9, width=0.6)
cbar.outline.set_linewidth(0.6)

# Labels
ax.set_xlabel(r"O$_2$ uptake rate (mmol gDW$^{-1}$ h$^{-1}$)", labelpad=6)
ax.set_ylabel(r"H$_2$ uptake rate (mmol gDW$^{-1}$ h$^{-1}$)", labelpad=6)

# Closed frame: keep all four spines, uniform weight
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_color("black")
ax.tick_params(direction="out", length=3.5, width=0.5, top=False, right=False)

annotation = (
    r"$\mathbf{Line\ of\ optimality}$" + "\n"
    rf"H$_2$ = {m_opt:.2f}$\cdot$O$_2$ {b_opt:+.2f}" + "\n\n"
    r"$\mathbf{Experimental\ fit}$" + "\n"
    rf"H$_2$ = {m_exp:.2f}$\cdot$O$_2$ {b_exp:+.2f}" + "\n"
    rf"$R^2$ = {r2_exp_fit:.2f}"
)

ax.text(
    0.7, 0.5, annotation,
    transform=ax.transAxes,
    va="top", ha="left",
    fontsize=8,
    linespacing=1.6,
    bbox=dict(
        boxstyle="round,pad=0.4",
        facecolor="white",
        edgecolor=(0, 0, 0, 0.5),
        linewidth=0.6,
        alpha=0.92,
    ),
)

ax.legend(
    frameon=True,
    framealpha=0.9,
    edgecolor=(0, 0, 0, 0.3),
    fancybox=False,
    loc="lower right",
    handlelength=1.8,
    borderpad=0.6,
)

plt.tight_layout()

# Save for paper
# plt.savefig("Results/PhPP/PhPP_heatmap_optimality_exp.pdf", bbox_inches="tight")
# plt.savefig("Results/PhPP/PhPP_heatmap_optimality_exp.png", bbox_inches="tight")

plt.show()