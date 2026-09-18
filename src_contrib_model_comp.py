import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Load the data
data_dir = "/Users/johanna/Uni/ws_25_26/ecohydrology/ecohydro_project/output"
kuehnhammer_fname = "contribution_data2.csv"
mixsire_fname = "summary_models.csv"

method2 = pd.read_csv(f"{data_dir}/{kuehnhammer_fname}")
method1 = pd.read_csv(f"{data_dir}/{mixsire_fname}")

# Process method2 data
data2 = method2[["tree_code", "contribution", "level", "Date", "species"]].copy()
data2 = data2[data2["Date"] != "2026-02-24"]
data2["tree_id"] = data2["tree_code"].str[1:2].astype(int)
data2["Method"] = "Kuehnhammer"
data2["SD"] = 0
data2["contribution"] = data2["contribution"] * 100
data2["contribution"] = np.abs(data2["contribution"])
data2 = data2.drop(columns=["tree_code"])
data2.loc[data2["species"] == "Scots Pine", "species"] = "Pine"

# Process method1 data
data1 = method1[["level", "treeid", "Mean", "SD", "species", "Date", "source"]].copy()
data1 = data1[(data1["Date"] != "2026-02-24") & (data1["source"] == "Label")]
data1["tree_id"] = data1["treeid"].str[1:2].astype(int)
data1["Method"] = "MixSIAR"
data1["Mean"] = data1["Mean"] * 100
data1["SD"] = data1["SD"] * 100
data1 = data1.rename(columns={"Mean": "contribution"})
data1 = data1.drop(columns=["source", "treeid"])

# Combine the data
plot_data = pd.concat([data1, data2], ignore_index=True)

# Ensure correct ordering of levels
plot_data["level"] = pd.Categorical(
    plot_data["level"],
    categories=["low", "mid", "high"],
    ordered=True
)
# Convert date for nicer facet labels
plot_data["Date"] = pd.to_datetime(plot_data["Date"]).dt.strftime("%Y-%m-%d")
# Create hierarchical hue label
plot_data["Method_tree"] = plot_data["Method"] + " | Tree " + plot_data["tree_id"].astype(str)

# color palettes
#mix_cols = {"low": "#7f5539", "mid": "#dda15e", "high": "#606c38"}
#k_cols = {"low": "#C2B3A6", "mid": "#EAD5BA", "high": "#B9BCA7"}

mix_cols = {"low": "#8B636C", "mid": "#CD8500", "high": "#CDAF95"}
k_cols = {"low": "#C1B1B5", "mid": "#E0C396", "high": "#E3D7CB"}

species = plot_data["species"].unique()
dates = sorted(plot_data["Date"].unique())
trees = sorted(plot_data["tree_id"].unique())

level_order = ["low", "mid", "high"]
methods = ["MixSIAR", "Kuehnhammer"]

plot_data["level"] = pd.Categorical(plot_data["level"], categories=level_order, ordered=True)

# spacing parameters
level_spacing = 0.28      # distance between level groups
method_spacing = 0.125     # distance between methods (smaller)
bar_width = 0.12

sns.set_theme(style="whitegrid", rc={"axes.grid": True,  # Ensure grid is enabled
                                     "grid.linestyle": "-",  # Default for both axes
                                     "axes.grid.axis": "y",  # Only show grid lines for the y-axis (horizontal lines)
                                     })

fig, axes = plt.subplots(len(species), len(dates), figsize=(10,5), sharey=True)

for i, sp in enumerate(species):
    for j, date in enumerate(dates):

        ax = axes[i, j]
        subset = plot_data[(plot_data.species == sp) & (plot_data.Date == date)]

        base_x = np.arange(len(trees))

        for li, level in enumerate(level_order):

            level_offset = (li - 1) * level_spacing

            for mi, method in enumerate(methods):

                method_offset = (mi - 0.5) * method_spacing

                xs = base_x + level_offset + method_offset

                values = []
                sds = []
                for t in trees:
                    row = subset[
                        (subset.tree_id == t) &
                        (subset.level == level) &
                        (subset.Method == method)
                        ]
                    v = row.contribution
                    sd = row.SD

                    values.append(v.values[0] if len(v) else 0)
                    sds.append(sd.values[0] if len(sd) else 0)

                color = mix_cols[level] if method == "MixSIAR" else k_cols[level]

                ax.bar(
                    xs,
                    values,
                    width=bar_width,
                    color=color,
                    linewidth=0.6,
                    label=f"{level} | {method}"
                )
                if method == "MixSIAR":
                    print(sds)
                    ax.errorbar(xs, values, yerr=sds, fmt=".", color="black", capsize=3)

        if i == len(species) - 1:
            ax.set_xticks(base_x)
            ax.set_xticklabels([f"Tree {t}" for t in trees])
        else:
            ax.set_xticks(base_x)
            ax.set_xticklabels([f"" for t in trees])

        ax.tick_params(
            axis="x",
            which="major",
            color="lightgrey",
            bottom=True,
            top=False,
            length=4,  # tick length
            width=1,
            direction="out"
        )

        if i == 0:
            ax.set_title(date, weight='bold', fontsize=14)

        if j == 0:
            ax.set_ylabel("Contribution [%]")

# Add species labels on the right side of each row
for i, sp in enumerate(species):
    # Get vertical center of the row
    pos = axes[i, -1].get_position()
    y_center = (pos.y0 + pos.y1) / 2

    fig.text(
        0.98, y_center + 0.03 if i == 0 else y_center, sp,
        rotation=-90,
        va="center",
        ha="right",
        weight="bold",
        fontsize=14
    )

# clean legend
handles, labels = ax.get_legend_handles_labels()
unique = dict(zip(labels, handles))
fig.legend(unique.values(), unique.keys(), title="Level | Method", loc="center",
           ncol=3, bbox_to_anchor=(0.5, 0.43))

plt.tight_layout(rect=[0,0,0.97,1])
sns.despine(left=True, right=True, top=True)

plt.show()

plt.savefig(f"{data_dir}/contr_model_comp_v4.png", dpi=700)


round(method1.groupby(["species", "level", "source"])["Mean"].mean() * 100, 2)
data2

sns.barplot(data=method1[(method1["Date"] == "2026-02-26") | (method1["Date"] == "2026-02-27")],
            x="source", y="Mean", hue="species")

# --- Only plot one day and species

fig, ax = plt.subplots(figsize=(5,4))

sp = "Douglas fir"
date = '2026-02-27'
j = 1

subset = plot_data[(plot_data.species == sp) & (plot_data.Date == date)]

base_x = np.arange(len(trees))

for li, level in enumerate(level_order):

    level_offset = (li - 1) * level_spacing

    method_offset = (0 - 0.5) * method_spacing

    xs = base_x + level_offset + method_offset

    values = []
    sds = []
    for t in trees:
        row = subset[
            (subset.tree_id == t) &
            (subset.level == level) &
            (subset.Method == 'MixSIAR')
            ]
        v = row.contribution
        sd = row.SD

        values.append(v.values[0] if len(v) else 0)
        sds.append(sd.values[0] if len(sd) else 0)

    color = mix_cols[level]

    ax.bar(
        xs,
        values,
        width=bar_width * 2.2,
        color=color,
        linewidth=0.6,
        label=f"{level}"
    )
    ax.errorbar(xs, values, yerr=sds, fmt=".", color="black", capsize=3)

ax.set_xticks(base_x)
ax.set_xticklabels([f"Tree {t}" for t in trees], fontsize=14)
ax.tick_params(axis='y', labelsize=13)

ax.tick_params(
    axis="x",
    which="major",
    color="lightgrey",
    bottom=True,
    top=False,
    length=4,  # tick length
    width=1,
    direction="out"
)
ax.set_title(f"{sp} - {date}", weight='bold', fontsize=14)
ax.set_ylabel("Label contribution [%]", fontsize=14)

# clean legend
handles, labels = ax.get_legend_handles_labels()
unique = dict(zip(labels, handles))
fig.legend(unique.values(), unique.keys(), title="", loc="center",
           ncol=3, bbox_to_anchor=(0.52, 0.05))

plt.tight_layout(rect=[0,0.05,0.97,1])
sns.despine(left=True, right=True, top=True)

plt.savefig(f"{data_dir}/contr_mixsiar.png", dpi=700)