# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.6
#   kernelspec:
#     display_name: cardio
#     language: python
#     name: python3
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from matplotlib import pyplot as plt
import PyComplexHeatmap as pch
import sys

sys.path.append('../../scripts')
import figure_functions as fig_func

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=120,
    facecolor="white",
    frameon=True,
)

plt.rcParams["axes.grid"] = False  # Disable grids for all plots
plt.rcParams["grid.color"] = "white"  # Ensure grids are fully invisible
plt.rcParams["grid.alpha"] = 0  # Remove any transparency effects
plt.grid(False)  # Turn off grids explicitly

# %matplotlib inline

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data/"
FIGURE_PATH = "/home/sychen9584/projects/cardio_paper/figures"

# %%
adata = sc.read_h5ad(os.path.join(DATA_PATH, "processed/scRNA_all.h5ad"))
adata.obs['cell_type_fine'] = adata.obs['cell_type_fine'].replace({
    'Granulocyte/Neutrophil': "Neutrophil",
    'Peri/Smooth Muscle': 'Smooth Muscle'
})

# %%
scenic_loom_input = os.path.join(DATA_PATH, "scenic/input.loom")
adj_output = os.path.join(DATA_PATH, "scenic/adj.tsv")

# %%
sc.pp.filter_genes(adata, min_cells=265) # expressed in at least 1% of the cells
adata.var.set_index("mouse_genesymbol", inplace=True)

# %%
# create basic row and column attributes for the loom file:
row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}
lp.create(scenic_loom_input, adata.X.transpose(), row_attrs, col_attrs)

# %% [markdown]
# # After running SCENIC on AWS, load in resulting regulon scores

# %%
# collect SCENIC AUCell output
lf = lp.connect(os.path.join(DATA_PATH, "scenic/output.loom"), mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)

# Extract gene names
genes = lf.ra["Gene"]
# Convert the binary regulon matrix to a DataFrame
regulons = pd.DataFrame(lf.ra["Regulons"], index=genes)
# Rename columns using TF names (from dtype field names)
regulons.columns = [col.replace("(+)", "") for col in lf.ra["Regulons"].dtype.names]

regulon_dict = {
    tf: regulons.index[regulons[tf] == 1].tolist()
    for tf in regulons.columns
}

lf.close()

# %%
auc_mtx['cell_type_fine'] = adata.obs['cell_type_fine'].copy()
regulon_celltype = auc_mtx.groupby('cell_type_fine').mean()

# %% [markdown]
# ## Figure 5A

# %%
tfs_to_plot = ['Bclaf1', "Ikzf2", "Ikzf3", "Bach2", "Irf9", "Irf7", "Zfp426", "Hey1", "Sox17", "Smad1",
               "Hltf", "Ets2", "Pparg", "Etv6", "Cux1", "Rara", "Nfkb2", "Zfp384", "Fosl1", "Pbx3", 
               "E2f1", "Maf", "Nfatc2", "Tcf7l2", "Zfp112", "Meis1", "Nr1h4", "Gli3", "Myc", "Jun", "Rela",
               "Mafg", "Hey2", "Zfp467", "Etv1", "Sox6", "Etv4", "Arntl", "Klf2", "Klf4"]

regulon_celltype.columns = regulon_celltype.columns.str.replace(r"\(\+\)", "", regex=True)
# Filter columns where the TF name appears before '('
filtered_columns = [col for col in regulon_celltype.columns if any(tf in col for tf in tfs_to_plot)]
regulon_celltype = regulon_celltype[filtered_columns]

# %%
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()

# %%
regulon_celltype_scaled = pd.DataFrame(scaler.fit_transform(regulon_celltype), columns=regulon_celltype.columns, index=regulon_celltype.index)

# %%
annot_df = pd.DataFrame(regulon_celltype_scaled.index)
annot_df["cell_type_fine"] = annot_df["cell_type_fine"].astype(str)
annot_df.index = annot_df["cell_type_fine"]
# Get colors
cell_type_colors, cell_type_fine_colors = fig_func.get_celltype_colors()  

# %%
with plt.rc_context({"figure.dpi": 150, "figure.frameon": True, 
                     "axes.labelsize": 10, "axes.titlesize": 10,"xtick.labelsize": 10, 
                    "ytick.labelsize": 10, "legend.fontsize": 8, "legend.title_fontsize": 8, "figure.titlesize": 16}):

    plt.figure(figsize=(6, 7))
    col_ha = pch.HeatmapAnnotation(
        Celltype=pch.anno_simple(annot_df.cell_type_fine, colors=cell_type_fine_colors, legend=True),  
        axis=1,  # Annotation on the x-axis (columns = cells)
        verbose=1, plot_legend=True, legend_kws={"fontsize": 8}
    )
    cm = pch.ClusterMapPlotter(
        data=regulon_celltype_scaled.T, 
        top_annotation=col_ha,
        row_cluster=True, col_cluster=False, row_dendrogram=False, col_dendrogram=False,
        vmin=-2, vmax=4, show_rownames=True, cmap="coolwarm", label="Regulon \n Activity")
    cm.ax.figure.savefig(os.path.join(FIGURE_PATH, "figure5/figure5a.png"), dpi=150, bbox_inches='tight', format='png')
    plt.show()

# %% [markdown]
# ## Figure 5B

# %%
adata.obs[['month', 'sample_num']] = adata.obs['sample'].str.extract(r'(m\d+)_(s\d+)')
auc_mtx['month'] = adata.obs['month'].copy()
regulon_celltype_month = auc_mtx.groupby(['cell_type_fine','month']).mean()

# %%
tfs_to_plot = ["Bclaf1", "Dbp", "Thra", "Hlf", "Rarb", "Klf13", "Stat1", "Irf9", "Irf7", "Elf2", "Fosl2", "Atf6"]
regulon_celltype_month.columns = regulon_celltype_month.columns.str.replace(r"\(\+\)", "", regex=True)
# Filter columns where the TF name appears before '('
filtered_columns = [col for col in regulon_celltype_month.columns if any(tf in col for tf in tfs_to_plot)]
regulon_celltype_month = regulon_celltype_month[filtered_columns]

# %%
regulon_celltype_month_scaled = pd.DataFrame(scaler.fit_transform(regulon_celltype_month), columns=regulon_celltype_month.columns, index=regulon_celltype_month.index)

# %%
# reorder months
month_order = ["m3", "m12", "m24"]
regulon_celltype_month_scaled = regulon_celltype_month_scaled.reindex(month_order, level=1)

# %%
annot_df = regulon_celltype_month_scaled.index.to_frame()
annot_df["cell_type_fine"] = annot_df["cell_type_fine"].astype(str)

# %%
with plt.rc_context({"figure.dpi": 150, "figure.frameon": True, 
                     "axes.labelsize": 10, "axes.titlesize": 10,"xtick.labelsize": 10, 
                    "ytick.labelsize": 10, "legend.fontsize": 8, "legend.title_fontsize": 8, "figure.titlesize": 16}):

    plt.figure(figsize=(5, 6))
    col_ha = pch.HeatmapAnnotation(
        Celltype=pch.anno_simple(annot_df.cell_type_fine, colors=cell_type_fine_colors, legend=True), 
        Timepoint=pch.anno_simple(annot_df.month, colors={"m3": "#f69697", "m12": "#9dd189", "m24": "#c6b0d6"}, legend = True),        
        axis=1,  # Annotation on the x-axis (columns = cells)
        verbose=1, plot_legend=True, legend_kws={"fontsize": 8}
    )
    cm = pch.ClusterMapPlotter(
        data=regulon_celltype_month_scaled.T, 
        top_annotation=col_ha,
        row_cluster=True, col_cluster=False, row_dendrogram=False, col_dendrogram=False,
        vmin=-2, vmax=4, show_rownames=True, cmap="coolwarm", label="Regulon \n Activity")
    cm.ax.figure.savefig(os.path.join(FIGURE_PATH, "figure5/figure5b.png"), dpi=150, bbox_inches='tight', format='png')
    plt.show()

# %% [markdown]
# ## Linear Mixed Model

# %%
import statsmodels.api as sm
import statsmodels.formula.api as smf

# %%
auc_mtx.columns = auc_mtx.columns.str.replace(r"\(\+\)", "", regex=True)
auc_mtx['sample_num'] = adata.obs['sample_num'].copy()

# %%
auc_long = auc_mtx.melt(id_vars=['cell_type_fine', 'month', 'sample_num'], var_name='TF', value_name='TF_activity')

# %%
auc_long['month'] = pd.Categorical(auc_long["month"], categories=["m3", "m12", "m24"], ordered=True)

# %%
import statsmodels.formula.api as smf

# Store results in a dictionary
results = []

for tf in auc_long["TF"].unique():
    sub_df = auc_long[auc_long["TF"] == tf]  # Filter data for the TF
    
    try:
        # Fit the LMM model
        model = smf.mixedlm("TF_activity ~ month", sub_df, groups=sub_df["cell_type_fine"], 
                            vc_formula={"sample_num": "0 + C(sample_num)"})
        result = model.fit()
        
        # Store results
        results.append({
            "TF": tf,
            "pvalue_12m": result.pvalues.get("month[T.m12]", float("nan")),  
            "pvalue_24m": result.pvalues.get("month[T.m24]", float("nan")), 
            "Converged": getattr(result, "converged", False)
        })
    
    except Exception as e:
        print(f"Error fitting LMM for {tf}: {e}")
        results.append({"TF": tf, "pvalue_12m": float("nan"), "pvalue_24m": float("nan"), "Converged": False})

# Convert to DataFrame
results_df = pd.DataFrame(results)

# %%
results_df['adj_pval_12m'] = sm.stats.multipletests(results_df['pvalue_12m'], method='fdr_bh')[1]
results_df['adj_pval_24m'] = sm.stats.multipletests(results_df['pvalue_24m'], method='fdr_bh')[1]

# %%
results_df.to_csv(os.path.join(DATA_PATH, "scenic/lmm_results.csv"), index=False)

# %%
auc_mtx.to_csv(os.path.join(DATA_PATH, "scenic/tf_auc_scores.csv"), index=True)

# %% [markdown]
# ## Figure 5C and D

# %%
auc_mtx = pd.read_csv(os.path.join(DATA_PATH, "scenic/tf_auc_scores.csv"), index_col=0)

# %%
fib_scores = auc_mtx[auc_mtx['cell_type_fine'].str.startswith("Fib")]
endo_scores = auc_mtx[auc_mtx['cell_type_fine'].str.startswith("Endo")]

# %%
fib_scores = fib_scores[['Stat1', "Irf7", "Atf4", "Yy1", "cell_type_fine", "month"]]
fib_scores_long = fib_scores.melt(id_vars=['cell_type_fine', 'month'], var_name='TF', value_name='TF_activity')

# %%
import json
with open(os.path.join(DATA_PATH, "scenic/regulons.json"), "w") as f:
    json.dump(regulon_dict, f, indent=4)

# %%
cell_type_colors, cell_type_fine_colors = fig_func.get_celltype_colors()

# %%
fig_func.tf_activity_lineplot(fib_scores_long, regulon_dict, 'Stat1', significance=[('m3', 'm12', 0.001), ('m3', 'm24', 0.001)],
                     legend=False, hue='cell_type_fine', palette=cell_type_fine_colors)
plt.show()

# %% [markdown]
# ### Statistic test for TF activity

# %%
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# %%
# Perform One-Way ANOVA
df = fib_scores_long[fib_scores_long['TF'] == "Yy1"]

model = ols('TF_activity ~ C(month)', data = df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print("One-Way ANOVA Results:\n", anova_table)
# Perform Tukey's HSD Test
tukey = pairwise_tukeyhsd(df['TF_activity'], df['month'], alpha=0.05)
print("\nTukey's HSD Test Results:\n", tukey)

# %%
endo_scores = endo_scores[["Hlf", "Irf2", "Dbp", "Elk4", "cell_type_fine", "month"]]
endo_scores_long = endo_scores.melt(id_vars=['cell_type_fine', 'month'], var_name='TF', value_name='TF_activity')

# %%
# Perform One-Way ANOVA
df = endo_scores_long[endo_scores_long['TF'] == "Hlf"]

model = ols('TF_activity ~ C(month)', data = df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print("One-Way ANOVA Results:\n", anova_table)
# Perform Tukey's HSD Test
tukey = pairwise_tukeyhsd(df['TF_activity'], df['month'], alpha=0.05)
print("\nTukey's HSD Test Results:\n", tukey)

# %%
fig_func.tf_activity_lineplot(endo_scores_long, regulon_dict, 'Hlf', significance=[('m3', 'm12', 0.001), ('m3', 'm24', 0.001)],
                     legend=False, hue='cell_type_fine', palette=cell_type_fine_colors)
plt.show()

# %% [markdown]
# ### GO analysis

# %%
from gprofiler import GProfiler

# %%
go_dict = {k: regulon_dict[k] for k in ["Stat1", "Irf7", "Atf4", "Yy1", "Hlf", "Atf6b", "Irf2", "Dbp", "Elk4"]}

# %%
gp = GProfiler(return_dataframe=True) #return pandas dataframe or plain python structures  )
go_df = gp.profile(organism='mmusculus', user_threshold=0.05, significance_threshold_method='fdr', sources=["GO:BP", "GO:CC", "KEGG"], query=go_dict)

# %%
go_df['GeneRatio'] = go_df['intersection_size']/go_df['query_size']

# %%
go_df.query('query=="Dbp"')

# %%
highlight_idx = [1, 4, 8, 46, 101,
                 10, 11, 13, 19, 56,
                 118, 168, 236, 326, 394,
                 2095, 2096, 2128, 2241, 2097,
                 1019, 1461, 1819, 2101, 2237,
                 1693, 2327, 2328, 2329, 2353, 
                 829, 1397, 2172, 2177, 2175, 
                 347, 572, 1319, 1210, 1410]
highlight_df = go_df.loc[highlight_idx, :]

# %%
highlight_df["logP"] = -np.log10(highlight_df["p_value"])
highlight_df = highlight_df.sort_values(['query', 'GeneRatio'], ascending=[True, False])

# %%
highlight_df['name'] = [
    "reg. of metabolic proc.","reg. of biosynthestic proc.","catabolic proc.","macromolecule catabolic proc.","protein ubiquitination",
    "cytoplasm", "organelle lumen","circadian rhythm", 'rhythmic process', "transcription regulator complex", 
    "organelle lumen", "transcription by RNA pol II", "cellcular localization", "chromatin organization", "cell division",
    "rhythmic process", "circdian rhythm", "circadian behavior", "lipid homeostasis", "re. to leptin",
    "reg. of biosynthestic proc.", "reg. of gene expression", "pos. reg. of growth", "transciption regulator complex", "cohesin complex",
    "re. to biotic stimulus", "defense re. to symbiont", "re. to virus", "defense re. to virus", "reg. of viral process",
    "defense re. to symbiont", "re. to virus", "defense re. to virus", "cytokine production", "re. to interferon-beta",
    "primary metabolic proc.", "protein metabolic proc.", "macromolecule modification", "lipoic acid metabolism", 'tRNA processing'
]

# %%
highlight_df

# %%
with plt.rc_context({"figure.figsize": (18, 4), "figure.dpi": 150, "figure.frameon": True, 
                     "axes.labelsize": 10, "axes.titlesize": 10, "xtick.labelsize": 10, 
                     "ytick.labelsize": 10, "legend.fontsize": 8, "legend.title_fontsize": 10, "figure.titlesize": 16}):

    fig, axes = plt.subplots(1, 4)
    # List of TFs to plot
    tf_list = ["Stat1", "Irf7", "Atf4", "Yy1"]

    # Loop through subplots and apply function
    for ax, tf in zip(axes, tf_list):
        fig_func.plot_go_dotplot(highlight_df, tf, ax=ax)  # Call function for each TF
        
    fig.subplots_adjust(wspace=1.5, left=0.05, right=0.85)
    
    fig.savefig(os.path.join(FIGURE_PATH, "figure5/figure5c_l.png"), dpi=300, bbox_inches='tight', format='png')
    plt.show()

# %%
with plt.rc_context({"figure.figsize": (18, 4), "figure.dpi": 150, "figure.frameon": True, 
                     "axes.labelsize": 10, "axes.titlesize": 10, "xtick.labelsize": 10, 
                     "ytick.labelsize": 10, "legend.fontsize": 8, "legend.title_fontsize": 10, "figure.titlesize": 16}):

    fig, axes = plt.subplots(1, 4)
    # List of TFs to plot
    tf_list = ["Hlf", "Irf2", "Dbp", "Elk4"]

    # Loop through subplots and apply function
    for ax, tf in zip(axes, tf_list):
        fig_func.plot_go_dotplot(highlight_df, tf, ax=ax)  # Call function for each TF
        
    fig.subplots_adjust(wspace=1.5, left=0.05, right=0.85)
    
    fig.savefig(os.path.join(FIGURE_PATH, "figure5/figure5d_l.png"), dpi=300, bbox_inches='tight', format='png')
    plt.show()

# %%
highlight_df.to_csv(os.path.join(DATA_PATH, "scenic/go_analysis_highlight.csv"), index=False)

# %%
