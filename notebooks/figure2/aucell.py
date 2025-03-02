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
import pandas as pd
import scanpy as sc
import decoupler as dc
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

import sys

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=120,
    facecolor="white",
    frameon=True,
)
# %matplotlib inline

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data/"
FIGURE_PATH = "/home/sychen9584/projects/cardio_paper/figures"

# %%
adata = sc.read_h5ad(os.path.join(DATA_PATH, "processed/scRNA_all.h5ad"))

# %%
# convert to human genes symbols by capitalizing the gene names
adata.var['mouse_genesymbol'] = adata.var_names
adata.var_names = adata.var_names.str.upper()

# %%
sasp_genes = [
  "Acvr1b", "Ang", "Angpt1", "Angptl4", "Areg", "Axl", "Bex3", 
  "Bmp2", "Bmp6", "C3", "Ccl1", "Ccl2", "Ccl20", "Ccl24", "Ccl26", 
  "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Cd55", "Cd9", "Csf1", 
  "Csf2", "Csf2rb", "Cst10", "Ctnnb1", "Ctsb", "Cxcl1", "Cxcl10", 
  "Cxcl12", "Cxcl16", "Cxcl2", "Cxcl3", "Cxcr2", "Dkk1", "Edn1", 
  "Egf", "Egfr", "Ereg", "Esm1", "Ets2", "Fas", "Fgf1", "Fgf2", 
  "Fgf7", "Gdf15", "Gem", "Gmfg", "Hgf", "Hmgb1", "Icam1", "Icam5", 
  "Igf1", "Igfbp1", "Igfbp2", "Igfbp3", "Igfbp4", "Igfbp5", "Igfbp6", 
  "Igfbp7", "Il10", "Il13", "Il15", "Il18", "Il1a", "Il1b", "Il2", 
  "Il6", "Il6st", "Il7", "Inha", "Iqgap2", "Itga2", "Itpka", "Jun", 
  "Kitl", "Lcp1", "Mif", "Mmp13", "Mmp10", "Mmp12", "Mmp13", 
  "Mmp14", "Mmp2", "Mmp3", "Mmp9", "Nap1l4", "Nrg1", "Pappa", 
  "Pecam1", "Pgf", "Pigf", "Plat", "Plau", "Plaur", "Ptbp1", "Ptger2", 
  "Ptges", "Rps6ka5", "Scamp4", "Selplg", "Sema3f", "Serpinb3a", 
  "Serpine1", "Serpine2", "Spp1", "Spx", "Timp2", "Tnf", "Tnfrsf11b", 
  "Tnfrsf1a", "Tnfrsf1b", "Tubgcp2", "Vegfa", "Vegfc", "Vgf", "Wnt16", "Wnt2"
]

sasp_genes = [g.upper() for g in sasp_genes if g.upper() in adata.var_names]

# %%
cardiac_fibrosis_genes = [
  "TGFB1", "MMP9", "SERPINE1", "FN1", "EDN1", "ACE2", "AGER", "COL1A1", 
  "GJA1", "GHSR", "FBN1", "MMP7", "PLEKHA7", "NFKBIA", "MFN2", "GHRL", 
  "DMD", "LGALS3", "AGT", "CCN2", "MIR21", "SMAD3", "ST2", "STAT3", 
  "TNF", "ACE", "SMAD2", "ACTB", "VPS51", "ANGPT2", "REN", "VEGFA", 
  "NLRP3", "MMP2", "FGF23", "TRPM7", "NR3C2", "MIR29A", "SPP1", "PLG", 
  "NFE2L2", "MIR155", "MIR143", "GAS5", "MAPK3", "IL6", "SIRT3", "SIRT1", 
  "CCL2", "SLC5A2", "TLR2", "BRD4", "PTEN", "TGM2", "FGF21", "TACR1", 
  "SOX9", "RNF19A", "MIR29B2", "MIR29B1", "HSD11B1", "GABPA", "NOX4", 
  "CYBB", "CYTL1", "DPP4", "APLN", "MAPK14", "PARP1", "CRK", "IL33", 
  "PLAU", "MIR378A", "AHSA1", "WWP2", "UTS2", "PIN1", "POSTN", "GRAP2", 
  "RNR2", "PPARD", "F2RL1", "AIMP2", "PDGFD", "DNMT3A", "EGFR", "POLDIP2", 
  "MAPK1", "NLRC5", "PPARG", "DECR1", "AGTR1"
]

cardiac_fibrosis_genes = [g for g in cardiac_fibrosis_genes if g in adata.var_names]

# %%
heart_failure_genes = [
  "PPARGC1A", "ADRB1", "GRK2", "ACE", "AGT", "AGTR1", "EDN1", "HIF1A", 
  "IL6", "NOS3", "NPPA", "NPPB", "ATP2A2", "TNF", "CCL2", "HMOX1", 
  "SOD2", "MSTN", "ADRB3", "HTR2B", "RAC1", "CRP", "CCN2", "EDNRA", 
  "ALB", "GCG", "NRG1", "IL1B", "NR3C2", "PIK3CG", "AVP", "PTH", 
  "REN", "VEGFA", "ADIPOQ", "GDF15", "SIRT1", "PPARG", "PPP1R1A", 
  "ADRA2C", "CYBB", "NFE2L2", "PTGS2", "NOX4", "SLC9A1", "TNFRSF1A", 
  "UCN2", "NOS2", "SERPINE1", "PON1", "PRL", "EDNRB", "VWF", "PEBP1", 
  "RETN", "TLR2", "CS", "CSF2", "EPHX2", "GSK3B", "IFNG", "CXCL8", 
  "NPR1", "OLR1", "POMC", "AVPR2", "RBP4", "XDH", "CSF3", "CFD", 
  "GATM", "ACACA", "ACADS", "SOD3", "TNNT2", "UCP1", "NRIP1", "CAT", 
  "HAND2", "ROCK2", "FBLN5", "DSTN", "CIDEA", "ADM", "ADRB2", "FASN", 
  "BAMBI", "NOX1", "GPX4", "CXCL2", "APCS", "HSPB1", "APOC1", "IGF1", 
  "INS", "ITGB1", "LGALS3", "MME", "MMP9", "ACLY", "ATP1A3", "ATP2A1", 
  "PCK1", "PDGFRA", "PDPK1", "GHRL", "PLAT", "PLN", "PRKAR2B", "MAP2K7", 
  "TBX20", "PTGS1", "ACE2", "RYR2", "SCN5A", "SOD1", "SOX4", "TP53", 
  "ELOVL6", "CSRP3", "CALR", "FIP1L1", "CASP3", "KAT8", "APLN", "PRKCA", 
  "SLC8A1", "IL10", "MTPN", "ERBB2", "IL17A", "MMP2", "POSTN", "SIRT3", 
  "NR3C1", "MYH6", "MYH7", "CASQ2", "HSF1", "MMP1", "PRKCB", "TIMP1", 
  "AGTR2", "ALOX15", "GNAQ", "NOS1", "NPY", "TNFRSF11B", "SCD", "BNIP3", 
  "ADAM17", "PARP1", "ENG", "HMGB1", "HP", "MAS1", "NCAM1", "OPA1", 
  "P2RX4", "FXYD1", "BCL2", "TH", "TIMP2", "UCP3", "CORIN", "CX3CR1", 
  "ECE1", "ERBB4", "HTR2A", "IL6R", "KCNE1", "KCNK2", "ARRB1", "PGF", 
  "PPP1CC", "BAX", "SRI", "TIMP4", "DDR1", "DYNLL1", "CASP12", "RAMP2", 
  "CHRM4", "CHRNA4", "CYP11B1", "DIO3", "ACAN", "DNMT1", "ENDOG", "TUBB", 
  "FKBP1B", "NPTXR", "FOS", "ALOX12", "FTH1", "BIRC5", "FAS", "FASLG", 
  "LIF", "MDH2", "MMP8", "MMP13", "MMP14", "MYOD1", "NGFR", "PLCD1", 
  "PLD2", "MED1", "RNLS", "GRK1", "SCN8A", "SELP", "UCP2", "FOSL1", 
  "OGT", "BECN1", "NOS1AP", "TOMM70", "CDKN2B-AS1", "SLC39A8", "TCF7L2", 
  "FTO", "ZGLP1", "C8orf37-AS1", "CERT1", "LINC01782", "LIPC-AS1", 
  "CXCR6", "CMTM7", "MAPK14", "CST3", "ADRA1A", "ADRA2B", "OTUD7A", 
  "DMD", "DNAH8", "DPP4", "DSPP", "APLNR", "CELSR2", "SIK3", "COPD", 
  "GCKR", "GH1", "GJA1", "AMPD1", "GLP1R", "ANGPT1", "ANGPT2", "GPR42", 
  "GRK5", "GRHL1", "IL1A", "ITPK1", "LCN2", "LIPC", "FADS1", "LMNA", 
  "MIR21", "MYBPC3", "ZFHX3", "ATP2B1", "MIR423", "POLK", "PIK3CA", 
  "PIK3CB", "PIK3CD", "PLCG1", "SEMA5B", "PPARA", "CHDH", "MAML3", 
  "ZCCHC8", "ANKS1B", "ACKR3", "CPNE5", "SUGP1", "REG1A", "RFC1", 
  "ACTB", "S100A1", "SLC5A2", "SLC6A8", "SSTR4", "ST2", "BRS3", 
  "TGFB1", "TLR4", "TM7SF2", "TNNI3", "TTN", "TTR", "VPS51", "BEST1", 
  "CA5A", "SLC30A3", "FSD1", "EHMT1", "CALCR", "NUBPL", "FGF23", 
  "C1orf21", "FSD1L", "ALDH1A2", "ZPR1", "BAZ1B", "LPAR2", "CYTH3", 
  "ADAMTSL1", "FADS2", "BAG3", "CELSR1"
]

heart_failure_genes = [g for g in heart_failure_genes if g in adata.var_names]

# %%
inflammation_genes = ["LACRT", "EGF", "CCK", "GALNS", "SCGB1A1", "AGTR1", "VIP", "IDO1", 
                        "CNTF", "EHF", "PON1", "SERPINE1", "C3AR1", "RETN", "PRKAA1", "DLL4",
                        "MUC1", "IFNL3", "AKT1", "SERPINF1", "PRSS3", "S100A12", "ELANE", 
                        "PRKCD", "TSLP", "TREML2", "COL2A1", "SULT2B1", "BIN1", "SERPING1", 
                        "HAMP", "CRYAB", "PADI1", "TFRC", "PDGFB", "PDGFA", "F11R", "PLA2G6",
                        "GATA3", "HIF1A", "C3", "C5", "PRDX5", "HMOX1", "TPSAB1", "CCL26",
                        "ANXA1", "CCL21", "CCL20", "FN1", "SIRT1", "CD200R1", "MMP10", "PIGF",
                        "ADA", "SERPINA3", "RET", "SERPINA1", "MEFV", "FPR2", "ADM", "ADIPOR1",
                        "ADIPOR2", "TIMP4", "TIMP1", "SERPINB4", "SERPINB3", "PPARA", "AHSG",
                        "STAT3", "NR1H4", "F2", "F7", "PPARG", "PLP1", "PLCB3", "APCS", "CCL13",
                        "GRN", "CCL11", "CDKN1A", "SETD7", "PTGS2", "TMSB4X", "NAMPT", "UCN3",
                        "GAST", "CCL19", "CCL18", "JAK2", "PTPN1", "CD163", "F11", "TFPI2", 
                        "DEFA3", "DEFA1", "AGT", "FOXP3", "CXCL10", "CXCL11", "CXCL12", 
                        "S100A7", "PF4", "SLC28A3", "SLC28A2", "FKBP4", "S100A9", "EPO", "CSF2",
                        "CSF1", "RASL11B", "AHR", "IL1RAP", "CLU", "TNF", "CXCL14", "CXCL16", 
                        "ZFP36", "NTF3", "CDH1", "NTF4", "SPP1", "SH2B3", "SLC15A1", "F2R",
                        "NFKB1", "PGF", "MYD88", "BIRC5", "LCN2", "PGR", "CPB2", "MAZ", 
                        "EIF4EBP1", "PHB", "RORC", "FGF1", "FASLG", "FGF2", "CYP2C19", 
                        "CYP19A1", "VTN", "FGF7", "TNFRSF9", "ITGA1", "HSP90AA1", "ANGPT2",
                        "ANGPT1", "BDNF", "TRPA1", "NRG1", "APOA5", "APOA4", "TNFRSF1A", 
                        "FCGR2A", "APOA1", "PYDC1", "PECAM1", "DSG2", "KCNK3", "NFE2L2", 
                        "GTF2A1", "KCNK9", "TREM1", "KNG1", "FCGR3A", "DUSP10", "TNFRSF17", 
                        "ACP5", "TNFRSF18", "VWF", "GZMA", "NGF", "KLF2", "ADCYAP1", "TXNIP", 
                        "CRH", "ANG", "BPI", "CRP", "FLT1", "WDR1", "PTEN", "ITGAL", "FCAR", 
                        "TNFSF13B", "THBD", "FCGR1A", "CCL5", "CCL4", "CCL3", "CCL2", "CCL1", 
                        "ATP7B", "MBL2", "APLN", "MAPK14", "IL23A", "LEP", "TP53", "CD200", 
                        "APP", "AQP7", "NOD1", "AQP8", "NOD2", "ITIH2", "NTN1", "AQP3", "AQP4",
                        "C4A", "STS", "C4B", "BSG", "FLT3LG", "EGR1", "TNFSF18", "TNFSF15", 
                        "TNFSF13", "TNFSF14", "PLA2G4D", "MIF", "IGF1", "TP73", "HNF1A", 
                        "SPRR1B", "PTGES", "CFTR", "SLC22A5", "IL1RN", "ICAM1", "CRHR2", 
                        "MT2A", "IL36A", "APOM", "TNFSF12", "TREX1", "IL12B", "TNFSF10", 
                        "IL12A", "APOD", "CTSB", "EDN1", "IL10", "IL15", "NCOA6", "TNFRSF12A", 
                        "IL13", "SFTPD", "VEGFC", "LIF", "LILRB1", "IL18", "CYP7B1", "VEGFA", 
                        "VEGFB", "IL1A", "TLR2", "TF", "IL1B", "BST1", "CEACAM1", "CD28", 
                        "TLR9", "TLR7", "TLR6", "TLR4", "TLR3", "CD40", "PTGDR2", "KLK1", 
                        "STC1", "TGFA", "TNFRSF11B", "ADRB3", "LIAS", "CD3E", "CX3CL1", 
                        "SOCS3", "SOCS1", "C1QBP", "CD38", "CCR8", "CCR7", "CTSG", "CCR6", 
                        "ASIC3", "ASIC1", "CAMP", "ASIC2", "ADIPOQ", "IL33", "IL32", "MME", 
                        "GCH1", "PARP1", "GAL", "IL17A", "NLRP10", "F2RL1", "CD48", "SAA1", 
                        "GHRL", "CD47", "CD40LG", "CD44", "IL22", "IL20", "SELPLG", "TNFRSF13C", 
                        "TNFRSF13B", "IL27", "CTNNB1", "IL25", "PTGER4", "PTN", "CXORF40A", 
                        "CST3", "CXCR3", "DMBT1", "ENPP7", "TGFB1", "CHRNB2", "PDGFRB", "IL4R", 
                        "TSC2", "TRPV1", "TYROBP", "VCAM1", "SELENOS", "FAS", "CD68", "PLAA", 
                        "DAP3", "RELA", "CHRNA4", "CD83", "MPO", "SCUBE1", "SCUBE2", "EFNB1", 
                        "TBX21", "NLRP3", "SMAD3", "HGF", "MIR22", "NAT1", "SELE", "SMAD7", 
                        "PROC", "SELP", "PROCR", "CARD17", "PRTN3", "C5AR1", "CALCA", "CACNA1C", 
                        "PLAU", "PLAT", "AGER", "NPPC", "TGM2", "NPPB", "CD99", "MMP1", 
                        "LIMK1", "MMP2", "MMP3", "HSPA4", "MMP9", "ADA2", "MIR217", "GH1", 
                        "BMP4", "IRF1", "IL6", "ELF3", "IFNG", "ALPL", "IRF7", "UCN", "CHI3L1", 
                        "DEFB4A", "CXCL8", "ABCB4", "LTB4R", "YWHAG", "ABCB1", "NCAM1", 
                        "MIR34A", "TNFAIP6", "TNFAIP3", "HMGB1", "CXCL2", "CXCL1", "OXT", 
                        "ADAR", "CXCL6", "HRH1", "BDKRB1", "BDKRB2", "MADCAM1", "CASP1", 
                        "TAC1", "YWHAH", "NOS3", "PLAUR", "S100B", "SOD1", "POMC", "CAT", 
                        "ACKR1", "MASP2", "ACKR2", "PTX3", "LTF", "CLEC12A", "THBS1", "NR3C2", 
                        "HSPD1", "NR3C1", "OASL", "CLN6", "DDT", "CNR2", "CNR1", "KYNU", "NOS2"]

inflammation_genes = [g for g in inflammation_genes if g in adata.var_names]

# %%
geneset_dict = {
    "SASP": sasp_genes,
    "Cardiac Fibrosis": cardiac_fibrosis_genes,
    "Heart Failure": heart_failure_genes,
    "Inflammation": inflammation_genes
}
geneset_df = pd.concat(pd.DataFrame({'geneset':k, 'gene':v}) for k, v in geneset_dict.items())
geneset_df = geneset_df.drop_duplicates(subset=["geneset", "gene"])

# %%
dc.run_aucell(mat=adata, net=geneset_df, source="geneset", target="gene", verbose=True, use_raw=False)

# %%
adata.obsm['aucell_estimate']

# %%
adata.obs['SASP_score'] = adata.obsm['aucell_estimate']["SASP"]

# %%
geneset_scores = adata.obsm['aucell_estimate']
geneset_scores['cell_type'] = adata.obs['cell_type']
geneset_scores['cell_type_fine'] = adata.obs['cell_type_fine']
geneset_scores[['month', 'sample']] = adata.obs['sample'].str.extract(r'(m\d+)_(s\d+)')

# %%
geneset_scores

# %%
geneset_scores.to_csv(os.path.join(DATA_PATH, "geneset_scores.csv"))

# %%
adata.write_h5ad(os.path.join(DATA_PATH, "processed/scRNA_all.h5ad"), compression="gzip")

# %% [markdown]
# ## Statistical significance

# %%
geneset_scores = pd.read_csv(os.path.join(DATA_PATH, "geneset_scores.csv"), index_col=0)

# %%
geneset_scores

# %% [markdown]
# ### SASP for Fib.1 Cxcl1

# %%
sasp_fib1 = geneset_scores[geneset_scores['cell_type_fine'] == 'Fib.1']
# Perform One-Way ANOVA
model = ols('SASP ~ C(month)', data=sasp_fib1).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print("One-Way ANOVA Results:\n", anova_table)
# Perform Tukey's HSD Test
tukey = pairwise_tukeyhsd(sasp_fib1['SASP'], sasp_fib1['month'], alpha=0.05)
print("\nTukey's HSD Test Results:\n", tukey)

# %% [markdown]
# ### SASP for Fib.6 Erbb4

# %%
sasp_fib6 = geneset_scores[geneset_scores['cell_type_fine'] == 'Fib.6']
# Perform One-Way ANOVA
model = ols('SASP ~ C(month)', data=sasp_fib1).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print("One-Way ANOVA Results:\n", anova_table)
# Perform Tukey's HSD Test
tukey = pairwise_tukeyhsd(sasp_fib6['SASP'], sasp_fib6['month'], alpha=0.05)
print("\nTukey's HSD Test Results:\n", tukey)

# %% [markdown]
# ### SASP for all clusters

# %%
# Perform One-Way ANOVA
model = ols('SASP ~ C(month)', data=geneset_scores).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print("One-Way ANOVA Results:\n", anova_table)
# Perform Tukey's HSD Test
tukey = pairwise_tukeyhsd(geneset_scores['SASP'], geneset_scores['month'], alpha=0.05)
print("\nTukey's HSD Test Results:\n", tukey)

# %% [markdown]
# ### SASP for MC.Cd209a

# %%
sasp_mc = geneset_scores[geneset_scores['cell_type_fine'] == 'MC.2']
# Perform One-Way ANOVA
model = ols('SASP ~ C(month)', data=sasp_mc).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print("One-Way ANOVA Results:\n", anova_table)
# Perform Tukey's HSD Test
tukey = pairwise_tukeyhsd(sasp_mc['SASP'], sasp_mc['month'], alpha=0.05)
print("\nTukey's HSD Test Results:\n", tukey)

# %% [markdown]
# ### SASP for neutrophil

# %%
sasp_neutrophil = geneset_scores[geneset_scores['cell_type_fine'] == 'Granulocyte/Neutrophil']
# Perform One-Way ANOVA
model = ols('SASP ~ C(month)', data=sasp_neutrophil).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print("One-Way ANOVA Results:\n", anova_table)
# Perform Tukey's HSD Test
tukey = pairwise_tukeyhsd(sasp_neutrophil['SASP'], sasp_neutrophil['month'], alpha=0.05)
print("\nTukey's HSD Test Results:\n", tukey)

# %% [markdown]
# ### Cardiac Fibrosis for all clusters

# %%
geneset_scores.rename(columns={"Cardiac Fibrosis": "Cardiac_Fibrosis", "Heart Failure": "Heart_Failure"}, inplace=True)

# %%
# Perform One-Way ANOVA
model = ols('Cardiac_Fibrosis ~ C(month)', data=geneset_scores).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print("One-Way ANOVA Results:\n", anova_table)
# Perform Tukey's HSD Test
tukey = pairwise_tukeyhsd(geneset_scores['Cardiac_Fibrosis'], geneset_scores['month'], alpha=0.05)
print("\nTukey's HSD Test Results:\n", tukey)

# %% [markdown]
# ### Heart Failure for all clusters

# %%
# Perform One-Way ANOVA
model = ols('Heart_Failure ~ C(month)', data=geneset_scores).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print("One-Way ANOVA Results:\n", anova_table)
# Perform Tukey's HSD Test
tukey = pairwise_tukeyhsd(geneset_scores['Heart_Failure'], geneset_scores['month'], alpha=0.05)
print("\nTukey's HSD Test Results:\n", tukey)

# %% [markdown]
# ### Inflammation for all clusters

# %%
# Perform One-Way ANOVA
model = ols('Inflammation ~ C(month)', data=geneset_scores).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print("One-Way ANOVA Results:\n", anova_table)
# Perform Tukey's HSD Test
tukey = pairwise_tukeyhsd(geneset_scores['Inflammation'], geneset_scores['month'], alpha=0.05)
print("\nTukey's HSD Test Results:\n", tukey)

# %%
