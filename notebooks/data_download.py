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

# %% [markdown]
# # Data Download from GEO

# %%
import GEOparse

# path to the raw data directory
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data/raw/"

# %%
gse = GEOparse.get_GEO(geo="GSE277702", destdir=DATA_PATH)

# %%
gse.show_metadata()

# %% [markdown]
# ## Download expression matrices

# %%
gse.download_supplementary_files(directory=DATA_PATH)

# %% [markdown]
# ## Wrangle directories and file names for easier use in the future

# %%
import os
import re

def mod_GEO_directory_name(orig_name: str):
    match = re.match(r".*__([0-9]+)_month_+(sc[A-Z]+)__sample_([0-9]+)", orig_name)
    if match:
        sc_exp = match.group(2)
        months = "m" + match.group(1)
        sample = "s" + match.group(3)
        return f"{sc_exp}_{months}_{sample}"
    
    return None # return None if the format doesn't match


# %%
for folder in os.listdir(DATA_PATH):
    folder_path = os.path.join(DATA_PATH, folder)
    if os.path.isdir(folder_path):
        new_name = mod_GEO_directory_name(folder)
        
        if new_name:
            new_folder_path = os.path.join(DATA_PATH, new_name)
            os.rename(folder_path, new_folder_path)
            print(f"Renamed: {folder} -> {new_name}")
        else:
            print(f"Skipped: {folder} (does not match expected pattern)")

# %%
# strip GSM prefix from individual data files within each sample directory
#prefix_pattern = r'^GSM\d+_\d+_(ATAC|GEX)_'
prefix_pattern = r'^features'

for folder in os.listdir(DATA_PATH):
    folder_path = os.path.join(DATA_PATH, folder)
    if os.path.isdir(folder_path):
        for filename in os.listdir(folder_path):
            new_name = re.sub(prefix_pattern, "barcodes", filename)
            
            if new_name != filename:
                old_path = os.path.join(folder_path, filename)
                new_path = os.path.join(folder_path, new_name)
                os.rename(old_path, new_path)
                print(f"Renamed: {filename} -> {new_name}")
            else:
                print(f"Skipped: {filename} (does not match pattern)")
