import tapioca as tp
import config as c



"""
Tapioca is executed when this file is run. Before running Tapioca, go to the config.py
file to set the input file and other parameters
"""
tp.run_tapioca(
input_file = c.input_file,
ref_channel = c.ref_channel,
pre_normalized = c.pre_normalized,
co_fractionation = c.co_fractionation,
conditions = c.conditions,
replicates = c.replicates,
tissue = c.tissue,
base_save_name = c.base_save_name,
full_model = c.full_model
)