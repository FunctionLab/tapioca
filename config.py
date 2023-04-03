"""
*************************************************************************************
To run Tapioca on your data, you will need to adjust the parameters within this section
to make them appropriate for your specific experiment
"""

"""
Here list the name of your input csv file. Be sure to include .csv at the end of the file name.
If your file is not normalized then it should go in the 'raw_input' folder. If it is pre normalized
it should be placed in the 'normalized_input' folder
"""
input_file = 'example.csv'

"""
If you have a reference channel, then set this to true so that it will be used in normalization
"""
ref_channel = True

"""
If your data is already normalized and in the correct format (see example file in normalized folder)
then you should set this to True. If you are using prenormalized data, make sure to place it in 
the 'normalized_input' folder
"""
pre_normalized = False

"""
If you are using co-fractionation data, set this variable to True. We recommend normalizing your
co-fractionation data yourself before putting it into Tapioca
"""
co_fractionation = False

"""
Here, list the names of the conditions used as they appear in the csv file you are running on.
For example, if you wrote 'ctrl' in your csv file, then you need to put 'ctrl', not 'Ctrl',
in this list.
"""
conditions = ['ctrl','experimental_condition_1','experimental_condition_2']

"""
Here, list the names of the replicates. As with conditions, this should exactly match the
names used in the csv file.
"""
replicates = ['replicate_1','replicate_2','replicate_3']

"""
Here you can choose to use the full model, or the base model. To use the base model set this variable to False.
If you choose to use the full model, you can choose to use or not use the tissue specific functional networks in
the next step.
"""
full_model = True

"""
Here list the address of the tissue functional network that you'd like to use. A full list
of tissue-specific functional networks can be found and downloaded at: https://hb.flatironinstitute.org/download
If you download a new functional network, be sure to place it in the 'tissue_functional_network'
folder. If you would like to not include a tissue specific functional network, set the variable
to ''
"""
tissue = 'global_top.gz'
# How to not use a tissue specific functional network
# tissue = ''

"""
List the base save name you would like to use. This will be combined with the condition and replicate
names to create the file name for the output (Tapioca predictions) file.
For example a base name of 'tapioca_predictions' with condition 'ctrl' and replicate 
'replicate_1' will result in a file name of 'tapioca_predictions_ctrl_replicate_1'.
"""
base_save_name = 'example_predictions'
