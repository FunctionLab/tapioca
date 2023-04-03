# Tapioca

---

## Contents

1. [What is Tapioca?](#what_is_tapioca)
2. [The Tapioca Website](#tapioca_website)
3. [Running Tapioca Locally](#local_usage)
4. [Citing Tapioca](#citation)


## <a name="what_is_tapioca"></a> What is Tapioca?
Tapioca, an integrative machine learning based computational pipeline for PPI networks in dynamics contexts. Tapioca integrates mass spectrometry (MS) curve data from thermal proximity coaggregation (TPCA), co-fractionation (CF-MS), or ion-based proteome-integrated solubility alteration (I-PISA), with protein physical properties, domains (PFAM, now [InterPro](https://www.ebi.ac.uk/interpro/)), and tissue-specific functional networks ([HumanBase](https://hb.flatironinstitute.org/)). Tapioca itself consists of eight distinct logistic regression models that take unique combinations of these features as inputs. The predictions of these eight models are combined into a final interaction score for each pair of proteins in a given dataset. 


## <a name="tapioca_website"></a> The Tapioca Website
The [Tapioca website](https://tapioca.princeton.edu/) contains Tapioca scores for protein-protein interactions for thousands of proteins detected throughout cell culture infections with herpes simplex virus 1 (HSV-1; Justice et al. 2021), human cytomegalovirus (HCMV; Hashimoto et al. 2020) and Kaposi’s sarcoma associated herpesvirus (KSHV; unpublished). Individual or groups of proteins can be searched within these datasets on the Tapioca website, and the data can be downloaded.


## <a name="local_usage"></a> Running Tapioca Locally
To run Tapioca on your computer you will need both this repository and [tissue specific functional networks](https://hb.flatironinstitute.org/download) from HumanBase (files are slightly too large to be uploaded to GitHub).
The following subsections describe how to install and set up Tapioca locally, and use it to predict protein-protein interactions.


### Installation and Set Up

To install Tapioca locally you should first clone this repository as follows.


Before cloning the repository, make sure you are inside the directory in which you would like the repository to be placed.

```
git clone https://github.com/FunctionLab/tapioca.git
```

To navigate into the Tapioca folder use the command:

```
cd tapioca 
```

We recommend managing Tapioca's dependencies with a [conda environment](https://www.anaconda.com/products/distribution) as follows:

```
conda env create -f tapioca_env.yml 
```

To activate the tapioca environment use the command:

```
conda activate Tapioca
```



### Downloading [Tissue Specific Functional Networks from HumanBase](https://hb.flatironinstitute.org/download)

Now you will need to download tissue specific functional networks [tissue specific functional networks](https://hb.flatironinstitute.org/download) from HumanBase. On the linked HumanBase webpage you will see a list of networks
for different tissues. Simply download click the download button in the Top edges column for which ever network(s) you wish to use. If you are unsure which one to use, we recommend starting with the Global network. Be sure to place the downloaded files in the /data/tissue_functional_networks/ directory so that Tapioca can access them. Tapioca is capable of running without these networks. If you'd like to do so, then you can skip this text.



### Usage overview

There are two main ways to run Tapioca. In either case, you will need to set a few variables such that Tapioca knows where your data is located, how to read it, and how you want Tapioca to run. These variables are:

#### input_file
This is the name of the csv file contain the data which Tapioca will run on. This cdv file should be located in the 'raw_input' folder. To see how to format your input file, please see the example.csv file in the raw_input folder.

#### ref_channel
If you included a reference sample, and which this sample to be used in normalization, please set this variable to True. If you have no idea what a reference sample is, you should probably set this variable to False.

#### pre_normalized
If you have already normalized your data, then this variable should be set to True, otherwise it should be False. If you are using co-fractionation data, we recommend normalizing your data first before inputing it into Tapioca. If this is set to True, be sure to place your data in the normalized_input folder, and see the example_normalized.csv file in the normalized_input folder for how to set up your file.

#### co_fractionation
If you are using co_fractionation data, this variable should be set to True.

#### conditions
This is a list of strings with the conditions you would like Tapioca to run on. The conditions you put in the list should be written exactly (including capitalization and spaces) as it is written in your input csv file. You do not need to list all the
conditions in your input csv file, only the ones you currently want Tapioca to run on.

#### replicates
This is a list of strings with the replicates you would like Tapioca to run on. The replicates you put in the list should be written exactly (including capitalization and spaces) as it is written in your input csv file. You do not need to list all the
replicates in your input csv file, only the ones you currently want Tapioca to run on.

#### full_model
If you would like to run the full Tapioca model (all 8 sub-models and weighted integration) then set this variable to True. If you would like to run only the base model, which considers only the data in your input file and nothing else, set this variable to False.

#### tissue
Specific the file name of the tissue specific functional network you would like Tapioca to use. Remember you must download these networks from HumanBase yourself (see above). If you would like to not use any networks, set this variable to ''. Tapioca can only use one network at a time. If you have conditions that should be run in different tissue contexts then run Tapioca multiple times with different tissues, including only the conditions relevant to that tissue (in the conditions variable).

#### base_save_name
Set the name that you would like to have your output (Tapioca predictions) named as. The predictions will appear in the predictions folder.


These variables need to be set in 1 of 2 places, depending on how you would like to run Tapioca. If you would like to run Tapioca by running the main.py file, these variables should be set in the config.py file. If you would like to run Tapioca using jupyter notebook, these variables can be set within the Tapioca_Notebook.ipynb file. To open this file, first run jupyter notebooks by typing the following command in your conda terminal:

```
jupyter notebook 
```

The jupyter notebook server should then open in your internet browser. Simply navigate to the tapioca folder and open the Tapioca_Notebook.ipynb file


## <a name="citation"></a> Citing Tapioca
Tapioca has recently been submitted for publication.
