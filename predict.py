import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import model as m
import scripts as scr


# This function computes the z scores, based on euclidean distance, for a set of curves
def create_z_score_dict(curve_dict, metric='euclidean', ex_distance=False, prot_mean_std_dict=False,
                        curve_points=[36.9, 40.2, 43.9, 46.6, 48.6, 52.7, 55.3, 58.5, 61.2, 64]):
    keys = curve_dict.keys()
    dataframe = pd.DataFrame(index=keys, columns=curve_points)
    curve_points_len = len(curve_points)
    for key in keys:
        values = curve_dict[key]
        if len(values) > curve_points_len:
            values = values[0:curve_points_len]
        while len(values) < curve_points_len:
            values = np.append(values, values[-1])
        dataframe.loc[key] = np.reshape(values, (curve_points_len))

    dataframe = dataframe.fillna(0)

    distances = pdist(dataframe.to_numpy(), metric=metric)
    dist_matrix = squareform(distances)
    distance_dataframe = pd.DataFrame(dist_matrix, index=keys, columns=keys)

    dist = distance_dataframe.to_numpy()
    keys = distance_dataframe.keys()
    dist_dict = {}
    for i in range(len(dist)):
        for j in range(len(dist)):
            if i == j:
                continue
            if (keys[i], keys[j]) in dist_dict:
                continue
            elif (keys[j], keys[i]) in dist_dict:
                continue

            if ex_distance:
                dist_dict[keys[i], keys[j]] = 1 / (dist[i][j] + 1)
            else:
                dist_dict[keys[i], keys[j]] = dist[i][j]

    if prot_mean_std_dict == False:
        del dist, keys, distance_dataframe
    else:
        prots = list(distance_dataframe.keys())
        prot_mean_std_dict = {}
        for prot in prots:
            val_list = list(distance_dataframe[prot])
            prot_mean_std_dict[prot] = (np.mean(val_list), np.std(val_list))
        del dist, keys, distance_dataframe

    z_score_dict = {}
    mean_dist = np.mean(list(dist_dict.values()))
    std_dist = np.std(list(dist_dict.values()))

    for pair, val in dist_dict.items():
        z_score_dict[pair] = (val - mean_dist) / std_dist

    if prot_mean_std_dict == False:
        return z_score_dict, dist_dict
    else:
        return z_score_dict, dist_dict, prot_mean_std_dict


# This function computes the z scores, based on euclidean distance, for multiple sets of curves
def prep_euc_dist_reps(curve_dict_list, ex_distance=True,
                       curve_points=[36.9, 40.2, 43.9, 46.6, 48.6, 52.7, 55.3, 58.5, 61.2,64]):
    z_score_dict_list = []
    for curve_dict in curve_dict_list:
        z_score_dict, dist_dict = create_z_score_dict(curve_dict, metric='euclidean', ex_distance=ex_distance,
                                                         curve_points=curve_points)
        del dist_dict

        z_score_dict_list.append(z_score_dict)

    return z_score_dict_list


# This function compares the predictions of different models
def compare_pair_dicts(pair_dict_1, pair_dict_2):
    score_1 = []
    score_2 = []

    for pair, val in pair_dict_1.items():
        if pair in pair_dict_2 or (pair[1], pair[0]) in pair_dict_2:
            if (pair[1], pair[0]) in pair_dict_2:
                pair_2 = (pair[1], pair[0])
            else:
                pair_2 = pair

            val_2 = pair_dict_2[pair_2]
            score_1.append(val)
            score_2.append(val_2)

    corr = np.corrcoef(x=score_1, y=score_2)[0][1]
    return corr


def model_predict(model_address, curve_dict,output_dir, save_name='tmp', model_type='B',batch_size=5000, shuffle=True,
                  prior_dict_address=None, mouse=False,
                  pfam=False, properties=False, x_=[36.9, 40.2, 43.9, 46.6, 48.6, 52.7, 55.3, 58.5, 61.2, 64]):

    #temp_dir = os.path.join(output_dir,"temp")
    if model_type.upper() == 'B':
        model = m.Base_Model(model_address=model_address)
        #print("Before predict Base model")
        model.predict(save_name,output_dir,curve_dict,batch_size=batch_size, shuffle=shuffle, x_=x_)

    else:
        model = m.Extended_Model(model_address=model_address)
        #print("Before predict Full model")
        if properties == True:
            if mouse == True:
                properties_dict_address_list = ['./data/model_features/mouse_strain_C57BL6J_master_properties_dict_1',
                                                './data/model_features/mouse_strain_C57BL6J_master_properties_dict_2']
            else:
                properties_dict_address_list = ['./data/model_features/master_properties_dict_1',
                                                './data/model_features/master_properties_dict_2',
                                                './data/model_features/master_properties_dict_3']
        else:
            properties_dict_address_list = None

        if pfam == True:
            if mouse == True:
                pfam_dict_address_list = ['./data/model_features/mouse_10090_pfam.tsv.gz']
            else:
                pfam_dict_address_list = ['./data/model_features/human_9606_pfam.tsv.gz',
                      './data/model_features/hsv1_10299_pfam.tsv.gz',
                      './data/model_features/hcmv_10360_pfam.tsv.gz',
                      './data/model_features/kshv_868565_pfam.tsv.gz']
        else:
            pfam_dict_address_list = None

        model.predict(save_name, curve_dict,output_dir,batch_size=batch_size,
                      prior_dict_address=prior_dict_address, pfam_dict_address_list=pfam_dict_address_list,
                      properties_dict_address=properties_dict_address_list, shuffle=shuffle, x_=x_)

    return save_name


# This function predicts protien-protein interactions
def tapioca_predict(base_address, tissue_address, curve_dict, curve_points, model_type, output_dir):

    print(f"Processing {base_address} {tissue_address} ")

    batch_size = 100000
    if int((len(curve_dict)*len(curve_dict)) / 2) < batch_size:
        batch_size = int((len(curve_dict)*len(curve_dict)) / 10)

    tapioca_base = './data/models/tapioca_base_submodel'
    tapioca_prop = './data/models/tapioca_prop_submodel'
    tapioca_pfam = './data/models/tapioca_pfam_submodel'
    tapioca_prior = './data/models/tapioca_prior_submodel'
    tapioca_prop_pfam = './data/models/tapioca_prop_pfam_submodel'
    tapioca_prior_pfam = './data/models/tapioca_prior_pfam_submodel'
    tapioca_prior_prop = './data/models/tapioca_prior_prop_submodel'
    tapioca_prior_prop_pfam = './data/models/tapioca_prior_prop_pfam_submodel'

    curve_points = [float(val) for val in curve_points]

 
    if model_type == 'B':
        models = './data/models/tapioca_base_submodel'
        output_pred_dir = output_dir
        print(f"models: {models}")
        print(f"output_pred_dir: {output_pred_dir}")
        output_filename = base_address+".tsv"
        model_predict(models, curve_dict, output_pred_dir,save_name=output_filename, model_type=model_type, batch_size=batch_size, shuffle=True,
                     pfam=True, properties=True, prior_dict_address=tissue_address, x_=curve_points)
        print(f"====FINAL output at {os.path.join(output_pred_dir,output_filename)}")
        return 'Done'


    #full model  model_type == 'F':
    #prediction will be in temp directory
    output_pred_dir = os.path.join(output_dir,"temp")
    temp_base_address = os.path.join(output_pred_dir,base_address)
    if tissue_address != '':
        pred_dict_address_list = [temp_base_address + '_base.tsv',
                                temp_base_address + '_prop.tsv',
                                temp_base_address + '_pfam.tsv',
                                temp_base_address + '_prior.tsv',
                                temp_base_address + '_prop_pfam.tsv',
                                temp_base_address + '_prior_pfam.tsv',
                                temp_base_address + '_prior_prop.tsv',
                                temp_base_address + '_prior_prop_pfam.tsv']


        models = [tapioca_base,tapioca_prop,tapioca_pfam,tapioca_prior,
                tapioca_prop_pfam,tapioca_prior_pfam,tapioca_prior_prop,
                tapioca_prior_prop_pfam]

    else:
        pred_dict_address_list = [temp_base_address + '_base.tsv',
                                temp_base_address + '_prop.tsv',
                                temp_base_address + '_pfam.tsv',
                                temp_base_address + '_prop_pfam.tsv']

        models = [tapioca_base, tapioca_prop, tapioca_pfam, tapioca_prop_pfam]

    

    print(f"models: {models}")
    print(f"output_pred_dir: {output_pred_dir}")
    model_predict(models, curve_dict, output_pred_dir,save_name=base_address, model_type=model_type, batch_size=batch_size, shuffle=True,
                     pfam=True, properties=True, prior_dict_address=tissue_address, x_=curve_points)

    #
    # Merge the results from temp predictions into final output directory
    #
    dict_list = []
    for pred_dict_address in pred_dict_address_list:
        dict_list.append(scr.tsv_to_dict(pred_dict_address))

    zs_dict = prep_euc_dist_reps([curve_dict], ex_distance=True, curve_points=curve_points)[0]
    corr_list = []
    for pred_dict in dict_list:
        corr_list.append(compare_pair_dicts(zs_dict, pred_dict))

    del zs_dict
    weighted_dict = {}
    for pair, val_0 in dict_list[0].items():
        val_list = [val_0]
        for pred_dict in dict_list[1:]:
            if pair in pred_dict:
                val_list.append(pred_dict[pair])
            elif (pair[1], pair[0]) in pred_dict:
                val_list.append(pred_dict[pair[1], pair[0]])
            else:
                val_list.append(0)

        weighted_val = 0
        for i, corr in enumerate(corr_list):
            weighted_val = weighted_val + corr * val_list[i]

        weighted_dict[pair] = weighted_val / np.sum(corr_list)

#    save_name = './predictions/' + base_address +".tsv"
    output = os.path.join(output_dir, base_address +".tsv")
    scr.save_dict_to_tsv_file(weighted_dict, output_filename=output)
    print(f"====FINAL output at {output}")

    for file in pred_dict_address_list:
        print("\tnot Removing "+file)
        #os.remove(file)

    return 'Done'




