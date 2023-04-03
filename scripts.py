import numpy as np
import pandas as pd
import math
from scipy.spatial.distance import pdist, squareform
import scipy.integrate as spint
import pickle


def convert_gene_dict(gene_dict, name_convention='g'):
    converted_gene_dict = {}
    for key, sub_dict in gene_dict.items():
        if name_convention == 'a':
            new_key = str(sub_dict['Accession'])
        elif name_convention == 'g':
            new_key = str(sub_dict['Name'])
        elif name_convention == 'name':
            new_key = str(key)
            key = str(sub_dict['Name'])
        elif name_convention == 'accession':
            new_key = str(key)
            key = str(sub_dict['Accession'])
        elif name_convention == 'na':
            new_key = sub_dict['Name']
            key = str(sub_dict['Accession'])
        elif name_convention == 'an':
            new_key = sub_dict['Accession']
            key = str(sub_dict['Name'])
        elif name_convention == 'ga':
            new_key = key
            key = str(sub_dict['Accession'])

        converted_gene_dict[str(new_key)] = str(key)

    return converted_gene_dict


def create_tissue_dict(address, name_convention='accession', relevant_prots=None):
    gene_dict_raw = load_object('./data/Biogrid/gene_dict_Homo_sapiens')
    gene_dict = convert_gene_dict(gene_dict_raw, name_convention=name_convention)
    del gene_dict_raw

    tissue_table = pd.read_table(address, header=0, names=['Gene1', 'Gene2', 'Value'])

    gene1s = tissue_table['Gene1']
    gene2s = tissue_table['Gene2']
    values = tissue_table['Value']
    tissue_dict = {}

    for i in range(len(gene1s)):
        gene1 = str(gene1s[i])
        gene2 = str(gene2s[i])
        val = values[i]
        if gene1 == gene2:
            continue
        if gene1 not in gene_dict or gene2 not in gene_dict:
            continue

        gene1 = gene_dict[gene1]
        gene2 = gene_dict[gene2]

        if gene1 == '-' or gene2 == '-':
            continue

        if (gene1, gene2) in tissue_dict or (gene2, gene1) in tissue_dict:
            continue

        if relevant_prots != None:
            if gene1 not in relevant_prots or gene2 not in relevant_prots:
                continue

        tissue_dict[gene1, gene2] = val

    return tissue_dict


def load_pfam_tsv(address):
    pfam = pd.read_table(address, low_memory=False, skiprows=2)
    column_headers = list(pfam.columns)
    fixed_column_headers = column_headers[0].split('> <')
    for i in range(len(fixed_column_headers)):
        fixed_column_headers[i] = fixed_column_headers[i].replace('#', '').replace('<', '').replace('>', '')

    pfam = pd.read_table(address, low_memory=False, skiprows=3)
    first_row = list(pfam.columns)

    pfam.columns = fixed_column_headers
    pfam.loc[len(pfam.index)] = first_row
    len(fixed_column_headers)

    return pfam


def process_pfam(address, max_evalue=None, loc_info=False):
    pfam = load_pfam_tsv(address)

    seq_ids = list(pfam['seq id'])
    hmm_accs = list(pfam['hmm acc'])
    hmm_names = list(pfam['hmm name'])
    types = list(pfam['type'])
    clans = list(pfam['clan'])
    evalues = list(pfam['E-value'])

    if loc_info == True:
        envelope_starts = list(pfam['envelope start'])
        envelope_ends = list(pfam['envelope end'])

    del pfam

    pfam_dict = {}

    for i in range(len(seq_ids)):
        seq_id = seq_ids[i]
        hmm_acc = hmm_accs[i]
        hmm_name = hmm_names[i]
        type_ = types[i]
        clan = clans[i]
        evalue = evalues[i]
        if loc_info == True:
            try:
                start =  envelope_starts[i]
                end =  envelope_ends[i]
            except:
                start = -1
                end = -1

        if max_evalue:
            if evalue > max_evalue:
                continue

        if seq_id not in pfam_dict:
            pfam_dict[seq_id] = {}

        if type_ not in pfam_dict[seq_id]:
            pfam_dict[seq_id][type_] = []

        if loc_info == True:
            pfam_dict[seq_id][type_].append((hmm_acc, hmm_name, clan, start, end))
        else:
            pfam_dict[seq_id][type_].append((hmm_acc, hmm_name, clan))

    return pfam_dict


def build_relevent_pfam_dict(pfam_addresses, loc_info=False):
    pfam_dicts = []
    for address in pfam_addresses:
        pfam_dict = process_pfam(address, loc_info=loc_info)
        pfam_dicts.append(pfam_dict)

    relevent_pfam_dict = {}


    for pfam in pfam_dicts:
        for gene in list(pfam.keys()):
            relevent_pfam_dict[gene] = pfam[gene]

    del pfam_dicts

    return relevent_pfam_dict


def create_domains_families_clans_lists(gene_pfam):
    gene_clans = []
    domains = []
    families = []
    clans = []

    if 'Domain' in gene_pfam:
        gene_domains_list = gene_pfam['Domain']
        gene_domains = [value[0] for value in gene_domains_list]
        for j, domain in enumerate(gene_domains):
            domains.append(domain)

        gene_clans = gene_clans + \
                     list(set([value[2] for value in gene_domains_list if value[2] != 'No_clan']))

    if 'Family' in gene_pfam:
        gene_families_list = gene_pfam['Family']
        gene_families = [value[0] for value in gene_families_list]
        for j, family in enumerate(gene_families):
            families.append(family)

        gene_clans = gene_clans + \
                     list(set([value[2] for value in gene_families_list if value[2] != 'No_clan']))

    for j, clan in enumerate(gene_clans):
        clans.append(clan)

    return domains, families, clans


def save_object(obj, save_name):
    f = open(save_name + ".pkl", "wb")
    pickle.dump(obj, f, -1)
    f.close()


def load_object(file_name):
    f = open(file_name + ".pkl", "rb")
    obj = pickle.load(f)
    f.close()
    return obj


def derivative_points(x, y, normalize=False):
    yd = []
    xd = []
    for i in range(len(x) - 1):
        slope = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
        yd.append(slope)
        xd.append(i)

    if normalize:
        yd = yd - np.min(yd)
        yd = yd / np.max(yd)

    xd = np.array(xd) / len(xd)
    return xd, yd


def calculate_weighted_distance(x, y1, y2, start_weight=0.5, weight_decay=0.5, return_derivative_list=False,
                                normalize=False):
    weight_list = [start_weight]
    distance_list = [calculate_parameter_euclidean_distance(y1, y2)]
    y1_list = [y1]
    y2_list = [y2]
    x1 = x
    x2 = x

    while len(x1) > 2:
        x1, y1 = derivative_points(x1, y1, normalize=normalize)
        x2, y2 = derivative_points(x2, y2, normalize=normalize)
        y1_list.append(y1)
        y2_list.append(y2)
        distance_list.append(calculate_parameter_euclidean_distance(y1, y2))
        weight_list.append(weight_list[-1] * weight_decay)

    distance = 0
    for i in range(len(distance_list)):
        distance = distance + distance_list[i] * weight_list[i]
    if return_derivative_list:
        return distance, y1_list, y2_list
    else:
        return distance


def calculate_parameter_euclidean_distance(paramters_1, parameters_2):
    euclidean_distance = 0
    for i, p1 in enumerate(paramters_1):
        distance = (paramters_1[i] - parameters_2[i]) ** 2
        euclidean_distance = euclidean_distance + distance

    euclidean_distance = math.sqrt(euclidean_distance)

    return euclidean_distance


def intergrate_curve(x, y, normalize=True):
    y_int = []
    for i in range(len(x) - 1):
        x_i = x[i:i + 2]
        y_i = y[i:i + 2]
        integral = spint.trapz(y_i, x=x_i)
        y_int.append(integral)

    if normalize:
        y_int = y_int - np.min(y_int)
        y_int = y_int / np.max(y_int)

    return x[0:-1], y_int


def calculate_integral_distance(x, y1, y2, start_weight=0.5, weight_decay=0.5, return_integral_list=False,
                                normalize=True):
    weight_list = [start_weight]
    distance_list = [calculate_parameter_euclidean_distance(y1, y2)]
    y1_list = [y1]
    y2_list = [y2]
    x1 = x
    x2 = x

    while len(x1) > 2:
        x1, y1 = intergrate_curve(x1, y1, normalize=normalize)
        x2, y2 = intergrate_curve(x2, y2, normalize=normalize)
        y1_list.append(y1)
        y2_list.append(y2)
        distance_list.append(calculate_parameter_euclidean_distance(y1, y2))
        weight_list.append(weight_list[-1] * weight_decay)

    distance = 0
    for i in range(len(distance_list)):
        distance = distance + distance_list[i] * weight_list[i]
    if return_integral_list:
        return distance, y1_list, y2_list
    else:
        return distance


def create_prot_mean_std_dict(curve_dict, metric='euclidean', ex_distance=False,
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


    prots = list(distance_dataframe.keys())
    prot_mean_std_dict = {}
    for prot in prots:
        val_list = list(distance_dataframe[prot])
        prot_mean_std_dict[prot] = (np.mean(val_list), np.std(val_list))
    del dist, keys, distance_dataframe

    return prot_mean_std_dict, dist_dict


# This function saves Tapioca predictions as a tsv
def network_dict_to_tsv_file(network_dict, savename='./data/test', n=2):
    f = open(savename + '.tsv', "w")
    for key, value in network_dict.items():
        line = ''
        skip_flag = False
        for i in range(n):
            gene = key[i]
            if gene == None:
                skip_flag = True
                break
            line = line + '\t' + gene
        if skip_flag:
            continue

        line = line +'\t'+ str(value) + '\n'
        f.write(line)

    f.close()


# This function reads Tapioca predictions
def tsv_to_dict(address,n=2):
    out_dict = {}
    with open(address, 'rt') as f:
        for line in f:
            try:
                gene_list = []
                if line[0] == '\t':
                    line=line[1:]
                value = float(line.split('\t')[n].replace('\n', '').replace('[', '').replace(']', ''))
                for i in range(n):
                    gene_list.append(line.split('\t')[i])
                out_dict[tuple(gene_list)] = float(value)
            except:
                continue

    return out_dict


def curve_dict_from_pd_master_curve_dict(master_curve_dict, replicate, condition=None, curve_points_len=10):
    curve_dict = {}
    for gene, rep_dict in master_curve_dict.items():
        if replicate in rep_dict:
            if condition == None:
                curve_point_dict = rep_dict[replicate]
            else:
                if condition in rep_dict[replicate]:
                    curve_point_dict = rep_dict[replicate][condition]
                else:
                    continue
            if len(curve_point_dict) > curve_points_len:
                curve = np.reshape(list(curve_point_dict.values())[0:-1], (curve_points_len, 1))
            else:
                curve = np.reshape(list(curve_point_dict.values()), (curve_points_len, 1))

            for i,val in enumerate(curve):
                if np.isnan(val):
                    if i == 0:
                        curve[i] = 1
                    elif i == len(curve)-1:
                        curve[i] = curve[i-1]*.9
                    else:
                        curve[i] = (curve[i - 1] + curve[i + 1]) / 2

            for i in range(len(curve)):
                if np.isnan(curve[i][0]):
                    if i ==0:
                        curve[i][0] = 1
                    elif i == len(curve)-1:
                        curve[i][0] = 0
                    else:
                        curve_i_plus_1 = curve[i+1]
                        curve_i_min_1 = curve[i-1]
                        if np.isnan(curve_i_plus_1) and i+2 < len(curve):
                            curve_i_plus_1 = curve[i+2]
                        if np.isnan(curve_i_plus_1) and i-2 > -1:
                            curve_i_min_1 = curve[i - 2]

                        curve[i][0] = (curve_i_min_1+curve_i_plus_1)/2

                        if np.isnan(curve[i][0]):
                            curve = None
                            break
            if type(curve) != type(None):
                if np.isnan(curve).sum() or np.isinf(curve).sum():
                    continue
                curve_dict[gene] = curve



    return curve_dict


def create_master_curve_dict(address,
                     curve_points=['36.9', '40.2', '43.9', '46.6', '48.6', '52.7','55.3', '58.5', '61.2', '64.0']):
    master_table = pd.read_csv(address)
    master_curve_dict = {}
    for i in range(len(master_table)):
        row = master_table.iloc[i]
        condition = str(row['condition'])
        rep = str(row['replicate'])
        acc = row['accession']
        temp_list = []
        for i in range(len(curve_points)):
            temp_list.append(row[curve_points[i]])

        if acc not in master_curve_dict:
            master_curve_dict[acc] = {}

        if rep not in master_curve_dict[acc]:
            master_curve_dict[acc][rep] = {}

        if condition not in master_curve_dict[acc][rep]:
            master_curve_dict[acc][rep][condition] = {}

        for i in range(len(curve_points)):
            master_curve_dict[acc][rep][condition][curve_points[i]] = temp_list[i]

    return master_curve_dict

def get_curve_points(address):
    table = pd.read_csv(address)
    columns = list(table.columns)
    curve_points = []
    for val in columns:
        if val.replace('.','',1).isdigit():
            curve_points.append(val)
    return curve_points