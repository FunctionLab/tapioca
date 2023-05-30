from sklearn.linear_model import SGDClassifier
from sklearn.calibration import CalibratedClassifierCV

import scripts as scr

import numpy as np
import itertools
import os

from scipy.spatial.distance import euclidean
from config import MODELS_FEAT_DIR

class Base_Data_Generator():

    def __init__(self, prot_ints_dict, pairs_list, stat_dict, batch_size=1, shuffle=True, prot_mean_std_dict = {},
                 x_={'default':[36.9, 40.2, 43.9, 46.6, 48.6, 52.7, 55.3, 58.5, 61.2, 64]}):
        self.base_init(prot_ints_dict, pairs_list, stat_dict, batch_size=batch_size,
                       shuffle=shuffle, x_=x_, prot_mean_std_dict=prot_mean_std_dict,)

    def base_init(self, prot_ints_dict, pairs_list, stat_dict, batch_size=1, shuffle=True,
                  prot_mean_std_dict = {},
                  x_={'default':[36.9, 40.2, 43.9, 46.6, 48.6, 52.7, 55.3, 58.5, 61.2, 64]}):


        self.prot_ints_dict = prot_ints_dict
        self.pairs_list = pairs_list
        self.stat_dict = stat_dict
        self.prot_mean_std_dict = prot_mean_std_dict
        self.dist_func = euclidean
        self.batch_size = batch_size
        self.shuffle = shuffle

        self.x_ = x_
        for key,item in self.x_.items():
            self.x_[key] = item - np.min(item)
            self.x_[key] = item / np.max(item)

        self.pair_dict = {}
        for pair in self.pairs_list:
            self.pair_dict[pair] = -1

        del self.pairs_list

        self.list_IDs = list(self.pair_dict.keys())
        self.on_epoch_end()
        self.create_deriv_dict()
        self.create_integral_dict()

    def create_deriv_dict(self):
        self.deriv_dict = {}

        for key, value_list in self.prot_ints_dict.items():

            for dataset_key in self.x_.keys():
                if dataset_key in key:
                    x_ = self.x_[dataset_key]
                    break
                else:
                    x_ = self.x_['default']

            temp_len = len(x_)
            if len(value_list) > temp_len:
                value_list = value_list[0:temp_len]
            while len(value_list) < temp_len:
                value_list = np.append(value_list, value_list[-1])
            value_list = np.reshape(value_list, (temp_len))

            distance, y1, y2 = scr.calculate_weighted_distance(
                x_, value_list, value_list, return_derivative_list=True, normalize=True)

            self.deriv_dict[key] = y1

    def create_integral_dict(self):
        self.integral_dict = {}
        for key, value_list in self.prot_ints_dict.items():
            for dataset_key in self.x_.keys():
                if dataset_key in key:
                    x_ = self.x_[dataset_key]
                    break
                else:
                    x_ = self.x_['default']

            temp_len = len(x_)
            if len(value_list) > temp_len:
                value_list = value_list[0:temp_len]
            while len(value_list) < temp_len:
                value_list = np.append(value_list, value_list[-1])
            value_list = np.reshape(value_list, (temp_len))

            distance, y1, y2 = scr.calculate_integral_distance(
                x_, value_list, value_list, return_integral_list=True, normalize=True)

            self.integral_dict[key] = y1

    def on_epoch_end(self):
        self.indexes = np.arange(len(self.list_IDs))
        # if self.shuffle == True:
        #     np.random.shuffle(self.indexes)

    def mean_euc_dist(self,y1, y2):
        cdist = []
        for i, val1 in enumerate(y1):
            val1 = float(val1)
            val2 = float(y2[i])
            if val1 == 0 and val2 == 0:
                cdist_ = 0
            else:
                cdist_ = (val1 - val2) ** 2
            cdist.append(cdist_)

        return (np.mean(cdist)) ** 0.5

    def base_data_generation(self, list_IDs_temp):
   
        batch_size = len(list_IDs_temp)

        cdist = np.empty((batch_size, 13))
        curve_corr = np.empty((batch_size, 7))
        rel_dist = np.empty((batch_size, 3))

        # Generate data
        for i, ID in enumerate(list_IDs_temp):
            gene1 = ID[0]
            gene2 = ID[1]
            y1 = self.deriv_dict[gene1]
            y2 = self.deriv_dict[gene2]
            yi1 = self.integral_dict[gene1]
            yi2 = self.integral_dict[gene2]
            # Store sample
            try:
                sample = gene1.split('*')[1]
                mean = self.stat_dict[sample]['mean']
                std = self.stat_dict[sample]['std']
            except:
                mean = self.stat_dict['mean']
                std = self.stat_dict['std']


            cdist[i, 0] = (self.dist_func(y1[0], y2[0]) - mean)/std
            cdist[i, 1] = euclidean(y1[1], y2[1])
            cdist[i, 2] = self.mean_euc_dist(y1[1], y2[1])
            cdist[i, 3] = euclidean(y1[2], y2[2])
            cdist[i, 4] = self.mean_euc_dist(y1[2], y2[2])
            cdist[i, 5] = euclidean(y1[3], y2[3])
            cdist[i, 6] = self.mean_euc_dist(y1[3], y2[3])
            cdist[i, 7] = euclidean(yi1[1], yi2[1])
            cdist[i, 8] = self.mean_euc_dist(yi1[1], yi2[1])
            cdist[i, 9] = euclidean(yi1[2], yi2[2])
            cdist[i, 10] = self.mean_euc_dist(yi1[2], yi2[2])
            cdist[i, 11] = euclidean(yi1[3], yi2[3])
            cdist[i, 12] = self.mean_euc_dist(yi1[3], yi2[3])


            curve_corr[i, 0] = np.corrcoef(y1[0], y2[0])[0][1]
            curve_corr[i, 1] = np.corrcoef(y1[1], y2[1])[0][1]
            curve_corr[i, 2] = np.corrcoef(y1[2], y2[2])[0][1]
            curve_corr[i, 3] = np.corrcoef(y1[3], y2[3])[0][1]
            curve_corr[i, 4] = np.corrcoef(yi1[1], yi2[1])[0][1]
            curve_corr[i, 5] = np.corrcoef(yi1[2], yi2[2])[0][1]
            curve_corr[i, 6] = np.corrcoef(yi1[3], yi2[3])[0][1]

            stats_1 = self.prot_mean_std_dict[gene1]
            stats_2 = self.prot_mean_std_dict[gene2]
            rel_1 = (euclidean(y1[0], y2[0]) - stats_1[0])/stats_1[1]
            rel_2 = (euclidean(y1[0], y2[0]) - stats_2[0])/stats_2[1]
            rel_dist[i,] = [np.min([rel_1,rel_2]), np.mean([rel_1,rel_2]), np.max([rel_1,rel_2])]

        return cdist, curve_corr, rel_dist,

    def data_generation(self, list_IDs_temp):
        out = self.base_data_generation(list_IDs_temp)
        X = np.concatenate(out, axis=-1)
        return X

    def __len__(self):
        if not self.batch_size:
            return 1

        return int(np.ceil(len(self.list_IDs) / self.batch_size))

    def __getitem__(self, index):
        if self.batch_size == None:
            list_IDs_temp =  self.list_IDs
        else:
            # Generate indexes of the batch
            indexes = self.indexes[index * self.batch_size:(index + 1) * self.batch_size]

            # Find list of IDs
            list_IDs_temp = [self.list_IDs[k] for k in indexes]

        # Generate data
        X = self.data_generation(list_IDs_temp)

        return X, list_IDs_temp


class Extended_Data_Generator(Base_Data_Generator):

    def __init__(self, prot_ints_dict, pairs_list, stat_dict, batch_size=1, shuffle=True, prot_mean_std_dict={},
                 x_={'default':[36.9, 40.2, 43.9, 46.6, 48.6, 52.7, 55.3, 58.5, 61.2, 64]}):
        super().__init__(prot_ints_dict, pairs_list, stat_dict, batch_size=batch_size,
                         shuffle=shuffle, prot_mean_std_dict=prot_mean_std_dict, x_=x_)

        self.base_init(prot_ints_dict, pairs_list, stat_dict,batch_size=batch_size, shuffle=shuffle, x_=x_,
                       prot_mean_std_dict=prot_mean_std_dict)

    def other_info_init(self, prior_dict=None, pfam_dict=None, properties_dict=None):

        self.prior_dict=prior_dict
        self.pfam_dict=pfam_dict
        self.properties_dict=properties_dict

        if self.pfam_dict != None:
            self.domain_co_occurance_dict = scr.load_object(MODELS_FEAT_DIR+'/human_herpesvirus_gs_pfam_domain_co_occurance_dict')
            self.family_co_occurance_dict = scr.load_object(MODELS_FEAT_DIR+'/human_herpesvirus_gs_pfam_family_co_occurance_dict')
            self.clan_co_occurance_dict = scr.load_object(MODELS_FEAT_DIR+'/human_herpesvirus_gs_pfam_clan_co_occurance_dict')

    def pfam_process(self, list1, list2, co_occurance_dict):
        shared_funcs = 0
        scores = []

        for l1 in list1:
            if l1 == '':
                continue

            l1 = l1.replace(' ', '')

            for l2 in list2:
                if l2 == '':
                    continue

                l2 = l2.replace(' ', '')

                if (l1, l2) in co_occurance_dict:
                    scores.append(co_occurance_dict[l1, l2])
                elif (l2, l1) in co_occurance_dict:
                    scores.append(co_occurance_dict[l2, l1])
                else:
                    scores.append(0)

                if l1 == l2:
                    shared_funcs = shared_funcs + 1

        if len(scores) > 0:
            output = [np.min(scores), np.mean(scores), np.max(scores), shared_funcs]
        else:
            output = [0, 0, 0, shared_funcs]

        return output

    def shared_attributes(self, attribute_list1, attribute_list2, norm=False):
        num_shared = 0
        for attr in attribute_list1:
            if attr in attribute_list2:
                num_shared = num_shared + 1

        if norm and max(len(attribute_list1), len(attribute_list2)) > 0:
            num_shared = num_shared / max(len(attribute_list1), len(attribute_list2))

        return num_shared

    def other_data_generation(self, list_IDs_temp):
        batch_size = len(list_IDs_temp)
        
        prior = np.empty((batch_size, 1))

        length = np.empty((batch_size, 2))
        gravy = np.empty((batch_size, 2))
        molecular_weight = np.empty((batch_size, 2))
        aromaticity = np.empty((batch_size, 2))
        instability_index = np.empty((batch_size, 2))
        isoelectric_point = np.empty((batch_size, 2))
        secondary_structure_fraction = np.empty((batch_size, 6))

        domains = np.empty((batch_size, 4))
        families = np.empty((batch_size, 4))
        clans = np.empty((batch_size, 4))

        # Generate data
        for i, ID in enumerate(list_IDs_temp):
            gene1 = ID[0]
            gene2 = ID[1]
            true_gene1 = gene1.split('*')[0].split('_')[0]
            true_gene2 = gene2.split('*')[0].split('_')[0]
            try:
                sample = gene1.split('*')[1]
            except:
                sample = 'predict'

            if self.prior_dict != None:
                tissue = self.prior_dict['info'][sample]
                tissue_dict = self.prior_dict[tissue]

                if (true_gene1, true_gene2) in tissue_dict:
                    prior[i,] = float(tissue_dict[true_gene1, true_gene2])
                elif (true_gene2, true_gene1) in tissue_dict:
                    prior[i,] = float(tissue_dict[true_gene2, true_gene1])
                else:
                    prior[i,] = 0.0

            if true_gene1 in self.pfam_dict and true_gene2 in self.pfam_dict:
                domains_1, families_1, clans_1 = scr.create_domains_families_clans_lists(self.pfam_dict[true_gene1])
                domains_2, families_2, clans_2 = scr.create_domains_families_clans_lists(self.pfam_dict[true_gene2])

                domains[i,] = self.pfam_process(domains_1, domains_2, self.domain_co_occurance_dict)
                families[i,] = self.pfam_process(families_1, families_2, self.family_co_occurance_dict)
                clans[i,] = self.pfam_process(clans_1, clans_2, self.clan_co_occurance_dict)

            else:
                domains[i,] = [0, 0, 0, 0]
                families[i,] = [0, 0, 0, 0]
                clans[i,] = [0, 0, 0, 0]

            if true_gene1 in self.properties_dict and true_gene2 in self.properties_dict:

                gene1_props = self.properties_dict[true_gene1]
                gene2_props = self.properties_dict[true_gene2]

                gene1_length = gene1_props['length'] * 0.001
                gene2_length = gene2_props['length'] * 0.001

                gene1_gravy = gene1_props['gravy'] * 1
                gene2_gravy = gene2_props['gravy'] * 1

                gene1_mw = gene1_props['molecular_weight'] * 0.0001
                gene2_mw = gene2_props['molecular_weight'] * 0.0001

                gene1_arom = gene1_props['aromaticity'] * 1
                gene2_arom = gene2_props['aromaticity'] * 1

                gene1_ii = gene1_props['instability_index'] * 0.01
                gene2_ii = gene2_props['instability_index'] * 0.01

                gene1_ip = gene1_props['isoelectric_point'] * 0.1
                gene2_ip = gene2_props['isoelectric_point'] * 0.1

                gene1_ssf = gene1_props['secondary_structure_fraction'] * 1
                gene2_ssf = gene2_props['secondary_structure_fraction'] * 1

                length[i,] = [(gene1_length + gene2_length) / 2.0, abs(gene1_length - gene2_length)]
                gravy[i,] = [(gene1_gravy + gene2_gravy) / 2.0, abs(gene1_gravy - gene2_gravy)]
                molecular_weight[i,] = [(gene1_mw + gene2_mw) / 2.0, abs(gene1_mw - gene2_mw)]
                aromaticity[i,] = [(gene1_arom + gene2_arom) / 2.0, abs(gene1_arom - gene2_arom)]
                instability_index[i,] = [(gene1_ii + gene2_ii) / 2.0, abs(gene1_ii - gene2_ii)]
                isoelectric_point[i,] = [(gene1_ip + gene2_ip) / 2.0, abs(gene1_ip - gene2_ip)]
                secondary_structure_fraction[i,] = list((gene1_ssf + gene2_ssf) / 2.0) + list(
                    abs(gene1_ssf - gene2_ssf))

            else:
                length[i,] = np.zeros(2)
                gravy[i,] = np.zeros(2)
                molecular_weight[i,] = np.zeros(2)
                aromaticity[i,] = np.zeros(2)
                instability_index[i,] = np.zeros(2)
                isoelectric_point[i,] = np.zeros(2)
                secondary_structure_fraction[i,] = np.zeros(6)

        return_dict = {}

        return_dict['prior'] = prior

        return_dict['length'] = length
        return_dict['gravy'] = gravy
        return_dict['molecular_weight'] = molecular_weight
        return_dict['aromaticity'] = aromaticity
        return_dict['instability_index'] = instability_index
        return_dict['isoelectric_point'] = isoelectric_point
        return_dict['secondary_structure_fraction'] = secondary_structure_fraction

        return_dict['domains'] = domains
        return_dict['families'] = families
        return_dict['clans'] = clans

        return return_dict

    def data_generation(self, list_IDs_temp):
        out = self.base_data_generation(list_IDs_temp)

        base = np.concatenate(out, axis=-1)

        other_info_dict = self.other_data_generation(list_IDs_temp)

        prior = other_info_dict['prior']

        length = other_info_dict['length']
        gravy = other_info_dict['gravy']
        molecular_weight = other_info_dict['molecular_weight']
        aromaticity = other_info_dict['aromaticity']
        instability_index = other_info_dict['instability_index']
        isoelectric_point = other_info_dict['isoelectric_point']
        secondary_structure_fraction = other_info_dict['secondary_structure_fraction']

        props = np.concatenate((length, gravy, molecular_weight, aromaticity,
                            instability_index, isoelectric_point, secondary_structure_fraction), axis=-1)



        domains = other_info_dict['domains']
        families = other_info_dict['families']
        clans = other_info_dict['clans']
        pfam = np.concatenate((domains, families, clans), axis=-1)

        if self.prior_dict != None:
            lr1 = base
            lr2 = np.concatenate((base, props), axis=-1)
            lr3 = np.concatenate((base, pfam), axis=-1)
            lr4 = np.concatenate((base, prior), axis=-1)
            lr5 = np.concatenate((base, props, pfam), axis=-1)
            lr6 = np.concatenate((base, prior, pfam), axis=-1)
            lr7 = np.concatenate((base, prior, props), axis=-1)
            lr8 = np.concatenate((base, prior, props, pfam), axis=-1)

            X = [lr1,lr2,lr3,lr4,lr5,lr6,lr7,lr8]
        else:
            lr1 = base
            lr2 = np.concatenate((base, props), axis=-1)
            lr3 = np.concatenate((base, pfam), axis=-1)
            lr4 = np.concatenate((base, props, pfam), axis=-1)

            X = [lr1, lr2, lr3, lr4]


        return X


class Base_Model():

    def __init__(self, model_address=None, loss='log', penalty='l2', max_iter=1000, n_jobs=-1, warm_start=True):
        if model_address!=None:
            self.load_model(model_address)
        else:
            self.create_model(loss=loss, penalty=penalty, max_iter=max_iter,
                              n_jobs=n_jobs, warm_start=warm_start,)

    def create_model(self,loss='log', penalty='l2', max_iter=200, n_jobs=-1,
                     warm_start=True, num_models=1):

        self.model_list = []
        for i in range(num_models):
            model = SGDClassifier(loss=loss, penalty=penalty, max_iter=max_iter,
                                           n_jobs=n_jobs, warm_start=warm_start)

            model = CalibratedClassifierCV(base_estimator=model, n_jobs=n_jobs, method='isotonic')

            self.model_list.append(model)

    def load_model(self, model_address):
        self.model_list = None
        self.model = None
        if type(model_address) == list:
            model_list = []
            for model in model_address:
                model_list.append(scr.load_object(model))
            self.model_list = model_list
        else:
            self.model = scr.load_object(model_address)

    def build_data_generator(self, prot_ints_dict, pairs_list, stat_dict,batch_size=1, shuffle=True,
                             prot_mean_std_dict={},
                             x_={'default':[36.9, 40.2, 43.9, 46.6, 48.6, 52.7, 55.3, 58.5, 61.2, 64]}):
        data_generator = Base_Data_Generator(prot_ints_dict, pairs_list, stat_dict,
                                                batch_size=batch_size, shuffle=shuffle, x_=x_,
                                                prot_mean_std_dict=prot_mean_std_dict,)

        return data_generator

    def predict(self, save_address,output_dir, prot_ints_dict, batch_size=1, shuffle=True,
                x_={'default':[36.9, 40.2, 43.9, 46.6, 48.6, 52.7, 55.3, 58.5, 61.2, 64]}):

        x_ = {'default': x_}

        prot_mean_std_dict, dist_dict = scr.create_prot_mean_std_dict(prot_ints_dict, ex_distance=False,
                                                                      metric='euclidean',
                                                                      curve_points=x_['default'])

        genes_list = list(prot_ints_dict.keys())
        pairs_list = list(itertools.combinations(genes_list, 2))


        stat_dict = {}
        stat_dict['mean'] = np.mean(list(dist_dict.values()))
        stat_dict['std'] = np.std(list(dist_dict.values()))

        self.data_generator = self.build_data_generator(prot_ints_dict.copy(), pairs_list,
                                                        stat_dict.copy(), batch_size=batch_size, shuffle=shuffle,
                                                        x_=x_, prot_mean_std_dict=prot_mean_std_dict.copy())

        del dist_dict, prot_ints_dict, prot_mean_std_dict

        num_batches = self.data_generator.__len__()
        self.data_generator.on_epoch_end()
        pred_dict = {}

        for i in range(num_batches):
            print('Batch: ' + str(i + 1) + '/' + str(num_batches))
            X, list_IDs_temp = self.data_generator.__getitem__(i)
            preds = self.model.predict_proba(X)

            for i in range(len(preds)):
                gene1_key = list_IDs_temp[i][0]
                gene2_key = list_IDs_temp[i][1]
                pred_dict[gene1_key, gene2_key] = preds[i][1]
        save_address = os.path.join(output_dir,save_address)
        scr.save_dict_to_tsv_file(pred_dict, output_filename=save_address)


class Extended_Model(Base_Model):
    def __init__(self, model_address=None, loss='log', penalty='l2', max_iter=1000, n_jobs=-1,
                     warm_start=True):

        super().__init__(model_address, loss, penalty, max_iter, n_jobs, warm_start)
        if model_address != None:
            self.load_model(model_address)
        else:
            self.create_model(loss=loss, penalty=penalty, max_iter=max_iter,
                              n_jobs=n_jobs, warm_start=warm_start)


    def build_data_generator(self, prot_ints_dict, pairs_list, stat_dict, prior_dict=None, pfam_dict=None,
                             properties_dict=None, batch_size=1, shuffle=True, prot_mean_std_dict={},
                             x_={'default':[36.9, 40.2, 43.9, 46.6, 48.6, 52.7, 55.3, 58.5, 61.2, 64]}):

        data_generator = Extended_Data_Generator(prot_ints_dict, pairs_list, stat_dict,
                                                batch_size=batch_size, shuffle=shuffle, x_=x_,
                                                prot_mean_std_dict=prot_mean_std_dict)

        data_generator.other_info_init(prior_dict=prior_dict, pfam_dict=pfam_dict, properties_dict=properties_dict,)

        return data_generator

    def fix_curve(self, curve, shape=(1, 10)):
        if len(curve) > 10:
            curve = curve[0:10]
        while len(curve) < 10:
            curve = np.append(curve, curve[-1])
        curve = np.reshape(curve, shape)

        return curve

    def network_dict_to_tsv_file(self, network_dict, savename='./data/test', n=2, mode='a'):

        f = open(savename + '.tsv', mode)
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

    def dict_cleaner(self, in_dict, relevant_keys):
        del_list = []
        for key in in_dict.keys():
            if key not in relevant_keys:
                del_list.append(key)

        for key in del_list:
            del in_dict[key]

        return in_dict

    def predict(self, save_address, prot_ints_dict, output_dir,
                prior_dict_address=None, go_dict_address=None, pfam_dict_address_list=None,
                properties_dict_address=None, batch_size=1, shuffle=True,
                x_={'default':[36.9, 40.2, 43.9, 46.6, 48.6, 52.7, 55.3, 58.5, 61.2, 64]}):

        x_ = {'default':x_}

        prot_mean_std_dict, dist_dict = scr.create_prot_mean_std_dict(prot_ints_dict, ex_distance=False,
                                                                     metric='euclidean', curve_points=x_['default'])
        genes_list = list(prot_ints_dict.keys())
        pairs_list = list(itertools.combinations(genes_list, 2))

        del genes_list

        prior_dict = None
        if prior_dict_address != '':
            relevant_prots = set()
            for prot in prot_ints_dict.keys():
                relevant_prots.add(prot)
            prior_dict = {}
            prior_dict['info'] = {}
            print(f"prior_dict_address: {prior_dict_address}, relevant_prots: {relevant_prots}")
            prior_dict['predict'] = scr.create_tissue_dict(prior_dict_address, relevant_prots=relevant_prots)
            prior_dict['info']['predict'] = 'predict'

        if type(properties_dict_address) == str():
            properties_dict = scr.load_object(properties_dict_address)
        else:
            properties_dict = {}
            for prop_address in properties_dict_address:
                properties_dict.update(scr.load_object(prop_address))

        relevant_prots = set()
        for prot in prot_ints_dict.keys():
            relevant_prots.add(prot)

        pfam_dict = scr.build_relevent_pfam_dict(pfam_dict_address_list)
        pfam_dict = self.dict_cleaner(pfam_dict, relevant_prots)

        stat_dict = {}
        stat_dict['mean'] = np.mean(list(dist_dict.values()))
        stat_dict['std'] = np.std(list(dist_dict.values()))


        self.data_generator = self.build_data_generator(prot_ints_dict.copy(), pairs_list, stat_dict.copy(),
                                                        prior_dict=prior_dict, pfam_dict=pfam_dict,
                                                        properties_dict=properties_dict, batch_size=batch_size,
                                                        shuffle=shuffle, x_=x_,
                                                        prot_mean_std_dict=prot_mean_std_dict.copy())


        del dist_dict, prot_mean_std_dict, prot_ints_dict

        num_batches = self.data_generator.__len__()
        self.data_generator.on_epoch_end()
        pred_dict = {}

        if prior_dict_address != '':
            save_add_on = ['_base', '_prop', '_pfam', '_prior', '_prop_pfam', '_prior_pfam', '_prior_prop',
                           '_prior_prop_pfam'
                           ]
        else:
            save_add_on = ['_base', '_prop', '_pfam', '_prop_pfam']

        #Creating temporary files
        save_address = os.path.join(output_dir,save_address)
        for i, model in enumerate(self.model_list):
            f = open(save_address + save_add_on[i] + '.tsv', 'w')
            f.close()


        for i in range(num_batches):
            print('Batch: ' + str(i + 1) + '/' + str(num_batches))
            X, list_IDs_temp = self.data_generator.__getitem__(i)

            if prior_dict_address != '':
                save_add_on = ['_base', '_prop', '_pfam', '_prior', '_prop_pfam', '_prior_pfam', '_prior_prop',
                               '_prior_prop_pfam'
                               ]
            else:
                save_add_on = ['_base', '_prop', '_pfam', '_prop_pfam']

            for j, model in enumerate(self.model_list):
                preds = model.predict_proba(X[j])
                pred_dict = {}

                for i in range(len(preds)):
                    gene1_key = list_IDs_temp[i][0]
                    gene2_key = list_IDs_temp[i][1]
                    pred_dict[gene1_key, gene2_key] = preds[i][1]
                #Appending to the aleady files
                self.network_dict_to_tsv_file(pred_dict, savename=save_address + save_add_on[j], mode='a')

        if self.model!= None:
            scr.network_dict_to_tsv_file(pred_dict, savename=save_address)

