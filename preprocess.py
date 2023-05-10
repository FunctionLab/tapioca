import preprocessing_helpers as ph
import numpy as np
import os


def process_in_file(input_filename,output_filename, ref_norm=True, logistic_curve_fit=True):
    # ------ load in data -------
    non_normed = ph.load_data.load_protein_data(input_filename)
    #in_file = in_file.replace('raw_input','normalized_input')
    if ref_norm:
        # ------- reference normalization -------

        ref_normed_data = ph.normalization.normalize_to_reference_channel(non_normed)
        cols = list(ref_normed_data.columns)
        del_cols = set()
        for col in cols:
            if 'unnamed' in col.lower():
                del_cols.add(col)

        for col in del_cols:
            del ref_normed_data[col]

        # ------- MOM normalization -------
        MOMy_normed_data = ph.normalization.MOMy_normalization(ref_normed_data, logistic_curve_fit=logistic_curve_fit)
        MOMy_normed_data.to_csv(output_filename)
        #normalized_file_address = in_file + '_normalized.csv'

    else:
        MOMy_normed_data = ph.normalization.MOMy_normalization(non_normed, logistic_curve_fit=logistic_curve_fit)
        MOMy_normed_data.to_csv(output_filename)
        #normalized_file_address = in_file + '_normalized.csv'


    #return normalized_file_address


def preprocess(input_file,output_dir,ref_norm,logistic_curve_fit):
    output_filename = os.path.join(output_dir,'input_normalized.csv')
    process_in_file(input_file,output_filename,ref_norm=ref_norm,
                                              logistic_curve_fit=logistic_curve_fit)
    return output_filename


def normalzie_cofrac_dict(cf_dict):
    del_list = set()
    for prot, rep_dict in cf_dict.items():
        for rep, cond_dict in rep_dict.items():
            for cond, val_dict in cond_dict.items():
                max_val = np.max(list(val_dict.values()))
                if max_val == 0:
                    del_list.add((prot, rep, cond))
                    continue
                for frac, val in val_dict.items():
                    cf_dict[prot][rep][cond][frac] = val / max_val

    for thing in del_list:
        del cf_dict[thing[0]][thing[1]][thing[2]]

    return cf_dict