import predict as p
import preprocess as pre
import scripts as scr
import os

def run_tapioca(input_file, ref_channel, pre_normalized, co_fractionation, tissue, base_save_name, full_model):
    if pre_normalized == False:
        if co_fractionation == False:
            logistic_curve_fit = True
            if os.path.exists('./raw_input/' + input_file):
                input_file = './raw_input/' + input_file
                normalized_file = pre.preprocess(input_file, ref_norm=ref_channel,
                                                 logistic_curve_fit=logistic_curve_fit)
            else:
                raise FileNotFoundError('Error, could not find the input file')
        else:
            logistic_curve_fit = False
            if os.path.exists('./raw_input/' + input_file):
                input_file = './raw_input/' + input_file
                normalized_file = pre.preprocess(input_file, ref_norm=ref_channel,
                                                 logistic_curve_fit=logistic_curve_fit)
            else:
                raise FileNotFoundError('Error, could not find the input file')
    else:
        if os.path.exists('./normalized_input/' + input_file):
            input_file = './normalized_input/' + input_file
            normalized_file = input_file
        else:
            raise FileNotFoundError('Error, could not find the input file')

    if full_model == False:
        model_type = 'B'
    else:
        model_type = 'F'

    curve_points = scr.get_curve_points(normalized_file)
    conditions = scr.get_conditions(normalized_file)
    replicates = scr.get_replicates(normalized_file)

    master_curve_dict = scr.create_master_curve_dict(normalized_file, curve_points=curve_points)

    if co_fractionation == True:
        master_curve_dict = pre.normalzie_cofrac_dict(master_curve_dict)

    tissue_address = './data/tissue_functional_networks/' + tissue

    for rep in replicates:
        for cond in conditions:
            rep_cond_base_save_name = base_save_name + '_' + cond + '_rep_' + rep
            curve_dict = scr.curve_dict_from_pd_master_curve_dict(master_curve_dict, rep, condition=cond,
                                                          curve_points_len=len(curve_points))
            p.tapioca_predict(rep_cond_base_save_name, tissue_address, curve_dict, curve_points, model_type)

