import predict as p
import preprocess as pre
import scripts as scr
import os

def  run_tapioca(input_file, ref_channel, pre_normalized, co_fractionation, tissue, base_save_name, full_model, output_dir):

    if not os.path.exists(input_file):
        raise FileNotFoundError(f'Error, could not find the input file: {input_file}')
    
    normalized_file = input_file
    if not pre_normalized:
        #normalize the input file
        logistic_curve_fit = not co_fractionation
        normalized_file = pre.preprocess(input_file, output_dir,ref_norm=ref_channel,
                                                 logistic_curve_fit=logistic_curve_fit)
    
        print(f"Normalized file: {normalized_file}")
    else:
        print(f"Input file already normalized: {normalized_file}")

    if full_model:
        model_type = 'F'
    else:
        model_type = 'B'

    curve_points = scr.get_curve_points(normalized_file)
    conditions = scr.get_conditions(normalized_file)
    replicates = scr.get_replicates(normalized_file)

    master_curve_dict = scr.create_master_curve_dict(normalized_file, curve_points=curve_points)

    if co_fractionation:
        master_curve_dict = pre.normalzie_cofrac_dict(master_curve_dict)

    for rep in replicates:
        for cond in conditions:
            rep_cond_base_save_name = base_save_name + '_' + cond + '_rep_' + rep
            curve_dict = scr.curve_dict_from_pd_master_curve_dict(master_curve_dict, rep, condition=cond,
                                                          curve_points_len=len(curve_points))
            p.tapioca_predict(rep_cond_base_save_name, tissue, curve_dict, curve_points, model_type,output_dir)
 

