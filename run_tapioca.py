import tapioca as tp
import os
import argparse
import datetime
import pandas as pd

# error messages
INVALID_FILETYPE_MSG = "\n Error: Invalid file format. %s must be a .csv file."
INVALID_PATH_MSG = "\n Error: Invalid file path/name. Path %s does not exist."
INVALID_CONTENT_MSG = "\n Error: Invalid data content/structure. %s must be reformatted to fit expected structure." \
                      " See example input file on the Tapioca GitHub (https://github.com/FunctionLab/tapioca)"

INVALID_NORM_CONTENT_MSG = "\n Error: Invalid data content/structure. %s must be reformatted to fit expected structure." \
                      " See example normlaized input file on the Tapioca GitHub" \
                      "(https://github.com/FunctionLab/tapioca)"


def validate_file(file_name, prenorm):
    #validate file name and path.
    if not valid_path(file_name):
        print(INVALID_PATH_MSG % (file_name))
        quit()
    elif not valid_filetype(file_name):
        print(INVALID_FILETYPE_MSG % (file_name))
        quit()
    elif not validate_file_structure(file_name,prenorm):
        if prenorm:
            print(INVALID_CONTENT_MSG % (file_name))
        else:
            print(INVALID_NORM_CONTENT_MSG % (file_name))
        quit()


def valid_filetype(file_name):
    # validate file type
    return file_name.endswith('.csv')


def valid_path(path):
    # validate file path
    print(path)
    print(os.path.exists(path))
    return os.path.exists(path)

def validate_file_structure(path,prenorm):
    # validate file contents
    table = pd.read_csv(path)
    cols = list(table.columns)
    pass_flag = True
    if prenorm:
        if 'accession' not in cols:
            print('\n Error: condition label not found in input file. Make sure condition is typed exactly as seen here and'
                  'that the condition label is in the first row of the first column of the file')
            pass_flag = False

        if 'condition' not in cols:
            print('\n Error: condition label not found in input file. Make sure condition is typed exactly as seen here and'
                  'that the condition label is in the first row of the second column of the file')
            pass_flag = False

        if 'replicate' not in cols:
            print('\n Error: replicate label not found in input file. Make sure replicate is typed exactly as seen here and'
                  'that the replicate label is in the first row of the second column of the file')
            pass_flag = False

        for col in cols[3:]:
            if not col.replace('.','').isnumeric():
                print('Error: all curve data points (ex. temperatures) should be a number, '+str(col)+'is not a number.'
                      )
                pass_flag = False
    else:
        if 'condition' not in cols:
            print('\n Error: condition label not found in input file. Make sure condition is typed exactly as seen here and'
                  'that the condition label is in the first row of the first column of the file')
            pass_flag = False
        if list(table['condition'])[1] != 'replicate':
            print('\n Error: replicate label not found in input file. Make sure replicate is typed exactly as seen here and'
                  'that the replicate label is in the third row of the first column of the file')
            pass_flag = False
        if list(table['condition'])[2] != 'accession':
            print('\n Error: accession label not found in input file. Make sure accession is typed exactly as seen here and'
                  'that the accession label is in the fourth row of the first column')
            pass_flag = False

    return pass_flag


def main():
    # create parser object
    parser = argparse.ArgumentParser(description="A Command Line Interface For Running Tapioca")

    # defining arguments for parser object
    parser.add_argument("-i", "--input", type=str,
                        metavar="raw_file", default=None,required=True,
                        help="The path to the input csv file")

    parser.add_argument("-b", "--output_base_name", type=str,
                        metavar="output_base_name", default="pred",
                        help="The base save name for the prediction files")

    parser.add_argument("-o", "--output", type=str,
                        metavar="output_dir_name", required=True,
                        help="The output directory where the predictions will be stored")

    parser.add_argument("-r", "--ref", type=int,
                        metavar="ref_channel", default=1,
                        help="0 to not perform Reference channel normalization. Default 1")

    parser.add_argument("-p", "--prenorm", type=int,
                        metavar="pre_normalized", default=0,
                        help="Set 1 when inputting pre-normalized data. Default 0")

    parser.add_argument("-c", "--cofrac", type=int,
                        metavar="co_fractionation", default=0,
                        help="Set 1 when inputting cofractionation data. Default 0")

    parser.add_argument("-f", "--fullmodel", type=int,
                        metavar="full_model", default=0,
                        help="Set 0 to use only the base submodel. Default 1")

    parser.add_argument("-t", "--tissue", type=str,
                        metavar="tissue", default='',
                        help="The path to the tissue-specific functional network you would like to use.")

    # parse the arguments from standard input
    args = parser.parse_args()
    print(args)
    args.ref = bool(args.ref)
    args.prenorm = bool(args.prenorm)
    args.cofrac = bool(args.cofrac)
    args.fullmodel = bool(args.fullmodel)

    # Validate the Input
    validate_file(args.input, args.prenorm)


    if not os.path.exists(args.output):
        try:
            os.mkdir(args.output) 
        except OSError as error: 
            parser.exit(1, message=f"Can't create output {args.output} . Check file permissions")

    temp_dir = os.path.join(args.output,"temp")
    if not os.path.exists(temp_dir):
        try:
            os.mkdir(temp_dir) 
        except OSError as error: 
            parser.exit(1, message=f"Can't create output temp directory >{temp_dir}<. Check file permissions")

    if args.tissue and not os.path.exists(args.tissue):
        parser.exit(1, message=f"The tissue network file: >{args.tissue}< doesn't exist")


    output_base_name = args.output_base_name
    # If savename is None then set a default name based on the time
    if not args.output_base_name:
        current_datatime = str(datetime.datetime.now()).replace('-', '').replace(':', '') \
                        .split('.')[0].replace(' ', '')
        output_base_name = current_datatime



    tp.run_tapioca(
        input_file=args.input,
        ref_channel=args.ref,
        pre_normalized=args.prenorm,
        co_fractionation=args.cofrac,
        tissue=args.tissue,        
        full_model=args.fullmodel,
        output_dir=args.output,
        base_save_name=output_base_name,

    )


if __name__ == "__main__":
    # calling the main function
    main()