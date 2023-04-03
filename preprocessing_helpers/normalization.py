from .fit_curves import dfg

def normalize_to_reference_channel(df, tol = 10**(-8)):
    
    '''
    
    Normalize data according to a pooled TMT reference channel from all samples. 
    After normalization, all reference channel values for a given protein should be the same.
    
    parameters:
        df: pd.DataFrame | (n x m) matrix where n is the num ber of observations (i.e. protein signal across conditions)
        and m is the number of temperatures plus one reference channel
        
        tol: float | tolerance for checking that final reference channel values across samples are equal
    
    returns:
        pd.DataFrame | df normalized to reference channel across samples
    
    '''

    # w/in plex normalization
    plex_ref_medians = df.groupby(['condition', 'replicate']).median()['REF']
    MOM = plex_ref_medians.median()
    plex_normed = df/MOM


    # across plex normalization
    across_plex_ref_normed = plex_normed['REF']/plex_normed['REF'].unstack(['condition', 'replicate']).mean(axis=1).loc[plex_normed.index.get_level_values('accession')].values
    ref_normed_data = plex_normed.apply(lambda x: x/across_plex_ref_normed)

    # check that all REF channels are equal across samples
    r = ref_normed_data['REF'].unstack(['condition', 'replicate'])
    r_sub = r.apply(lambda x: x - ref_normed_data['REF'].unstack(['condition', 'replicate']).mean(axis=1))

    if not all(r.where(r.isnull(), r_sub).stack(['condition', 'replicate'])<tol):
        raise ValueError('Normalization to reference channel has failed . . . all reference channels across samples are not equal')

    return ref_normed_data.drop('REF', axis=1)
        
def MOMy_normalization(df,logistic_curve_fit=True):
    
    '''
    
    Normalize data using the median of the means method. Calculate per-sample medians, the median of \
    per-sample medians, and fit the median of the medians (MOM) array using the log-logistic function. \
    Ultimately, this curve is used to adjust the per-sample medians and the input matrix itself.
    
    parameters:
        df: pd.DataFrame | (n x m) matrix where n is the number of observations (i.e. protein signal across conditions)
        and m is the number of temperatures).
        
    returns
        pd.DataFrame | MOMy normalized df
        
    Note: if provided a dataframe with a "REF" column (i.e. a reference channel, this function will throw a ValueError)
        
    '''
    
    check_needs_normalized(df)

    # divide by lowest temperature column (i.e. normalize so that the values are relative to the "most soluble" temperature)
    df_normed = df.apply(lambda x: x/df[df.columns.min()])
    
    # get the median of each temperature column, grouping conditions and replicates together
    medians_by_temp = df_normed.groupby(['condition', 'replicate']).median()
    
    # now get the median of the medians
    MOMs = medians_by_temp.median()

    if logistic_curve_fit == True:
        # fit the MOM curve
        residual = dfg(MOMs.index, MOMs.values)

        # adjust the medians of each sample by the residual
        medians_adj = medians_by_temp.apply(lambda x: x-residual, axis=1)

    # adjust the normalized df by the adjusted medians
    return df_normed.groupby(['condition', 'replicate']).apply(lambda x: x-medians_adj.loc[x.name, :])

def check_needs_normalized(df):
    
    '''
    
    Make sure that the dataframe does not include a "REF" (i.e. reference channel) column.
    
    parameters:
        df: pd.DataFrame
    
    returns:
        None
        
    '''
    
    try:
        df.columns = df.columns.astype(float)
    
    except TypeError:

        issue_cols = []

        for col in df.columns:
            try:
                float(col)
            except ValueError:
                issue_cols.append(col)

        raise ValueError("Dataframe columns {} cannot be coerced to type 'float'... did you forget to normalize based on a reference channel?".format(issue_cols))

