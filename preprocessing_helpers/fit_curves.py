import pandas as pd
import numpy as np
from scipy import optimize as opt

def log_logistic_func(x, a, b, c):
    
    '''
    
    Log-logistic function definition: 
    
    parameters:
        x: array | points at which to evaluate the function
        a: float | function parameter (also the curve inflection point (e.g. the Tm)) (between T(min) and T(max))
        b: float | function parameter (b > 0)
        c: float | function parameter (0 < c < 1)
        
    returns:
        r: array | log-logistic function evaluated at x given parameters a, b, and c 
    
    '''
    
    r = c + (1-c)/(1+np.exp(b*(np.log(x)-np.log(a))))
    
    return r 
    
def fit_func(func, x_data, y_data, **kwargs):
    
    '''
    
    Leverages the Scipy optimizer to optmize func over its parameters relative to y_data given x_data.
    
    parameters:
        func: function
        x_data: np.array | input data to func
        y_data: np.array | array against which to optimize the function parameters
        
    returns:
        opt_params: np.array | optimized function parameters in array form; ordered according to input order in the func
    
    '''
   
    opt_params  = opt.curve_fit(func, x_data, y_data, **kwargs)[0]
    
    return opt_params

def dfg(x, y):
    
    '''
    
    Fit the log-logistic function parameters for x given y.
    
    parameters:
        x: np.array
        y: np.array
        
    returns:
        df_g: pd.Series | log-logistic curve evaluated at x when parameters are optimized against y
    
    '''
    
    a, b, c = fit_func(log_logistic_func, x, y, maxfev=5000)
    df_g = pd.Series([log_logistic_func(x, a, b, c) for x in x], index = x)

    return df_g
