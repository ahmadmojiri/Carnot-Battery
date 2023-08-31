#This file runs the annual optimisation of the storage system operation.

from projdirs import modeldir, optdir
from subprocess import check_output
from package.make_dzn_file_with_variable_fcas import make_dzn_file_with_variable_fcas
from package.get_NEM_time_series import get_time_series
import matplotlib.pyplot as mpl
from numpy import array


def optimise(simparams):
    """
simparams is a dictionary with the parameters for dzn file."""
    
    curdir =  optdir + "arbitrage\\"
        
    make_dzn_file_with_variable_fcas(**simparams)
    
    output = str(check_output(["c:\\Program Files\\minizinc\\minizinc", "--soln-sep", '""', "--search-complete-msg", '""', "--solver", "OSICBC", curdir + "arbitrage_plus_variable_fcas.mzn", curdir + "arbitrage_plus_variable_fcas.dzn"]))
    
    output = output.replace('\\r\\n', '')
    #output = output.replace('\r', '')
    output = output.replace('[', '')
    output = output.replace(']', '')
    output = output.replace('b', '')
    output = output.replace('"', '')
   
    
    Pin, Penergy, Pfcas, Q, obj = output.split(';')
    
    Q = array(Q[:-6].split(',')).astype(float)
    Pin = array(Pin[2:].split(',')).astype(float)
    Penergy = array(Penergy.split(',')).astype(float)
    Pfcas = array(Pfcas.split(',')).astype(float)
    Obj = float(obj.replace("'",''))
     
    results = {}
    
    results['Q'] = Q
    results['Pin'] = Pin
    results['Penergy'] = Penergy
    results['Pfcas'] = Pfcas
    results['obj'] = Obj
    
    return results

def load_simparams():
    
    simparams = {}
    simparams['eta_in'] = 0.99 #charging efficiency [nondimensional]
    simparams['eta_out'] = 0.4 #discharge efficiency [nondimensional]
    simparams['self_disch'] = 0 # self-discharge during each interval
    simparams['Q_min_frac'] = 0 # miminum charge level
    simparams['Q_init_frac'] = 0.5 #
    simparams['Pin_max'] = 50 #maximum charging power [?]
    simparams['Pin_min_frac'] = 0 #?
    simparams['Pout_max'] = 50 #maximum (electric?) discharge power [?]
    simparams['Pout_min_frac'] = 0 #?
    simparams['SH'] = 8 #?
    simparams['dt'] = 0.5 #simulation time interval? [h?]
    
    return simparams


if __name__ == "__main__":

    start_year, start_month, start_day, start_hour = 2011, 3, 25, 0
    end_year, end_month, end_day, end_hour = 2013, 11, 3, 0
    state = "SA"

    alldata = True
    
    c = get_time_series(state, start_year, start_month, start_day, start_hour,
                        end_year, end_month, end_day, end_hour)
    
    simparams = load_simparams()
    simparams['N'] = len(c)
    simparams['c'] = c
    
    results = optimise(simparams)
    
    fig = mpl.figure()
    ax = fig.add_subplot(111)
    ax.plot(results['Pout'])
    ax.plot(-results['Pin'])
