#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 21:06:01 2020

@author: jeff
"""

from .sim_classes import PSF_plant
from subprocess import check_output
from projdirs import optdir
from package.make_dzn_file_PV_plant import make_dzn_file_PV_plant
from numpy import array
import sys
import pandas as pd
import matplotlib.pyplot as plt

pd.plotting.register_matplotlib_converters()

windows = sys.platform == 'win32'

def load_simparams(c, P_in_panels):
    """This function sets the default storage parameters"""    
    
    simparams = {}

    simparams['c'] = c
    simparams['N'] = len(c)

    simparams['dt'] = 0.5 #simulation time interval [h]
    simparams['eta_rh'] = 0.99 #charging efficiency [nondimensional]
    simparams['eta_pc'] = 0.4 #discharge efficiency [nondimensional]
    simparams['eta_inv'] = 0.96
    simparams['eta_ts'] = 0.99**(1./24)
    simparams['Q_min_frac'] = 0.1 #Minimum charge in tank
    simparams['Q_init_frac'] = 0.5 #
    simparams['P_ts_in_max'] = 150 #Maximum electrical power into the tank elements
    simparams['P_ts_out_max'] = 50 #Maximum generation by power cycle
    simparams['P_inv_in_max'] = PSF_plant.power_AC / 1000 / simparams['eta_inv']
    simparams['P_in_panels'] = P_in_panels[:len(c)]     
    simparams['SH'] = 10 #Storage hours
    

    return simparams

def optimise(simparams):
    """simparams is a dictionary containing the storage
    parameters for the dzn file."""    
    curdir =  optdir + "arbitrage%s" %['/', '\\'][windows]
    mzdir = ['/home/jeff/local/software/MiniZinc/bin/',
             'C:\\Program Files\\minizinc\\'][windows]
    
    make_dzn_file_PV_plant(**simparams)

    output = str(check_output([mzdir + 'minizinc', "--soln-sep", '""', 
                               "--search-complete-msg", '""', "--solver", 
                               "COIN-BC", curdir + "arbitrage_PV_plant.mzn", 
                               curdir + "arbitrage_PV_plant.dzn"]))

    output = output.replace("""\\r\\n""\\r\\n""\\r\\n\'""", '')
    output = output.replace("b\'[", '')
    output = output.replace('[', '')
    output = output.replace(']', '')
    
    if not windows: output = output.replace('\\n""\\n""\\n\'','')

    P_in_grid, P_out_pass, P_inv_in, P_out_grid, Q = output.split(';')
    
    P_in_grid = array(P_in_grid.split(',')).astype(float)
    P_out_pass = array(P_out_pass.split(',')).astype(float)
    P_inv_in = array(P_inv_in.split(',')).astype(float)
    P_out_grid = array(P_out_grid.split(',')).astype(float)
    Q = array(Q.split(',')).astype(float)
    obj = sum((P_out_grid-P_in_grid) * simparams['c']) * simparams['dt']
   
    results = {}
    results['Q'] = Q
    results['P_in_panels'] = array(simparams['P_in_panels'])
    results['P_in_grid'] = P_in_grid
    results['P_out_grid'] = P_out_grid
    results['P_out_pass'] = P_out_pass
    results['P_inv_in'] = P_inv_in
    results['obj'] = obj
    
        #calculate the equivalent electrical energy leaving the thermal battery
    results['P_ts_out'] = results['P_out_grid'] - results['P_inv_in'] * PSF_plant.eta_inv  #(MW-e)
    
    #calculate the electrical energy being used to charge the thermal battery
    results['P_ts_in'] = results['P_in_panels'] + results['P_in_grid'] - \
        results['P_out_pass'] - results['P_inv_in']
        
    results['P_ts_in_panels'] = results['P_ts_in'] - results['P_in_grid']
    
    P_ts_in_panels_cur = results['P_in_panels'] - simparams['P_inv_in_max']
    P_ts_in_panels_cur[P_ts_in_panels_cur < 0] = 0
    results['P_ts_in_panels_cur'] = P_ts_in_panels_cur
    
    results['P_ts_in_panels_vol'] = results['P_ts_in_panels'] - P_ts_in_panels_cur
    results['P_ts_in_grid'] = results['P_in_grid'] 
    results['timeindex'] = simparams['P_in_panels'].index.to_timestamp()
    results['c'] = simparams['c']
    
    return results

def plot_results(results):
    
    #plt.rcParams['figure.subplot.bottom'] = 0.14
        
    fig1 = plt.figure()#figsize = (160./25.4, 160./25.4))
    ax1 = fig1.add_subplot(311)
    ax2 = fig1.add_subplot(312)
    ax3 = fig1.add_subplot(313)
    
    index = results['timeindex']
    
    #just to be clear that none of the
    # grid purchased electricity should be going to the inverter, and only to storage
    
    # bar1 = ax1.bar(rrp.index, P_ts_out, simparams['dt']/24, -P_ts_out)[0]
    # bar2 = ax1.bar(rrp.index, P_ts_in_panels, simparams['dt']/24, 0)[0]
    # bar3 = ax1.bar(rrp.index, P_ts_in_grid, simparams['dt']/24, P_ts_in_panels)[0]
    
    # fig1.legend([bar1, bar2, bar3], 
    #            ['P_ts_out', 'P_ts_in_panels', 'P_ts_in_grid'], loc = 9, ncol = 5)
    
    line1 = ax1.plot(index, -results['P_ts_out'], label = 'P_ts_out')[0]
    line2 = ax1.plot(index, results['P_ts_in_panels_cur'], label = 'P_ts_in_panels_cur')[0]
    line3 = ax1.plot(index, results['P_ts_in_panels_vol'], label = 'P_ts_in_panels_vol')[0]
    line4 = ax1.plot(index, results['P_ts_in_grid'], label = 'P_ts_in_grid')[0]
     
    fig1.legend([line1, line2, line3, line4], 
            ['P_ts_out', 'P_ts_in_panels_cur', 'P_ts_in_panels_vol', 
             'P_ts_in_grid'], loc = 9, ncol = 4)
    
    ax1.grid(True)
    ax1.set_xlabel('')
    ax1.set_ylabel('Power (MW-e)')
    ax1.axis([index[0], index[-1], -60, 60])
    
    ax2.plot(index, results['Q'][:-1])
    ax2.grid(True)
    ax2.set_xlabel('')
    ax2.set_ylabel('Q (MWh-th)')
    ax2.set_xlim([index[0], index[-1]])
    
    ax3.plot(index, results['c'])
    ax3.axis([index[0], index[-1],0,200])
    ax3.grid(True)
    ax3.set_xlabel('Date')
    ax3.set_ylabel('Spot Price\n(\$/MWh)')
    
    plt.show()