# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 22:14:09 2022

@author: Ahmad Mojiri
"""


from projdirs import optdir
from subprocess import check_output
from projdirs import datadir
from datetime import timedelta
import numpy as np
import pandas as pd
import random, math, os, pdb



def make_dzn_file_TB(ID, N, dt, eta_in, eta_out, Q_max, Q_min, Q_init,
                     Pin_max, Pin_min, Pout_max, Pout_min,
                     Pout, rrp, r_reg_rrp, l_reg_rrp):
    """
    Create a text file for minizinc optimisation of thermal battery
    """
    
    
    string = """
N = %i;

dt = %.2f;      %%time difference between sample points (also, interval length) (h)

eta_in = %.2f;                      %%efficiency of electricity to heat
eta_out = %.2f;                     %%efficiency of heat to electricity

Q_max = %.2f;    %%maximum stored energy (MWh-th)
Q_min = %.2f;               %%minimum stored energy (MWh-th)
Q_init = %.2f;              %%initial stored energy (MWh-th)

Pin_max = %.2f;                     %%maximum power to storage (MW-th)
Pin_min = %.2f;           %%minimum power to storage (MW-th)

Pout_max = %.2f;                    %%maximum power to grid (MW-th)
Pout_min = %.2f;         %%minimum power to grid (MW-th)

Pout = %s;

rrp = %s;

r_reg_rrp = %s;

l_reg_rrp = %s;

 """ %(N, dt, eta_in, eta_out, Q_max, Q_min, Q_init,
                     Pin_max, Pin_min, Pout_max, Pout_min,
                     str(Pout), str(rrp), str(r_reg_rrp), str(l_reg_rrp))
    with open(optdir + "thermal battery\\ThermalBattery_%s.dzn"%(str(ID)),
              "w") as text_file:
        text_file.write(string)


def optimise_TB(ID, simparams):
    """
    
    Parameters: simparams includes a range of inputs that are defined by 
    load_simparams_TB() function and can be adjusted individually.
    ----------
    simparams : TYPE
        This function optimises thermal battery schdule for minimizing
        the cost of charging including raise and lower regulation FCAS
        services.

    Returns results['Q']: the state of charge
            results['Pin']: input power in MW_e
            results['Pout']: output power in MW_th
            results['obj']: cost of charging over the optimisation cycle in
                            AUD; this cost shoould be then modified by applying
                            the random FCAS dispatch signals
            results['Praise']: Power allocated and enabled for raise regulation
            results['Plower']: Power allocated and enabled for lower regualtion
    -------
    
    """

    curdir =  optdir + "thermal battery\\"
    mzdir = 'C:\\Program Files\\minizinc\\'
    make_dzn_file_TB(ID, **simparams)
    
    output = str(check_output([mzdir + 'minizinc', "--soln-sep", '""',
                               "--search-complete-msg", '""', "--solver", "COIN-BC",
                               curdir + "ThermalBattery.mzn",
                               curdir + "ThermalBattery_%s.dzn"%(str(ID))]))
    
    
    output = output.replace("""\\r\\n""\\r\\n""\\r\\n\'""", '')
    output = output.replace("b\'[", '')
    output = output.replace('[', '')
    output = output.replace(']', '')
    
    Pin, Pout, Praise, Plower, Q, obj = output.split(';')
    
    Q = np.array(Q.split(',')).astype(float)
    Pin = np.array(Pin.split(',')).astype(float)
    Pout = np.array(Pout.split(',')).astype(float)
    Praise = np.array(Praise.split(',')).astype(float)
    Plower = np.array(Plower.split(',')).astype(float)
    
    
    results = {}
    results['Q'] = Q
    results['Pin'] = Pin
    results['Pout'] = Pout
    results['obj'] = obj
    results['Praise'] = Praise
    results['Plower'] = Plower
    
    if os.path.exists(curdir + "ThermalBattery_%s.dzn"%(str(ID))):
        os.remove(curdir + "ThermalBattery_%s.dzn"%(str(ID)))
    else:
        print("The file does not exist")
        
    return(results)







def revenue_TB(prices, start, period, horizon, ID):
    """
    This optimisation code is based on a constant load of 1 MW.

    Parameters
    ----------
    start : DateTime object
        SPECIFIES THE STARTING TIME OF SIMULATION.
    prices : DataFrame
        DATAFRAME CONTAINING THE PRICE TIMESERIES OF ENERGY AND FCAS FOR A YEAR
        STARTING AT start
    period : float
        THE NUMBER OF DAYS THAT THE CALCULATIONS SHOULD RUN.
    horizon: integer
        THE NUMBER OF DAYS COMPRISING THE FORECAST HORIZON

    Returns
    -------
    A datafram containing details of the battery operation

    """
    

    simparams = {}
    
    simparams['dt'] = 5/60
    dt = simparams['dt']
    
    simparams['eta_in'] = 1;
    simparams['eta_out'] = 1;
    
    simparams['Pin_min'] = 0 #MW
    simparams['Pin_max'] = 10 #MW
    
    simparams['Pout_min'] = 0 #MW
    simparams['Pout_max'] = 1 #MW
    
    
    revenue = 0
    Pnet = 0
    Q=[]
    
    Results = pd.DataFrame(columns=['date_time', 'rrp', 'raise_rrp',
                                    'lower_rrp', 'Pin', 'Praise',
                                    'Plower', 'FCAS_service', 'Pnet',
                                    'dQ', 'Q', 'revenue'])
    
    for t, start in enumerate(pd.date_range(start,
                                            start+timedelta(days=period),
                                            freq='5Min', closed='left')):
        print(start, ID, period, horizon)
        end = start + timedelta(days=horizon)
        price_slice = prices.loc[start:end]
    
        simparams['rrp'] = price_slice.RRP.tolist()
        simparams['r_reg_rrp'] = price_slice.RAISEREGRRP.tolist()
        simparams['l_reg_rrp'] = price_slice.LOWERREGRRP.tolist()
    
        simparams['N'] = len(simparams['rrp'])
        simparams['Pout'] = (simparams['Pout_max']*np.ones(simparams['N'])).tolist()
    
        Q_max = (10-simparams['Pin_max'] * simparams['dt'])
    
        simparams['Q_min'] =  0 
        simparams['Q_max'] = math.floor(100 * Q_max)/100 #MWh
    
        if t==0:
            simparams['Q_init'] = simparams['Q_min'] #MWh
        else:
            simparams['Q_init'] = Q[1]
    
        optimised = optimise_TB(ID, simparams)
    
        Q = optimised['Q']
    
        FCAS_service = random.sample(['raise','idle','lower'],1)[0]
    
        if (FCAS_service == 'raise'):
            if Q[1]>=optimised['Praise'][0]*dt:
                Q[1] = Q[1] - optimised['Praise'][0]*dt
                dQ = Pnet*dt
                Pnet = optimised['Pin'][0] - optimised['Praise'][0]
                revenue = (
                        optimised['Praise'][0]*dt*simparams['r_reg_rrp'][0] +
                        optimised['Plower'][0]*dt*simparams['l_reg_rrp'][0] -            
                        Pnet*dt*simparams['rrp'][0]
                            )
            else:
                Pnet = optimised['Pin'][0]
                Q[1] = Q[1]
                dQ = Pnet*dt
                revenue = (
                        optimised['Plower'][0]*dt*simparams['l_reg_rrp'][0] -            
                        Pnet*dt*simparams['rrp'][0]
                            )
    
    
        elif (FCAS_service == 'lower'):
            if Q_max - Q[1] >= optimised['Plower'][0]*dt:
                Q[1] = Q[1] + optimised['Plower'][0]*dt
                Pnet = optimised['Pin'][0] + optimised['Plower'][0]
                dQ = Pnet*dt
                revenue = (
                        optimised['Praise'][0]*dt*simparams['r_reg_rrp'][0] +
                        optimised['Plower'][0]*dt*simparams['l_reg_rrp'][0] -            
                        Pnet*dt*simparams['rrp'][0]
                            )
            else:
                Pnet = optimised['Pin'][0]
                Q[1] = Q[1]
                revenue = (
                        optimised['Praise'][0]*dt*simparams['r_reg_rrp'][0] -            
                        Pnet*dt*simparams['rrp'][0]
                            )
    
        elif FCAS_service == 'idle':
            Pnet = optimised['Pin'][0]
            dQ = Pnet*dt
            Q[1] = Q[1]
            revenue = (
                optimised['Praise'][0]*dt*simparams['r_reg_rrp'][0] +
                optimised['Plower'][0]*dt*simparams['l_reg_rrp'][0] -            
                Pnet*dt*simparams['rrp'][0]
                        )
    
        Results = Results.append({'date_time': start,
                                  'rrp': simparams['rrp'][0],
                                  'raise_rrp': simparams['r_reg_rrp'][0],
                                  'lower_rrp': simparams['l_reg_rrp'][0],
                                  'Pin': optimised['Pin'][0],
                                  'Praise': optimised['Praise'][0],
                                  'Plower': optimised['Plower'][0],
                                  'Pnet': Pnet,
                                  'FCAS_service': FCAS_service,
                                  'revenue': revenue,
                                  'dQ': dQ,
                                  'Q': Q[1]},
                                ignore_index=True)        
        
    return(Results)
    # Time the code                
    # %load_ext line_profiler
    # %lprun -f test test()
                
                                                  
    # note: if at any interval, the starting Q is lower than what is required to meet the load, 
    # the battery will not be enabled for raise service in that interval.