#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 15:04:20 2019

@author: jeff
"""

#This file runs the annual optimisation of the storage system operation.

from subprocess import check_output
from projdirs import datadir, optdir
from package.make_dzn_file import make_dzn_file, make_dzn_file_parallel
from numpy import array
import package.get_NEM_data as gd
from datetime import datetime, time
import package.sql_manager as sm
from calendar import monthrange
import numpy as np
import pandas as pd
import pdb
import sys
import matplotlib.pyplot as mpl

win32 = sys.platform == 'win32'

def load_simparams():
    """This function sets the default storage parameters"""    
    simparams = {}
    simparams['eta_in'] = 0.99 #charging efficiency [nondimensional]
    simparams['eta_out'] = 0.4 #discharge efficiency [nondimensional]
    simparams['Q_min_frac'] = 0 # miminum charging allowed?
    simparams['Q_init_frac'] = 0.0 #
    simparams['Pin_max'] = 50 #maximum charging power [?]
    simparams['Pin_min_frac'] = 0 #?
    simparams['Pout_max'] = 50 #maximum (electric?) discharge power [?]
    simparams['Pout_min_frac'] = 0 #?
    simparams['SH'] = 10 #?
    simparams['dt'] = 0.5 #simulation time interval? [h?]
    simparams['loss'] = 0 
    return simparams


def optimise(simparams):
    """simparams is a dictionary containing the storage
    parameters for the dzn file."""    
    curdir =  optdir + "arbitrage%s" %['/', '\\'][win32]
    mzdir = ['/home/jeff/local/software/MiniZinc/bin/',
             'C:\\Program Files\\minizinc\\'][win32]
    make_dzn_file(**simparams)

    output = str(check_output([mzdir + 'minizinc', "--soln-sep", '""', "--search-complete-msg", '""', "--solver", "COIN-BC", curdir + "arbitrage.mzn", curdir + "arbitrage.dzn"]))

    if win32:
        output = output.replace("b\'[", '')
        output = output.replace("""\\r\\n""\\r\\n""\\r\\n\'""", '')
    else:
        output = output.replace("b'[", '')
        output = output.replace("""\\n""\\n""\\n\'""", '')

    output = output.replace('[', '')
    output = output.replace(']', '')

    Pin, Pout, Q = output.split(';')
    
    Pin = array(Pin.split(',')).astype(float)
    Pout = array(Pout.split(',')).astype(float)
    Q = array(Q.split(',')).astype(float)
    obj = sum((Pout-Pin) * simparams['c']) * simparams['dt']
   
    results = {}
    results['Q'] = Q
    results['Pin'] = Pin
    results['Pout'] = Pout
    results['obj'] = obj
    
    return results


def optimise_parallel(ID, simparams):
    """ simparams is a dictionary with the parameters for dzn file.
    Here ID is a unique string to distinguish between models for parallel processing.
    Without this ID, the parallel processing messes up the minizinc data files."""

    curdir =  optdir + "arbitrage%s" %['/', '\\'][win32]
    mzdir = ['/home/jeff/local/software/MiniZinc/bin/',
             'C:\\Program Files\\minizinc\\'][win32]

    make_dzn_file_parallel(ID, **simparams)
    
    output = str(check_output([mzdir + "minizinc", "--soln-sep", '""', "--search-complete-msg", '""', "--solver", "COIN-BC", curdir + "arbitrage.mzn",curdir + "arbitrage_%s.dzn"%(ID)]))

    if win32:
        output = output.replace("b\'[", '')
        output = output.replace("""\\r\\n""\\r\\n""\\r\\n\'""", '')
    else:
        output = output.replace("b'[", '')
        output = output.replace("""\\n""\\n""\\n\'""", '')

    output = output.replace('[', '')
    output = output.replace(']', '')

    Pin, Pout, Q = output.split(';')
    Q = array(Q.split(',')).astype(float)
    Pin = array(Pin.split(',')).astype(float)
    Pout = array(Pout.split(',')).astype(float)
    obj = sum((Pout-Pin) * simparams['c']) * simparams['dt']

    results = {}
    results['Q'] = Q
    results['Pin'] = Pin
    results['Pout'] = Pout
    results['obj'] = obj
    return results

def calc_rolling_storage(Year, State='VIC', SH=10, RTE=40, loss=0, price_cap=14500):
    """
    This function calculates the storage value based on forecast prices
    generated at the begining of each dispatch interval.
    It automatically writes the data into the database file.
    """
    # Create the database file to write the storage values
    db = 'storage_value_rolling_forecast_%s.db'%(State)
    table = 'storage_value'
    cols = ['state CHAR(3)', 'date_time CHAR', 'sh INT',
            'rte INT', 'cap INT', 'loss INT', 'window INT',
            'Pin float', 'Pout float', 'Q_start float',
            'Q_end float', 'obj float']
    sm.create_table(table, db, cols)
    cols_index = ['state', 'date_time', 'rte', 'cap', 'sh', 'loss']
    sm.create_unique_index(db,table, cols_index, 'idx')

    rrp_yearly = gd.load_rrp_cal(Year, State)
    for Month in np.arange(1,13,1):
        forecasts_monthly = gd.load_rrp_forecast_monthly(State, Year, Month)
        
        time_period = pd.date_range(datetime(year=Year, month=Month, day=1, hour=0, minute=0),
                  datetime(year=Year, month=Month, day=monthrange(Year,Month)[1]), freq='30min')
        df= pd.DataFrame()
        #

        # determine whether you want to cap the spot price market or not:
        capping = False #or Fals
        simparams = load_simparams()
        Q_init=0
        obj = 0
        for Time in time_period:
            print(State, Time, 'RTE=',float(RTE)/100, 'SH=',SH)

            # Load the forecast data
            prices_forecast  = forecasts_monthly[(forecasts_monthly['generated']==Time)]
            intervals = len(prices_forecast)
            c=prices_forecast['price'].tolist()

            if not c:
                continue
            simparams['c'] = c
            simparams['N'] = len(c)
            simparams['SH'] = SH
            simparams['eta_out'] = float(RTE)/100
            simparams['loss'] = float(loss)/100
            simparams['Q_init_frac'] = Q_init/(SH*simparams['Pout_max']/simparams['eta_out'])

            results = optimise_parallel(str(Year),simparams)
            Pin = results['Pin'][0]
            Pout = results['Pout'][0]
            
            price_spot = rrp_yearly[Time]
            if capping:
                price_spot=min(price_cap, price_spot)
            obj = (results['Pout'][0]-results['Pin'][0]) * price_spot * simparams['dt']
            Q_start = Q_init
            Q_end = Q_init + (Pin*simparams['eta_in'] - Pout/simparams['eta_out'])* simparams['dt']
            Q_init = Q_end   

            df = df.append(   pd.DataFrame(  [
                        [State,str(Time), str(SH), str(RTE),
                        str(price_cap),str(loss),str(intervals),
                        Pin, Pout, Q_start, Q_end,obj]
                                                ] ) )
        cols_recording = ['state', 'date_time', 'sh',
                              'rte', 'cap', 'loss', 'window',
                              'Pin', 'Pout', 'Q_start',
                              'Q_end', 'obj']
        sm.replace_into_db(df, db, table, cols_recording)
            
    return(Time)

def calc_S_with_daily_forecast(Year, state, Time=time(4,0,0), RTE=90, SH=10,
                               Loss=0, cap=14500, write=True):
    """This function calculates the battery performance for a year under
    optimisation with daily forecast prices for the next trading day.
    'state': e.g. 'NSW
    'Year': e.g. 2018
    'time': the time at which the forecast data is generated.
    if 'write=True', the function writes the data into the database that is defined by 'db'
    """
    db='storage_value_daily_forecast_%s.db'%(state)
    cols_recording = ['state', 'date_time', 'sh', 'rte', 'loss','cap', 'Pin','Pout', 'obj_forecast', 'obj_actual']
        
    DATA = pd.DataFrame(columns = cols_recording )
    simparams = load_simparams()
    simparams['SH'] = SH
    simparams['eta_out'] = float(RTE)/100
    simparams['loss'] = float(Loss)/100
    Q_init=0            
    rrp_yearly = gd.load_rrp_cal(Year, state).append(
        gd.load_rrp_actual(datetime(year=Year+1, month=1, day=1),state)['spot_price'])
    rrp_yearly[rrp_yearly>cap]=cap
    
    for Month in np.arange(1,13,1):
        forecast_prices = gd.load_rrp_24hr_forecast(state,Year,Month,Time)
        if Month<12:
              Date_Range = pd.date_range(
                    datetime(year=Year, month=Month, day=1,hour=4),
                    datetime(year=Year, month=Month+1, day=1))
        else:
              Date_Range = pd.date_range(
                    datetime(year=Year, month=Month, day=1,hour=4),
                    datetime(year=Year+1, month=Month, day=1))
              
        for Date in Date_Range:
            # Load the forecast data
            c_forecast = (
                forecast_prices[forecast_prices['generated'].dt.date == Date.date()]
                )
            time_stamps = c_forecast['date_time']
            if c_forecast.empty:
                  continue
            simparams['c'] = c_forecast['price'].tolist()
            simparams['N'] = len(simparams['c'])
            simparams['Q_init_frac'] = Q_init/(SH*simparams['Pout_max']/simparams['eta_out'])

            results = optimise_parallel(str(Year),simparams)
            obj_forecast = simparams['dt'] * (results['Pout']-results['Pin']) *\
                             c_forecast['price']
            obj_actual = simparams['dt'] * (results['Pout']-results['Pin']) *\
                             rrp_yearly[time_stamps].tolist()
            
#            pdb.set_trace()

            time_stamp_to_record = pd.date_range(time_stamps.tolist()[0], periods=48,
                                                 freq='0.5H')
            n = len(time_stamp_to_record)
#            pdb.set_trace()
            DATA = DATA.append(pd.DataFrame({'state': state,
                                             'date_time':time_stamp_to_record,
                                             'sh': SH,
                                             'rte': RTE,
                                             'loss': Loss,
                                             'cap': cap,
                                             'Pin': results['Pin'][0:n],
                                             'Pout': results['Pout'][0:n],                                           
                                             'obj_forecast':obj_forecast[0:n],
                                             'obj_actual': obj_actual[0:n]
                                            }), ignore_index=True)
            Q_init = results['Q'][-1]
        print('Month %d completed!'%(Month))
        
    if write:
        DATA=DATA[cols_recording]
        DATA['date_time'] = DATA['date_time'].astype(str)
        sm.replace_into_db(DATA, db, 'storage_value', cols_recording,
            create_unique_idx=True,
            idx_cols=cols_recording[0:-4])

    return(DATA)                

def calc_S_perfect_foresight(State, Year, RTE, SH, Loss, Cap,
                             write=False, db='storage_value_arbitrage.db',
                             table='storage_value'):
    """This function calculates the storage value with perfect foresight.
        write: specifies whether to write the data into the database or not.
        db: the database to write
        table: the table in he db to write
    """

    simparams = load_simparams()
      
    cols = ['state CHAR(3)', 'year INT', 'sh INT',  'rte INT', 'cap INT',
              'loss INT', 'obj float']
    idx = ['state', 'year', 'rte', 'cap', 'sh', 'loss']
      
    print(State, Year, float(RTE)/100, SH)
      
    c = gd.load_rrp_cal(Year, State)
    c[c>Cap]=Cap
    c = c.tolist()
      
    simparams['c'] = c
    simparams['N'] = len(c)
    simparams['SH'] = SH
    simparams['eta_out'] = float(RTE)/100
    simparams['loss'] = float(Loss)/100
      
    results = optimise_parallel(str(Year), simparams)
      
      
    # register the generated data into the database
    if write:
        DATA = pd.DataFrame({'state': [State],
                           'year':[int(Year)],
                           'sh': [int(SH)],
                           'rte': [int(RTE)],
                           'cap': [Cap],
                           'loss': [Loss],
                           'obj': [results['obj']]
                           })
        sm.create_table(table, db, cols, create_unique_idx=True, idx_cols=idx)
        sm.replace_into_db(DATA, db,table, cols)
    return(results)


def test_parallel(a,b,c):
    print(a+b+c)
    return(a+b+c)
    

if __name__ == "__main__":

    state = "SA"
    year = 2014
    window = 31

    Date = pd.date_range(datetime(year=year, month=1, day=1),
                         datetime(year=year, month=12, day=31))[0]

    rrp = pd.DataFrame()
    for Day in pd.date_range(Date, periods=window):
        rrp = rrp.append(pd.DataFrame(gd.load_rrp_actual(Day, state)['spot_price']))
    
    c=rrp['spot_price'].tolist()

    simparams = load_simparams()

    simparams['c'] = c
    simparams['N'] = len(c)
    
    results = optimise_parallel(100, simparams)