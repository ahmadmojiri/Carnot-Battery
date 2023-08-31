#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 15:04:20 2019

@author: jeff
"""

#This file runs the annual optimisation of the storage system operation.

from subprocess import check_output
from projdirs import datadir
from projdirs import optdir
from package.make_dzn_file import make_dzn_file,\
                                  make_dzn_file_parallel,\
                                  make_dzn_file_PJM_parallel,\
                                  make_dzn_file_capcon,\
                                  make_dzn_file_ramped,\
                                  make_dzn_file_ECPP
from numpy import array
import package.get_NEM_data as gd
from datetime import datetime, time, timedelta
import package.sql_manager as sm
import numpy as np
import pandas as pd
import sys, pdb,os

windows = sys.platform == 'win32'

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
    curdir =  optdir + "arbitrage%s" %['/', '\\'][windows]
    mzdir = ['/home/jeff/local/software/MiniZinc/bin/',
             'C:\\Program Files\\minizinc\\'][windows]
    make_dzn_file(**simparams)

    output = str(check_output([mzdir + 'minizinc', "--soln-sep", '""', "--search-complete-msg", '""', "--solver", "COIN-BC", curdir + "arbitrage.mzn", curdir + "arbitrage.dzn"]))

    output = output.replace("""\\r\\n""\\r\\n""\\r\\n\'""", '')
    output = output.replace("b\'[", '')
    output = output.replace('[', '')
    output = output.replace(']', '')
    
    if not windows: output = output.replace('\\n""\\n""\\n\'','')

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



def optimise_ramped(ID,simparams):
    """simparams is a dictionary containing the storage
    parameters for the dzn file."""    
    curdir =  optdir + "arbitrage%s" %['/', '\\'][windows]
    mzdir = ['/home/jeff/local/software/MiniZinc/bin/',
             'C:\\Program Files\\minizinc\\'][windows]
    make_dzn_file_ramped(ID,**simparams)

    output = str(check_output([mzdir + 'minizinc', "--soln-sep", '""',
                               "--search-complete-msg", '""', "--solver",
                               "COIN-BC", curdir + "arbitrage_NEM_with_ramp.mzn",
                               curdir + "arbitrage_ramped_%s.dzn"%(ID)   ]))

    output = output.replace("""\\r\\n""\\r\\n""\\r\\n\'""", '')
    output = output.replace("b\'[", '')
    output = output.replace('[', '')
    output = output.replace(']', '')
    
    if not windows: output = output.replace('\\n""\\n""\\n\'','')

    Pin, Pout, Q = output.split(';')
    
    Q = array(Q.split(',')).astype(float)
    Pin = array(Pin.split(',')).astype(float)
    Pout = array(Pout.split(',')).astype(float)
#    obj = sum((Pout-Pin) * simparams['c']) * simparams['dt']
   
    results = {}
    results['Q'] = Q
    results['Pin'] = Pin
    results['Pout'] = Pout
    print('%s completed!'%(ID))
    if os.path.exists(optdir + "arbitrage%sarbitrage_ramped_%s.dzn"%(['/', '\\'][windows],ID)):
        os.remove(optdir + "arbitrage%sarbitrage_ramped_%s.dzn"%(['/', '\\'][windows],ID))
    else:
        print("The file does not exist")      
    return results




def optimise_parallel(ID, simparams):
    """ simparams is a dictionary with the parameters for dzn file.
    Here ID is a unique string to distinguish between models for parallel processing.
    Without this ID, the parallel processing messes up the minizinc data files."""
    curdir =  optdir + "arbitrage%s" %['/', '\\'][windows]

    mzdir = ['/home/jeff/local/software/MiniZinc/bin/',
             'C:\\Program Files\\minizinc\\'][windows]

    make_dzn_file_parallel(ID, **simparams)
    
    output = str(check_output([mzdir + 'minizinc', "--soln-sep", '""', 
                               "--search-complete-msg", '""', "--solver", 
                               "COIN-BC", curdir + "arbitrage.mzn", 
                               curdir + "arbitrage_%s.dzn"%(ID)]))
    
    output = output.replace("""\\r\\n""\\r\\n""\\r\\n\'""", '')
    output = output.replace("b\'[", '')
    output = output.replace('[', '')
    output = output.replace(']', '')
    
    if not windows: output = output.replace('\\n""\\n""\\n\'','')

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
    if os.path.exists(optdir + "arbitrage%sarbitrage_%s.dzn"%(['/', '\\'][windows],ID)):
        os.remove(optdir + "arbitrage%sarbitrage_%s.dzn"%(['/', '\\'][windows],ID))
    else:
        print("The file does not exist")
    return results





def optimise_ramped_parallel(Year,Region,SH,RTE,Ramp_rate,Pin_max,Pout_max):
      """
      This function calculates the SV for a storage system with a linear
      discharge ramp rate. Here, ramp up rate=ramp down rate.
      Rampe rate: % of full power rating per minute
      """
      c = gd.load_rrp_cal(Year, Region)
      c = c.resample('5T').pad()
      simparams = load_simparams()
      simparams['eta_out'] = float(RTE)/100
      simparams['R']=Pout_max*Ramp_rate*60 #MW/hr
      simparams['c']=c.tolist()
      simparams['N']=len(simparams['c'])
      simparams['dt']=0.5/6
      simparams['Pin_max']=Pin_max
      simparams['Pout_max']=Pout_max
      
      ID = '%s_%d_%d_%d' %(Region, Year, RTE, int(100*Ramp_rate))
      results = optimise_ramped(ID,simparams)
#      pdb.set_trace()
      
      Results = pd.DataFrame( {'Pout':results['Pout'],
                               'Pin':results['Pin'],
                               'Q':results['Q']} )
      
      Data = pd.DataFrame()
      Data['date_time']= c.index
      Data['price']=c.tolist()
      Data['Pout_mean']=Results.Pout.rolling(2).mean()[1:].reset_index(drop=True)
      Data['Pin_mean']=Results.Pin.rolling(2).mean()[1:].reset_index(drop=True)
      
                           
      SV = simparams['dt']*( Data.price*(Data.Pout_mean-Data.Pin_mean)  ).sum()
      return({'Region': Region,
            'Year': Year,
            'RTE': RTE,
            'SH':SH,
            'Ramp_rate': Ramp_rate,
            'Pin_max':Pin_max,
            'Pout_max': Pout_max,
            'SV': SV
                  })
      
      
      
      


def optimise_PJM_parallel(Zone, Year, Pin_max, Pout_max, SH, RTE, Loss,
                          market='DA+RT', details= 'off'):
    """
    simparams is a dictionary with the parameters for dzn file.
    Here ID is a unique string to distinguish between models for parallel processing.
    Without this ID, the parallel processing messes up the minizinc data files.
    This function optimises the storage for wither of the following cases:
          1- 'DA': arbitrage only in DA market
          2- 'RT': arbitrage only in RT market
          3- 'DA+RT': arbotrage simultaneuosly in DA and RT markets
    """
    curdir =  optdir + "arbitrage%s" %['/', '\\'][windows]

    mzdir = ['/home/jeff/local/software/MiniZinc/bin/',
             'C:\\Program Files\\minizinc\\'][windows]
    lmps_da = gd.load_lmp_zone_year(Zone, Year, lmp_type='da', Year_type='cal')
    lmps_rt = gd.load_lmp_zone_year(Zone, Year, lmp_type='rt', Year_type='cal')
    time_stamp_to_record = lmps_da.index.tolist()
    simparams = load_simparams()
    LMP_da = lmps_da.tolist()
    LMP_rt = lmps_rt.tolist()
    simparams['Pout_max'] = Pout_max
    simparams['Pin_max'] = Pin_max
    simparams['c_da'] = LMP_da
    simparams['c_rt'] = LMP_rt
    simparams['N'] = len(LMP_da)
    simparams['SH'] = SH
    simparams['eta_out'] = float(RTE)/100
    simparams['loss'] = float(Loss)/100
    simparams['dt'] = 1    
    
    if market=='DA+RT':
          simparams['Pda_frac'] = 1
          simparams['Prt_frac'] = 1
    elif market=='DA':
          simparams['Pda_frac'] = 1
          simparams['Prt_frac'] = 0
    elif market=='RT':
          simparams['Pda_frac'] = 0
          simparams['Prt_frac'] = 1
          
   
    ID = str(Year)+Zone.replace('/','') +str(RTE) + str(Loss) + str(SH)         
    make_dzn_file_PJM_parallel(ID, **simparams)
    
    output = str(check_output([mzdir + 'minizinc', "--soln-sep", '""', 
                               "--search-complete-msg", '""', "--solver", 
                               "COIN-BC", curdir + "arbitrage_PJM.mzn", 
                               curdir + "arbitrage_PJM_%s.dzn"%(ID)]))
    
    output = output.replace("""\\r\\n""\\r\\n""\\r\\n\'""", '')
    output = output.replace("b\'[", '')
    output = output.replace('[', '')
    output = output.replace(']', '')  
    output = output.replace("b'", '')
    if not windows: output = output.replace('\\n""\\n""\\n\'','')

    obj, Pin_da, Pin_rt, Pout_da, Pout_rt, Q = output.split(';')
    obj = array(obj).astype(float)
    Q = array(Q.split(',')).astype(float)
    Pin_da = array(Pin_da.split(',')).astype(float)
    Pin_rt = array(Pin_rt.split(',')).astype(float)
    Pout_da = array(Pout_da.split(',')).astype(float)
    Pout_rt = array(Pout_rt.split(',')).astype(float)
    #obj = sum((Pout-Pin) * simparams['c']) * simparams['dt']

    results = {}
    results['Q'] = Q
    results['Pin_da'] = Pin_da
    results['Pin_rt'] = Pin_rt
    results['Pout_da'] = Pout_da
    results['Pout_rt'] = Pout_rt
    results['lmp_da'] = LMP_da
    results['lmp_rt'] = LMP_rt
    results['obj'] = obj.item(0)
    results['date_time'] = time_stamp_to_record
    print(ID, ' ', results['obj'])
    if os.path.exists(optdir + "arbitrage%sarbitrage_PJM_%s.dzn"%(['/', '\\'][windows],ID)):
        os.remove(optdir + "arbitrage%sarbitrage_PJM_%s.dzn"%(['/', '\\'][windows],ID))
    else:
        print("The file does not exist")
    Output = pd.DataFrame({'zone': [Zone],
                           'year': [Year],
                           'rte': [RTE],
                           'pin_max': [Pin_max],
                           'pout_max': [Pout_max],
                           'sh': [SH],
                           'loss': [Loss],
                           'arbi_type':[market],                              
                           'obj': [results['obj']]})
    
    return [Output,(Output,results)][details=='on']
    



def optimise_PJM_parallel_Ec(Zone, Year, Pin_max, Pout_max, Ec, RTE, Loss,
                          market='DA+RT', details= 'off'):
    """
    This function calculates the value of stprage with constant storage capacity.
    The SH here varies with Ec and Pout.
    simparams is a dictionary with the parameters for dzn file.
    Here ID is a unique string to distinguish between models for parallel processing.
    Without this ID, the parallel processing messes up the minizinc data files.
    This function optimises the storage for wither of the following cases:
          1- 'DA': arbitrage only in DA market
          2- 'RT': arbitrage only in RT market
          3- 'DA+RT': arbotrage simultaneuosly in DA and RT markets
    """
    SH = Ec*RTE/(100*Pout_max)
    Output = optimise_PJM_parallel(Zone, Year, Pin_max, Pout_max, SH, RTE,
                                   Loss, market='DA+RT', details= 'off')
    Output['Ec']=Ec
    return Output





def calc_SV_capcon(State, Year, RTE, Loss, Pin_max, Pout_max, SH,
                   Strike_prices, cap_prices):
    """ 
    Calcuate the SV for cap contracts in the NEM.
    Opimises the SV with an imposed contraint forcing the storage to
    necessarily discharge at inervals with RRP greater that the strike price.
    Without this ID, the parallel processing messes up the minizinc data files.
    Strike_prices: the strike price of the cap contract in each quarter [s1,s2,s3,s4]
    cap_price: the price of the contract ($/MWh) in each quarter [c1,c2,c3,c4]
    c=0 means that cap contract has not been sold for that quarter
    """
    curdir =  optdir + "arbitrage%s" %['/', '\\'][windows]

    mzdir = ['/home/jeff/local/software/MiniZinc/bin/',
             'C:\\Program Files\\minizinc\\'][windows]
    
    
    Data = gd.load_rrp_fy(Year,State).reset_index()
    Data['pout_constr']=0
    
    simparams = load_simparams()
    simparams['loss'] = Loss
    simparams['Pin_max'] = Pin_max
    simparams['Pout_max'] = Pout_max
    simparams['SH'] = SH
    simparams['eta_out'] = RTE/100
    
    
      
    for C,q in enumerate(cap_prices):
          if C>0:
                Data.loc[(Data['date_time'].dt.quarter==q+1)&
                        (Data['spot_price']>Strike_prices[q])
                        ,['spot_price','pout_constr']]=[Strike_prices[q],1]
    
    Data = Data.set_index('date_time').sort_index()
    Data=Data.reset_index()
    #pdb.set_trace()            
    simparams['c'] = Data['spot_price'].tolist()
    simparams['Pout_constr']= (Data[Data['pout_constr']>0].index+1).tolist() 
    
    ID = State+str(Year)+str(RTE)
    make_dzn_file_capcon(ID, **simparams)
    
    output = str(check_output([mzdir + 'minizinc', "--soln-sep", '""', 
                               "--search-complete-msg", '""', "--solver", 
                               "COIN-BC", curdir + "arbitrage_NEM_capcon.mzn", 
                               curdir + "arbitrage_NEM_capcon_%s.dzn"%(ID)]))
    #pdb.set_trace()
    output = output.replace("""\\r\\n""\\r\\n""\\r\\n\'""", '')
    output = output.replace("b\'[", '')
    output = output.replace('[', '')
    output = output.replace(']', '')  
    output = output.replace("b'", '')
    if not windows: output = output.replace('\\n""\\n""\\n\'','')

    obj, Pin, Pout, Q, Pout_constr = output.split(';')
    obj = array(obj).astype(float)
    Q = array(Q.split(',')).astype(float)
    Pin = array(Pin.split(',')).astype(float)
    Pout = array(Pout.split(',')).astype(float)
    Pout_constr = array(Pout_constr.split(',')).astype(int)
    
    results = {}
    results['Q'] = Q
    results['Pin'] = Pin
    results['Pout'] = Pout
    results['Pout_constr'] = Pout_constr
    results['obj'] = obj
    if os.path.exists(optdir + "arbitrage%sarbitrage_NEM_capcon_%s.dzn"%(['/', '\\'][windows],ID)):
        os.remove(optdir + "arbitrage%sarbitrage_NEM_capcon_%s.dzn"%(['/', '\\'][windows],ID))
    else:
        print("The file does not exist")
    return results







def calc_rolling_storage(Year, State='VIC', SH=10, RTE=40, loss=0, price_cap=14500):
    """
    This function calculates the storage value based on forecast prices
    generated at the begining of each trading interval.
    It automatically writes the data into the database file.
    """
    # Create the database file to write the storage values
    db = 'storage_value_rolling_forecast_%s.db'%(State)
    table = 'storage_value'
    cols = ['state CHAR(3)', 'date_time CHAR', 'sh INT',
            'rte INT', 'cap INT', 'loss INT', 'window INT',
            'Pin float', 'Pout float', 'RRP float', 'Q_start float',
            'Q_end float', 'obj float']
    sm.create_table(table, db, cols)
    cols_index = ['state', 'date_time', 'rte', 'cap', 'sh', 'loss']
    sm.create_unique_index(db,table, cols_index, 'idx')

    rrp_yearly = gd.load_rrp_cal(Year, State).append(
        gd.load_rrp_actual(datetime(year=Year+1, month=1, day=1),State)['spot_price'])
    rrp_yearly[rrp_yearly>price_cap]=price_cap
    
    for Month in np.arange(1,13,1):
        forecasts_monthly = gd.load_rrp_forecast_monthly(State, Year, Month)
        if Month<12:
              time_period = pd.date_range(
                          datetime(year=Year, month=Month, day=1, hour=0, minute=0),
                          datetime(year=Year, month=Month+1, day=1),
                          # datetime(year=Year, month=Month, day=2),
                          freq='30min')
        else:
              time_period = pd.date_range(
                          datetime(year=Year, month=Month, day=1, hour=0, minute=0),
                          datetime(year=Year+1, month=1, day=1),
                          freq='30min')
                    
        df= pd.DataFrame()
        
        # determine whether you want to cap the spot price market or not:
        simparams = load_simparams()
        Q_init=0
        obj = 0
        for Time in time_period:
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

            results = optimise_parallel(State+str(Year)+str(RTE),simparams)
            Pin = results['Pin'][0]
            Pout = results['Pout'][0]
            
            price_spot = rrp_yearly[Time]
            obj = (results['Pout'][0]-results['Pin'][0]) * price_spot * simparams['dt']
            Q_start = Q_init
            Q_end = Q_init + (Pin*simparams['eta_in'] - Pout/simparams['eta_out'])* simparams['dt']
            Q_init = Q_end   

            df = df.append(
                        pd.DataFrame(
                                    [[State,str(Time), SH, RTE,
                                      price_cap, loss, intervals,
                                      Pin, Pout,price_spot, Q_start,
                                      Q_end, obj]]
                                    )
                        )
            print(State, Time, 'RTE=',float(RTE)/100, 'SH=',SH)
            
        cols_recording = ['state', 'date_time', 'sh',
                          'rte', 'cap', 'loss', 'window',
                          'Pin', 'Pout', 'RRP', 'Q_start',
                          'Q_end', 'obj']
        sm.replace_into_db(df, db, table, cols_recording)
    return(Time)



def calc_S_with_daily_forecast(Year, state, Time=time(4,0,0), RTE=90, SH=10,
                               Loss=0, cap=14500, Horizon=48, write=True ,
                               #DB='SV_variable_forecast_horizon_detailed_%s.db',
                               DB='SV_variable_forecast_horizon_%s.db',
                               detailed=False):
    """This function calculates the battery performance for a year under
    optimisation with daily forecast prices for the next trading day.
    'state': e.g. 'NSW
    'Year': e.g. 2018
    'time': the time at which the forecast data is generated.
    if 'write=True', the function writes the data into the database that is defined by 'db'
    """
    db=DB%(state)
    cols_recording = ['state', 'forecast_time', 'forecast_depth', 'date_time',
                      'sh', 'rte', 'loss','cap', 'Pin','Pout', 'Q', 'RRP',
                      'obj_forecast', 'obj_actual']
        
    DATA = pd.DataFrame(columns = cols_recording)
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
                    datetime(year=Year, month=Month, day=1, hour=Time.hour, minute=Time.minute),
                    datetime(year=Year, month=Month+1, day=1))
        else:
              Date_Range = pd.date_range(
                    datetime(year=Year, month=Month, day=1, hour=Time.hour, minute=Time.minute),
                    datetime(year=Year, month=12, day=31))
              
        for Date in Date_Range:
            # Load the forecast data
            c_forecast = (
                forecast_prices[forecast_prices.generated == Date]
                ).iloc[0:Horizon]
            time_stamps = c_forecast.date_time
            #pdb.set_trace()
            if c_forecast.empty:
                  #pdb.set_trace()
                  continue
            simparams['c'] = c_forecast.price.tolist()
            simparams['N'] = len(simparams['c'])
            simparams['Q_init_frac'] = Q_init/(SH*simparams['Pout_max']/simparams['eta_out'])
            print(Date)
            results = optimise_parallel(str(Year)+str(state)+str(RTE)+\
                                        str(Horizon),simparams)
            obj_forecast = simparams['dt'] * (results['Pout']-results['Pin']) *\
                             c_forecast['price']
            rrp = rrp_yearly[time_stamps].tolist()
            obj_actual = simparams['dt'] * (results['Pout']-results['Pin']) *\
                             rrp
            
            time_stamp_to_record = pd.date_range(time_stamps.tolist()[0], periods=48,
                                                 freq='0.5H')
            n = len(time_stamp_to_record)
            DATA = DATA.append(pd.DataFrame({'state': state,
                                             'forecast_time': str(Time),
                                             'forecast_depth': len(c_forecast),
                                             'date_time':time_stamp_to_record,
                                             'sh': SH,
                                             'rte': RTE,
                                             'loss': Loss,
                                             'cap': cap,
                                             'Pin': results['Pin'][0:n],
                                             'Pout': results['Pout'][0:n],
                                             'Q': results['Q'][1:n+1],
                                             'RRP': rrp[0:n],
                                             'obj_forecast':obj_forecast[0:n],
                                             'obj_actual': obj_actual[0:n]
                                            }), ignore_index=True)
            #pdb.set_trace()
            Q_init = results['Q'][n]
        print('state:%s month=%d-%d RTE=%d completed!'%(state, Month,Year,RTE))
        
    if write:
          if detailed:
                DATA=DATA[cols_recording]
                DATA['date_time'] = DATA.date_time.astype(str)
                sm.replace_into_db(DATA, db, 'storage_value', cols_recording,
                                    create_unique_idx=True,
                                    idx_cols=cols_recording[0:-6])
          else:
                DATA_simple = pd.DataFrame({
                                        'region': [state],
                                        'forecast_time': str(Time),
                                        'forecast_depth': Horizon,
                                        'date_time':Year,
                                        'sh': SH,
                                        'rte': RTE,
                                        'loss': Loss,
                                        'cap': cap,
                                        'obj_forecast':DATA.obj_forecast.sum(),
                                        'obj_actual': DATA.obj_actual.sum()
                                              }
                                           )
                cols=DATA_simple.columns.tolist()
                sm.replace_into_db(DATA_simple,db,'storage_value',cols,
                                   create_unique_idx=True,idx_cols=cols[0:8])
    return(DATA)







def calc_SV_cal_daily_forecast(Year, state, Time=time(4,0,0), Onset=time(23,0,0),
                               RTE=90, SH=10, Loss=0, cap=14500, Horizon=48,
                               write=True,
                               DB='SV_CalDay_forecast_time_%s.db',
                               detailed=False):
    """This function calculates the battery performance for a year under
    optimisation with forecast prices. Day is a calander day. i.e. the storage
    cycle onset is always 00:00.
    'state': e.g. 'NSW
    'Year': e.g. 2018
    'Time': the time at which the forecast data is generated.
    if 'write=True', the function writes the data into the database that is defined by 'db'
    """
    db=DB%(state)
    cols_recording = ['state', 'forecast_time', 'forecast_depth', 'date_time',
                      'sh', 'rte', 'loss','cap', 'Pin','Pout', 'RRP',
                      'obj_forecast', 'obj_actual']
        
    DATA = pd.DataFrame(columns = cols_recording)
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
                    datetime(year=Year, month=Month, day=1, hour=Time.hour, minute=Time.minute),
                    datetime(year=Year, month=Month+1, day=1))
        else:
              Date_Range = pd.date_range(
                    datetime(year=Year, month=Month, day=1, hour=Time.hour, minute=Time.minute),
                    datetime(year=Year, month=12, day=31))
              
        for Date in Date_Range:
            # Load the forecast data
            c_forecast =(forecast_prices[forecast_prices.generated == Date])
            ONSET = Date.replace(hour=Onset.hour,minute=Onset.minute)
            c_forecast = c_forecast[
                        (c_forecast.date_time>ONSET)&
                        (c_forecast.date_time<=(ONSET+timedelta(days=1)))
                                    ]
             
            time_stamps = c_forecast.date_time
            if c_forecast.empty:
                  #pdb.set_trace()
                  continue
            simparams['c'] = c_forecast.price.tolist()
            simparams['N'] = len(simparams['c'])
            simparams['Q_init_frac'] = Q_init/(SH*simparams['Pout_max']/simparams['eta_out'])
            print(Date)
            results = optimise_parallel(str(Year)+str(Time).replace(':','_')+str(state)+str(RTE)+\
                                        str(Horizon),simparams)
            obj_forecast = simparams['dt'] * (results['Pout']-results['Pin']) *\
                             c_forecast['price']
            rrp = rrp_yearly[time_stamps].tolist()
            obj_actual = simparams['dt'] * (results['Pout']-results['Pin']) *\
                             rrp
            
            time_stamp_to_record = pd.date_range(time_stamps.tolist()[0], periods=48,
                                                 freq='0.5H')
            n = len(time_stamp_to_record)
            DATA = DATA.append(pd.DataFrame({'state': state,
                                             'forecast_time': str(Time),
                                             'forecast_depth': len(c_forecast),
                                             'date_time':time_stamp_to_record,
                                             'onset': Onset,
                                             'sh': SH,
                                             'rte': RTE,
                                             'loss': Loss,
                                             'cap': cap,
                                             'Pin': results['Pin'][0:n],
                                             'Pout': results['Pout'][0:n],
                                             'RRP': rrp[0:n],
                                             'obj_forecast':obj_forecast[0:n],
                                             'obj_actual': obj_actual[0:n]
                                            }), sort=False, ignore_index=True)
            Q_init = results['Q'][-1]
        print('state:%s month=%d-%d RTE=%d completed!'%(state, Month,Year,RTE))
        
    if write:
          if detailed:
                DATA=DATA[cols_recording]
                DATA['date_time'] = DATA.date_time.astype(str)
                sm.replace_into_db(DATA, db, 'storage_value', cols_recording,
                                    create_unique_idx=True,
                                    idx_cols=cols_recording[0:-4])
          else:
                DATA_simple = pd.DataFrame({
                                        'region': [state],
                                        'forecast_time': str(Time),
                                        'forecast_depth': Horizon,
                                        'date_time':Year,
                                        'onset': str(Onset),
                                        'sh': SH,
                                        'rte': RTE,
                                        'loss': Loss,
                                        'cap': cap,
                                        'obj_forecast':DATA.obj_forecast.sum(),
                                        'obj_actual': DATA.obj_actual.sum()
                                              }
                                           )
                cols=DATA_simple.columns.tolist()
                sm.replace_into_db(DATA_simple,db,'storage_value',cols,
                                   create_unique_idx=True,idx_cols=cols[0:9])
    return(DATA)
    
    
    
    




def calc_S_perfect_foresight(State, Year, RTE, SH, Loss, Cap,
                             write=False, db='storage_value_perfect',
                             table='storage_value'):
    """This function calculates the storage value with perfect foresight.
        write: specifies whether to write the data into the database or not.
        db: the database to write
        table: the table in the db to write
    """

    simparams = load_simparams()
      
    cols_recording = ['state', 'date_time', 'sh', 'rte',
                      'loss','cap', 'Pin','Pout','RRP']
    idx = ['state', 'date_time', 'rte', 'cap', 'sh', 'loss']
      
    print(State, Year, float(RTE)/100, SH)
      
    c = gd.load_rrp_cal(Year, State)
    time_stamp_to_record = c.index.astype(str)
    c[c>Cap]=Cap
    c = c.tolist()
      
    simparams['c'] = c
    simparams['N'] = len(c)
    simparams['SH'] = SH
    simparams['eta_out'] = float(RTE)/100
    simparams['loss'] = float(Loss)/100
      
    results = optimise_parallel(str(Year)+State+str(RTE)+str(SH), simparams)

    # register the generated data into the database
    if write:
          n = len(time_stamp_to_record)
          DATA = pd.DataFrame({'state': State,
                               'date_time':time_stamp_to_record,
                               'sh': SH,
                               'rte': RTE,
                               'loss': Loss,
                               'cap': Cap,
                               'Pin': results['Pin'][0:n],
                               'Pout': results['Pout'][0:n],
                               'RRP': c[0:n]
                               })
          DB=db+'_%s.db'%(State)
          sm.create_table(table, DB, cols_recording,
                          create_unique_idx=True, idx_cols=idx)
          sm.replace_into_db(DATA, DB,table, cols_recording)
          input('hello ...')
    return(results)




    
def calc_simple_SV_PF(State, Year, RTEs, SHs, Losses, Caps, Pin_maxs, Pout_maxs,
                             write=False, db='storage_value_PF_simple',
                             table='storage_value'):
    """This function calculates the storage value with perfect foresight.
        It does not record the half-hourly data.
        write: specifies whether to write the data into the database or not.
        db: the database to write
        table: the table in he db to write
    """

    simparams = load_simparams()
      
    cols_recording = ['state', 'date_time', 'sh', 'rte',
                      'loss','cap', 'Pin_max','Pout_max', 'SV']
    idx = ['state', 'date_time', 'rte', 'cap', 'sh', 'loss',
           'Pin_max', 'Pout_max']
    DATA = pd.DataFrame(columns=cols_recording)
    rrps = gd.load_rrp_cal(Year, State)
    for RTE in RTEs:
          for SH in SHs:
                for Loss in Losses:
                      for Cap in Caps:
                          for Pin_max in Pin_maxs:
                              for Pout_max in Pout_maxs:
                                  print(State, Year, 'RTE=%d'%(RTE),
                                        'SH=%d'%(SH), 'Pin_max=%d'%(Pin_max),
                                        'Pout_max=%d'%(Pout_max))
                                  c = rrps
                                  c[c>Cap]=Cap
                                  c = c.tolist()
                                  simparams['Pout_max']=Pout_max
                                  simparams['Pin_max']=Pin_max
                                  simparams['c'] = c
                                  simparams['N'] = len(c)
                                  simparams['SH'] = SH
                                  simparams['eta_out'] = float(RTE)/100
                                  simparams['loss'] = float(Loss)/100
                                  results = optimise_parallel(
                                              State+str(Year)+str(RTE)+
                                              str(SH).replace('.','_'),
                                              simparams)
                                  SV =  sum(0.5*(results['Pout']-results['Pin']) *c)
                                  # register the generated data into the database
                                  if write:
                                      DATA = DATA.append({'state': State,
                                                           'date_time':Year,
                                                           'sh': SH,
                                                           'rte': RTE,
                                                           'loss': Loss,
                                                           'cap': Cap,
                                                           'SV': SV,
                                                           'Pin_max':Pin_max,
                                                           'Pout_max':Pout_max
                                                           }, ignore_index=True)
    DB=db+'_%s.db'%(State)
    sm.create_table(table, DB, cols_recording,
                    create_unique_idx=True, idx_cols=idx)
    sm.replace_into_db(DATA, DB,table, cols_recording)
    string = ('%s %d completed!'%(State, Year))
    return([string,results][write==False])





def calc_simple_SV_PF_Eout(State, Year, RTEs, SHs, Losses, Caps, Pin_maxs, Pout_maxs,
                             write=False, db='storage_value_PF_with_Eout',
                             table='storage_value'):
    """This function calculates the storage value with perfect foresight.
        It does not record the half-hourly data.
        It stores the energy the storage generates.
        write: specifies whether to write the data into the database or not.
        db: the database to write
        table: the table in he db to write
    """

    simparams = load_simparams()
      
    cols_recording = ['state', 'date_time', 'sh', 'rte',
                      'loss','cap', 'Pin_max','Pout_max', 'SV', 'Eout']
    idx = ['state', 'date_time', 'rte', 'cap', 'sh', 'loss',
           'Pin_max', 'Pout_max']
    DATA = pd.DataFrame(columns=cols_recording)
    rrps = gd.load_rrp_cal(Year, State)
    for RTE in RTEs:
          for SH in SHs:
                for Loss in Losses:
                      for Cap in Caps:
                          for Pin_max in Pin_maxs:
                              for Pout_max in Pout_maxs:
                                  print(State, Year, 'RTE=%d'%(RTE),
                                        'SH=%d'%(SH), 'Pin_max=%d'%(Pin_max),
                                        'Pout_max=%d'%(Pout_max))
                                  c = rrps
                                  c[c>Cap]=Cap
                                  c = c.tolist()
                                  simparams['Pout_max']=Pout_max
                                  simparams['Pin_max']=Pin_max
                                  simparams['c'] = c
                                  simparams['N'] = len(c)
                                  simparams['SH'] = SH
                                  simparams['eta_out'] = float(RTE)/100
                                  simparams['loss'] = float(Loss)/100
                                  results = optimise_parallel(
                                              State+str(Year)+str(RTE)+
                                              str(SH).replace('.','_'),
                                              simparams)
                                  SV =  sum(0.5*(results['Pout']-results['Pin']) *c)
                                  E_out = sum(0.5*(results['Pout']))
                                 # register the generated data into the database
                                  if write:
                                      DATA = DATA.append({'state': State,
                                                           'date_time':Year,
                                                           'sh': SH,
                                                           'rte': RTE,
                                                           'loss': Loss,
                                                           'cap': Cap,
                                                           'SV': SV,
                                                           'Eout':E_out,
                                                           'Pin_max':Pin_max,
                                                           'Pout_max':Pout_max
                                                           }, ignore_index=True)
    DB=db+'_%s.db'%(State)
    sm.create_table(table, DB, cols_recording,
                    create_unique_idx=True, idx_cols=idx)
    sm.replace_into_db(DATA, DB,table, cols_recording)
    string = ('%s %d completed!'%(State, Year))
    return([string,results][write==False])
    
    





    

def calc_simple_SV_PF_PJM_da_rt(zone_list, Year, RTEs, SHs, Losses, Pin_maxs, Pout_maxs,
                             write=False, db='storage_value_PF_simple_PJM_dart.db',
                             table='storage_value'):
    """This function calculates the storage value with perfect foresight for
       day ahead and real time markets in the PJM simultaneously.
        It does not record the hourly LMP values.
        write: specifies whether to write the data into the database or not.
        db: the database to write
        table: the table in he db to write
    """

    # simparams = load_simparams()
    simparams = {}
    cols_recording = ['zone', 'date_time', 'sh', 'rte',
                      'loss', 'Pin_max','Pout_max', 'SV']
    idx = ['zone', 'date_time', 'rte', 'sh', 'loss',
           'Pin_max', 'Pout_max']
    DATA = pd.DataFrame(columns=cols_recording)
    # lmps_da = gd.load_lmp_cal(Year, lmp_type='da')
    # lmps_rt = gd.load_lmp_cal(Year, lmp_type='rt')
    for Zone in zone_list:
          for RTE in RTEs:
                for SH in SHs:
                      for Loss in Losses:
                          for Pin_max in Pin_maxs:
                              for Pout_max in Pout_maxs:
                                  print(Zone, Year, 'RTE=%d'%(RTE),
                                        'SH=%d'%(SH), 'Pin_max=%d'%(Pin_max),
                                        'Pout_max=%d'%(Pout_max))
                                  # LMP_da = lmps_da[lmps_da['pnode_name']==Zone].\
                                  # set_index('datetime_beginning_ept')\
                                  # ['total_lmp_da'].sort_index()
                                  # LMP_rt = lmps_rt[lmps_rt['pnode_name']==Zone].\
                                  # set_index('datetime_beginning_ept')\
                                  # ['total_lmp_rt'].sort_index()
                                  #pdb.set_trace()
                                  # LMP_da = LMP_da.tolist()
                                  # LMP_rt = LMP_rt.tolist()
                                  simparams['Zone']=Zone
                                  simparams['Pout_max']=Pout_max
                                  simparams['Pin_max']=Pin_max
                                   # simparams['c_da'] = LMP_da
                                   # simparams['c_rt'] = LMP_rt
                                   # simparams['N'] = len(LMP_da)
                                  simparams['SH'] = SH
                                  simparams['Year'] = Year
                                  simparams['RTE'] = RTE 
                                  # simparams['eta_out'] = float(RTE)/100
                                  simparams['Loss'] = float(Loss)/100
                                  # simparams['dt'] = 1
                                  
                                  results = optimise_PJM_parallel(**simparams)
                                  SV =  results['obj']
                                  # register the generated data into the database
                                  if write:
                                      DATA = DATA.append({'zone': Zone,
                                                          'date_time':Year,
                                                          'sh': SH,
                                                          'rte': RTE,
                                                          'loss': Loss,
                                                          'SV': float(SV),
                                                          'Pin_max':Pin_max,
                                                          'Pout_max':Pout_max
                                                          }, ignore_index=True)
    sm.create_table(table, db, cols_recording,
                    create_unique_idx=True, idx_cols=idx)
    sm.replace_into_db(DATA, db,table, cols_recording)
    string = ('%s %d completed!'%(Zone, Year))
    return(string)
    
    


def calc_simple_SV_PF_const_Ec(State, Year, RTEs, Ecs, Losses, Caps, Pin_maxs, Pout_maxs,
                             write=False, db='storage_value_PF_simple_const_Ec',
                             table='storage_value'):
    """This function calculates the storage value with perfect foresight
        for constant storage capacity (Ec).
        It does not record the half-hourly data.
        write: specifies whether to write the data into the database or not.
        db: the database to write
        table: the table in db to write
    """

    simparams = load_simparams()
      
    cols_recording = ['state', 'date_time', 'Ec', 'rte',
                      'loss','cap', 'Pin_max','Pout_max', 'SV']
    idx = ['state', 'date_time', 'rte', 'cap', 'Ec', 'loss',
           'Pin_max', 'Pout_max']
    DATA = pd.DataFrame(columns=cols_recording)
    rrps = gd.load_rrp_cal(Year, State)
    for RTE in RTEs:
          for Ec in Ecs:
                for Loss in Losses:
                      for Cap in Caps:
                          for Pin_max in Pin_maxs:
                              for Pout_max in Pout_maxs:
                                  print(State, Year, 'RTE=%d'%(RTE),
                                        'Ec=%d'%(Ec), 'Pin_max=%d'%(Pin_max),
                                        'Pout_max=%d'%(Pout_max))
                                  c = rrps
                                  c[c>Cap]=Cap
                                  c = c.tolist()
                                  simparams['Pout_max']=Pout_max
                                  simparams['Pin_max']=Pin_max
                                  simparams['c'] = c
                                  simparams['N'] = len(c)
                                  simparams['SH'] = Ec/Pout_max
                                  simparams['eta_out'] = float(RTE)/100
                                  simparams['loss'] = float(Loss)/100
                                  results = optimise_parallel(State+str(Year)+str(Pin_max)+str(Pout_max), simparams)
                                  SV =  sum(0.5*(results['Pout']-results['Pin']) *c)
                                  print('Pin=',Pin_max,'Pout=',Pout_max, 'SV=',SV)
                                 # register the generated data into the database
                                  if write:
                                      DATA = DATA.append({'state': State,
                                                           'date_time':Year,
                                                           'Ec': Ec,
                                                           'rte': RTE,
                                                           'loss': Loss,
                                                           'cap': Cap,
                                                           'SV': SV,
                                                           'Pin_max':Pin_max,
                                                           'Pout_max':Pout_max
                                                           }, ignore_index=True)
    DB=db+'_%s.db'%(State)
    sm.create_table(table, DB, cols_recording,
                    create_unique_idx=True, idx_cols=idx)
    sm.replace_into_db(DATA, DB,table, cols_recording)
    string = ('%s %d completed!'%(State, Year))
    return(string)



def calc_simple_SV_PF_const_Ec_comparison(State, Year, RTE, Ec, Loss, 
                                          Cap, Pin_max, Pout_max):
    """This function calculates the storage value with perfect foresight
        for constant storage capacity (Ec).
        It does not record the half-hourly data.
        write: specifies whether to write the data into the database or not.
        db: the database to write
        table: the table in db to write
    """

    simparams = load_simparams()
      
    rrps = gd.load_rrp_cal(Year, State)

    print(State, Year, 'RTE=%d'%(RTE),
          'Ec=%d'%(Ec), 'Pin_max=%d'%(Pin_max),
          'Pout_max=%d'%(Pout_max))
    c = rrps
    c[c>Cap]=Cap
    c = c.tolist()
    simparams['Pout_max']=Pout_max
    simparams['Pin_max']=Pin_max
    simparams['c'] = c
    simparams['N'] = len(c)
    simparams['SH'] = Ec/Pout_max
    simparams['eta_out'] = float(RTE)/100
    simparams['loss'] = float(Loss)/100
    results = optimise_parallel(str(Year), simparams)
    SV =  sum(0.5*(results['Pout']-results['Pin']) *c)
                                 # register the generated data into the database

    print('%s %d completed!'%(State, Year))
    return(SV)



def calc_SV_rolling_horizon_NEM(Year, State, SH, RTE, Loss, cap, Time_gen, Window,
                            write=True):
      """
      This function calculates and records detailed storage performance data
      for cases that the storage has access to the exact RRPs into the future.
      The depth of foresight into the future is determined by 'Window'.
      The foresight starts at 'Time_gen' and is extended once a day at 'Time_gen'.
      Time_gen is the start of the foresight cycle e.g. time(4,0,0).
      """
      table = 'storage_value'
      cols = ['state CHAR(3)', 'date_time CHAR', 'sh INT',
        'rte INT', 'cap INT', 'loss INT', 'Time_gen CHAR' , 'window INT',
        'SV float']
      db='storage_value_short_foresight_%s.db'%(State)
      sm.create_table(table, db, cols)
      
      rrp = gd.load_rrp_cal(Year, State).append(
                        gd.load_rrp_cal(Year+1, State))
      rrp[rrp>cap]=cap
      cols_recording = ['state', 'date_time', 'sh', 'rte', 'loss','cap',
                        'Time_gen', 'window', 'SV']
      cols_DATA = ['state', 'date_time', 'sh', 'rte', 'loss','cap',
                   'Time_gen', 'window', 'revenue']
        
      DATA = pd.DataFrame(columns = cols_DATA )
    
      time_period = pd.date_range(
                          datetime(year=Year, month=1, day=1,
                                   hour=Time_gen.hour, minute=Time_gen.minute),
                          datetime(year=Year, month=12, day=28, hour=0, minute=0),
                          freq='1D')
      
      simparams = load_simparams()
      Q_init=0
      for Time in time_period:
            time_stamps = pd.date_range(Time,Time+timedelta(days=Window),
                                        freq='0.5H')  
            price_list  = rrp.loc[time_stamps]
            c=price_list.tolist()

            if not c:
                continue
            simparams['c'] = c
            simparams['dt'] = 0.5
            simparams['N'] = len(c)
            simparams['SH'] = SH
            simparams['eta_out'] = float(RTE)/100
            simparams['loss'] = float(Loss)/100
            simparams['Q_init_frac'] = Q_init/(SH*simparams['Pout_max']/simparams['eta_out'])
            ID = State+str(Year)+str(Window)+(Time).strftime('%Y%m%d%H%M')+str(RTE)
            results = optimise_parallel(ID,simparams)
            Pin = results['Pin'][0:48]
            Pout = results['Pout'][0:48]
            Q = results['Q'][0:49]
            
            revenue = (Pout - Pin) * c[0:48] * simparams['dt']
            Q_init = Q[-1]

            time_stamp_to_record = time_stamps[0:48]
            DATA = DATA.append(pd.DataFrame({'state': State,
                                             'date_time':time_stamp_to_record,
                                             'sh': SH,
                                             'rte': RTE,
                                             'loss': Loss,
                                             'cap': cap,
                                             'Time_gen': str(Time_gen),
                                             'window': Window,
                                             'revenue': revenue
                                            }), ignore_index=True)
            
            print(State, Time, 'RTE=',float(RTE)/100, 'SH=',SH)
      if write:
        data_to_record= pd.DataFrame({'state': State,
                                      'date_time':Year,
                                      'sh': SH,
                                      'rte': RTE,
                                      'loss': Loss,
                                      'cap': cap,
                                      'Time_gen': str(Time_gen),
                                      'window': Window,
                                      'SV': DATA.revenue.sum()
                                      }, index=[0])
        # pdb.set_trace()
        data_to_record = data_to_record[cols_recording]
        sm.replace_into_db(data_to_record, db, 'storage_value', cols_recording,
            create_unique_idx=True,
            idx_cols=cols_recording[0:-1])
      print('Completed at ', pd.Timestamp.now())  
      return()






def calc_SV_rolling_horizon_PJM(Year, Region, SH, RTE, Loss, Time_gen, Window,
                            write=True):
      """
      This function calculates and records detailed storage performance data
      for cases that the storage has access to the exact RT and DA prices 
      into the future.
      The depth of foresight into the future is determined by 'Window'.
      The foresight starts at 'Time_gen' and is extended once a day at 'Time_gen'.
      Time_gen is the start of the foresight cycle e.g. time(4,0,0) .
      """
      curdir =  optdir + "arbitrage%s" %['/', '\\'][windows]
      mzdir = ['/home/jeff/local/software/MiniZinc/bin/',
             'C:\\Program Files\\minizinc\\'][windows]
      db = 'storage_value_short_foresight_%s.db'%(Region)
      cols_recording = ['region', 'date_time', 'sh', 'rte', 'loss',
                        'Time_gen', 'window', 'SV']
      cols_DATA = ['region', 'date_time', 'sh', 'rte', 'loss',
                   'Time_gen', 'window', 'revenue']
        
      DATA = pd.DataFrame(columns = cols_DATA )
    
      time_period = pd.date_range(
                          datetime(year=Year, month=1, day=1,
                                   hour=Time_gen.hour, minute=Time_gen.minute),
                          datetime(year=Year, month=12, day=30, hour=0, minute=0),
                          freq='1D')
      
      lmps_da = gd.load_lmp_zone_year(Region, Year, lmp_type='da', Year_type='cal').\
              append(gd.load_lmp_zone_year(Region, Year+1, lmp_type='da', Year_type='cal'))  
    
      lmps_rt = gd.load_lmp_zone_year(Region, Year, lmp_type='rt', Year_type='cal').\
              append(gd.load_lmp_zone_year(Region, Year+1, lmp_type='rt', Year_type='cal'))  
    
      
      simparams = load_simparams()
      
      
      Q_init_frac=0
      for Time in time_period:
            time_stamps = pd.date_range(Time,Time+timedelta(days=Window),
                                        freq='1H')
            LMP_da = lmps_da.reindex(time_stamps).fillna(method='ffill')
            LMP_rt = lmps_rt.reindex(time_stamps).fillna(method='ffill')
            
            simparams['c_da'] = LMP_da.tolist()
            simparams['c_rt'] = LMP_rt.tolist()
            simparams['N'] = len(LMP_da)
            simparams['SH'] = SH
            simparams['eta_out'] = float(RTE)/100
            simparams['loss'] = float(Loss)/100
            simparams['dt'] = 1
            simparams['Pout_max'] = 50
            simparams['Pin_max'] = 50
            # we assume market='DA+RT'
            simparams['Pda_frac'] = 1
            simparams['Prt_frac'] = 1
            simparams['Q_init_frac'] = Q_init_frac
            
            
            
            ID = str(Year)+Region.replace('/','') +str(RTE) + str(Loss) +\
                    str(SH) + str(Window)
                    
            make_dzn_file_PJM_parallel(ID, **simparams)
          
            output = str(check_output([mzdir + 'minizinc', "--soln-sep", '""', 
                                     "--search-complete-msg", '""', "--solver", 
                                     "COIN-BC", curdir + "arbitrage_PJM.mzn", 
                                     curdir + "arbitrage_PJM_%s.dzn"%(ID)]))
          
            output = output.replace("""\\r\\n""\\r\\n""\\r\\n\'""", '')
            output = output.replace("b\'[", '')
            output = output.replace('[', '')
            output = output.replace(']', '')  
            output = output.replace("b'", '')
            if not windows: output = output.replace('\\n""\\n""\\n\'','')
      
            obj, Pin_da, Pin_rt, Pout_da, Pout_rt, Q = output.split(';')
            obj = array(obj).astype(float)
            Q = array(Q.split(',')).astype(float)
            Pin_da = array(Pin_da.split(',')).astype(float)
            Pin_rt = array(Pin_rt.split(',')).astype(float)
            Pout_da = array(Pout_da.split(',')).astype(float)
            Pout_rt = array(Pout_rt.split(',')).astype(float)
          
            results = {}
            results['Q'] = Q
            results['Pin_da'] = Pin_da
            results['Pin_rt'] = Pin_rt
            results['Pout_da'] = Pout_da
            results['Pout_rt'] = Pout_rt
            results['lmp_da'] = LMP_da
            results['lmp_rt'] = LMP_rt
            results['obj'] = obj.item(0)
            results['date_time'] = time_stamps
#            print(ID, ' ', results['obj'])
            
            Pin_da = results['Pin_da'][0:24]
            Pout_da = results['Pout_da'][0:24]
            Pin_rt = results['Pin_rt'][0:24]
            Pout_rt = results['Pout_rt'][0:24]
            Q = results['Q'][0:25]
            revenue = (Pout_da - Pin_da) * LMP_da[0:24] * simparams['dt']+\
                      (Pout_rt - Pin_rt) * LMP_rt[0:24] * simparams['dt']
            Q_init_frac = Q[-1]/(50*SH*100/RTE)

            time_stamp_to_record = time_stamps[0:24]
            DATA = DATA.append(pd.DataFrame({'region': Region,
                                             'date_time':time_stamp_to_record,
                                             'sh': SH,
                                             'rte': RTE,
                                             'loss': Loss,
                                             'Time_gen': str(Time_gen),
                                             'window': Window,
                                             'revenue': revenue
                                            }), ignore_index=True)
            
            print(Region, Time, 'RTE=',float(RTE)/100, 'SH=',SH, 'Window=',
                  Window)
      if write:
        data_to_record= pd.DataFrame({'region': Region,
                                      'date_time':Year,
                                      'sh': SH,
                                      'rte': RTE,
                                      'loss': Loss,
                                      'Time_gen': str(Time_gen),
                                      'window': Window,
                                      'SV': DATA.revenue.sum()
                                      }, index=[0])
        data_to_record = data_to_record[cols_recording]
        sm.replace_into_db(data_to_record, db, 'storage_value', cols_recording,
            create_unique_idx=True,
            idx_cols=cols_recording[0:-1])
      print('Completed at ', pd.Timestamp.now())  
      return()
      
      
      




def calc_SV_naive_forecast(Year, State, SH, RTE, Loss, cap, Time_gen, write=True):
      """
      This function calculates and records detailed storage performance data
      for cases that the storage naive forecast in which the RRP of day n
      is assumed to be the same as day n-1.
      The foresight starts at 'Time_gen' and is extended once a day at 'Time_gen'.
      Time_gen is the start of the foresight cycle e.g. time(4,0,0) .
      """
      db = 'storage_value_naive_foresight_%s.db'%(State)
      table = 'storage_value'
      cols = ['state CHAR(3)', 'date_time CHAR', 'sh INT',
        'rte INT', 'cap INT', 'loss INT', 'Time_gen CHAR',
        'SV float']

      sm.create_table(table, db, cols)
      
      rrp = gd.load_rrp_cal(Year, State).append(
                        gd.load_rrp_cal(Year+1, State))
      rrp[rrp>cap]=cap
      
      cols_recording = ['state', 'date_time', 'sh', 'rte', 'loss','cap',
                        'Time_gen', 'SV']
      cols_DATA = ['state', 'date_time', 'sh', 'rte', 'loss','cap',
                   'Time_gen', 'revenue']
        
      DATA = pd.DataFrame(columns = cols_DATA )
    
      time_period = pd.date_range(
                          datetime(year=Year, month=1, day=1,
                                   hour=Time_gen.hour, minute=Time_gen.minute),
                          datetime(year=Year+1, month=1, day=1, hour=0, minute=0),
                          freq='1D')
      simparams = load_simparams()
      Q_init=0
      for Time in time_period[1:]:
            time_stamps_today = pd.date_range(Time,Time+timedelta(days=1),
                                        freq='0.5H')  
            time_stamps_yesterday = pd.date_range(Time-timedelta(days=1),Time,
                                        freq='0.5H')  
            rrp_fc = rrp.loc[time_stamps_yesterday].tolist()
            rrp_ac = rrp.loc[time_stamps_today].tolist()
            
            if not rrp_fc:
                continue
            simparams['c'] = rrp_fc
            simparams['N'] = len(rrp_fc)
            simparams['SH'] = SH
            simparams['eta_out'] = float(RTE)/100
            simparams['loss'] = float(Loss)/100
            simparams['Q_init_frac'] = Q_init/(SH*simparams['Pout_max']/simparams['eta_out'])

            results = optimise_parallel(str(Year)+State+str(RTE),simparams)
            Pin = results['Pin'][0:48]
            Pout = results['Pout'][0:48]
            Q = results['Q'][0:49]
            
            revenue = (Pout - Pin) * rrp_ac[0:48] * simparams['dt']
            Q_init = Q[-1]

            time_stamp_to_record = time_stamps_today[0:48]
            DATA = DATA.append(pd.DataFrame({'state': State,
                                             'date_time':time_stamp_to_record,
                                             'sh': SH,
                                             'rte': RTE,
                                             'loss': Loss,
                                             'cap': cap,
                                             'Time_gen': str(Time_gen),
                                             'revenue': revenue
                                            }), ignore_index=True)
            
            print(State, Time, 'RTE=',float(RTE)/100, 'SH=',SH)
      if write:
        data_to_record= pd.DataFrame({'state': State,
                                      'date_time':Year,
                                      'sh': SH,
                                      'rte': RTE,
                                      'loss': Loss,
                                      'cap': cap,
                                      'Time_gen': str(Time_gen),
                                      'SV': DATA['revenue'].sum()
                                      }, index=[0])
        data_to_record = data_to_record[cols_recording]
        sm.replace_into_db(data_to_record, db, 'storage_value', cols_recording,
            create_unique_idx=True,
            idx_cols=cols_recording[0:-1])
      print('Completed at ', pd.Timestamp.now())  
      return()




def calc_SV_synthetic_forecast(Year, state, Time=time(0,0,0), Onset=time(0,0,0),
                               RTE=90, SH=10, Loss=0, Pin_max=50, Pout_max=50, Horizon=48,
                               STD=0.1, N=1,
                               write=True,
                               DB='SV_daily_synthetic_forecast_%s.db',
                               detailed=False):
    """This function calculates the battery performance for a year under
    optimisation with synthetic forecast prices. 
    'state': e.g. 'NSW
    'Year': e.g. 2018
    'Time': the time at which the forecast data is generated.
    if 'write=True', the function writes the data into the database that is defined by 'db'
    """
    db=DB%(state)
    cols_recording = ['N','state', 'forecast_time', 'forecast_depth', 'forecast_std',
                      'date_time','sh', 'rte', 'loss', 'Pin','Pout','Q_init', 'RRP',
                      'obj_actual']
        
    DATA = pd.DataFrame(columns = cols_recording)
    simparams = load_simparams()
    simparams['SH'] = SH
    simparams['Pin_max'] = Pin_max
    simparams['Pout_max'] = Pout_max
    simparams['eta_out'] = float(RTE)/100
    simparams['loss'] = float(Loss)/100
    
    rrp_yearly = gd.load_rrp_cal(Year, state).append(
             gd.load_rrp_cal(Year+1, state)[0:480])
    
    
    Q_init=0
    forecast = pd.DataFrame(rrp_yearly)
    forecast['MEAN'] = 1
    forecast['STD'] = STD
    forecast['random'] = np.random.normal(loc=forecast.MEAN, scale=forecast.STD)
    forecast['price'] = forecast.spot_price*forecast.random
    
    Dates = pd.date_range(pd.datetime(Year,1,1), pd.datetime(Year+1,1,1))[0:-1]
    
    for DAY in Dates:
          ONSET = DAY.replace(hour=Onset.hour,minute=Onset.minute)
          c_forecast = forecast[
                        (forecast.index>=ONSET)&
                        (forecast.index<(ONSET+timedelta(hours = Horizon*0.5)))
                                    ]
             
          time_stamps = c_forecast.index
          time_stamp_to_record = time_stamps[0:49]
          n = len(time_stamp_to_record)
          
          if c_forecast.empty:
                #pdb.set_trace()
                continue
          
          simparams['c'] = c_forecast.price.tolist()
          simparams['N'] = len(simparams['c'])
          simparams['Q_init_frac'] = Q_init/(SH*simparams['Pout_max']/simparams['eta_out'])
          
          print(DAY)
          results = optimise_parallel(str(Year)+str(Time).replace(':','_')+str(state)+str(RTE)+\
                                      str(Horizon)+str(int(STD*100))+str(N),simparams)
          
          
          rrp = rrp_yearly[time_stamps].tolist()[0:n]
          obj_actual = simparams['dt'] * (
                      (results['Pout'][0:n]-results['Pin'][0:n]) * rrp 
                                              )
            
          
          DATA = DATA.append(pd.DataFrame({'N':N,
                                           'state': state,
                                           'forecast_time': str(Time),
                                           'forecast_depth': len(c_forecast),
                                           'forecast_std': STD,
                                           'date_time':time_stamp_to_record,
                                           'onset': Onset,
                                           'sh': SH,
                                           'rte': RTE,
                                           'loss': Loss,
                                           'Pin': results['Pin'][0:n],
                                           'Pout': results['Pout'][0:n],
                                           'Q_init': results['Q'][0:n],
                                           'RRP': rrp[0:n],
                                           'obj_actual': obj_actual[0:n]
                                           }), sort=False, ignore_index=True)


          Q_init = results['Q'][n]
                      
          
    DATA_simple = pd.DataFrame({ 'N':N, 
                                'region': [state],
                                'forecast_time': str(Time),
                                'forecast_depth': Horizon,
                                'forecast_std': STD,
                                'date_time':Year,
                                'onset': str(Onset),
                                'sh': SH,
                                'rte': RTE,
                                'loss': Loss,
                                'Pin_max':Pin_max,
                                'Pout_max':Pout_max,
                                'obj_actual': DATA.obj_actual.sum()
                                        }
                                     )    
        
    
    if detailed:
          DATA=DATA[cols_recording]
          DATA['date_time'] = DATA.date_time.astype(str)
          output = DATA
          if write:
                sm.replace_into_db(DATA, db, 'storage_value', cols_recording,
                                    create_unique_idx=True,
                                    idx_cols=cols_recording[0:12])
    else:
          output =DATA_simple 
          if write:
                cols=DATA_simple.columns.tolist()
                sm.replace_into_db(DATA_simple,db,'storage_value',cols,
                                   create_unique_idx=True,idx_cols=cols[0:12])
                
    return(output)
    
    
    




def correct_rolling_SV(State, Year, RTE):
    SH=10
    loss=0
    price_cap=14500
    db = 'storage_value_rolling_forecast_%s.db'%(State)
    table = 'storage_value'
    cols = [x[0] for x in sm.list_columns(db,table)]
    data = sm.get_data(cols,table,db)
    data['date_time'] = pd.to_datetime(data['date_time'])
    data.set_index('date_time',inplace=True)
      
    DATA = data[ data['rte']==RTE]
    DATA = DATA.loc[pd.Timestamp(year=Year, month=12, day=1, hour=0, minute=0)]
    
#    pdb.set_trace()
    Q_init = DATA['Q_end']
      
    rrp_yearly = gd.load_rrp_cal(Year, State).append(
                gd.load_rrp_actual(datetime(year=Year+1, month=1, day=1),State)['spot_price'])
          
    for Month in np.arange(12,13,1):
          forecasts_monthly = gd.load_rrp_forecast_monthly(State, Year, Month)
          if Month<12:
                time_period = pd.date_range(
                            datetime(year=Year, month=Month, day=1, hour=0, minute=0),
                            datetime(year=Year, month=Month+1, day=1),
                            freq='30min')
          else:
                time_period = pd.date_range(
                            datetime(year=Year, month=Month, day=1, hour=0, minute=0),
                            datetime(year=Year+1, month=1, day=1),
                            freq='30min')
                          
          df= pd.DataFrame()
      
          # determine whether you want to cap the spot price market or not:
          simparams = load_simparams()
          obj = 0
          for Time in time_period:
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
                obj = (results['Pout'][0]-results['Pin'][0]) * price_spot * simparams['dt']
                Q_start = Q_init
                Q_end = Q_init + (Pin*simparams['eta_in'] - Pout/simparams['eta_out'])* simparams['dt']
                Q_init = Q_end
                df = df.append(
                              pd.DataFrame(
                                          [[State,str(Time), SH, RTE,
                                            price_cap, loss, intervals,
                                            Pin, Pout,price_spot, Q_start,
                                            Q_end, obj]]
                                          )
                              )
                print(State, Time, 'RTE=',float(RTE)/100, 'SH=',SH)
                  
    cols_recording = ['state', 'date_time', 'sh',
                              'rte', 'cap', 'loss', 'window',
                              'Pin', 'Pout', 'RRP', 'Q_start',
                              'Q_end', 'obj']
    sm.replace_into_db(df, db, table, cols_recording)
    return()



if __name__ == "__main__":
    
    state = "SA"
    year = 2014
    window = 3


    date_curtailment = pd.date_range(start = str(year), end = str(year + 1), 
                                     freq = 'H')
    
    date_rrp = pd.date_range(start = str(year), end = str(year + 1))
        
    rrp = pd.DataFrame()
    for Day in date_rrp:
        rrp = rrp.append(pd.DataFrame(gd.load_rrp_actual(Day, state)['spot_price']))
    
    c=rrp['spot_price'].tolist()

    simparams = load_simparams()

    simparams['c'] = c
    simparams['N'] = len(c)
    
    results = optimise(simparams)






#######################################################
#####  Start of electrifying coal-fired power plants #####
#######################################################
    
    
def optimise_ECPP(ID,simparams):
    """simparams is a dictionary containing the storage
    parameters for the dzn file."""    
    curdir =  optdir + "arbitrage%s" %['/', '\\'][windows]
    mzdir = ['/home/jeff/local/software/MiniZinc/bin/',
             'C:\\Program Files\\minizinc\\'][windows]
    make_dzn_file_ECPP(ID,**simparams)

    output = str(check_output([mzdir + 'minizinc', "--soln-sep", '""',
                               "--search-complete-msg", '""', "--solver",
                               "COIN-BC", curdir + "arbitrage_ECPP.mzn",
                               curdir + "arbitrage_ECPP_%s.dzn"%(ID)   ]))

    output = output.replace("""\\r\\n""\\r\\n""\\r\\n\'""", '')
    output = output.replace("b\'[", '')
    output = output.replace('[', '')
    output = output.replace(']', '')
    
    if not windows: output = output.replace('\\n""\\n""\\n\'','')

    Pin, Pout, Q, Obj = output.split(';')
    
    Q = array(Q.split(',')).astype(float)
    Pin = array(Pin.split(',')).astype(float)
    Pout = array(Pout.split(',')).astype(float)
    Obj = array(Obj).astype(float)
   
    results = {}
    results['Q'] = Q
    results['Pin'] = Pin
    results['Pout'] = Pout
    results['obj'] = Obj
    print('%s completed!'%(ID))
    if os.path.exists(optdir + "arbitrage%sarbitrage_ECPP_%s.dzn"%(['/', '\\'][windows],ID)):
        os.remove(optdir + "arbitrage%sarbitrage_ECPP_%s.dzn"%(['/', '\\'][windows],ID))
    else:
        print("The file does not exist")      
    return results


   
def optimise_ECPP_parallel(Year,Region,SH,RTE,Ramp_rate,Pin_max,Pout_max,
                           detailed=False):
      """
      This function calculates the revenue for a storage system with a linear
      discharge ramp rate. Here, ramp up rate=ramp down rate.
      Rampe rate: fraction of full power rating per minute
      """
      c = gd.load_rrp_cal(Year, Region)
      simparams = load_simparams()
      simparams['eta_out'] = float(RTE)/100
      simparams['Pin_max'] = Pin_max
      simparams['Pout_max'] = Pout_max
      simparams['R']=Pout_max*Ramp_rate*60 #MW/hr
      simparams['c']=c.tolist()
      simparams['N']=len(simparams['c'])
      simparams['dt']=0.5
      ID = '%s_%d_%d_%d' %(Region, Year, RTE, int(100*Ramp_rate))
      results = optimise_ECPP(ID,simparams)

      if detailed:
            return({'Region': Region,
                    'Year': Year,
                    'RTE': RTE,
                    'SH':SH,
                    'Ramp_rate': Ramp_rate,
                    'Pin_max':Pin_max,
                    'Pout_max': Pout_max,
                    'Revenue': results['obj'],
                    'Pout':results['Pout'],
                    'Pin':results['Pin'],
                    'Q':results['Q']})
      else:
            return({'Region': Region,
                    'Year': Year,
                    'RTE': RTE,
                    'SH':SH,
                    'Ramp_rate': Ramp_rate,
                    'Pin_max':Pin_max,
                    'Pout_max': Pout_max,
                    'Revenue': results['obj']
                        })
                      
                      
                      
                      
                      
                                          
