
# coding: utf-8

# In[ ]:

#!/usr/bin/env python3


# load the FCAS and energy data from the databases   

import sqlite3, os
from projdirs import datadir
from datetime import datetime, timedelta, date, time
import datetime as dt
import pandas as pd
import package.sql_manager as sm
import numpy as np
import sys
import pdb

arbdir = 'arbitrage%s' %['/', '\\'][sys.platform == 'win32']

def load_rrp_cal(Y, region):
    """
    Y: calendar year between 2009 and 2018
    region: state name on the NEM one of ['QLD', 'TAS', 'VIC', 'SA', 'NSW']
    """
    db = 'spot_price.db'
    table = 'actual_price'
    cols = ['state', 'date_time', 'price']
    price_spot = sm.get_data(cols, table, db)
    
    price_spot.columns = ['region', 'date_time', 'spot_price'] # define the data columns 
    price_spot['date_time'] = pd.to_datetime(price_spot['date_time'])
    data = price_spot[(price_spot['date_time'].dt.year == Y) &
                      (price_spot['region'] == region)]
    
    # the following lines were editted to fix the transition to 5 min intervals
    # in the NEM
    data = data.set_index('date_time').sort_index().resample('30min').mean()
    
    return(data.spot_price)


def load_rrp_month(Y, M, region):
    """
    Y: calendar year between 2009 and 2018
    M: calendar month between 1 and 12
    region: state name on the NEM one of ['QLD', 'TAS', 'VIC', 'SA', 'NSW']
    """
    db = 'spot_price.db'
    table = 'actual_price'
    cols = ['state', 'date_time', 'price']
    price_spot = sm.get_data(cols, table, db)
    
    price_spot.columns = ['region', 'date_time', 'spot_price'] # define the data columns 
    price_spot['date_time'] = pd.to_datetime(price_spot['date_time'])
    data = price_spot[(price_spot['date_time'].dt.year == Y) &
                      (price_spot['date_time'].dt.month == M) & 
                      (price_spot['region'] == region)]
    
    return(data.set_index('date_time').sort_idex.spot_price)
    
    
    
def load_rrp_cal_stamped(Y, region):
    """
    Y: calendar year between 2009 and 2018
    region: state name on the NEM one of ['QLD', 'TAS', 'VIC', 'SA', 'NSW']
    """
    # choose the state and date here
    state = "'" + region + "'"
    year = Y

    # connect to spot price database and load the data
    connection = sqlite3.connect(datadir + arbdir + 'spot_price.db')
    cursor = connection.cursor()
    cursor.execute('SELECT * FROM actual_price WHERE state =' + state)
    connection.close
    
    price_spot = pd.DataFrame([row for row in cursor.fetchall()])#.drop(0, axis=1)
    price_spot.columns = ['region', 'date_time', 'spot_price'] # define the data columns 
    price_spot['date_time'] = pd.to_datetime(price_spot['date_time'])
    data = price_spot[price_spot['date_time'].dt.year == year]
    
    return(data)



def load_rrp_fy(FY, State):
   
    """
    Loads the spot price in a region for a financial year.
    FY: financial year, defined as the year that financial year ends at.
    e.g. for financial year 2018-19, FY=2019
    State is one of ['QLD', 'TAS', 'VIC', 'SA', 'NSW']
    """
    DATA = load_rrp_cal(FY-1,State).append(load_rrp_cal(FY,State))
    return(DATA.loc['%d-07-01'%(FY-1):'%d-06-30'%(FY)])
    
   

def load_rrp_forecast(date, region):
    """
    This function loads the forecast predispatch prices
    for a caertain day. This data is published just before
    12:00 midnight for the next trading day.
    'date' must be in the form of pd.to_datetime('2009/03/02')
    'region' should be one of the followings:
    'SA', 'VIC', 'NSW', 'QLD', 'TAS'
    """
    Region = "'" + region + "'"
    connection = sqlite3.connect(datadir + arbdir + 'predispatch.db')
    cursor = connection.cursor()
    cursor.execute('SELECT * FROM forecast_price WHERE state =' + Region)
    price_pre = pd.DataFrame([row for row in cursor.fetchall()]).drop(0, axis=1)
    connection.close
    price_pre.columns = ['region', 'date_time', 'forecast_price']# define the columns of data
    price_pre['date_time'] = pd.to_datetime(price_pre['date_time'])
    price_pre = price_pre.set_index('date_time').sort_index()
    return (price_pre.loc[date: date+timedelta( hours = 23, minutes = 30)])


def load_rrp_24hr_forecast(state, Year, Month, forecast_time):
    """
    This function loads the forecast predispatch prices
    for a caertain year.
    'Year' is an integer 
    'state' should be one of the followings:
    'SA', 'VIC', 'NSW', 'QLD', 'TAS'
    'Month' is an anteger
    'forecast_time': is the time at which the forecast data is generated at. 
    """
    db = '\\Predispatch rolling\\predispatch_rolling_%s.db' %(date(Year, Month,1).strftime("%Y%m")) 
    data = sm.get_data('N/A', 'N/A', db)
    data['generated'] = pd.to_datetime(data.generated)
    data['date_time'] = pd.to_datetime(data.date_time)
    data = data.sort_values(['generated','date_time'])
    return(data[(data.state==state) & (data.generated.dt.time==
                forecast_time)])

def load_rrp_forecast_rolling(state, timestamp):
    """ 'state': should be one of the five states in the NEM.
    'timestamp': should be a panda Timestamp indicating the point in time
    that the predispatch forecast prices have been generated by NEMDE.
    This function loads all the forecast prices for the intervals
    in that NEMDE run.
    """

    year =  timestamp.year
    month = timestamp.month
    db_file = datadir + 'arbitrage\\Predispatch rolling\\predispatch_rolling_%s.db' %(date(year, month,1).strftime("%Y%m")) 
    if not(os.path.exists(db_file)):
        return('nill', 'Database was not found!')
    else:    
        connection = sqlite3.connect(db_file)
        cursor = connection.cursor()
        cursor.execute(" select * from forecast_price_rolling")
        data = pd.DataFrame(cursor.fetchall())
        data.columns=['state','generated','date_time','price']
        data['generated'] = pd.to_datetime(data['generated'])
        data['date_time'] = pd.to_datetime(data['date_time'])
        data = data.sort_values(['generated','date_time'])
        data_selected = data[(data['generated']==timestamp) & (data['state']==state)]
        return (len(data_selected), data_selected)
    
                                 
def load_rrp_forecast_monthly(State, year, month):
    """This function load the forecast prices of a full month in a year
    into a dataframe. """
    db = '\\Predispatch rolling\\predispatch_rolling_%s.db' %(date(year, month,1).strftime("%Y%m")) 
    data = sm.get_data('N/A', 'N/A', db)
    data['generated'] = pd.to_datetime(data['generated'])
    data['date_time'] = pd.to_datetime(data['date_time'])
    data = data[data['state']==State].sort_values(['generated','date_time'])
    return(data)
    



def load_rrp_actual(date, region):
    """
    This function loads the dispatch prices
    for a caertain day.
    'date' must be in the form of pd.to_datetime('2009/03/02')
    'region' should be one of the followings:
    'SA', 'VIC', 'NSW', 'QLD', 'TAS'
    """
    
    price_spot = sm.get_data(cols='N/A', table='N/A',db='spot_price.db')
    price_spot.columns = ['region', 'date_time', 'spot_price'] # define the data columns 
    price_spot['date_time'] = pd.to_datetime(price_spot['date_time'])
    price_spot = price_spot.set_index('date_time').sort_index()
    
    return(price_spot.loc[date: date+timedelta( hours = 23, minutes = 30)])


def load_actual_rrp_at(region, timestamp):
    """
    This function loads the dispatch price
    for a caertain trading interval.
    'timestamp': must be a Timestamp
    'region':  must be one of the followings:
    'SA', 'VIC', 'NSW', 'QLD', 'TAS'
    """
    Region = "'" + region + "'"
    connection = sqlite3.connect(datadir + arbdir + 'spot_price.db')
    cursor = connection.cursor()
    cursor.execute('SELECT * FROM actual_price WHERE state =' + Region)
    price_spot = pd.DataFrame([row for row in cursor.fetchall()]).drop(0, axis=1)
    price_spot.columns = ['region', 'date_time', 'spot_price'] # define the data columns 
    price_spot['date_time'] = pd.to_datetime(price_spot['date_time'])
    price_spot = price_spot.set_index('date_time').sort_index()
    connection.close
    
    return(price_spot.loc[timestamp]['spot_price'])

   

def load_fcas_rrp(Year, region, fcas_type='all'):
    
    """
    this function loads the price of fcas service
    Year: the year ending at.
    region is one of ['QLD', 'TAS', 'VIC', 'SA', 'NSW']
    fcas_type: leave empty to get all fcas prices, or choose one of
    the following: 'RAISE6SEC',	'RAISE60SEC',	'RAISE5MIN',
	'RAISEREG',	'LOWER6SEC',	'LOWER60SEC',	'LOWER5MIN',
	'LOWERREG'
    """
    db = 'FCAS.db'
    FCAS = sm.get_data('N/A', 'N/A', db)
    
    FCAS.SETTLEMENTDATE = pd.to_datetime(FCAS.SETTLEMENTDATE)
    result = FCAS[(FCAS.REGIONID==region)&
                  (FCAS.SETTLEMENTDATE.dt.year==Year)]
    if fcas_type=='all':
        Result = result
    else:
        Result = result[['SETTLEMENTDATE',fcas_type+'RRP']]
    return Result
 
    
    


def load_fcas_disp(FY, region, fcas_market):
    """
    this function loads the enabled fcas service
    FY: financial year, defined as the year that financial year ends at.
    e.g. for financial year 2018-19, FY=2019
    region is one of ['QLD', 'TAS', 'VIC', 'SA', 'NSW']
    """
    state = "'" + region + "'"
    connection = sqlite3.connect(datadir + arbdir + 'FCAS.db')
    cursor = connection.cursor()
    cursor.execute('SELECT * FROM FCAS_dispatch WHERE REGIONID =' + state)
    data = pd.DataFrame([row for row in cursor.fetchall()])
    connection.close
    
    if data.empty:
        print('data is not available!')
    else:
        data.drop(0, axis=1, inplace=True)
        data.columns = ['date_time', 'region', 'total_demand', 'lower_5min_disp', 'lower_60sec_disp', 'lower_6sec_disp',
                        'raise_5min_disp','raise_60sec_disp', 'raise_6sec_disp', 'lower_reg_disp', 'raise_reg_disp']
##        data[list(data.columns)[2:]] = data[data.applymap(np.isreal)][list(data.columns)[2:]].fillna(0)
        data['date_time'] = pd.to_datetime(data['date_time'])
        data.set_index('date_time', inplace=True)
        data.drop_duplicates(keep=False, inplace=True)
    return(data[pd.Timestamp(FY-1,8,1): pd.Timestamp(FY,8,1)][fcas_market+'_disp'])



def load_lmp_cal(Y, lmp_type='da'):
    """
    This function loads the day ahead or reat ime LMP data for a specific year
    for all the zones across the PJM.
    Y: calendar year between 2010 and 2019
    """
    db = 'spot_price_PJM_%s.db'%(lmp_type)
    table = sm.list_tables(db)[0]
    cols = [col[0] for col in sm.list_columns(db,table)]
    price_spot = sm.get_data(cols, table, db)
    price_spot['datetime_beginning_ept'] = pd.to_datetime(price_spot['datetime_beginning_ept'])
    price_spot = price_spot.sort_values('datetime_beginning_ept')
    data = price_spot[(price_spot['datetime_beginning_ept'].dt.year == Y)]
    return(data)
    
def load_lmp_zone(Zone, lmp_type='da'):
    """
    This function loads the real time or day ahead LMP data for a specific zone
    in the PJM for the period of 2010 until 2019.
    lmp_type: 'da' or 'rt'
    """
    db = 'spot_price_PJM_%s.db'%(lmp_type)
    table = sm.list_tables(db)[0]
    cols = [col[0] for col in sm.list_columns(db,table)]
    price_spot = sm.get_data(cols, table, db)
    price_spot['datetime_beginning_ept'] = pd.to_datetime(price_spot['datetime_beginning_ept'])
    price_spot = price_spot.sort_values('datetime_beginning_ept')
    data = price_spot[(price_spot['pnode_name'] == Zone)].\
          sort_values('datetime_beginning_ept').\
          set_index('datetime_beginning_ept')['total_lmp_%s'%(lmp_type)]
    
    return(data)
    
def load_lmp_zone_year(Zone, Year, lmp_type='da', Year_type='cal'):
    """
    This function loads the day ahead or real time LMP data for a specific zone
    in the PJM for a specific year.
    lmp_type: 'da' or 'rt'
    Year_type = 'cal': calendar year or 'fy': financial year
    Convention: FY ends in fy
    """
    data = load_lmp_zone(Zone, lmp_type)
    if Year_type=='cal':
          data = data[data.index.year==Year]
    elif Year_type=='fy':
          data = data.loc['%d-10-01'%(Year-1):'%d-09-30'%(Year)]
                
    return(data)



def get_capa_price(Year,Zone):
    """
    This functions loads the cpacity price $/MW for 
    a Zone in PJM for a specific year
    """
    db = 'capacity_price_PJM.db'
    data = sm.get_data('N/A','N/A', db)
    return(data.loc[data['zone']==Zone, 'FY%ito%i'%(Year,Year-2000+1)].values[0])
    
    

def get_NEM_cap_price(Region, Year, Averaging='annually'):
      """
      This function loads all the cap prices for a NEM region in a year
      """
      data = sm.get_data('N/A', 'N/A','cap_prices_NEM.db')
      data.year=data.year.astype(int)
      Data = data[(data.averaging==Averaging)&
                  (data.region==Region)&
                  (data.year==Year)]
      return(Data) 



def convert_curr(Amount, Year, exchange='AUD-to-USD'):
      """
      This function converts between AUD and USD
      Year: the year of exchange
      """
      exc_rate = pd.DataFrame({'Year':np.arange(2010,2023,1),
                         'AUD/USD':[0.92, 1.03, 1.04, 0.97, 0.9,
                                    0.75, 0.74, 0.77, 0.75, 0.7,
                                    0.69, 0.75, 0.69]
                         }).set_index('Year')
      if exchange== 'AUD-to-USD':
            return(exc_rate.loc[Year].values[0]*Amount)
      else:
            return(Amount/exc_rate.loc[Year])
            
            
def time_value(Amount, Year, Year_ref, infl_rate):
      value = Amount*(1+infl_rate)**(Year_ref-Year)
      return(value)
            
            
            
            

