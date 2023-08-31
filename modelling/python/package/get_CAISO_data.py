
import sqlite3, os
from projdirs import datadir
from datetime import datetime, timedelta, date, time
import datetime as dt
import pandas as pd
import package.sql_manager as sm
import numpy as np
import sys

arbdir = 'arbitrage%s' %['/', '\\'][sys.platform == 'win32']

def load_lmp_cal(Y, hub='TH_NP15'):
    """
    Y: calendar year between 2013 and 2020
    hub: hub name in CAISO; one of ['TH_NP15', 'TH_SP15', 'TH_ZP26'];
    """
    db = 'CAISO_LMP.db'
    table = 'actual_price'
    cols = ['date_time', 'price', 'hub']
    LMP = sm.get_data(cols, table, db)
    
    
    LMP['date_time'] = pd.to_datetime(LMP['date_time'])
    data = LMP[(LMP['date_time'].dt.year == Y) &
                      (LMP['hub'] == hub)]
    
    
    
    return(data.set_index('date_time').sort_index().price)



















            

            
            
            
            

