# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 14:40:27 2020

@author: Ahmad Mojiri
"""


import os
from projdirs import datadir
from datetime import datetime, timedelta, date, time
import datetime as dt
import pandas as pd
import package.sql_manager as sm
import numpy as np
import sys
import package.get_NEM_data as gd

arbdir = 'arbitrage%s' %['/', '\\'][sys.platform == 'win32']

def count_price_spikes(Year, State, Strike_price):
      """
      Counts the number of reading intervals in each day with
      spot price greated that the strike price
      """
      RRP = gd.load_rrp_cal(Year,State).reset_index()
      days = pd.date_range('%d-01-01'%(Year),'%d-12-31'%(Year))
      daily_spikes = pd.DataFrame(columns=['date', 'spikes#'])
      for Day in days:
          daily_rrp = RRP[(RRP['date_time'].dt.date==Day.date())]
          Spikes_count = (daily_rrp['spot_price']>=Strike_price).sum()
          new_row = pd.DataFrame(
                      {'date': [Day.date()],
                       'spikes#': [Spikes_count]})
          daily_spikes = daily_spikes.append(new_row,ignore_index=True)
      return(daily_spikes)
      