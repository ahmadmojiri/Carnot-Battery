
# coding: utf-8

# In[ ]:

#!/usr/bin/env python3

"""
Created on Mon Jun 24 15:51:21 2019

@author: jeff
"""

import pandas as pd
import numpy as np
from projdirs import datadir

def is_leap_year(year):
    """
Extra days occur in years which are multiples of four with the exception of 
centennial years not divisible by 400"""
    return int(year % 4 == 0) - int(year % 100 == 0) + int(year % 400 == 0)

def get_yearly_data(state, year, startidx = 0 , endidx = 17520, alldata = False):

    filename = datadir + "NEM\\" + state + '\\' + year + '\\' + year + "_RRP.csv"
    
    c = np.loadtxt(filename, delimiter = ',')
    
    if alldata == True:
        c = c[0:len(c)]
    else:
        c = c[startidx:endidx]  
    
    return c 

def get_time_indices(year, start_month, start_day, start_hour, 
                   end_month, end_day, end_hour):
    """
INPUTS
year: integer (year integer must be in full, eg. 2010, not 10)
start_month, start_day, start_hour: integers 
end_month, end_day, end_hour: integers

Note: possible hours for the day are defined in the range [0,23]

Note: when the end index that this function returns is used in the 'get_yearly_data' function, 
the timestep corresponding to end_hour will not be included in the time series returned.
This is because the spot price for the 2nd-last timestep applies up until the end_hour

RETURNS
start_idx: integer - the start index for the year
end_idx: integer - the end index for the year
"""
    leapyear = is_leap_year(year)

    montharray = [31,28+leapyear,31,30,31,30,31,31,30,31,30,31]

    start_idx = 48 * (int(np.sum(montharray[:start_month-1])) + start_day - 1) + 2 * start_hour
    end_idx = 48 * (int(np.sum(montharray[:end_month-1])) + end_day - 1) + 2 * end_hour + 2  

    return start_idx, end_idx

def get_time_series(state, start_year, start_month, start_day, start_hour, 
                        end_year, end_month, end_day, end_hour):

    """
INPUTS 
state: string ['NSW', 'TAS', ... etc]
start_year: integer - This is the year expressed in YYYY format
start_month: The number of the start calendar month
start_day: The number of the start day within the calendar month
start_hour: The number of the start hour within the day in the range [0-23]
end variables follow the same convention

Note: leap years are accounted for

OUTPUT
c: numpy array - This is the time series to be used in the arbitrage optimisation
"""

    c = np.array([])
    
    years = np.arange(start_year, end_year + 1)
    
    for i in range(len(years)):
        
        if i != 0:
            start_month_now, start_day_now, start_hour_now = 1, 1, 0
        else: 
            start_month_now, start_day_now, start_hour_now = start_month, start_day, start_hour
        
        if i != (len(years) - 1): 
            end_month_now, end_day_now, end_hour_now = 12, 31, 23
        else: 
            end_month_now, end_day_now, end_hour_now = end_month, end_day, end_hour
        
        start_idx, end_idx = get_time_indices(years[i], start_month_now, start_day_now, start_hour_now,
                                   end_month_now, end_day_now, end_hour_now)
        
        print(start_idx, end_idx)
        
        c = np.append(c, get_yearly_data(state, str(years[i]), start_idx, end_idx), axis = 0)
        
    return c