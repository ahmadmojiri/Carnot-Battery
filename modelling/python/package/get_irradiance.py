#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 21:17:38 2020

@author: jeff
"""

#from solar_model import calc_angles
import numpy as np
import pandas as pd
from projdirs import datadir
import matplotlib.pyplot as mpl
from datetime import datetime
# from .sim_classes import plantlist

def load_solar_grid(ts, datatype = 'DNI'):
    
    year = str(ts.year)

    ts_UTC = ts.tz_convert('UTC')
    
    sourcedata = '/home/jeff/local/data/gridded/'
        
    grid = np.loadtxt(sourcedata + datatype + "/" + year + \
    "/solar_%s_%s%02d%02d_%02dUT.txt" %(datatype.lower(), ts_UTC.year, 
                                        ts_UTC.month, ts_UTC.day, ts_UTC.hour), \
                                        skiprows = 6)
     
    print(ts_UTC.year, ts_UTC.month, ts_UTC.day, ts_UTC.hour)
    
    return grid


def calc_grid_idx(lat, long):
    "Note: latitude and longitude are in degrees."
    
    N_rows =  679
    cellsize = 0.0499954
    x_ll = 112.025005
    y_ll = -43.9750009

    x_idx = int(np.round((long - x_ll) / cellsize))
    y_idx = N_rows - int(np.round((lat - y_ll) / cellsize))
 
    return x_idx, y_idx


def plot_grid(grid, ax = None, vmin = 0, vmax = 1000):
    
    if ax == None: 
        fig = mpl.figure()
        ax = fig.add_subplot(111)
    ax.imshow(grid, cmap = 'inferno', vmin = vmin, vmax = vmax)
    return ax

def plot_grid_loc(xidx, yidx, ax):
    ax.plot(xidx, yidx, 'x', color = 'lime', ms = 12, mew = 2)

def get_grid_irradiance(grid, xidx, yidx):
    
    irradiance = grid[yidx, xidx]
        
    return irradiance


def calc_aperture_sun_cosine(tilt, azi_coll, alt, azi_sun):
    """
This function calculates the cosine angle between the solar vector and
the surface of the PV panel. See:

http://powerfromthesun.net/Book/chapter04/chapter04.html

Equation 4.3.
    
ARGUMENTS
NAME        VAR       TYPE        FORMAT

tilt        beta      float      radians
azi_coll    gamma     float      radians
alt         alpha     float      radians
azi_sun     A         float      radians
"""

    return np.sin(alt) * np.cos(tilt) + \
           np.cos(alt) * np.sin(tilt) * np.cos(azi_coll - azi_sun)

def extract_gridded_irrd(time, coordlist = [[-35.2834600, 149.1280700]],
                         datatype = 'DNI'):
    
    #asks for the list of coordinates because once you've got the gridded data
    #loaded you may as well extract as many points from it as possible.
    
    grid = load_solar_grid(time, datatype)
    
    irrd_list = []
    
    for lat, lon in coordlist:
        irrd_list += [get_grid_irradiance(grid, *calc_grid_idx(lat, lon))]
     
    return irrd_list


def extract_gridded_irrd_ts(ts, coordlist = [[-35.2834600, 149.1280700]], 
                            datatype = 'DNI', plantnames = None):
    
    if plantnames: 
        columns = plantnames    
    else: columns = range(len(coordlist))
    
    ts_hourly = pd.date_range(ts[0].floor(freq = 'H'), 
                              ts[-1].ceil(freq = 'H'), 
                              freq = 'H')
    
    irrd = pd.DataFrame(index = ts_hourly, columns = columns, dtype = float)
    
    for time in range(len(irrd.index)):
        irrd.iloc[time,:] = extract_gridded_irrd(irrd.index[time], 
                                                 coordlist, datatype)
        print(irrd.index[time])
        
    return irrd


def extract_gd_and_save(ts, plantlist, datatype = 'DNI', debug = True):
    """
Extracts the gridded data and saves it as pandas Series.
When 'debug' is True, the DataFrame is save to a temporary folder.  
"""   

    numplants = len(plantlist)
    filename = datadir + "irradiance/grid_derived/" + ['','temp/'][debug] + "%s_2015.pkl" %(datatype)

    plantnames = [i for i in map(getattr, plantlist, ['name'] * numplants)]
    latitudes = [i for i in map(getattr, plantlist, ['latitude'] * numplants)]
    longitudes = [i for i in map(getattr, plantlist, ['longitude'] * numplants)]
    
    coordlist = np.array([latitudes, longitudes]).T
    
    irrd = extract_gridded_irrd_ts(ts, coordlist, datatype, plantnames)
    
    irrd.to_pickle(filename)


    
def load_grid_irradiance_data(datatype='DNI', debug = True):
    
    filename = datadir + "irradiance/grid_derived/" + ['','temp/'][debug] + "%s_2015.pkl" %(datatype)
    irrd = pd.read_pickle(filename)
    
    return irrd

def load_minute_irradiance_data(dr, siteid = '072150', datatype = 'DNI'):
    
    dr_month = pd.date_range(dr[0], dr[-1], freq = 'MS')
    
    ts = pd.Series()
    
    dtype_idx = {'GHI': 7,
                 'DNI':12,
                 'DHI':17}
    
    for i, month in enumerate(dr_month):
        filename = datadir + "irradiance/minute/" + "sl_%s_%s_%02d.txt" %(siteid, month.year, month.month)
        data = np.loadtxt(filename, delimiter = ',', skiprows = 1, dtype = object)
        dates = data[:,2:7].astype(int)
        indices = [pd.Timestamp(datetime(*i)) for i in dates.astype(int)]
        irrd = data[:,dtype_idx[datatype]]
        irrd[irrd == '       '] = np.nan
        irrd = irrd.astype(float)
        series = pd.Series(data = irrd, index = indices)
        ts = ts.append(series)
    
    ts.index = ts.index.tz_localize('Australia/Sydney')

    ts = ts.loc[dr[0]:dr[-1]]
    
    return(ts)

def resample_irradiance_data(ts, irrd):
    
    samples = pd.DataFrame(index = ts, columns = irrd.columns)
    irrd = irrd.assign(interp = 'no', dtype = str)
    samples = samples.assign(interp = 'yes', dtype = str)
    
    irrd_new = irrd.append(samples)
    irrd_new = irrd_new.sort_index()
    
    irrd_new.iloc[:,:-1] = irrd_new.iloc[:,:-1].interpolate('time')
    irrd_new = irrd_new[irrd_new['interp'] == 'yes']

    return irrd_new
    

def make_annual_dr(year = '2015', freq = 'H', tz = 'Australia/Sydney'):
    return pd.date_range(start = year, end = str(int(year)+1), 
                         freq = 'H', tz = tz)[:-1]

def make_daily_dr(year = '2015', day = 1, freq = 'H', tz = 'Australia/Sydney'):
    starttime = pd.Timestamp('2015') + pd.Timedelta(days = day - 1)
    endtime = pd.Timestamp('2015') + pd.Timedelta(days = day)
    return pd.date_range(start = starttime, end = endtime, freq = freq, tz = tz)[:-1]
