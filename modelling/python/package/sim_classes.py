#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 23:03:39 2020

@author: jeff

This module defines classes for simulations of PV plants.

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as mpl
import requests
from projdirs import datadir

#Use the basic solar model that is defined in the solar_model module
# from .solar_model import calc_azi_alt
#==OR==
#use pysolar
# from pysolar.solar import get_position

pd.plotting.register_matplotlib_converters()

# from .get_irradiance import load_grid_irradiance_data

# class Sun():
#     def __init__(self, lat, lon, ts_local, elev = 0):
#         self.lat = lat
#         self.lat_rad = lat * np.pi / 180.
#         self.lon = lon
#         self.ts_local = ts_local
#         self.ts_UTC = ts_local.tz_convert('UTC')        
#         self.elev = elev       
#         self.calc_sunpos()
#         self.load_irradiance()
     
                       
#     def calc_sunpos(self, mode = 1):
        
#         if mode == 0:
#             #This model does not have any information about the equation of time.
#             self.azis, self.alts = calc_azi_alt(self.lat_rad, 
#                                                 np.array(self.ts_local.dayofyear), 
#                                                 np.array(self.ts_local.hour))
            
 
#         if mode == 1:
#             print("Calculating solar angles ...")
#             self.azis, self.alts = np.array([[i,j] for i, j in map(get_position, 
#                                         [self.lat] * len(self.ts_local), 
#                                         [self.lon] * len(self.ts_local), 
#                                         self.ts_local.to_pydatetime(), 
#                                         [self.elev] * len(self.ts_local))]).T
        
            
#     def load_irradiance(self):
#         """Currently this function assumes that it is an annual simulation. 
#         Maybe in the future it will be worth defining the simulation to work
#         for an arbitrary time series."""
        
#         ts_hours = pd.date_range(self.ts[0], self.ts[-1], freq = 'H')
        
#         self.DNI = load_grid_irradiance_data(year = self.ts_local.year[0], 
#                                         datatype = 'DNI')
#         # self.DNI = np.roll(self.DNI, 1)
        
#         self.GHI = load_grid_irradiance_data(year = self.ts_local.year[0],
#                                         datatype = 'GHI')
#         # self.GHI = np.roll(self.GHI, 1)
        
#         self.DNI[np.where(self.DNI < 0)] = 0
#         self.GHI[np.where(self.GHI < 0)] = 0
        
#         #create a array with a default value of zero for the 
#         #horizontal surface solar cosine efficiency
#         #self.eta_cos_horiz = np.zeros(len(self.alts))
        
#         #create an array of inidices where the solar altitude is greater than 0
#         #sunup_idx = np.where(self.alts > 0)
        
#         #calculate horizontal cosine factor at times when the sun is up
#         #self.eta_cos_horiz[sunup_idx] = np.sin(self.alts[sunup_idx] * np.pi/180)
#         self.eta_cos_horiz = np.sin(self.alts * np.pi/180)
        
#         #calculate the diffuse solar irradiance
#         self.DSI = self.GHI - self.DNI * self.eta_cos_horiz


#     def plot_series(self, plottype = 'DNI', plotstr = '', ax = None):
        
#         if ax == None:
#             fig = mpl.figure()
#             ax = fig.add_subplot(111)
        
#         if plotstr != '':
#             plotstr = "'%s' ," %plotstr
        
#         exec("ax.plot(self.ts_local, self.%s, %slabel = '%s')" \
#              %(plottype, plotstr, plottype))
        

#     def plot_irrd_all(self, ax = None):
        
#         if ax == None:
#             fig = mpl.figure()
#             ax = fig.add_subplot(111)
        
#         self.plot_series('DNI', ax=ax)
#         self.plot_series('GHI', ax=ax)
#         self.plot_series('DSI', ax=ax)
                
#         fig.legend()

       
#     def plot_sun_angles(self, ax = None):
        
#         if ax == None:
#             fig = mpl.figure()
#             ax = fig.add_subplot(111)
            
#         self.plot_series('alts', ax=ax)
#         self.plot_series('azis', ax=ax)
        
#     def plot_compare(self, type1 = 'DNI', type2 = 'alts', ax = None):

#         if ax == None:
# #See http://web.archive.org/web/20120618121009/http://notes.brooks.nu/2008/03/plotting-on-left-and-right-axis-simulateously-using-matplotlib-and-numpy
#             fig = mpl.figure()
#             axL = fig.add_subplot(111)
#         else: 
#             axL = ax
#             fig = ax.figure
            
#         axR = fig.add_subplot(111, sharex=axL, frameon = False)
#         axR.yaxis.tick_right()
#         axR.yaxis.set_label_position("right")
        
#         self.plot_series(type1, ax=axL)
#         self.plot_series(type2, ax=axR, plotstr = 'r')
        
#         fig.legend()
        
        
# class PlantPerf(): 
#     def __init__(self, plant, sun):
#         self.plant = plant
#         self.sun = sun
#         self.calc_eta_cos()
#         self.calc_panel_irrd()

#     def calc_panel_irrd(self):

#         #calculate the direct panel irradiance
#         self.DPI = self.sun.DNI * self.eta_cos

#         #calculate the total panel irradiance
#         self.TPI = self.sun.DSI + self.DPI

#     def calc_eta_cos(self):
#         """
#     This function calculates the cosine angle between the solar vector and
#     the surface of the PV panel. See:
    
#     http://powerfromthesun.net/Book/chapter04/chapter04.html Equation 4.3.
        
#     ARGUMENTS
#     NAME        VAR       TYPE        FORMAT
    
#     tilt        beta      float      radians
#     azi_coll    gamma     float      radians
#     alt         alpha     float      radians
#     azi_sun     A         float      radians
#     """
    
#         self.eta_cos = np.sin(self.sun.alts) * np.cos(self.plant.tilt) + \
#                np.cos(self.sun.alts) * np.sin(self.plant.tilt) * \
#                    np.cos(self.plant.azimuth - self.sun.azis)





#2020-01-18 Panel orientations for NSP are known, and the PSF azimuth is known,
#otherwise the values have been assumed. Still need to work out which ones are
#east-west facing.

import sqlite3

connection = sqlite3.connect(datadir + 'plants/' + 'plants.db')

plantcols = np.array(['region', 'latitude', 'longitude', 'power_AC', 'elev', 
                 'tilt', 'azimuth', 'oversize', 'design', 'losses', 'eta_inv', 
                 'C_pv_cap', 'C_pv_OM'])

plantrows = np.array(['PSF', 'GrSF', 'DSH', 'GuSF', 'NSP', 'MSF', 'BHP'])

datatypes = np.array(['TEXT', 'REAL', 'REAL','REAL', 'REAL', 'REAL', 'REAL', 
                      'REAL', 'REAL', 'REAL', 'REAL', 'REAL', 'REAL'])

typedict = dict(zip(plantcols, datatypes))

plantlib = np.array([
['PSF', 'NSW', -33.1129309, 148.0742085, 55, 0, np.pi/6, 0, 1.2, 1., 14, 0.96, 1111e3, 20e3],
['GrSF', 'NSW', -34.3174099, 146.1185502, 30, 0, np.pi/6, 0, 1.2, 1., 14, 0.96, 1111e3, 20e3],
['DSH', 'NSW', -32.2728558, 148.6476476, 25, 0, np.pi/6, 0, 1.2, 1., 14, 0.96, 1111e3, 20e3],
['GuSF', 'NSW', -34.6094078, 149.4730104, 10, 0, np.pi/6, 0, 1.2, 1., 14, 0.96, 1111e3, 20e3],
['NSP', 'NSW', -31.5660058, 147.0642427, 102, 0, 25*np.pi/180, 0, 1.2, 1., 14, 0.96, 1111e3, 20e3],
['MSF', 'NSW', -29.5697345, 149.8646149, 56, 0, np.pi/6, 0, 1.2, 1., 14, 0.96, 1111e3, 20e3],
['BHP', 'NSW', -31.9882521, 141.3896637, 53, 0, np.pi/6, 0, 1.2, 1., 14, 0.96, 1111e3, 20e3]])

df = pd.DataFrame(plantlib[:,1:], index = plantlib[:,0], columns = plantcols)

df.to_sql('PV', con = connection, if_exists = 'replace', dtype = typedict)


#REFERENCES
#PSF: https://arena.gov.au/projects/parkes-solar-farm/
#GrSF: https://arena.gov.au/projects/griffith-solar-farm/
#DSH: https://arena.gov.au/projects/dubbo-solar-hub/
#GuSF: https://arena.gov.au/projects/gullen-solar-farm/*
#NSP: https://www.agl.com.au/about-agl/how-we-source-energy/nyngan-solar-plant*
#MSF: https://arena.gov.au/projects/moree-solar-farm/
#BHP: https://www.agl.com.au/about-agl/how-we-source-energy/broken-hill-solar-plant