# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 14:48:30 2016

@author: jecu01
"""

import pylab as pl

import numpy as np

epsilon = 1e-8

def hour_angle(hour):
    return 15./180 * pl.pi * (hour - 12)
    
def declination(day):
    return 23.45 / 180 * pl.pi * pl.sin(2 * pl.pi/365 * (284 + day))
    
def zenith_angle(lat, ha, dec):
    return pl.arccos(pl.cos(lat) * pl.cos(dec) * pl.cos(ha) + pl.sin(lat) * pl.sin(dec))

def alt_angle(zen):
    return pl.pi/2 - zen
    
def azi_angle(zen, lat, dec, ha):
    ha[ha == 0] += epsilon
    return np.sign(ha) * abs(np.arccos(np.round((np.cos(zen) * np.sin(lat) - \
                                                 np.sin(dec))/(np.sin(zen) * \
                                                 np.cos(lat)), 10)))

def calc_azi_alt(lat, day, hour):
    ha = hour_angle(hour)
    dec = declination(day)
    
    zen = zenith_angle(lat, ha, dec)
    alt = alt_angle(zen)
    azi = azi_angle(zen, lat, dec, ha)
    
    return azi, alt

