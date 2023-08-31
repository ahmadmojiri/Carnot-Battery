#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:53:54 2019

@author: jeff
"""

#generating the fake data to test arbitrage.mzn

import numpy as np
import matplotlib.pyplot as mpl

def make_sin_price(numdays, dt, c_change, c_avg, period):
    
    numhrs = numdays * 24
    
    timesteps = int(numhrs / dt)
    
    t = np.linspace(0, numhrs-dt, timesteps)
    
    c = -np.sin(2 * np.pi / period * t) * c_change + c_avg
    
    fig = mpl.figure()
    
    ax = fig.add_subplot(111)
    ax.plot(t,c)
    
    c_str = '['
    
    for i in range(len(c)-1):
        c_str += '%.2f, ' %c[i]
        
    c_str += '%.2f]' %c[-1]
    
    return c_str


if __name__ == "__main__":
    
    numdays = 7                 #days
    dt = 0.5                    #hrs
    c_change = 20               #$ (plus-minus )
    c_avg = 40                  #$
    period = 24                 #hrs
    
    make_sin_price(numdays, dt, c_change, c_avg, period)