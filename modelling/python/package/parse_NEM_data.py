#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 13:55:45 2019

@author: jeff
"""

#load a year of data from the AEMO CSV files

import numpy as np
from projdirs import datadir

def make_yearly_data_file(state, year):

    curdir = datadir + "NEM\\" + state + '\\' + year + "\\" 
    
    RRP = []
    
    for i in range(1,13):
    
        RRP += list(np.loadtxt(curdir + "PRICE_AND_DEMAND_%s%02d_%s1.csv" %(year,i,state), 
                               delimiter = ',', skiprows = 1, usecols = [3]))
        
    RRP = np.array(RRP)
    
    np.savetxt(curdir + year + "_RRP.csv", RRP)
    

if __name__ == "__main__":

    state = "NSW"
    year = "2018"
