#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tuesday Oct 8 11:07:53 2019

@author: Ahmad
"""

from projdirs import optdir
import sys
import numpy as np

win = sys.platform == 'win32'
connector = ['/','\\'][win]

#make a dzn file for the minizinc optimisation for when the thermal battery
#is collocated with a PV plant and accepts any curtailed energy from it.

def make_dzn_file_PV_plant_SH_PB(N, dt, eta_rh, eta_pc, eta_inv, eta_ts, 
                                 Q_min_frac, Q_init_frac, 
                                 P_inv_in_max, c_pb, c_ts, c_rh, c, P_in_panels, 
                                 Q_ts_max = None, P_ts_in_max = None, 
                                 P_ts_out_max=None):

    P_in_panels = np.array2string(np.array(P_in_panels), 
                            threshold = len(P_in_panels) + 1,
                            separator = ',', 
                            formatter = {'float_kind' : lambda x: "%.3f" %x})
    
    if Q_ts_max != None: 
        Q_ts_max_str = "Q_ts_max = %.8f;" %Q_ts_max
    else: Q_ts_max_str = ''
    
    if P_ts_in_max != None:
        P_ts_in_max_str = "P_ts_in_max = %.8f;" %P_ts_in_max
    else: P_ts_in_max_str = ''
    
    if P_ts_out_max != None: 
        P_ts_out_max_str = "P_ts_out_max = %.8f;" %P_ts_out_max
    else: P_ts_out_max_str = ''
    
    string = """
N = %i;

dt = %.2f;      %%time difference between sample points (also, interval length) (h)

eta_rh = %.2f;                      %%efficiency of electricity to heat
eta_pc = %.2f;                     %%efficiency of heat to electricity
eta_inv = %.2f;                     %%design inverter efficiency from PVWatts
eta_ts = %.2f;                      %%hourly thermal energy storage efficiency

%s
Q_ts_min = %.2f * Q_ts_max;               %%minimum stored energy (MWh-th)
Q_ts_init = %.2f * Q_ts_max;              %%initial stored energy (MWh-th)

%s
%s
P_inv_in_max = %.2f;                %%minimum power to storage (MW-e)
                         
c_pb = %i;
c_ts = %i;
c_rh = %i;
c = %s;                              %%Electricity spot price
P_in_panels = %s;                    %%DC flow from the PV plant

""" %(N, dt, eta_rh, eta_pc, eta_inv, eta_ts, Q_ts_max_str, Q_min_frac, 
    Q_init_frac, P_ts_in_max_str, P_ts_out_max_str, P_inv_in_max, c_pb, c_ts, c_rh,
    str(c), P_in_panels)

    with open(optdir + "arbitrage%sarbitrage_PV_plant_SH_PB.dzn" %connector, "w") as text_file:
        text_file.write(string)