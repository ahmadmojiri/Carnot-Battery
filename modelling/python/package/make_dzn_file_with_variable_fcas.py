#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tuesday Oct 8 11:07:53 2019

@author: Ahmad
"""

from projdirs import optdir

#make a dzn file for the minizinc optimisation for when the thermal battery
#is locat

def make_dzn_file_with_variable_fcas(N, dt, eta_in, eta_out, self_disch,
                                     Q_min_frac, Q_init_frac, Pin_max,
                                     Pin_min_frac,Pout_max, Pout_min_frac,
                                    SH, rrp, fcas_rrp, fcas_disp):
    
    string = """
N = %i;

dt = %.2f;      %%time difference between sample points (also, interval length) (h)

eta_in = %.2f;                      %%efficiency of electricity to heat
eta_out = %.2f;                     %%efficiency of heat to electricity
self_discharge = %0.2f;             %%self_discharge losses per day

Q_max = SH * Pout_max / eta_out;    %%maximum stored energy (MWh-th)
Q_min = %.2f * Q_max;               %%minimum stored energy (MWh-th)
Q_init = %.2f * Q_max;              %%initial stored energy (MWh-th)

Pin_max = %.2f;                     %%maximum power to storage (MW-th)
Pin_min = %.2f * Pin_max;           %%minimum power to storage (MW-th)

Pout_max = %.2f;                    %%maximum power to grid (MW-th)
Pout_min = %.2f * Pout_max;         %%minimum power to grid (MW-th)


SH = %.2f;                          %%storage hours (h)

rrp = %s;
fcas_rrp = %s;
fcas_disp = %s;
""" %(N, dt, eta_in, eta_out, self_disch, Q_min_frac, Q_init_frac,
      Pin_max, Pin_min_frac, Pout_max, Pout_min_frac, SH, str(rrp),
      str(fcas_rrp), str(fcas_disp))

    with open(optdir + "arbitrage\\arbitrage_plus_variable_fcas.dzn", "w") as text_file:
        text_file.write(string)

