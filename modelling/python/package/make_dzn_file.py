#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 14:19:53 2019

@author: jeff
"""

from projdirs import optdir
import sys
#make a dzn file for the minizinc simulation.

win32 = sys.platform == 'win32'


########################################
def make_dzn_file(N, dt, eta_in, eta_out, Q_min_frac, Q_init_frac, Pin_max, Pin_min_frac,
                  Pout_max, Pout_min_frac, SH, loss, c):

    cstring = "["
    
    for i in range(len(c)-1):
        cstring += '%.2f, ' %c[i]
        
    cstring += '%.2f]' %c[-1]
    
    string = """
N = %i;

dt = %.2f;      %%time difference between sample points (also, interval length) (h)

eta_in = %.2f;                      %%efficiency of electricity to heat
eta_out = %.2f;                     %%efficiency of heat to electricity

Q_max = SH * Pout_max / eta_out;    %%maximum stored energy (MWh-th)
Q_min = %.2f * Q_max;               %%minimum stored energy (MWh-th)
Q_init = %.9f * Q_max;              %%initial stored energy (MWh-th)

Pin_max = %.2f;                     %%maximum power to storage (MW-th)
Pin_min = %.2f * Pin_max;           %%minimum power to storage (MW-th)

Pout_max = %.2f;                    %%maximum power to grid (MW-th)
Pout_min = %.2f * Pout_max;         %%minimum power to grid (MW-th)

SH = %.2f;                          %%storage hours (h)
loss = %.2f;                        %% loss percent per day
c = %s
""" %(N, dt, eta_in, eta_out, Q_min_frac, Q_init_frac, Pin_max, Pin_min_frac,
      Pout_max, Pout_min_frac, SH, loss, cstring)

    with open(optdir + "arbitrage%sarbitrage.dzn" %['/', '\\'][win32], "w") as text_file:
        text_file.write(string)
        



#########################################   
def make_dzn_file_ramped(ID, N, dt, R, eta_in, eta_out, Q_min_frac, Q_init_frac, Pin_max, Pin_min_frac,
                  Pout_max, Pout_min_frac, SH, loss, c):

    cstring = "["
    
    for i in range(len(c)-1):
        cstring += '%.3f, ' %c[i]
        
    cstring += '%.3f]' %c[-1]
    
    string = """
N = %i;
dt = %.3f;      %%time difference between sample points (also, interval length) (h)
R = %.3f;        %%rampe rate MW/hr

eta_in = %.3f;                      %%efficiency of electricity to heat
eta_out = %.3f;                     %%efficiency of heat to electricity

Q_max = SH * Pout_max / eta_out;    %%maximum stored energy (MWh-th)
Q_min = %.3f * Q_max;               %%minimum stored energy (MWh-th)
Q_init = %.9f * Q_max;              %%initial stored energy (MWh-th)

Pin_max = %.3f;                     %%maximum power to storage (MW-th)
Pin_min = %.3f * Pin_max;           %%minimum power to storage (MW-th)

Pout_max = %.3f;                    %%maximum power to grid (MW-th)
Pout_min = %.3f * Pout_max;         %%minimum power to grid (MW-th)

SH = %.3f;                          %%storage hours (h)
loss = %.3f;                        %% loss percent per day
c = %s
""" %(N, dt, R, eta_in, eta_out, Q_min_frac, Q_init_frac, Pin_max, Pin_min_frac,
      Pout_max, Pout_min_frac, SH, loss, cstring)

    with open(optdir + "arbitrage%sarbitrage_ramped_%s.dzn" 
              %(['/', '\\'][win32],ID), "w") as text_file: text_file.write(string)

##########################################


def make_dzn_file_parallel(ID, N, dt, eta_in, eta_out, Q_min_frac, Q_init_frac,
                           Pin_max, Pin_min_frac, Pout_max, Pout_min_frac, SH,
                           loss, c):

    cstring = "["
    
    for i in range(len(c)-1):
        cstring += '%.2f, ' %c[i]
        
    cstring += '%.2f]' %c[-1]
    
    string = """
N = %i;

dt = %.2f;      %%time difference between sample points (also, interval length) (h)

eta_in = %.2f;                      %%efficiency of electricity to heat
eta_out = %.2f;                     %%efficiency of heat to electricity

Q_max = SH * Pout_max / eta_out;    %%maximum stored energy (MWh-th)
Q_min = %.2f * Q_max;               %%minimum stored energy (MWh-th)
Q_init = %.9f * Q_max;              %%initial stored energy (MWh-th)

Pin_max = %.2f;                     %%maximum power to storage (MW-th)
Pin_min = %.2f * Pin_max;           %%minimum power to storage (MW-th)

Pout_max = %.2f;                    %%maximum power to grid (MW-th)
Pout_min = %.2f * Pout_max;         %%minimum power to grid (MW-th)

SH = %.2f;                          %%storage hours (h)
loss = %.2f;                        %% loss percent per day
c = %s
""" %(N, dt, eta_in, eta_out, Q_min_frac, Q_init_frac, Pin_max, Pin_min_frac,
      Pout_max, Pout_min_frac, SH, loss, cstring)

    #print("Saving dzn file ...")
    
    with open(optdir + "arbitrage%sarbitrage_%s.dzn"%(['/', '\\'][win32],ID), "w") as text_file:
        text_file.write(string)






def make_dzn_file_PJM_parallel(ID, N, dt, eta_in, eta_out, Q_min_frac,
                               Q_init_frac, Pin_max, Pin_min_frac,
                               Pout_max, Pout_min_frac, Pda_frac, Prt_frac, SH,
                               loss, c_da, c_rt):
      """
      This function creates the dzn file for optimising the storage
      simultaneously for the day ahead an real time markets of the PJM
      """
      
      c_da_string = str(c_da)
      c_rt_string = str(c_rt)
    
      string = """
N = %i;

dt = %.2f;      %%time difference between sample points (also, interval length) (h)

eta_in = %.2f;                      %% efficiency of electricity to heat
eta_out = %.2f;                     %% efficiency of heat to electricity

Q_max = SH * Pout_max / eta_out;    %% maximum stored energy (MWh-th)
Q_min = %.2f * Q_max;               %% minimum stored energy (MWh-th)
Q_init = %.9f * Q_max;              %% initial stored energy (MWh-th)

Pin_max = %.2f;                     %% maximum power to storage (MW-th)
Pin_min = %.2f * Pin_max;           %% minimum power to storage (MW-th)

Pout_max = %.2f;                    %% maximum power to grid (MW-th)
Pout_min = %.2f * Pout_max;         %% minimum power to grid (MW-th)

Pda_frac = %.2f;                    %% maximum power into/from DA market
Prt_frac = %.2f;                    %% minimum power into/from RT market

SH = %.2f;                          %% storage hours (h)
loss = %.2f;                        %% loss percent per day
c_da = %s;                           
c_rt = %s                           
""" %(N, dt, eta_in, eta_out, Q_min_frac, Q_init_frac, Pin_max, Pin_min_frac,
      Pout_max, Pout_min_frac, Pda_frac, Prt_frac, SH, loss, c_da_string, c_rt_string)

    
      with open(optdir + "arbitrage%sarbitrage_PJM_%s.dzn"%(['/', '\\'][win32],
                                                            ID), "w") as text_file:
            text_file.write(string)
            




def make_dzn_file_PJM_capacity(ID, N, dt, eta_in, eta_out, q_max, q_min,
                               q_constraint, Pin_max, Pin_min,
                               SH, loss, Pout_rt_string, c_rt):
      """
      This function creates the dzn file for optimising the storage
      for the capacity market in the PJM
      """
      
      c_rt_string = str(c_rt)
    
      string = """
N = %i;

dt = %.2f;      %%time difference between sample points (also, interval length) (h)

eta_in = %.2f;                      %%efficiency of electricity to heat
eta_out = %.2f;                     %%efficiency of heat to electricity

Q_max = %.2f ;                      %%maximum stored energy (MWh-th)
Q_min = %.2f ;                      %%minimum stored energy (MWh-th)

Pin_max = %.2f;                     %%maximum power to storage (MW-th)
Pin_min = %.2f;           %%minimum power to storage (MW-th)

SH = %.2f;                          %%storage hours (h)
loss = %.2f;                        %% loss percent per day

Q_constraint = %s ;               %% start of day interval number
Pout_rt = %s;                        %%intervals in which emergency occurs
c_rt = %s;                           
""" %(N, dt, eta_in, eta_out, q_max, q_min, 
      Pin_max, Pin_min, SH, loss, q_constraint, Pout_rt_string,
      c_rt_string)

    
      with open(optdir + "arbitrage%sarbitrage_PJM_capacity_%s.dzn"%(['/', '\\'][win32],
                                                            ID), "w") as text_file:
            text_file.write(string)
            
            


def make_dzn_file_capcon(ID, dt, eta_in, eta_out, Q_min_frac,
                         Q_init_frac, Pin_max, Pin_min_frac, Pout_max,
                         Pout_min_frac, SH, loss, Pout_constr, c):
      """
      This function creates the dzn file for optimising the storage
      for the capacity constract market in the NEM
      """
      N = len(c)
      S = len(Pout_constr)
      c_string = str(c)
      Pout_constr_string = str(Pout_constr)
    
      string = """
N = %i;
S = %i;                             %% number of intervals with prices>strike price

dt = %.2f;      %%time difference between sample points (also, interval length) (h)

eta_in = %.2f;                      %%efficiency of electricity to heat
eta_out = %.2f;                     %%efficiency of heat to electricity

Q_init = %.2f ;                     %%maximum initial stored energy (0-1)
Q_max = SH * Pout_max / eta_out;    %%maximum stored energy (MWh-th)
Q_min = %.2f ;                      %%minimum stored energy (MW-th)

Pin_max = %.2f;                     %%maximum power to storage (MW-e)
Pin_min = %.2f;                     %%minimum power to storage (MW-e)

Pout_max = %.2f;                    %%maximum discharge power (MW-e)
Pout_min = %.2f;                    %%minimum discharge (MW-e)

SH = %.2f;                          %%storage hours (h)
loss = %.2f;                        %% loss percent per day

pout_constraint = %s;               %% intervals in which storage must discharge
c = %s;                           
""" %(N, S, dt, eta_in, eta_out, Q_init_frac, Q_min_frac, 
      Pin_max, Pin_min_frac, Pout_max, Pout_min_frac, SH, loss,
      Pout_constr_string, c_string)

    
      with open(optdir + "arbitrage%sarbitrage_NEM_capcon_%s.dzn"%(['/', '\\'][win32],
                                                            ID), "w") as text_file:
            text_file.write(string)




#######################################################
#####  Start of electrifying coal-fired power plants #####
#######################################################
            
            
def make_dzn_file_ECPP(ID, N, dt, R, eta_in, eta_out, Q_min_frac, Q_init_frac,
                       Pin_max, Pin_min_frac, Pout_max, Pout_min_frac,
                       SH, loss, c):

    cstring = "["
    
    for i in range(len(c)-1):
        cstring += '%.3f, ' %c[i]
        
    cstring += '%.3f]' %c[-1]
    
    string = """
N = %i;
dt = %.3f;      %%time difference between sample points (also, interval length) (h)
R = %.3f;        %%rampe rate MW/hr

eta_in = %.3f;                      %%efficiency of electricity to heat
eta_out = %.3f;                     %%efficiency of heat to electricity

Q_max = SH * Pout_max / eta_out;    %%maximum stored energy (MWh-th)
Q_min = %.3f * Q_max;               %%minimum stored energy (MWh-th)
Q_init = %.9f * Q_max;              %%initial stored energy (MWh-th)

Pin_max = %.3f;                     %%maximum power to storage (MW-th)
Pin_min = %.3f * Pin_max;           %%minimum power to storage (MW-th)

Pout_max = %.3f;                    %%maximum power to grid (MW-th)
Pout_min = %.3f * Pout_max;         %%minimum power to grid (MW-th)

SH = %.3f;                          %%storage hours (h)
loss = %.3f;                        %% loss percent per day
c = %s
""" %(N, dt, R, eta_in, eta_out, Q_min_frac, Q_init_frac, Pin_max, Pin_min_frac,
      Pout_max, Pout_min_frac, SH, loss, cstring)

    with open(optdir + "arbitrage%sarbitrage_ECPP_%s.dzn" 
              %(['/', '\\'][win32],ID), "w") as text_file: text_file.write(string)