#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 21:06:01 2020

@author: jeff
"""

from subprocess import check_output
from projdirs import optdir, datadir
from . import sql_manager as sm
from numpy import array, ceil
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import requests
from scipy.optimize import minimize, minimize_scalar
from copy import copy

pd.plotting.register_matplotlib_converters()

win = sys.platform == 'win32'
connector = ['/','\\'][win]

class Simulation():
    """This is the calss that brings all the variables together and creates
    the minzinc file for optimisation, returning results that can be plotted
    and investigated."""
    
    ##Things to do:
    # * Allow for all the time series to be resampled at interval 'dt', with 
    #   a specific interpolation technique.
    
    def __init__(self, Time, Financials, Market, Storage, PowerCycle, 
                 Plant = None, Inverter = None, opt_os = 0, opt_rh = 0,
                 opt_st = 0, opt_pc = 0, region = None, 
                 write = True, verbose = True, purchase = True):
        self.Time = Time
        self.Financials = Financials
        self.Market = Market
        self.Plant = Plant
        self.Storage = Storage
        self.Inverter = Inverter
        self.PowerCycle = PowerCycle
        self.region = region
        self.optimised = False
        self.write = write
        self.verbose = verbose
        self.purchase = purchase
        if not opt_os + opt_rh + opt_st + opt_pc: self.optimisation = False
        self._set_optimisation_bools(opt_os, opt_rh, opt_st, opt_pc)
        
        opt_bools = np.array([opt_os, opt_rh, opt_st, opt_pc])
        self.opt_int = opt_bools.dot(2**np.arange(opt_bools.size)[::-1])        
        
        self.i = 0

        self._calc_SH()
        self._calc_CH()
        self._calc_SPH()
        self._calc_CPH()

        if Plant:
            if self.region == None:
                self.region = Plant.region
            self.P_inv_in_max = Plant.power_AC / Inverter.eta_inv
            
        if self.region == None: self.region = 'NSW'

        self.c = Market.load_rrp(Time.start_time, Time.end_time, self.region)
        self.c.index = self.c.index + pd.Timedelta(hours = self.Time.dt/2)
        
        self.Storage.set_eta_st(self.PowerCycle.eta_pc)
        self._resample_time_series()
        self._calc_cost_factors()
               
        #from NREL annual technology baseline spreadsheet https://atb.nrel.gov/
        #$1111/kW installed plus $20/(kW-yr) O&M cost
    
    def __repr__(self):
        
        repr_str = ''
        
        if self.Plant:
            repr_str += "Plant oversizing factor = %.2f" %self.Plant.oversize + '\n'
            
        if self.Storage.E_st_max > 0 and self.PowerCycle.P_pc_out_max > 0:
            repr_str += "Storage hours = %.1f" %(self.Storage.E_st_max / \
                  self.PowerCycle.P_pc_out_max) + '\n'
        else: repr_str += "Storage hours = 0 \n" 
        repr_str += "Heater power in = %.1f" %self.Storage.P_st_in_max + '\n'
        repr_str += "Power cycle size = %.1f" %self.PowerCycle.P_pc_out_max + '\n'
        repr_str += "Revenue = %i" %self.revenue + '\n' 
        if self.Plant:
            repr_str += "PV plant cost = %.2f" %-self.C_pv_time + '\n'
        repr_str += "Thermal heaters' cost = %.2f" %-self.C_rh_time + '\n'
        repr_str += "Thermal energy storage cost = %.2f" %-self.C_st_time + '\n'
        repr_str += "Power cycle cost = %.2f" %-self.C_pc_time + '\n'
        repr_str += 'profit = %i\n' %self.profit
        repr_str += 'objective function string = %s\n' %self.obj_string
        repr_str += 'Internal rate of return = %.4f' %self.IRR
           
        return(repr_str)

    def set_acc(self, acc):
        
        self.Financials.acc = acc
        self._calc_cost_factors()
     
    # def _write_data(self)
    def _calc_SH(self):
        self.SH = self.Storage.E_st_max / self.PowerCycle.P_pc_out_max
    
    def _calc_CH(self):
        self.CH = self.Storage.E_st_max / (self.PowerCycle.eta_pc * self.Storage.P_st_in_max)
        
    def _calc_SPH(self):
        self.SPH = self.PowerCycle.P_pc_out_max / self.Storage.E_st_max
    
    def _calc_CPH(self):
        self.CPH = self.PowerCycle.eta_pc * self.Storage.P_st_in_max / self.Storage.E_st_max 
    
    def _calc_cost_factors(self):
        
        self.c_rh = self.Financials.calc_annual_cost(self.Storage.C_rh_cap) * \
            self.Financials.USD2AUD * self.Time.yearfrac
        self.c_st = (self.Financials.calc_annual_cost(self.Storage.C_st_cap) + \
                    self.Storage.C_st_OM) * self.Financials.USD2AUD * \
                    self.Time.yearfrac
        self.c_pc = (self.Financials.calc_annual_cost(self.PowerCycle.C_pc_cap) + \
                    self.PowerCycle.C_pc_OM) * self.Financials.USD2AUD * \
                    self.Time.yearfrac
        if self.Plant:
            self.c_pv = (self.Financials.calc_annual_cost(self.Plant.C_pv_cap) + \
            self.Plant.C_pv_OM) * self.Financials.USD2AUD * \
            self.Time.yearfrac

    def _calc_total_costs(self):       
        self.C_rh_cap_total = self.Storage.C_rh_cap_total * self.Financials.USD2AUD
        self.C_st_cap_total = self.Storage.C_st_cap_total * self.Financials.USD2AUD
        self.C_pc_cap_total = self.PowerCycle.C_pc_cap_total * self.Financials.USD2AUD
        if self.Plant: 
            self.C_pv_cap_total = self.Plant.C_pv_cap_total * \
                self.Financials.USD2AUD

        self.C_st_OM_total = self.Storage.C_st_OM_total * self.Financials.USD2AUD
        self.C_pc_OM_total = self.PowerCycle.C_pc_OM_total * self.Financials.USD2AUD
        if self.Plant:
            self.C_pv_OM_total = self.Plant.C_pv_OM_total * self.Financials.USD2AUD
            
        #calculate capital cost and OM of total system
        self.CAPEX = self.C_pc_cap_total + self.C_st_cap_total + \
            self.C_rh_cap_total
        self.OM = self.C_pc_OM_total + self.C_st_OM_total
        if self.Plant:
            self.CAPEX += self.C_pv_cap_total
            self.OM += self.C_pv_OM_total 

    def _set_optimisation_bools(self, opt_os, opt_rh, opt_st, opt_pc):
        
        self.opt_os = opt_os
        self.opt_rh = opt_rh
        self.opt_st = opt_st
        self.opt_pc = opt_pc

    def _resample_time_series(self):
        ##This might need testing to see whether the resulting DC time series has
        #resampling the DC time series from the PV plant to be every half hour.
        #This function might be able to be made more general in case there
        #are other time series that need to redone
        ##WARNING: 
        ##
        ##Duration no more than one year total, though arbitrary start date OK. 
        ##
        if self.Plant:
            DC = self.Plant.DC
            
            period_range = pd.period_range(start = self.Time.start_time, 
                                           end = self.Time.end_time, 
                                           freq = 'H')[:-1]
            
            has_leap_year = len(period_range[(period_range.day == 29) & \
                                             (period_range.month == 2)]) > 0
            
            if has_leap_year:
                #just add another 31st Dec - I know it's dodgy, but won't matter 
                DC = np.append(DC, DC[-24:])
            
            period_range_ref = pd.period_range(start = self.Time.start_time.year,
                                               end = self.Time.start_time.year + 1, 
                                               freq = 'H')[:-1]
            
            start_idx = np.where(period_range_ref == period_range[0])[0]
            
            DC = np.roll(DC, -start_idx)
            
            DC = DC[:len(period_range)]
            
            self.P_in_plant = pd.Series(DC, index = period_range)
            self.P_in_plant = self.P_in_plant.resample('%iT' %int(self.Time.dt*60)).pad()
            self.P_in_plant.index = self.Time.int_index
            
        
        if not 'period' in str(self.c.index.dtype):
            dt_c = pd.Timedelta(self.c.index[1] - self.c.index[0]).seconds / 60
            self.c.index = self.c.index.to_period('%iT' %dt_c)
        self.c = self.c.resample('%iT' %int(self.Time.dt*60)).pad()
        self.c.index = self.c.index.to_timestamp()

        
    def _make_mzn_file(self):         
        
        # st_optstrs = ['float: c_st;', 'var ', ' - c_st * E_st_max / eta_st', 
        #                     ' ++ [";"] ++  [show(E_st_max)]']
        # if not self.opt_st: st_optstrs = [''] * len(st_optstrs)
            
        
        # rh_optstrs = ['float: c_rh;', 'var ', ' - c_rh * P_st_in_max', 
        #                     ' ++ [";"] ++ [show(P_st_in_max)]']
        # if not self.opt_rh: rh_optstrs = [''] * len(rh_optstrs)
        
        # pc_optstrs = ['float: c_pc;', 'var ', ' - c_pc * P_pc_out_max', 
        #                   ' ++ [";"] ++ [show(P_pc_out_max)]']
        # if not self.opt_pc: pc_optstrs = [''] * len(pc_optstrs)
        
        
        #Delete these and uncomment code above to include sizing in minizinc 
        #optimisation
        st_optstrs = [''] * 4
        rh_optstrs = [''] * 4
        pc_optstrs = [''] * 4
        
        plant_strs = ['float: eta_inv;',
                      'float: P_inv_in_max;  %%maximum power to the inverter (MW-e)]',
                      'array[tint] of float: P_in_plant;  %%DC electricity generated by the PV plant',
                      'array[tint] of var float: P_out_pass;  %%power lost due to full tank',
                      'array[tint] of var float: P_inv_in;     %%power into the inverter (MW-e)',
                      ' - 1e5 * P_out_pass[i]',
                      'constraint forall(i in tint)(P_out_pass[i] >= 0);',
                      'constraint forall(i in tint)(P_inv_in[i] >= 0);',
                      ' + P_inv_in_max',
                      'constraint forall(i in tint)(P_out_pass[i] <= P_inv_in_max);',
                      'constraint forall(i in tint)(P_inv_in[i] <= P_inv_in_max);',
                      'constraint forall(i in tint)(P_inv_in[i] <= P_in_plant[i]);',
                      'P_out_pass[i] + P_inv_in[i] - P_in_plant[i] + ',
                      ' + P_inv_in[i] * eta_inv',
                      ' ++ [";"] ++ [show(P_out_pass)] ++ [";"] ++ [show(P_inv_in)]']
        
        purch_str = ['0','P_st_in_max'][self.purchase]
        
        if not self.Plant: plant_strs = [''] * len(plant_strs)
        
        self._mzn_string = """
int: N;

set of int: tsamp = 1..N+1;    %%number of sample points
set of int: tint = 1..N;       %%number of intervals

float: dt;                %%time difference between sample points (also, interval length) (h)
float: eta_rh;            %%efficiency of resistance heater in tank
float: eta_st;            %%efficiency of storage (for thermal system this is the power cycle efficiency)
%s
%s
%s
%s
%s
%sfloat: E_st_max;          %%maximum stored energy (MWh-th)
%sfloat: P_st_in_max;       %%maximum power to storage (MW-e)
%sfloat: P_pc_out_max;

var float: E_st_min;          %%minimum stored energy (MWh-th)
var float: E_st_init;         %%initial stored energy (MWh-th)

array[tint] of float: c;                 %%energy spot price
%s

%%flows into and out of the system (plus P_in_plant)

array[tint] of var float: P_in_grid;    %%power purchased from the electricity grid (MW-e)
array[tint] of var float: P_out_grid;   %%total power generated (MW-e)
%s

%%flows within the system
array[tint] of var float: P_st_in;      %%power into the thermal energy storage (MW-e)
array[tint] of var float: P_pc_out;     %%power out of the thermal energy storage tank (MW-e)
%s
array[tsamp] of var float: E;        %%Energy in storage (MWh-th)

var float: obj;

obj = sum(i in tint)(c[i] * (P_out_grid[i] - P_in_grid[i]) * dt%s)%s%s%s;

%%=== CONSTRAINTS ===
constraint E[1] = E_st_init;
constraint E[N] = E[1];

%%"minimum" constraints
constraint forall(i in tint)(P_in_grid[i] >= 0);
constraint forall(i in tint)(P_out_grid[i] >= 0);
%s
constraint forall(i in tint)(P_st_in[i] >= 0);
constraint forall(i in tint)(P_pc_out[i] >= 0);
%s
constraint forall(i in tsamp)(E[i] >= E_st_min);
%%constraint P_st_in_max >= 0;
%%constraint P_pc_out_max >= 0;
%%constraint E_st_max >= 0;

%%"maximum" constraints
constraint forall(i in tint)(P_in_grid[i] <= %s);
constraint forall(i in tint)(P_out_grid[i] <= P_pc_out_max%s);
%s
constraint forall(i in tint)(P_st_in[i] <= P_st_in_max);
constraint forall(i in tint)(P_pc_out[i] <= P_pc_out_max);
%s
%s
constraint forall(i in tsamp)(E[i] <= E_st_max);
constraint P_st_in_max < 1e4;
constraint P_pc_out_max < 1e4;
constraint E_st_max < 1e5;

%%energy balances
constraint forall(i in tint)(P_in_grid[i] = %sP_st_in[i]);
constraint forall(i in tint)(P_out_grid[i] = P_pc_out[i]%s);
constraint forall(i in tint)(E[i+1] = E[i] + (P_st_in[i] * eta_rh * eta_st - P_pc_out[i]) * dt);    

solve maximize obj;

output [show(obj)] ++ [";"] ++ [show(P_in_grid)] ++ [";"] ++ [show(P_out_grid)] ++ [";"] ++ [show(P_pc_out)] ++ [";"] ++ [show(E)]%s%s%s%s;
        """%(plant_strs[0], st_optstrs[0], rh_optstrs[0], pc_optstrs[0], plant_strs[1], st_optstrs[1],
         rh_optstrs[1], pc_optstrs[1], plant_strs[2], plant_strs[3], plant_strs[4],
         plant_strs[5], st_optstrs[2], rh_optstrs[2], pc_optstrs[2], plant_strs[6], plant_strs[7],
         purch_str, plant_strs[8], plant_strs[9], plant_strs[10], plant_strs[11], plant_strs[12], 
         plant_strs[13], plant_strs[14], st_optstrs[3], rh_optstrs[3], pc_optstrs[3])

        self.obj_string = "obj = sum(i in tint)(c[i] * (P_out_grid[i] - \
P_in_grid[i]) * dt)%s%s%s;" %(st_optstrs[2], rh_optstrs[2], pc_optstrs[2])
    
        with open(optdir + "arbitrage%s" %connector[win] + \
                  'general_optimisation.mzn', 'w') as text_file:
            text_file.write(self._mzn_string)

    def _make_dzn_file(self):
        if self.Plant: 
            P_in_plant_str = 'P_in_plant = ' + np.array2string(np.array(self.P_in_plant), 
                                    threshold = len(self.P_in_plant) + 1,
                                    separator = ',', 
                                    formatter = {'float_kind' : lambda x: "%.3f" %x})
            eta_inv_str = 'eta_inv = %.2f;' %self.Inverter.eta_inv
            P_inv_in_max_str = 'P_inv_in_max = %.2f;' %self.P_inv_in_max ##move to inverter?
        else: 
            P_in_plant_str = ''
            eta_inv_str = ''
            P_inv_in_max_str = ''
        
        # if not self.opt_st: 
        E_st_max_str = "E_st_max = %.8f;" %self.Storage.E_st_max
        c_st_str = ''
        # else: 
            # E_st_max_str = ''
            # c_st_str = 'c_st = %.2f;' %self.c_st
        
        # if not self.opt_rh: 
        P_st_in_max_str = "P_st_in_max = %.8f;" %self.Storage.P_st_in_max
        c_rh_str = ''
        # else: 
            # P_st_in_max_str = ''
            # c_rh_str = 'c_rh = %.2f;' %self.c_rh
        
        # if not self.opt_pc: 
        P_pc_out_max_str = "P_pc_out_max = %.8f;" %self.PowerCycle.P_pc_out_max
        c_pc_str = ''
        # else: 
            # P_pc_out_max_str = ''
            # c_pc_str = 'c_pc = %.2f;' %self.c_pc
        
                  
        self._dzn_string = """
N = %i;

dt = %.2f;      %%time difference between sample points (also, interval length) (h)

eta_rh = %.2f;                      %%efficiency of electricity to heat
eta_st = %.2f;                     %%efficiency of heat to electricity
%s

%s
E_st_min = 0;
E_st_init = %.2f * E_st_max;              %%initial stored energy (MWh-th)

%s
%s
%s               
                         
%s
%s
%s


c = %s;                              %%Electricity spot price

%s                                  %%DC flow from the PV plant

            """ %(self.Time.N, self.Time.dt, self.Storage.eta_rh, 
                  self.Storage.eta_st, eta_inv_str, E_st_max_str, 
                  self.Storage.E_st_init_frac, P_st_in_max_str, P_pc_out_max_str, 
                  P_inv_in_max_str, c_pc_str, c_st_str, c_rh_str, 
                  str(list(self.c)), P_in_plant_str)
            
        with open(optdir + "arbitrage%s" %connector[win] + \
                  'general_optimisation.dzn', 'w') as text_file:
            text_file.write(self._dzn_string)
    
    def _run_minizinc(self):
        
        mzdir = ['/home/jeff/local/software/MiniZinc/bin/',
                 'C:\\Program Files\\minizinc\\'][win]
        
        curdir = optdir + "arbitrage%s" %connector[win]
        
        output = str(check_output([mzdir + 'minizinc', "--soln-sep", '""', 
                                   "--search-complete-msg", '""', "--solver", 
                                   "COIN-BC", curdir + "general_optimisation.mzn", 
                                   curdir + "general_optimisation.dzn"]))

        self.output = copy(output)

        output = output.replace("""\\r\\n""\\r\\n""\\r\\n\'""", '')
        output = output.replace("b\'", '')
        output = output.replace('[', '')
        output = output.replace(']', '')
        
        if not win: output = output.replace('\\n""\\n""\\n\'','')
        
        self._minizinc_output = output
        

    def _process_mzn_output_str(self):
        
        base_vars_str = ', '.join(['self.Results.obj', 'self.Results.P_in_grid', 'self.Results.P_out_grid', 'self.Results.P_pc_out', 'self.Results.E']) + ', '
        pv_vars_str = ['', ', '.join(['self.Results.P_out_pass', 'self.Results.P_inv_in']) + ', '][bool(self.Plant)]
        #opt_st_str = 'self.Storage.E_st_max, '
        #opt_rh_str = 'self.Storage.P_st_in_max, '
        #opt_pc_str = 'self.PowerCycle.P_pc_out_max, '
        
        vars_str = ', '.join([base_vars_str, pv_vars_str])
        #, opt_st_str, opt_rh_str, opt_pc_str])
        
        vars_str = vars_str.replace(' ,', '')
        
        exec_str = vars_str + "= self._minizinc_output.split(';')"
        exec_str = exec_str.replace(', =', ' =')
        
        exec(exec_str)

        self.Results.P_in_grid = pd.Series(self.Results.P_in_grid.split(','), 
                                    index = self.Time.int_index, dtype = float)

        self.Results.P_out_grid = pd.Series(self.Results.P_out_grid.split(','),
                                    index = self.Time.int_index, dtype = float)

        self.Results.P_pc_out = pd.Series(self.Results.P_pc_out.split(','), 
                                    index = self.Time.int_index, dtype = float)

        self.Results.E = pd.Series(self.Results.E.split(','),
                                   index = self.Time.samp_index, dtype = float)
        
        if self.Plant:
            
            self.Results.P_out_pass = pd.Series(
                self.Results.P_out_pass.split(','), 
                index = self.Time.int_index, dtype = float)
            
            self.Results.P_inv_in = pd.Series(
                self.Results.P_inv_in.split(','),
                index = self.Time.int_index, dtype = float)
            
            self.Results.P_inv_out = self.Results.P_inv_in * \
                                     self.Inverter.eta_inv
            
            self.Results.P_st_in = self.P_in_plant + self.Results.P_in_grid - \
                self.Results.P_out_pass - self.Results.P_inv_in
            
            self.Results.P_st_in_plant = self.Results.P_st_in - \
                self.Results.P_in_grid
            
            self.Results.P_st_in_plant_cur = \
                pd.Series(np.zeros(len(self.Time.int_index)), 
                                   index = self.Time.int_index)
            
            t_curtail = np.round(self.Results.P_inv_in, 2) == \
                        np.round(self.P_inv_in_max, 2)
            
            self.Results.P_st_in_plant_cur[t_curtail] = \
                self.Results.P_st_in_plant[t_curtail]
            
            self.Results.P_st_in_plant_vol = self.Results.P_st_in_plant - \
                self.Results.P_st_in_plant_cur
                   
            self.Results.P_st_out = self.Results.P_out_grid - \
                                    self.Results.P_inv_out 


        else: 
            self.Results.P_st_in = self.Results.P_in_grid
            self.Results.P_st_out = self.Results.P_out_grid
    
        if self.opt_st: 
            self.Storage.set_E_st_max(float(self.Storage.E_st_max))
            self._calc_SH()
        if self.opt_pc: self.PowerCycle.set_P_pc_out_max(float(self.PowerCycle.P_pc_out_max))
        if self.opt_rh: self.Storage.set_P_st_in_max(float(self.Storage.P_st_in_max))
        

    def _sim_details_str(self):
        
        sim_intro_str = "=== SIMULATION DETAILS ===\n\n"
        if self.Plant:
            plant_name_str = "Plant name: " + self.Plant.name + '\n'
            plant_region_str = "State: " + self.Plant.region + '\n'
            plant_os_str = "Oversizing factor: " + str(self.Plant.oversize) + '\n'
        else: 
            plant_name_str = ''
            plant_region_str = self.region + '\n'
            plant_os_str = ''
            
        eta_st_str = "Storage efficiency: %.2f\n" %self.Storage.eta_st
        start_time_string = "Simulation start time: %s\n" %str(self.Time.start_time)
        end_time_string = "Simulation end time: %s\n\n\n" %str(self.Time.end_time)
        
        opt_intro_str = "=== OPTIMISATION DETAILS ===\n\n"
        if self.Plant: 
            opt_os_str = "Optimising plant oversizing: %s\n" %['no', 'yes'][self.opt_os]
        else:
            opt_os_str = ""
        opt_rh_str = "Optimising tank charge rate: %s\n" %['no', 'yes'][self.opt_rh]
        opt_st_str = "Optimising tank size: %s\n" %['no', 'yes'][self.opt_st]
        opt_pc_str = "Optimising power cycle size: %s\n" %['no','yes'][self.opt_pc]
        
        details_str = sim_intro_str + plant_name_str + plant_region_str + \
            plant_os_str + eta_st_str + start_time_string + end_time_string + \
            opt_intro_str + opt_os_str + opt_rh_str + opt_st_str + opt_pc_str 
            
        if self.verbose: print(details_str)

    def _calc_costs_profit(self):
        #The suffix 'time' denotes that this is the fraction of the total
        #captial cost that is paid over the time interval of the simulation
        if self.Plant:
            self.C_pv_time = self.c_pv * self.Plant.power_AC * \
                self.Plant.oversize
        self.C_st_time = self.c_st * self.Storage.E_st_max / self.Storage.eta_st
        self.C_rh_time = self.c_rh * self.Storage.P_st_in_max
        self.C_pc_time = self.c_pc * self.PowerCycle.P_pc_out_max
    
        self.revenue = sum(self.c * (self.Results.P_out_grid - \
                                 self.Results.P_in_grid)) * self.Time.dt
    
        self.annual_revenue = self.revenue / self.Time.yearfrac 
    
        self.profit = self.revenue - self.C_st_time - self.C_rh_time - \
            self.C_pc_time

        if self.Plant: 
            self.profit -= self.C_pv_time
            
    def _calc_IRR(self):
           
        self.IRR_list = [-self.CAPEX] + \
            [self.annual_revenue - self.OM] * \
            self.Financials.lifetime
                         
        self.IRR = np.irr(self.IRR_list)
   
    def _optimise_base(self, write = True, assign = True, 
                      db='storage_value_PV_perfect', table='storage_value'):

        self._sim_details_str()
        self._make_mzn_file()
        self._make_dzn_file()
        self._run_minizinc()
        self.Results = Results()
        self._process_mzn_output_str()
        self._calc_costs_profit()
        self._calc_total_costs()
        self._calc_IRR()
        
        if self.verbose: print(self)

    def _optimise_os(self, os):
        
        self.Plant.set_oversize(os)
        self._resample_time_series()
        self._optimise_base()
                
        return -self.IRR  
    
    def _optimise_general(self, varlist):

        print("ITERATION NUMBER: %i" %self.i)
        print("OPTIMISING FOR STORAGE SIZE %.2f MWh-e" %self.Storage.E_st_max)

        eps = 1e-4
        
        index = 0
        
        if 'os' in self.x0_dict.keys():
            varlist[0] = abs(varlist[0])
            if varlist[0] < 0.01: varlist[0] = 0.01
            if (abs(varlist[0] / self.Plant.oversize - 1)) > eps:
                self.Plant.set_oversize(varlist[0])
            index += 1

        if 'P_st_in_max' in self.x0_dict.keys():
            varlist[index] = abs(varlist[index])
            self.Storage.set_P_st_in_max(varlist[index])
            index += 1

        if 'E_st_max' in self.x0_dict.keys():
            varlist[index] = abs(varlist[index])
            self.Storage.set_E_st_max(varlist[index])
            self._calc_SH()
            index += 1

        if 'P_pc_out_max' in self.x0_dict.keys():
            varlist[index] = abs(varlist[index])
            self.PowerCycle.set_P_pc_out_max(varlist[index])
        
        if self.verbose: print(varlist)
        
        self._resample_time_series()
        self._optimise_base()
        
        self.i += 1
        
        return -self.IRR  
                            
    
    def optimise(self):
        
        if self.opt_int != 0:
            self.x0_dict = {}
            
            # dict({'P_st_in_max':np.nan, 'P_pc_out_max':np.nan, 
            #                 'E_st_max':np.nan, 'os':np.nan})
        
            bounds_dict = {}
            
            # dict({'P_st_in_max':np.nan, 'P_pc_out_max':np.nan, 
            #                     'E_st_max':np.nan, 'os':np.nan})
            
        
            if self.opt_os == 1:
                self.x0_dict['os'] = self.Plant.oversize
                bounds_dict['os'] = [1e-4,3]
            if self.opt_rh == 1:
                self.x0_dict['P_st_in_max'] = self.Storage.P_st_in_max
                bounds_dict['P_st_in_max'] = [0,1000]
            if self.opt_st == 1:
                self.x0_dict['E_st_max'] = self.Storage.E_st_max
                bounds_dict['E_st_max'] = [0,10000]
            if self.opt_pc == 1:
                self.x0_dict['P_pc_out_max'] = self.PowerCycle.P_pc_out_max
                bounds_dict['P_pc_out_max'] = [0,1000]
        

            x0 = list(self.x0_dict.values())
            bounds = list(bounds_dict.values())


            minimize(self._optimise_general, x0 = x0, bounds = bounds, 
                     method = 'Nelder-Mead', options = {'xatol':1})

    
        # elif self.opt_os and bool(self.Plant):
        #     print("Optimising plant oversizing.")
        #     self.Plant.set_oversize(minimize_scalar(self._optimise_os, 
        #                                             bounds = (0.5,3), 
        #                                             method = 'bounded', 
        #                                             tol = 0.001)['x'])
        
        else: 
            self._optimise_base()
            self._calc_total_costs()
            self._calc_IRR()
        
        self.optimised = True
    
    def sensitivity(self):
        "Performs a sensitivity analysis across the four sizing variables"
            
        opt_bools = [copy(self.opt_os), copy(self.opt_rh), 
                     copy(self.opt_st), copy(self.opt_pc)]
        
        size_params = []
        
        size_params += [copy(self.Storage.E_st_max),
                        copy(self.Storage.P_st_in_max),
                        copy(self.PowerCycle.P_pc_out_max)]
        
        if self.Plant: size_params += [copy(self.Plant.oversize)]
        
        self._set_optimisation_bools(0, 0, 0, 0)

        #In case the user to have defined their own sensitivity range beforehand
        if not getattr(self, 'Sens', False):
            self.Sens = Sensitivity()
            
        E_st_max = self.Sens.range * self.Storage.E_st_max
        P_st_in_max = self.Sens.range * self.Storage.P_st_in_max
        P_pc_out_max = self.Sens.range * self.PowerCycle.P_pc_out_max
    
        self.Sens.E_st_max_profit = []
        self.Sens.E_st_max_IRR = []
        self.Sens.P_st_in_max_profit = []
        self.Sens.P_st_in_max_IRR = []
        self.Sens.P_pc_out_max_profit = []
        self.Sens.P_pc_out_max_IRR = []
        
        for stor_size in E_st_max:
            self.Storage.set_E_st_max(stor_size)
            self._optimise_base()
            self.Sens.E_st_max_profit.append(self.profit)
            self.Sens.E_st_max_IRR.append(self.IRR)
        self.Storage.set_E_st_max(size_params[0])
        
        for heater_size in P_st_in_max:
            self.Storage.set_P_st_in_max(heater_size)
            self._optimise_base()
            self.Sens.P_st_in_max_profit.append(self.profit)
            self.Sens.P_st_in_max_IRR.append(self.IRR)
        self.Storage.set_P_st_in_max(size_params[1])
        
        for powercycle_size in P_pc_out_max:
            self.PowerCycle.set_P_pc_out_max(powercycle_size)
            self._optimise_base()
            self.Sens.P_pc_out_max_profit.append(self.profit)
            self.Sens.P_pc_out_max_IRR.append(self.IRR)            
        self.PowerCycle.set_P_pc_out_max(size_params[2])
        
        if self.Plant:
            os = self.Sens.range * self.Plant.oversize
            self.Sens.os_profit = []
            self.Sens.os_IRR = []

            for oversize in os:
                self._optimise_os(oversize)
                self.Sens.os_profit.append(self.profit)
                self.Sens.os_IRR.append(self.IRR)
            self.Plant.set_oversize(size_params[3])
            self._resample_time_series()

        self._optimise_base()
        self._set_optimisation_bools(*opt_bools)
    

    def plot_results(self):
        
        fig1 = plt.figure()
        ax11 = fig1.add_subplot(311)
        ax12 = fig1.add_subplot(312)
        ax13 = fig1.add_subplot(313)
        
        #just to be clear that none of the
        # grid purchased electricity should be going to the inverter, and only to storage
        
        # bar1 = ax1.bar(rrp.index, P_st_out, simparams['dt']/24, -P_st_out)[0]
        # bar2 = ax1.bar(rrp.index, P_st_in_plant, simparams['dt']/24, 0)[0]
        # bar3 = ax1.bar(rrp.index, P_st_in_grid, simparams['dt']/24, P_st_in_plant)[0]
        
        # fig1.legend([bar1, bar2, bar3], 
        #            ['P_st_out', 'P_st_in_plant', 'P_st_in_grid'], loc = 9, ncol = 5)
        
        labellist = ['P_st_out', 'P_in_grid']
        
        line1 = ax11.plot(-self.Results.P_st_out, label = labellist[0])[0]
        line2 = ax11.plot(self.Results.P_in_grid, label = labellist[1])[0]
        linelist = [line1, line2]
        

        if self.Plant:
            labellist += ['P_out_pass', 'P_st_in_plant_cur', 'P_st_in_plant_vol']
            line3 = ax11.plot(self.Results.P_out_pass, label = labellist[2])[0] 
            line4 = ax11.plot(self.Results.P_st_in_plant_cur, label = labellist[3])[0]
            line5 = ax11.plot(self.Results.P_st_in_plant_vol, label = labellist[4])[0]
            linelist += [line3, line4, line5]
        
        fig1.legend(linelist, labellist, loc = 9, ncol = len(linelist))
        
        ax11.grid(True)
        ax11.set_xlabel('')
        ax11.set_ylabel('Power (MW-e)')
        ax11.axis([self.Time.int_index[0], 
                  self.Time.int_index[-1],
                  -1.1 * self.PowerCycle.P_pc_out_max, 
                   1.1 * self.Storage.P_st_in_max])
        
        ax12.plot(self.Results.E)
        ax12.set_ylim([0, 1.1*max(self.Results.E)])
        ax12.grid(True)
        ax12.set_xlabel('')
        ax12.set_ylabel('E (MWh-e)')
        ax12.axis([self.Time.int_index[0], 
                  self.Time.int_index[-1],
                  -0.1 * self.Storage.E_st_max, 
                  1.1 * self.Storage.E_st_max])
                  
        ax13.plot(self.c)
        ax13.axis([self.Time.samp_index[0], 
                  self.Time.samp_index[-1],
                  0, 1.1*max(self.c)])
        ax13.grid(True)
        ax13.set_xlabel('Date')
        ax13.set_ylabel('Spot Price\n(\$/MWh)')
        
        
        fig2 = plt.figure()
        ax21 = fig2.add_subplot(211)
        ax21.step(self.Time.int_index, self.Results.P_in_grid, label = 'purchased')
        if self.Plant:
            ax21.step(self.Time.int_index, self.Results.P_st_in_plant_cur, label = 'curtailed')
            ax21.step(self.Time.int_index, self.Results.P_st_in_plant_vol, label = 'voluntary')
            ax21.step(self.Time.int_index, self.Results.P_out_pass, label = 'bypass')
        ax21.legend()
        
        if self.Plant:
            ax22 = fig2.add_subplot(212)
            ax22.step(self.Time.int_index, self.P_in_plant, label = 'plant in')
            ax22.step(self.Time.int_index, self.Results.P_inv_in, label = 'inverter in')
            ax22.step(self.Time.int_index, self.Results.P_st_in_plant_cur, label = 'curtailed')
            ax22.step(self.Time.int_index, self.Results.P_st_in_plant_vol, label = 'voluntary')
            ax22.legend()
        
        plt.show()


class Sensitivity:
    def __init__(self, vmin = 0.5, vmax = 1.5, interval = 0.1):
        self.range = np.arange(vmin, vmax, interval)
    
    def plot(self):
        
        fig = plt.figure()
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        if getattr(self, 'os_profit', False):
            ax1.plot(self.range, self.os_profit, label = 'oversize')
            ax2.plot(self.range, self.os_IRR, label = 'oversize')
            
        ax1.plot(self.range, self.E_st_max_profit, label = 'storage size')
        ax2.plot(self.range, self.E_st_max_IRR, label = 'storage size')
        ax1.plot(self.range, self.P_st_in_max_profit, label = 'charging capacity')
        ax2.plot(self.range, self.P_st_in_max_IRR, label = 'charging capacity')
        ax1.plot(self.range, self.P_pc_out_max_profit, label = 'generation capacity')
        ax2.plot(self.range, self.P_pc_out_max_IRR, label = 'generation capacity')
        
        plt.legend()
        
        

class Results():
    pass

class Inverter():
    def __init__(self, eta_inv):
        self.eta_inv = eta_inv

class Market():
    def __init__(self, market = 'rrp'):
        "Market can be straight arbitrage, fcas etc."
        ##still need to set this up so that it gets the prices for a given
        #market. 'market' variable does nothing at the moment. 
        self.market = market
    
    def load_rrp(self, start_time, end_time, region):
        """
        start_time: time stamp for the start of the interval being simulated
        end_time: time stamp for the end of the interval being simulated
        region: state name on the NEM one of ['QLD', 'TAS', 'VIC', 'SA', 'NSW']
        """
        start_stamp = pd.Timestamp(start_time)
        end_stamp = pd.Timestamp(end_time)
        
        if region in ['QLD', 'TAS', 'VIC', 'SA', 'NSW']: 
            db = 'spot_price.db'
            cols = ['state', 'date_time', 'price']
        elif region in ['SYS', 'SE1', 'SE2', 'SE3', 'SE4', 'FI', 'DK1', 'DK2',\
                      'Oslo', 'Kr.sand', 'Bergen', 'Molde', 'Tr.heim', \
                      'TromsÃ¸', 'EE', 'ELE', 'LV', 'LT']: 
            db = 'spot_price_np.db'
            cols = ['region', 'date_time', 'spot_price']
        else: 
            print('region not found')
            return
        
        table = 'actual_price'
        
        price_spot = sm.get_data(cols, table, db)
        
        price_spot.columns = ['region', 'date_time', 'spot_price'] # define the data columns 
        price_spot['date_time'] = pd.to_datetime(price_spot['date_time'])
        
        data = price_spot[(price_spot['region'] == region) & 
                          (price_spot['date_time'] >= start_stamp) & 
                          (price_spot['date_time'] < end_stamp)]
               
        return data.set_index('date_time')['spot_price'].astype(float).interpolate()

class PV_plant:
    def __init__(self, name):
        
        db_dir = datadir + 'plants%s' %connector
        db = 'plants.db'
        #load dataframe with plant params
        var_df = dict(sm.get_table(db, 'PV', db_dir).loc[name])   
        
        self.name = name            #short name of plant
        self._load_variables(**var_df)
        self._calc_power_DC()
        self.C_pv_cap_total = self.C_pv_cap * self.power_AC * self.oversize
        self.C_pv_OM_total = self.C_pv_OM * self.power_AC * self.oversize
        
    def _load_variables(self, region, latitude, longitude, power_AC, elev, tilt,
                       azimuth, oversize, design, losses, eta_inv, C_pv_cap,
                       C_pv_OM):

        self.region = region        #Australian state in which the plant is located
        self.latitude = latitude    #[degrees] site coordinate
        self.longitude = longitude  #[degrees] site coordinate
        self.elev = elev            #[m] elevation of sites
        self.power_AC = power_AC    #[MW] design output AC
        self.tilt = tilt            #[rad] panel tilt
        self.azimuth = azimuth      #[rad] panel azimuth angle  
        self.oversize = oversize    #[kW/m^2] (derived from values for PSF, assumed for others)
        self.design = design        #[kW/m^2] (assumed)
        self.losses = losses
        self.eta_inv = eta_inv      #[-] The efficiency of the inverter
        self.C_pv_cap = C_pv_cap
        self.C_pv_OM = C_pv_OM   

    def set_oversize(self, oversize):
        print('Setting plant oversize %.2f' %oversize)
        self.oversize = oversize
        self._calc_power_DC()
        self.C_pv_cap_total = self.C_pv_cap * self.power_AC * self.oversize
        self.C_pv_OM_total = self.C_pv_OM * self.power_AC * self.oversize

    def set_C_pv_OM(self, C_pv_OM):
        self.C_pv_OM = C_pv_OM
        self.C_pv_OM_total = self.C_pv_OM * self.power_AC * self.oversize

    def set_C_pv_cap(self, C_pv_cap):
        self.C_pv_cap = C_pv_cap
        self.C_pv_cap_total = self.C_pv_cap * self.power_AC * self.oversize

    def _calc_power_DC(self):
        self.power_DC = self.power_AC * self.oversize / self.eta_inv  #[kW] design output DC
        self._get_PVWatts()
        

    def _get_PVWatts(self):
        #info on API: http://developer.nrel.gov/docs/solar/pvwatts/v6/
        
        api_key = "7C11bMIUFytK1Y0S06jj8ByfUCDGqlEYYARFq7fF"
        website = "https://developer.nrel.gov/api/pvwatts/v6.json?"
        
        req = website + "api_key=%s&" %api_key
        req = req + "lat=%s&" %self.latitude
        req = req + "lon=%s&" %self.longitude
        req = req + "system_capacity=%s&" %(self.power_DC * 1e3) # (input in kW)  
        req = req + "azimuth=%s&" %np.rad2deg(self.azimuth)
        req = req + "tilt=%s&" %np.rad2deg(self.tilt)
        req = req + "dc_ac_ratio=%s&" %self.oversize
        req = req + "array_type=2&"     #single-axis tracking
        req = req + "module_type=0&"
        req = req + "losses=%s&" %self.losses 
        req = req + "dataset=intl&"
        req = req + "timeframe=hourly&"
        req = req + "inv_eff=%s" %(self.eta_inv * 100)
        
        self.request = requests.get(req)
        
        self.data_dict = self.request.json()
        
        self.outputs = self.data_dict['outputs']
        
        self.station_info = self.data_dict['station_info']
        
        self.DNI = np.array(self.outputs['dn']) / 1000      #[kW/m^2]
        self.DC = np.array(self.outputs['dc']) / 1e6       #[MW] (output in W)
        self.AC = np.array(self.outputs['ac']) / 1e6       #[MW]
        self._calc_eta_inv_rt()
        self._calc_DC_to_AC()
        self._calc_DC_curtailed()
        

    def _calc_eta_inv_rt(self):
        P_DC0 = self.power_AC / self.eta_inv   #[kW]
        zeta = self.DC / P_DC0                           #[-]
        zeta[zeta > 1] = 1
        self.eta_inv_rt = self.eta_inv / 0.9637 * \
            (-0.0162 * zeta - 0.0059 / zeta + 0.9858)
    
    def _calc_DC_to_AC(self):
        #Here it's calculated how much DC power was required to generate
        #the AC power that the plant produced.
        self.DC_to_AC = self.AC / self.eta_inv_rt
    
    def _calc_DC_curtailed(self):   
        
        self.DC_curtailed = self.DC - self.DC_to_AC
        #Next, clean up small errors between the DC_to_AC that I calc'd and that 
        #from the AC and DC time series returned by PVWatts when there is no
        #curtailment taking place. Some error will still be present during
        #times of curtailment.
        self.DC_curtailed[self.DC_curtailed < 0] = 0 
        
class Time():
    def __init__(self, start_time, end_time, dt = 0.5):
        self.start_time = pd.Timestamp(start_time)
        self.end_time = pd.Timestamp(end_time)
        self.dt = dt
        self.interval = pd.Timedelta(self.end_time - self.start_time).days * 24
        self.N = int(self.interval / dt)
        self.samp_index = pd.date_range(start = self.start_time, 
                                           end = self.end_time, 
                                           freq = '%iT' %int(60*dt))
##      This should be reviewed - currently defined as halfway in time bins
        self.int_index = self.samp_index[:-1] + pd.Timedelta(hours = dt / 2)
        
        self.leapyear = len(self.int_index[(self.int_index.day == 29) & \
                                           (self.int_index.month == 2)]) > 0

        self.yearfrac = self.interval / ((365+self.leapyear) * 24)
        

class PowerCycle():
    #C_pc comes from turchi-techreport-2013 power block cost of 
    #$1000/kW-e generation capacity plus $350/KW-e balance of plant 
    #(steam generator)
    def __init__(self, P_pc_out_max, eta_pc = 0.4, 
                 C_pc_cap = 1350e3, C_pc_OM = 0, eta_pc_type = 'fixed',
                 C_pc_cap_type = 'fixed'):
        
        """
Captial costs of power cycles and the efficiency of the power cycle can 
be set to obey power laws and such that 
eta_st_type in ['fixed', 'lovegrove']
C_pc_cap_type in ['fixed', 'lovegrove']

C_pc_cap has the dual purpose of setting the fixed price or setting the
reference price for the 'lovegrove' mode of valuing the power cycle.

eta_pc has the dual purpose of setting the fixed efficiency or setting the
reference efficiency for the 'lovegrove' mode of setting the power cycle
efficiency."""
        
        self.C_pc_OM = C_pc_OM  #This might be more of a fixed cost ...
        self.eta_pc_type = eta_pc_type
        self.C_pc_cap_type = C_pc_cap_type
        self.C_pc_cap = C_pc_cap
        self.set_P_pc_out_max(P_pc_out_max)
        self._set_eta_pc(eta_pc)

        self.C_pc_OM_total = P_pc_out_max * C_pc_OM
        
    def set_P_pc_out_max(self, P_pc_out_max):
        print('Setting power cycle size to %.2f' %P_pc_out_max)
        self.P_pc_out_max = P_pc_out_max
        self._set_C_pc_cap_total()
        if self.eta_pc_type == 'lovegrove': self._set_eta_pc(False)

    def _set_C_pc_cap_total(self):
        if self.C_pc_cap_type == 'fixed':
            self.C_pc_cap_total = self.C_pc_cap * self.P_pc_out_max
        elif self.C_pc_cap_type == 'lovegrove':
            self.C_pc_cap_total = 100. * self.C_pc_cap * \
                (self.P_pc_out_max / 100.) ** 0.7       ## This needs reviewing

    def _set_eta_pc(self, eta_pc):
        if self.eta_pc_type == 'fixed':
            self.eta_pc = eta_pc
        elif self.eta_pc_type == 'lovegrove':
            A = 0.5856736175970182
            B = 0.05961115093356057
            self.eta_pc = 0.42 * (1 - A * np.exp(-self.P_pc_out_max * B))
            
class Financials():
    def __init__(self, lifetime = 20, period = 1./12, acc = 0.067, 
                 USD2AUD = 1.5):
        """All input prices are in USD
        Note that the default acc is from NREL's ATB spreadsheet"""
        self.USD2AUD = USD2AUD     
        self.period = period
        self.lifetime = lifetime
        self.periods = int(lifetime / period)
        self.acc = acc
        self.pcc = 1-(1-acc)**period     #periodic cost of capital

    def calc_annual_cost(self, capital_cost):
        """
    This function calculates the annual cost which is appropriate as the
    cost coefficient to be used in the the MiniZinc optimisations.
        
    Cost: capital cost of the system ($)
    lifetime: expected lifetime of the system (years)
    period: compounding period (years)
    acc: annual cost of capital as fraction between 0 and 1.
    """
        periodic_cost = capital_cost * (self.pcc * (1 + self.pcc)**self.periods) / \
            (((1 + self.pcc)**self.periods) - 1)
        annual_cost = periodic_cost / self.period
    
        return annual_cost  
    
    def calc_capital_cost(self, annual_cost):
        periodic_cost = self.period * annual_cost
        
        return periodic_cost * (((1 + self.pcc)**self.periods) - 1) / \
            (self.pcc * (1 + self.pcc)**self.periods)
    
class Storage():
    def __init__(self, P_st_in_max = 50, E_st_max = 200., E_st_init_frac = 0, 
                 C_st_cap = 25e3, C_st_OM = 0, C_rh_single = 100., C_rh_OM = 0, 
                 eta_rh = 1., V_rh = 240, I = 25, eta_st = None,
                 C_st_cap_type = 'fixed'):
        
        ##This need to be set up to take tank size OR cost of tank and
        #heat size OR max power into tank. The respective former options allow
        #the simulation to find an optimum size for each while the latter
        #fixed these sizes.
        #Typical values for the costs of the thermal storage and the
        #resistance heater can be C_st = 30e3, C_rh = 16.7e3

        ##Not sure it's worth including a loss but can try later
        #self.eta_st = 0.99**(1./24)
        
        """
Captial costs of power cycles and the efficiency of the power cycle can 
be set to obey power laws and such that 
C_st_cap_type in ['fixed', 'lovegrove']"""

        self.I = I
        self.V_rh = V_rh
        self.C_rh_single = C_rh_single
        self.E_st_init_frac = E_st_init_frac
        self.C_st_cap_type = C_st_cap_type
        self.C_st_OM = C_st_OM
        self.C_st_cap = C_st_cap
        self.eta_rh = eta_rh
        self.eta_st = None
        self.calc_rh_cost()
        self.C_rh_OM = C_rh_OM
        self.set_P_st_in_max(P_st_in_max)
        self.set_E_st_max(E_st_max)
        
    def set_eta_st(self, eta_st):
        self.eta_st = eta_st
        self._set_C_st_cap_OM()

    def set_E_st_max(self, E_st_max):
        print('Setting storage size %.2f' %E_st_max)
        self.E_st_max = E_st_max
        if self.eta_st != None: self._set_C_st_cap_OM()

    def _set_C_st_cap_OM(self):
        if self.C_st_cap_type == 'lovegrove':
            self.C_st_cap_total = 1429. * self.C_st_cap * \
                (self.E_st_max / self.eta_st / 1429.) ** 0.8
        self.C_st_cap_total = self.C_st_cap * self.E_st_max / self.eta_st
        self.C_st_OM_total = self.C_st_OM * self.E_st_max / self.eta_st
        
    def set_P_st_in_max(self, P_st_in_max):
        print('Setting storage charge rate %.2f' %P_st_in_max)
        self.P_st_in_max = P_st_in_max 
        self.C_rh_cap_total = self.C_rh_cap * P_st_in_max
        self.C_rh_OM_total = self.C_rh_OM * P_st_in_max
        
    def calc_rh_cost(self):
        #7/3/2020
        #This product: https://www.alibaba.com/product-detail/10kw-Immersion-Heater-400v_62345866250.html?spm=a2700.details.deiletai6.9.6998U3UAU3UABf
        #costs 100 USD per elements. At 240V, this translates to 6kW, or 167 units per
        #megawatt, which is $16,700 per megawatt. Annualised assuming 30 year plant
        #lifetime and 7.25% cost of money with monthly repayments, this equals about: 
        #25 amps of current calculated from P=VI with V = 400 and P = 10kW 
 
        P_rh_single = self.I * self.V_rh  
        
        N_rh = ceil(1e6 / P_rh_single)          #number per megawatt capacitiy
        self.C_rh_cap = self.C_rh_single * N_rh         #$ per megawatt
        
    
def write_summary_data(sim):

    optimise = sim.opt_int != 0
    
    if sim.Time.yearfrac != 1.0 and sim.Time.start_time.month != 1 and \
        sim.Time.start_time.day != 1:
        print('not saving partial years to database!')
        return
    
    print('Saving data to database ...')
    
    table = ['summary','optimised'][optimise]
    
    print(table)
    
    db='storage_value_PF_PVplant_summary.db'
    cols_recording = ['region', 'E_st_max', 'P_st_in_max',
                      'P_pc_out_max', 'oversize', 'eta_st', 'eta_inv', 
                      'year', 'P_in_grid_total', 'P_out_grid_total', 
                      'P_in_plant_total', 'P_out_pass_total', 
                      'P_inv_in_total', 'P_st_in_plant_cur_total',
                      'P_st_in_plant_vol_total', 'C_in_grid_total', 
                      'C_out_grid_total', 'C_out_pass_total',
                      'C_inv_out_total', 'C_st_in_plant_cur_total',
                      'C_st_in_plant_vol_total']
    
    if optimise: 
        cols_recording += ['opt_os', 'opt_rh', 'opt_st', 'opt_pc', 
                           'C_pv_cap', 'C_pv_OM', 'C_st_cap', 'C_st_OM',
                           'C_rh_cap', 'C_rh_OM', 'C_pc_cap', 'C_pc_OM']

        if sim.Plant:
            C_pv_cap = [np.nan, sim.Plant.C_pv_cap][sim.opt_os]
            C_pv_OM = [np.nan, sim.Plant.C_pv_cap][sim.opt_os]
        else: 
            C_pv_cap = np.nan
            C_pv_OM = np.nan
        
        C_st_cap = [np.nan, sim.Storage.C_st_cap][sim.opt_st]
        C_st_OM = [np.nan, sim.Storage.C_st_OM][sim.opt_st]
        C_rh_cap = [np.nan, sim.Storage.C_rh_cap][sim.opt_rh]
        C_rh_OM = [np.nan, sim.Storage.C_rh_OM][sim.opt_rh]
        C_pc_cap = [np.nan, sim.PowerCycle.C_pc_cap][sim.opt_pc]
        C_pc_OM = [np.nan, sim.PowerCycle.C_pc_OM][sim.opt_pc]
        
    
    P_in_grid_total = sum(sim.Time.dt * sim.Results.P_in_grid.values)
    P_out_grid_total = sum(sim.Time.dt * sim.Results.P_out_grid.values)
    C_in_grid_total = sum((sim.Time.dt * sim.Results.P_in_grid * sim.c).values)
    C_out_grid_total = sum((sim.Time.dt * sim.Results.P_out_grid * sim.c).values)        
    
    if sim.Plant:
        P_out_pass_total = sum(sim.Time.dt * sim.Results.P_out_pass.values)
        C_out_pass_total = sum((sim.Time.dt * sim.Results.P_out_pass * sim.c).values)
        P_in_plant_total = sum(sim.Time.dt * sim.P_in_plant.values)
        P_inv_in_total = sum(sim.Time.dt * sim.Results.P_inv_in.values)
        P_st_in_plant_cur_total = sum(sim.Time.dt * \
                                  sim.Results.P_st_in_plant_cur.values)
        P_st_in_plant_vol_total = sum(sim.Time.dt * \
                                  sim.Results.P_st_in_plant_vol.values)

        C_inv_out_total = sum((sim.Time.dt * sim.Results.P_inv_in * \
                               sim.Inverter.eta_inv * sim.c).values)
        C_st_in_plant_cur_total = sum((sim.Time.dt * sim.Results.P_st_in_plant_cur * \
                                       sim.Inverter.eta_inv * sim.c).values)
        C_st_in_plant_vol_total = sum((sim.Time.dt * sim.Results.P_st_in_plant_vol * \
                                       sim.Inverter.eta_inv * sim.c).values)
    else:
        P_out_pass_total = np.nan
        C_out_pass_total = np.nan
        P_in_plant_total = np.nan
        P_inv_in_total = np.nan
        P_st_in_plant_cur_total = np.nan
        P_st_in_plant_vol_total = np.nan
        C_inv_out_total = np.nan
        C_st_in_plant_cur_total = np.nan
        C_st_in_plant_vol_total = np.nan
                      
    idx = ['region', 'E_st_max', 'P_st_in_max','P_pc_out_max', 
           'oversize', 'eta_st', 'eta_inv', 'year']
    
    DATAdict = {'region' : sim.region, 
                'E_st_max': sim.Storage.E_st_max, 
                'P_st_in_max': sim.Storage.P_st_in_max,
                'P_pc_out_max': sim.PowerCycle.P_pc_out_max, 
                'oversize': 0,
                'eta_st': sim.Storage.eta_st, 
                'eta_inv': sim.Inverter.eta_inv, 
                'year': sim.Time.int_index[0].year, 
                'P_in_grid_total': P_in_grid_total,
                'P_out_grid_total': P_out_grid_total,
                'P_in_plant_total': P_in_plant_total,
                'P_out_pass_total': P_out_pass_total,
                'P_inv_in_total': P_inv_in_total,
                'P_st_in_plant_cur_total': P_st_in_plant_cur_total,
                'P_st_in_plant_vol_total': P_st_in_plant_vol_total,
                'C_in_grid_total': C_in_grid_total,
                'C_out_grid_total': C_out_grid_total,
                'C_out_pass_total': C_out_pass_total,
                'C_inv_out_total': C_inv_out_total,
                'C_st_in_plant_cur_total': C_st_in_plant_cur_total,
                'C_st_in_plant_vol_total': C_st_in_plant_vol_total}

    if optimise: 
        idx += ['opt_os', 'opt_rh', 'opt_st', 'opt_pc']
        DATAdict.update({'opt_os': sim.opt_os, 
                         'opt_rh': sim.opt_rh, 
                         'opt_st': sim.opt_st,
                         'opt_pc': sim.opt_pc,
                         'C_pv_cap': C_pv_cap,
                         'C_pv_OM': C_pv_OM,
                         'C_st_cap': C_st_cap,
                         'C_st_OM': C_st_OM, 
                         'C_rh_cap': C_rh_cap,
                         'C_rh_OM': C_rh_OM,
                         'C_pc_cap': C_pc_cap,
                         'C_pc_OM': C_pc_OM})
    
    DATA = pd.DataFrame(columns=cols_recording)
    DATA = DATA.append(DATAdict, ignore_index = True)

    if sim.Plant: DATA['oversize'] = sim.Plant.oversize
    sm.create_table(table, db, cols_recording,
                    create_unique_idx=True, idx_cols=idx)
    sm.replace_into_db(DATA, db, table, cols_recording)     


def write_normalised_data(sim):
    
    if sim.Plant: 
        print('plant not yet supported for dimensionless database')
        return
    
    optimise = sim.opt_int != 0
    
    if sim.Time.yearfrac != 1.0 and sim.Time.start_time.month != 1 and \
        sim.Time.start_time.day != 1:
        print('not saving partial years to database!')
        return
    
    print('Saving data to database ...')
    
    table = ['summary','optimised'][optimise]
    
    print(table)
    
    db='storage_value_PF_PVplant_normalised2.db'
    cols_recording = ['region', 'year', 'SPH', 'CPH', 'eta_pc', 
                      'P_in_grid_total', 'P_out_grid_total', 
                      'C_in_grid_total', 'C_out_grid_total']
    
    if optimise: 
        cols_recording += ['opt_rh', 'opt_pc', 
                           'C_rh_cap', 'C_rh_OM', 
                           'C_pc_cap', 'C_pc_OM']
        
        C_rh_cap = [np.nan, sim.Storage.C_rh_cap][sim.opt_rh]
        C_rh_OM = [np.nan, sim.Storage.C_rh_OM][sim.opt_rh]
        C_pc_cap = [np.nan, sim.PowerCycle.C_pc_cap][sim.opt_pc]
        C_pc_OM = [np.nan, sim.PowerCycle.C_pc_OM][sim.opt_pc]
        
    
    P_in_grid_total = sum(sim.Time.dt * sim.Results.P_in_grid.values)
    P_out_grid_total = sum(sim.Time.dt * sim.Results.P_out_grid.values)
    C_in_grid_total = sum((sim.Time.dt * sim.Results.P_in_grid * sim.c).values)
    C_out_grid_total = sum((sim.Time.dt * sim.Results.P_out_grid * sim.c).values)        
    
    idx = ['region', 'year', 'SPH', 'CPH', 'eta_pc']
    
    DATAdict = {'region' : sim.region, 
                'SPH': sim.SPH,
                'CPH': sim.CPH,
                'year': sim.Time.int_index[0].year, 
                'eta_pc': sim.PowerCycle.eta_pc,
                'P_in_grid_total': P_in_grid_total,
                'P_out_grid_total': P_out_grid_total,
                'C_in_grid_total': C_in_grid_total,
                'C_out_grid_total': C_out_grid_total}

    if optimise: 
        idx += ['opt_rh', 'opt_pc']
        DATAdict.update({'opt_rh': sim.opt_rh, 
                         'opt_pc': sim.opt_pc,
                         'C_rh_cap': C_rh_cap,
                         'C_rh_OM': C_rh_OM,
                         'C_pc_cap': C_pc_cap,
                         'C_pc_OM': C_pc_OM})
    
    DATA = pd.DataFrame(columns=cols_recording)
    DATA = DATA.append(DATAdict, ignore_index = True)

    if sim.Plant: DATA['oversize'] = sim.Plant.oversize
    sm.create_table(table, db, cols_recording,
                    create_unique_idx=True, idx_cols=idx)
    sm.replace_into_db(DATA, db, table, cols_recording) 

def write_schedule_data(sim):

    if sim.Time.yearfrac != 1.0 and sim.Time.start_time.month != 1 and \
        sim.Time.start_time.day != 1:
        print('not saving partial years to database!')
        return
    
    if sim.optimisation:
        print('not saving optimisation runs to database!')
        return
    print('Saving data to database ...')
    
    table = 'schedule'
    db='storage_value_PF_PVplant_schedules.db'
    cols_recording = [ 'region', 'E_st_max', 'P_st_in_max',
                      'P_pc_out_max', 'oversize', 'eta_st', 'eta_inv', 
                      'date_time', 'P_in_grid', 'P_out_grid', 'P_in_plant',  
                      'P_out_pass', 'P_inv_in', 'c']
    idx = ['region', 'E_st_max', 'P_st_in_max','P_pc_out_max', 
           'oversize', 'eta_st', 'eta_inv', 'date_time']
    
    DATA = pd.DataFrame({'region' : sim.region, 
                         'E_st_max': sim.Storage.E_st_max, 
                         'P_st_in_max': sim.Storage.P_st_in_max,
                         'P_pc_out_max': sim.PowerCycle.P_pc_out_max, 
                         'oversize': None,
                         'eta_st': sim.Storage.eta_st, 
                         'eta_inv': sim.Inverter.eta_inv, 
                         'date_time': sim.Time.int_index.astype(str), 
                         'P_in_grid': sim.Results.P_in_grid.values, 
                         'P_out_grid': sim.Results.P_out_grid.values, 
                         'P_in_plant': sim.P_in_plant.values,  
                         'P_out_pass': sim.Results.P_out_pass.values, 
                         'P_inv_in': sim.Results.P_inv_in.values, 
                         'c': sim.c.values})

    if sim.Plant: DATA['oversize'] = sim.Plant.oversize
    sm.create_table(table, db, cols_recording,
                    create_unique_idx=True, idx_cols=idx)
    sm.replace_into_db(DATA, db, table, cols_recording)
