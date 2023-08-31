# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 15:15:42 2017

@author: jecu01
"""

from parse_SP import parse_NEM_data
from parse_DNI import get_minute_solar
from run_minizinc import run_minizinc
from projdirs import mzn_model_dir, gms_model_dir, temp_dir, results_dir
from utilities import savefile, openfile
import pylab as pl
import sys
import cPickle as pickle
from os import remove
from copy import copy
from datetime import datetime, timedelta

sys.dont_write_bytecode = True

class SimInput():
    
    def __init__(self, program = 0, model = 0, solver = 0, design_DNI = 1000, 
                 DNI_location = 'Wagga', power_region = 'NSW', timestep_mins = 30, 
                 start_date = '20130901', end_date = '20130907', SM = 2, SH = 8, 
                 P_max = 50e6, eta_rcv = 0.9, eta_cycle = 0.4, eta_tank = 0.95, 
                 P_0 = 0, Qdot_re0 = 0, Qdot_se0 = 0, Qdot_rs0 = 0, Qdot_pass0 = 0, 
                 Q_s0_frac = 0.5, P_min_frac = 0.4, Q_smin_frac = 45./700, N_on = 0, 
                 N_off = 0, R_dis = 35e6, R_chg = 80e6, T_down = 2, T_up = 2, U_0 = 0,
                 robust = False, G_err_type = 0, G_err_scale = 0.2, numscenarios = 1, 
                 P_scenario = [1], gamma = 1):
              
        """
G_err_type      {0, 1, 2}    0: constant   1: fractional   2: altitude-related
G_err_scale     scaling the magnitude of the error
"""
              
        self.model_dir = [mzn_model_dir, gms_model_dir][program]            
        
        self.set_solver(solver)
        self.set_model(model)
        self.set_program(program)
        
        self.DNI_location = DNI_location
        self.power_region = power_region
        self.timestep_mins = timestep_mins
        
        self.start_date = datetime(int(start_date[:4]), int(start_date[4:6]), int(start_date[6:8]))
        self.end_date = datetime(int(end_date[:4]), int(end_date[4:6]), int(end_date[6:8])) + timedelta(days = 1)
        
        self.load_DNI_data()
        self.load_SP_data()
        self.design_DNI = design_DNI
        self.SM = SM
        self.SH = SH
        
        self.Q_s0_frac = Q_s0_frac
        self.Q_smin_frac = Q_smin_frac
        self.P_min_frac = P_min_frac
        self.G_err_scale = G_err_scale
        
        self.calc_duration()        
        self.calc_N_timesteps()
        self.P_max = P_max
        self.eta_rcv = eta_rcv
        self.eta_cycle = eta_cycle
        self.eta_tank = eta_tank
        self.P_0 = P_0
        self.Qdot_re0 = Qdot_re0
        self.Qdot_se0 = Qdot_se0
        self.Qdot_rs0 = Qdot_rs0
        self.Qdot_pass0 = Qdot_pass0
        self.robust = robust
        self.calc_Q_smax() 
        self.calc_Q_smin()
        self.calc_Q_s0()
        
        if self.model == 'dominguez':
            self.N_on = N_on
            self.N_off = N_off
            self.R_dis = R_dis
            self.R_chg = R_chg
            self.T_down = T_down
            self.T_up = T_up
            self.U_0 = U_0
            self.calc_P_min()            
            
            if self.robust:
                self.numscenarios = numscenarios
                self.P_scenario = P_scenario
                self.calc_G_err()
                #ALGORITHM FOR CONSTANT ERROR (with tapering for low DNI)
#                self.G_err = pl.array([G_err_scale * self.design_DNI] * self.N_timesteps)    #should be updates for different error distributions in a separate function
#                self.G_err[pl.where(self.G==0)] = 0
#                low_DNI_values = pl.where(self.G <= G_err_scale * design_DNI * 2)
#                self.G_err[low_DNI_values] = G_err_scale * self.G[low_DNI_values]
                self.gamma = gamma
            
        self.calc_A_field()
        
        self.modelled = False

    def calc_days(self):
        pass

    def set_solver(self, value):
        self.solver = ['SCIP', 'CPLEX'][value]
        self.modelled = False
    
    def set_model(self, value):    
        self.model = ['cirocco', 'dominguez'][value]
        self.modelled = False

    def set_program(self, value):
        self.program = ['minizinc', 'GAMS'][value]
        self.modelled = False

    def calc_G_err(self):
        self.G_err = pl.zeros(self.N_timesteps)
        self.G_err[~pl.isnan(self.G)] = self.G[~pl.isnan(self.G)] * self.G_err_scale
    
    def calc_duration(self):
        "duration in hours"
        self.duration = ((self.end_date - self.start_date).days) * 24 

    def calc_Q_smax(self):
        self.Q_smax = self.SH * self.P_max / (self.eta_cycle * self.eta_tank)
        self.modelled = False
        
    def calc_Q_smin(self):
        self.Q_smin = self.Q_smax * self.Q_smin_frac
        
    def calc_A_field(self):
        self.A_field = self.P_max * self.SM / (self.design_DNI * self.eta_cycle * self.eta_rcv)
        self.modelled = False
        
    def calc_N_timesteps(self):
        self.N_timesteps = int(self.duration * 60 / self.timestep_mins)
        self.modelled = False
        
    def calc_Q_s0(self):
        self.Q_s0 = self.Q_smax * self.Q_s0_frac
        self.modelled = False
          
    def set_Q_smax(self, value, dependent):
        print """The variable that was changed to accord with this change was %s.
You can select from 'SH', 'P_max', 'eta_cycle', 'eta_tank'.""" %dependent
        self.Q_smax = value
        if dependent == 'SH':
            SH = self.Q_smax * self.eta_cycle * self.eta_tank / self.P_max
            self.set_SH(SH)
            self.calc_Q_smin()      #so that is scales with the new value of Q_smax according to Q_smin_frac (not fully integrated)
        if dependent == 'P_max':
            P_max = self.Q_smax * self.eta_cycle * self.eta_tank / self.SH
            self.set_P_max(P_max)
            self.calc_Q_smin()
        if dependent == 'eta_cycle':
            eta_cycle = self.SH * self.P_max / (self.eta_tank * self.Q_smax)
            self.set_eta_cycle(eta_cycle)
            self.calc_Q_smin()
        if dependent == 'eta_tank':
            eta_tank = self.SH * self.P_max / (self.eta_cycle * self.Q_smax)    
            self.set_eta_tank(eta_tank)
            self.calc_Q_smin()
        self.modelled = False

    def set_A_field(self, value, dependent = 'P_max'):
        print """The variable that was changed to accord with this change was %s.
You can select from 'P_max', 'SM', 'design_DNI', 'eta_cycle', 'eta_rcv'.""" %dependent
        self.A_field = value
        if dependent == 'P_max':
            P_max = self.A_field * self.design_DNI * self.eta_cycle * self.eta_rcv / self.SM
            self.set_P_max(P_max)
        if dependent == 'SM':
            SM = self.A_field * self.design_DNI * self.eta_cycle * self.eta_rcv / self.P_max
            self.set_SM(SM)
        if dependent == 'design_DNI':
            design_DNI =  self.P_max * self.SM / (self.A_field * self.eta_cycle * self.eta_rcv)           
            self.set_design_DNI(design_DNI)
        if dependent == 'eta_cycle':
            eta_cycle = self.P_max * self.SM / (self.design_DNI * self.A_field * self.eta_rcv)
            self.set_eta_cycle(eta_cycle)
        if dependent == 'eta_rcv':
            eta_rcv = self.P_max * self.SM / (self.design_DNI * self.eta_cycle * self.A_field)
            self.set_eta_rcv(eta_rcv)
        self.modelled = False

    def set_P_min(self, value, dependent = 'P_min_frac'):
        print """The variable that was changed to accord with this change was %s.
You can select from 'P_min_frac' and 'P_max'.""" %dependent
        self.P_min = value        
        if dependent == 'P_min_frac':
            P_min_frac = self.P_min / self.P_max
            self.set_P_min_frac(P_min_frac)
        if dependent == 'P_max':
            P_max = self.P_min / self.P_min_frac
            self.set_P_max(P_max)
        self.modelled = False
    
    def set_N_timesteps(self, value, dependent = 'duration'):
        print """The variable that was changed to accord with this change was %s.
You can select from 'duration' and 'timestep'.""" %dependent
        self.N_timesteps = value        
        if dependent == 'duration':
            duration = self.N_timesteps / self.timestep
            self.set_duration(duration)
        if dependent == 'timestep':
            timestep = self.N_timesteps / self.duration
            self.set_timestep(timestep)
        self.modelled = False

    def set_Q_s0(self, value, dependent = 'Q_s0_frac'):
        print """The variable that was changed to accord with this change was %s.
You can select from 'Q_smax' and 'Q_s0_frac'.""" %dependent
        self.Q_s0 = value
        if dependent == 'Q_s0_frac':
            Q_s0_frac = self.Q_s0 / self.Q_smax
            self.set_Q_s0_frac(Q_s0_frac)
        if dependent == 'Q_smax':
            Q_smax = self.Q_s0 / self.Q_s0_frac
            subdependent = raw_input("""To change Q_smax you need to 
choose a dependent variable from the following: 'P_max', 'SM', 'design_DNI', 'eta_cycle', 'eta_rcv': """)
            self.set_Q_smax(Q_smax, subdependent)
        self.modelled = False

    def load_DNI_data(self):
        self.time, self.G = get_minute_solar(self.DNI_location, self.start_date, self.end_date, self.timestep_mins)
        "Start dates are inclusive."
        self.modelled = False

    def load_SP_data(self):
        self.SP = parse_NEM_data(self.power_region, self.start_date, self.end_date, self.timestep_mins)   
        ###FUNCTION NEEDS TO BE UPDATED TO COVER MULTIPLE MONTHS
        self.modelled = False

    def set_eta_rcv(self, value, debug = False):
        self.eta_rcv = value
        if debug: A_field_temp = copy(self.A_field)        
        self.calc_A_field()
        if debug: print A_field_temp - self.A_field        
        self.modelled = False
                                
    def set_eta_cycle(self, value):
        self.eta_cycle(value)
        self.calc_Q_smax()
        self.calc_A_field()
        self.modelled = False

    def set_G_err_scale(self, value):
        self.G_err_scale = value
        self.calc_G_err()
        self.modelled = False

    def set_gamma(self, value):
        if self.robust == False: print "OK, but this model is not robust."
        self.gamma = value
        
    def set_eta_tank(self, value):
        self.eta_tank = value
        self.modelled = False
        
    def set_SH(self, value):
        self.SH = value
        self.calc_Q_smax()
        self.calc_Q_s0()
        self.modelled = False
        
    def set_P_max(self, value):
        self.P_max
        if self.model == 'dominguez':
            self.calc_P_min()
        self.modelled = False
        
    def set_SM(self, value):    #set dependents afterwards..
        self.SM = value
        self.calc_A_field()
        self.modelled = False
        
    def set_design_DNI(self, value):
        self.design_DNI = value
        self.modelled = False

    def set_DNI_location(self, value):
        self.DNI_location = value 
        self.load_DNI_data()
        self.modelled = False
        
    def set_month(self, value):
        self.month = value
        self.load_DNI_data()
        self.modelled = False
        
    def set_year(self, value):
        self.year = value
        self.modelled = False

    def set_timestep(self, value):
        self.timestep = value
        self.calc_N_timesteps()
        self.modelled = False
        
    def set_power_region(self, value):
        """value is string corresponding to the energy market region. This part
        of the program needs some work because it only works for NEM time series."""        
        self.power_region = value        
        self.modelled = False
    
    def set_duration(self, value):
        self.duration = value
        self.calc_N_timesteps()
        self.modelled = False
        
    def set_Q_s0_frac(self, value):
        self.Q_s0_frac = value
        self.modelled = False
    
    def calc_P_min(self):
        "Does not apply to cirocco model"
        self.P_min = self.P_min_frac * self.P_max
        self.modelled = False

    def set_P_min_frac(self, value):
        "Does not apply to the cirocco model"        
        self.P_min_frac = value
        self.calc_P_min()
        self.modelled = False
             
    def save(self, filename):
        with open(filename, 'w') as f:
            pickle.dump(self, f)
                     

class Simulation():
    def __init__(self, SimInput):

        """
solver = {0,1}     %{SCIP, CPLEX}    
program = {0,1}    %{minizinc, GAMS}
model = {0,1}      %{cirocco, dominguez}
robust = {0,1}     %{non-robust, robust}
"""   #solver type does not currently apply to the execution of the minizinc models
        self.SimInput = SimInput
        self.robuststr = ['', '_robust'][self.SimInput.robust]
                    
    def _make_dzn_file(self):
        
        with open(self.SimInput.model_dir + self.SimInput.model + self.robuststr + '_GENERIC.dzn', 'r') as f: 
                sim_file = ''.join(i for i in f)

        DNI_str = '['    
        DNI_err_str = '['
        SP_str = '[|'        
        GAMMA_str = '['
        Pi_str = '['        
        
        for i in self.SimInput.G[:self.SimInput.N_timesteps]:
            if pl.isnan(i): 
                DNI_str += '0.0, '
            else: DNI_str += str(i) + ', '
        DNI_str = DNI_str[:-2] + ']'
    
        for i in self.SimInput.SP[:self.SimInput.N_timesteps]:
            SP_str += str(i) + ', '
        SP_str = SP_str[:-2] + '|]'    #currently only suitable for one scenario
        
        for i in self.SimInput.G_err[:self.SimInput.N_timesteps]:
            DNI_err_str += str(i) + ', '
        DNI_err_str = DNI_err_str[:-2] + ']'
        
        for i in range(self.SimInput.N_timesteps):
            GAMMA_str += str(self.SimInput.gamma) + ', '
        GAMMA_str = GAMMA_str[:-2] + ']'

        for i in self.SimInput.P_scenario:
            Pi_str += '%f' %i + ', '
        Pi_str = Pi_str[:-2] + ']'
                
        sim_file = sim_file.replace('TIMESTEPS', str(self.SimInput.N_timesteps))
        sim_file = sim_file.replace('MAXTURB', str(self.SimInput.P_max))
        sim_file = sim_file.replace('MAXSTORE', str(self.SimInput.Q_smax))
        sim_file = sim_file.replace('MINSTORE', str(self.SimInput.Q_smin))
        sim_file = sim_file.replace('ETA_RCV', str(self.SimInput.eta_rcv))
        sim_file = sim_file.replace('ETA_CYCLE', str(self.SimInput.eta_cycle))
        sim_file = sim_file.replace('ETA_TANK', str(self.SimInput.eta_tank))
        sim_file = sim_file.replace('A_FIELD', str(self.SimInput.A_field))
        sim_file = sim_file.replace('POWER0', str(self.SimInput.P_0))
        sim_file = sim_file.replace('RCVELEC0', str(self.SimInput.Qdot_re0))
        sim_file = sim_file.replace('STOREELEC0', str(self.SimInput.Qdot_se0))
        sim_file = sim_file.replace('RCVSTORE0', str(self.SimInput.Qdot_rs0))
        sim_file = sim_file.replace('PASS0', str(self.SimInput.Qdot_pass0))
        sim_file = sim_file.replace('TANK0', str(self.SimInput.Q_s0))
        if self.SimInput.model == 'dominguez':        
            sim_file = sim_file.replace('MINTURB', str(self.SimInput.P_min))
            sim_file = sim_file.replace('ONOFF0', ['0', '1'][self.SimInput.U_0])
            sim_file = sim_file.replace('N_ON', str(self.SimInput.N_on))            
            sim_file = sim_file.replace('N_OFF', str(self.SimInput.N_off))
            sim_file = sim_file.replace('DISCHARGE', str(self.SimInput.R_dis))
            sim_file = sim_file.replace('CHARGE', str(self.SimInput.R_chg))
            sim_file = sim_file.replace('MINDOWNTIME', str(self.SimInput.T_down))
            sim_file = sim_file.replace('MINUPTIME', str(self.SimInput.T_up))
            if self.SimInput.robust:
                sim_file = sim_file.replace('NUMSCENARIOS', str(self.SimInput.numscenarios))
                sim_file = sim_file.replace('P_SCENARIO', str(Pi_str))
                sim_file = sim_file.replace('GAMMALIST', str(GAMMA_str))
                sim_file = sim_file.replace('DNIERR', str(DNI_err_str))
            
            
        sim_file = sim_file.replace('DNILIST', str(DNI_str))
        sim_file = sim_file.replace('SPLIST', str(SP_str))
        
        self.sim_file = sim_file
        
        with open(temp_dir + self.SimInput.model + self.robuststr + '_current.dzn', 'w') as f:
            f.write(self.sim_file)

    def save_dzn_file(self, filename):
        with open(filename, 'w') as f:
            f.write(self.sim_file)
            
    def del_dzn_file(self):
        "just for some housekeeping should I need it later"
        remove(temp_dir + self.SimInput.model + self.robuststr + '_current.dzn')

    def run_optimization(self):
        "Eventually will allow to run GAMS and minizinc and choose the solver"
        
        if self.SimInput.program == 'minizinc':
            self._make_dzn_file()
            model_results = run_minizinc(self.SimInput.model_dir + self.SimInput.model + self.robuststr + "_model.mzn",
                                         temp_dir + self.SimInput.model + self.robuststr + "_current.dzn")
            
            #append the time series information to the results from the model
            model_results = pl.array([list(self.SimInput.time)] + list(model_results))
            
        elif self.SimInput.program == 'GAMS':
            print "GAMS is not yet implemented"
    
        self.Results = Results(*model_results)
        self.SimInput.modelled = True
        
        
class Results():
    def __init__(self, timestamps, Q_s, Qdot_re, Qdot_se, Qdot_rs, Qdot_pass, P, G, SP, U, search_complete):
        self.timestamps = pl.array(timestamps)
        self.Q_s = Q_s[1:] / 1e6
        self.Qdot_re = Qdot_re[1:] / 1e6
        self.Qdot_se = Qdot_se[1:] / 1e6
        self.Qdot_rs = Qdot_rs[1:] / 1e6
        self.Qdot_pass = Qdot_pass[1:] / 1e6
        self.P = P[1:] / 1e6
        self.G = G
        self.SP = SP
        self.U = U[1:]
        self.search_complete = search_complete
        self.revenue = sum(self.P * self.SP) * 0.5          #0.5 value should reflect time-steps in hours
    
    def save(self, filename = 'results.dat'):
        with open(filename, 'w') as f:
            pickle.dump(self, f)
         
    def save_raw(self, filename = 'results_raw.dat'):
        results_list = []
        results_list.append(pl.arange(0,168,0.5))
        results_list.append(self.Q_s)
        results_list.append(self.Qdot_re)
        results_list.append(self.Qdot_se)
        results_list.append(self.Qdot_rs)
        results_list.append(self.Qdot_pass)
        results_list.append(self.P)
        results_list.append(self.G)
        results_list.append(self.SP)
        results_list.append(self.U)
        results_list.append(self.revenue)        
        
        with open(filename, 'w') as f:
            pickle.dump(results_list, f)
            
         
    def bar_plot(self, ax = None):
        
        c1 = 'r'
        c2 = 'g'
        c3 = 'b'
        
        if ax == None:         
            fig = pl.figure(figsize = (160 / 25.4, 120 / 25.4))
            ax = fig.add_subplot(111)

        self.timesteps = pl.arange(0,168,0.5)       ###temporary solution

        bar1 = ax.bar(self.timesteps, self.Qdot_re, 0.5, 0, color = c1, linewidth = 0)    
        bar2 = ax.bar(self.timesteps, self.Qdot_se, 0.5, self.Qdot_re, color = c2, linewidth = 0)    
        bar3 = ax.bar(self.timesteps, self.Qdot_rs, 0.5, self.Qdot_re + self.Qdot_se, color = c3, linewidth = 0)     
        

        ax.grid(True)

        ax.set_xticks(pl.linspace(0,168,8))    
    
        ax.axis([0,168,0,220])

        ax.set_ylabel('$\dot{Q}$ (MW$_\mathrm{th}$)')    

if __name__ == "__main__":
    siminput = SimInput(start_date = '20130901', end_date = '20130901', robust = True, model = 1)

    sim = Simulation(siminput)
    sim.run_optimization()
    