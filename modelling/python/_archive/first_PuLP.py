#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 19:25:47 2020

@author: jeff
"""

#First scheduling problem solved using PuLP

from package.projdirs import datadir
import numpy as np
import pulp as plp

prob = plp.LpProblem("test1", plp.LpMinimize)

# indices
tsamp = 
                                 


# Objective
prob += x + 4*y + 9*z

# Constraints
prob += x+y <= 5
prob += x+z >= 10
prob += -y+z == 7

COIN().solve(prob)

# Solution
for v in prob.variables():
    print(v.name, "=", v.varValue)

print("objective =", value(prob.objective))
