#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 11:22:00 2020

@author: jeff
"""

#This example was cut and pasted from https://pypi.org/project/PuLP/1.1/
#It is an example of the PuLP linear programming package for python and I
#well test it to see whether it capable of solving the same problems that 
#we have been using MiniZinc to solve. It would be nice to not have to call
#MiniZinc but rather just be able to call the solvers from within Python.

from pulp import *

prob = LpProblem("test1", LpMinimize)

# Variables
x = LpVariable("x", 0, 4)
y = LpVariable("y", -1, 1)
z = LpVariable("z", 0)

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