# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding:utf-8 -*-
"""
Created on Thu Jan 24 2019
@author: Pierre-Jean Barjhoux
"""

from numpy import array
from categorical_data import build_gamma
from bilevel import bilevel_solver

# bound on vertical displacement
disp_bound = 1. # mm

# loads
Fx = 0. # N
Fy = 200e3 # N
f = array([Fx, Fy])

# build categorical design space
gamma = build_gamma()

# numerical options for lower optimization (refer to scipy SLSQP doc.)
opt = {'ftol': 1e-6,  # tolerancy on weight
       'eps': 1e-5,   # finite differences step size
       "disp": False} # True to display details on lower level optimizations

# Bi-level parameters
c0 = array([1,2,3])             # initial categorical design variables
a0 = array([1999.,1999.,1999.]) # initial guess of areas (lower level optimization)
l_b = [100.,]*3                 # lower bounds on areas
u_b = [2000.,]*3                # upper bounds on areas

# call solver
a, c, w = bilevel_solver(c0, a0, f, opt, l_b, u_b, disp_bound, gamma)

