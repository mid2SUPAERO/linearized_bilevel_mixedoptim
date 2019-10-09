# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding:utf-8 -*-
"""
Created on Thu Jan 24 2019
@author: Pierre-Jean Barjhoux
"""
from numpy import array

def build_gamma():
    '''
    Builds categorical data associated to choices in \Gamma = {1,2,3}
    -> possible to add categorical variables
    :return: a dictionary with data associated to each choice
    '''
    # init categorical space
    gamma = {}

    # catalog 1 : material AL2139 + I stiffener
    gamma[1] = {}
    gamma[1]['density'] = 2.8e-6          # kg/mmm^3
    gamma[1]['young_modulus'] = 7.1e+4    # MPa
    gamma[1]['allow_comp'] = 2.00e+2 # MPa
    gamma[1]['allow_tens'] = 1.5e+2  # MPa
    gamma[1]['poisson_coeff'] = 0.3
    gamma[1]['stiffener'] = 'I', [5.,50.0,40.0]

    # catalog 2 : material AL2024 + I stiffener
    gamma[2] = {}
    gamma[2]['density'] = 2.77e-6         # kg/mmm^3
    gamma[2]['young_modulus'] = 7.4e+4    # MPa
    gamma[2]['allow_comp'] = 2.10e+2 # MPa
    gamma[2]['allow_tens'] = 1.6e+2  # MPa
    gamma[2]['poisson_coeff'] = 0.33
    gamma[2]['stiffener'] = 'I', [5.,50.0,40.0]

    # catalog 3 : material TA6V + I stiffener
    gamma[3] = {}
    gamma[3]['density'] = 4.43e-6         # kg/mmm^3
    gamma[3]['young_modulus'] = 11.0e+4   # MPa
    gamma[3]['allow_comp'] = 8.6e+2  # MPa
    gamma[3]['allow_tens'] = 11.e+2  # MPa
    gamma[3]['poisson_coeff'] = 0.33
    gamma[3]['stiffener'] = 'I', [5.,50.0,40.0]

    return gamma

def build_properties(c, gamma):
    '''
    Builds properties given a configuration c in gamma
    :param c: categorical design variables
    :param gamma: categorical design space
    :return: dict of properties
    '''
    prop = {}
    prop['densities'] = array([gamma[c_i]['density'] for c_i in c])
    prop['young_moduli'] = array([gamma[c_i]['young_modulus'] for c_i in c])
    prop['allow_comp'] = array([gamma[c_i]['allow_comp'] for c_i in c])
    prop['allow_tens'] = array([gamma[c_i]['allow_tens'] for c_i in c])
    prop['poisson_coeffs'] = array([gamma[c_i]['poisson_coeff'] for c_i in c])
    prop['stiffeners'] = [gamma[c_i]['stiffener'] for c_i in c]
    return prop