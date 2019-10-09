# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding:utf-8 -*-
"""
Created on Thu Jan 24 2019
@author: Pierre-Jean Barjhoux
"""
from scipy.optimize import minimize
from numpy import empty, array, append, concatenate
from categorical_data import build_properties
from truss_3b import compute_stress_3btruss_analytical, compute_weight, bckl_euler_I, bckl_local_I, l
from categorical_data import build_gamma

# retrieve categorical data
gamma = build_gamma()

def stress_constr(stress, a, c):
    '''
    Computes stress constraints
    :param stress:
    :param a: areas
    :param c: categorical variables
    :return: a matrix (n,m) of constraints
             with n=3, m=4 (tension, compression, euler and local buckling)
    '''
    # retrieve material + stiffeners definitions
    prop = build_properties(c, gamma)
    stiffeners = prop['stiffeners']
    allow_comp = prop['allow_comp']
    allow_tens = prop['allow_tens']
    young_moduli = prop['young_moduli']
    poisson_coeffs = prop['poisson_coeffs']

    # init stress matrix s
    s = empty((3,4))
    # constraint on tension allowable
    s[:,0] = stress - allow_tens
    # constraint on compression allowable
    s[:,1] = - stress - allow_comp
    # euler and local buckling constraints
    bckl_euler = array([])
    bckl_local = array([])
    for i, ((type_i, geom_i), a_i) in enumerate(zip(stiffeners,a)): # for each element
        E_i = young_moduli[i]
        nu_i = poisson_coeffs[i]
        a_i = a[i]
        l_i = l[i]
        # specific formulation of buckling wrt stiffener definition
        if type_i == 'I':
            bckl_euler_i = bckl_euler_I(E_i, a_i, geom_i, l_i)
            bckl_local_i = bckl_local_I(E_i, nu_i, geom_i)
        elif type_i == 'T':
            msg = 'No buckling limit stress implementation for stifferner %s'%type_i
            raise NotImplementedError(msg)
        else:
            raise NotImplementedError
        bckl_euler = append(bckl_euler, bckl_euler_i)
        bckl_local = append(bckl_local, bckl_local_i)
    # euler buckling constraints
    s[:,2] = - stress - bckl_euler
    # local buckling constraints
    s[:,3] = - stress - bckl_local
    # non negative constraints for scipy SLSQP
    return -s

def displactements_constr(u, disp_bound):
    '''
    Computes the constraint on displacements
    :param u: displacements
    :param disp_bound: bound on displacements
    :return: difference between vertical displacements on the free node
             and a user defined bound value
    '''
    dof_bound = 1 # python index : bound on 'y' axis, 3rd node
    delta = array([disp_bound - u[dof_bound]])
    return delta

def build_constraints(c, f, disp_bound):
    '''
    Builds the method that evaluates all the constraints
    at given loads, categorical variables, and bound on displacements
    :param c: categorical variables
    :param f: external loads
    :param disp_bound: bound on displacements
    :return: python method that computes constraints
    '''
    def cons(a):
        '''
        Method that evaluates all the constraints in the optimizer.
        :param a:
        :return: stress and displacements constraints (as vector)
        '''
        u, stress = compute_stress_3btruss_analytical(a, c, f)
        s = stress_constr(stress, a, c)
        delta = displactements_constr(u, disp_bound)
        return concatenate((s.flatten(), delta))
    return cons

def build_objective(c):
    '''
    Builds optimization objective
    :param c: categorical variables
    :return: python method that computes weight
    '''
    def obj(a):
        '''
        Compute objective function
        :param a: areas
        :return: weight (scalar)
        '''
        return compute_weight(a, c)
    return obj

def build_psi(f, a0, opt, l_b, u_b, disp_bound):
    '''
    Build lower optimization problem
    :param f: external loads
    :param a0: initial areas
    :param opt: solver options
    :param l_b: lower bounds
    :param u_b: upper bounds
    :param disp_bound: bounds on displacements
    :return: python method that solves the lower problem, as Psi
    '''
    # formatting bounds for SLSQP
    bnds = zip(l_b, u_b)

    # build lower problem, with categorical variables as parameters
    def psi(c):
        # build constraints evaluation function (displacements and stress)
        cons = build_constraints(c, f, disp_bound)
        # build objective function
        obj = build_objective(c)
        # call to scipy solver
        slsqp_cons = {'type' : 'ineq',
                      'fun'  : cons}
        res = minimize(obj, a0, tol = 1e-8, method='SLSQP', bounds=bnds,
                       constraints = slsqp_cons, options=opt)
        return res.x, res.fun

    return psi

