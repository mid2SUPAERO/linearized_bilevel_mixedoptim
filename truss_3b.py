# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding:utf-8 -*-
"""
Created on Thu Jan 24 2019
@author: Pierre-Jean Barjhoux
"""
from numpy import array, dot, sum, empty, append
from math import sqrt, pi
from categorical_data import build_properties, build_gamma
from numpy.linalg import solve

# lengths
ul = 1000. # unitary length (mm)
l = array([sqrt(2)*ul, ul, sqrt(2)*ul])
gamma = build_gamma()

def compute_weight(a, c):
    '''
    Computes weight
    :param a: areas
    :param c: categorical variables
    :return: weight
    '''
    prop = build_properties(c, gamma)
    rho = prop['densities']
    return sum(rho*a*l)

def compute_glob_stiffness(a, c):
    '''
    Computes global stiffness matrix K where clampled nodes are removed
    (partial K)
    :param a: areas
    :param c: categorical variables
    :return: global stiffness matrix
    '''
    prop = build_properties(c, gamma)
    E = prop['young_moduli']
    K = empty((2,2))
    K[0,0] = E[0]*a[0]/(2*l[0]) + E[2]*a[2]/(2*l[2])
    K[0,1] = E[0]*a[0]/(2*l[0]) - E[2]*a[2]/(2*l[2])
    K[1,0] = E[0]*a[0]/(2*l[0]) - E[2]*a[2]/(2*l[2])
    K[1,1] = E[0]*a[0]/(2*l[0]) + E[1]*a[1]/l[1] + E[2]*a[2]/(2*l[2])
    return K

def compute_stress_3btruss_analytical(a, c, f):
    '''
    Computes stresses in elements, and displacements of the free node
    :param a: areas
    :param c: categorical variables
    :param f: external loads
    :return: displacements of free node, stress on the 3 bars
    '''
    # computes displaceemnts of the free node
    K = compute_glob_stiffness(a,c)
    u_freenode = solve(K, f)
    u = append(u_freenode, array([0,0]))

    # computes stresses
    prop = build_properties(c, gamma)
    E = prop['young_moduli']
    stress = empty(3)
    k1e = array([[-sqrt(2)/2, -sqrt(2)/2, 0., 0.],[0., 0., -sqrt(2)/2, -sqrt(2)/2]])
    stress[0] = (E[0]/l[0])*dot(dot(array([-1,1]),k1e),u)
    k2e = array([[0., -1., 0., 0.],[0., 0., 0, -1.]])
    stress[1] = (E[1]/l[1])*dot(dot(array([-1,1]),k2e),u)
    k3e = array([[sqrt(2)/2, -sqrt(2)/2, 0., 0.],[0., 0., sqrt(2)/2, -sqrt(2)/2]])
    stress[2] = (E[2]/l[2])*dot(dot(array([-1,1]),k3e),u)
    return u_freenode, stress

def bckl_local_I(E_i, nu_i, geom_i):
    '''
    Computes local buckling for I-stiffener
    :param E_i: young modulus
    :param nu_i: poisson coefficient
    :param geom_i: cross section of current stiffener
    :return: local buckling limit stress
    '''
    bckl_local = 4. * pi**2 * E_i / \
                (12. * (1. - nu_i**2)) * (geom_i[0] / geom_i[1])**2
    return bckl_local

def bckl_euler_I(E_i, a_i, geom_i, l_i):
    '''
    Computes Euler buckling for I-stiffener
    :param E_i: young modulus
    :param a_i: area
    :param geom_i: cross section of current stiffener
    :param l_i: length
    :return: euler buckling limit stress
    '''
    # quadratic moment of reference cross-section
    quad_moment_0 = (geom_i[2] * (geom_i[1] + 2. * geom_i[0]) ** 3 / 12.) - \
          (geom_i[2] - geom_i[0]) * geom_i[1] ** 3 / 12.
    aref =  geom_i[0] * geom_i[1] + 2 * geom_i[0] * geom_i[2]
    # scaling of reference quadratic moment
    quad_moment = quad_moment_0 * (a_i/aref)**2
    # Euler buckling computation
    bckl_euler = pi ** 2 * E_i * quad_moment / \
                 (a_i * l_i ** 2)
    return bckl_euler