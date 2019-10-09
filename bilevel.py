# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding:utf-8 -*-
"""
Created on Thu Jan 24 2019
@author: Pierre-Jean Barjhoux
"""

from lower_pb import build_psi
from numpy import empty, argmin, all, copy
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

def termination_criterion(w_hist, eps):
    '''
    Proposed termination criterion :
    stops when stationarity on weight is detected
    -> could be replaced by another criterion
    :param w_hist: optimal weight history
    :return: bool
    '''
    if len(w_hist) == 1:
        # if first iteration, continue
        return True
    else:
        if abs(w_hist[-1] - w_hist[-2]) > eps:
            # if no stationarity on weight, continue
            return True
        else:
            # otherwise stop
            return False

def display(w, c, a, iteration = None):
    '''
    Display key values
    '''
    iter_msg = ""
    if iteration is not None:
        iter_msg = "Iteration (%i) | " % iteration
    logger.info(iter_msg + "w* = %.2f (kg)"%w)
    logger.info(iter_msg + "c* = %s"%str(c))
    logger.info(iter_msg + "a* = %s (mm) \n"%str(a))

def bilevel_solver(c0, a0, f, opt, l_b, u_b, disp_bound, gamma):
    '''
    Bi-level Solver
    :param c0: initial categorical variables (parameters)
    :param a0: initial areas (design variables)
    :param f: external loads
    :param opt: scipy solver options
    :param l_b: lower bounds on areas
    :param u_b: upper bounds on areas
    :param disp_bound: bound on displacement
    :param gamma: categorical database
    :return: optimal areas, choices, weight
    '''
    # build lower level optimization problem
    psi = build_psi(f, a0, opt, l_b, u_b, disp_bound)

    # initialization
    w_hist = [psi(c0)[1]] # optimal weight history
    c_hist = [c0]      # categorical variables history
    W = empty((3,3))   # temporary matrix of optimal weights per iter
    eps = 1e-3         # epsilon for termination criteria
    k = 1              # iteration number
    while termination_criterion(w_hist, eps):
        logger.info("Starts iteration %i"%k)
        # enumeration of elements
        for elt in xrange(3):
            c = copy(c_hist[-1])
            # enumeration of choices
            for choice in gamma:
                c[elt] = choice
                if all(c == c_hist[-1]):
                    # avoid optimization already performed at previous iteration
                    W[elt, choice-1] = w_hist[-1]
                else:
                    # lower level optimization
                    W[elt, choice-1] = psi(c)[1]
       	logger.info("Iteration (k) | \n W = "+str(W))
        # for each element, take the choice that lead to optimal weight
        c = argmin(W, axis=1) + 1 # warning : python index starts by 1
        # lower level optimization at new categorical solution
        a, w = psi(c)
        # decrease trategy if needed
        if w > w_hist[-1]:
            raise NotImplementedError('Implement your decrease strategy.')
        # display
        display(w, c, a, k)
        # store results, increment k
        c_hist.append(c)
        w_hist.append(w)
        k += 1

    logger.info('End of Bi-level :')
    display(w, c, a)

    return a, c, w