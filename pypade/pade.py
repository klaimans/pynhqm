#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def schlessinger_point_method(x,y):
    '''
    schlessinger_point_method
    -------------------------
    Author = Shachar Klaiman
    Date   = 24/05/2017
    
    L. Schlessinger "Use of Analyticity in the Calculation of Nonrelativistic Scattering Aplitudes",
    Phys. Rev. 167 (5), 1411 (1968)
    
    This function implements Eq 2.10 from the above paper (see also the appendix).
    
    input : 
    --------
    x - vector with the points (abscissas) on which we have the function values.
    y - vector with the function values at the points in x
    '''

    N = len(y)
    a = []
    a.append((y[0]/y[1]-1)/(x[1]-x[0]))
    # Loop to construct the continue fraction 
    for l in range(1,N-1):
        CF = a[0] * (x[l+1] - x[0]) / (1 - (y[0] / y[l+1]))
        for lp in range(1,l):
            CF = a[lp]*(x[l+1] - x[lp])/(1 + CF)
        CF += 1
        a.append(CF / (x[l] - x[l+1]))
        
    return a
    
def padepoint(x,points,values):
    '''
    Pade approximation using the Schlessinger Point method
    ------------------------------------------------------
    Author = Shachar Klaiman
    Date   = 24/05/2017
    '''
    a = schlessinger_point_method(points,values)
    N = len(a) + 1
    CF = a[N-2] * (x - points[N-2])
    for l in range(N-3,-1,-1):
        CF = a[l] * (x - points[l]) / (1 + CF)
    F = values[0] / (1 + CF)
    return F

if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    points = [1,4,9,16]
    values = [1,2,3,4]
    fitted = padepoint(np.asarray([1,4,9,16,25]),points,values)
    plt.plot(points,values,'.-')
    plt.plot([1,4,9,16,25],fitted,'o')