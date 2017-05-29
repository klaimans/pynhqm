#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 29 16:39:49 2017

@author: Shachar Klaiman
"""
from pyeig.fbr import eig
from pypade.pade import padepoint
import numpy as np
import matplotlib.pyplot as plt


fVx = lambda x: (x**2/2 - 0.8)*np.exp(-0.1*x**2)
NB = 20
space = 1
basis = "Sine"
Niter = 50
Lvec = np.linspace(15,30,Niter)
evMAT = np.zeros([NB,len(Lvec)])
for i,L in enumerate(Lvec):
    print("Iteration # {}".format(i))
    eigen = eig(NB,fVx,basis_type=basis,L=L,space=space)
    evMAT[:,i] = eigen.ev

plt.plot(Lvec,evMAT.transpose())
plt.ylim([np.min(eigen.Vx),np.max(eigen.Vx)+1])
plt.show()
#plt.plot(Lvec[0:-1],np.diff(evMAT.transpose(),axis=0)[:,0:10])
#plt.ylim([-0.1,0])
#plt.show()

#eigen = eig(num_basis,fVx,basis_type="Sine",L=30,space=1)
#ev = eigen.ev
#EFC = eigen.EFC
#Vx = eigen.Vx
#x = eigen.fbr.Tev
#plt.plot(x,Vx)
#plt.plot([-3,3],[ev,ev])
#plt.ylim([np.min(Vx),np.max(Vx)+1])
#plt.show()