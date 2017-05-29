#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 29 13:41:20 2017

@author: Shachar Klaiman
"""
import numpy as np

class fbr:
    def __init__(self,num_basis,basis_type = "Sine",**kwargs):
        self.basis_type = basis_type
        self.nb = num_basis
        try:
            self.ndvr = kwargs["num_dvr"]
        except:
            self.ndvr = self.nb
        if basis_type == "Sine":
            self.L = kwargs["L"]
            try:
                self.space = kwargs["space"]
            except:
                self.space = 1
            try:
                self.x0 = kwargs["x0"]
            except:
                if self.space == 1:
                    self.x0 = 0
                elif self.space == 3:
                    self.x0 = self.L / 2
                else:
                    raise ValueError("Space argument can only be 1 or 3")
       
            self.basisfun = lambda n,x: np.sqrt(2.0/self.L) \
                            * np.sin(n*np.pi*(x + self.x0)/self.L)
            self.dx2 = np.diag(np.arange(1,self.nb+1)**2
                                   *np.pi**2/(self.L**2))
            self.Tm, self.Tev = self.diag_x()
        else:
            raise ValueError("Unknown basis type")
    def diag_x(self):
        if self.basis_type == "Sine":
            XMAT=np.zeros([self.ndvr,self.ndvr])
            for s in range(0,self.ndvr):
                XMAT[s,s] = self.x0
                for sp in range(s+1,self.ndvr):
                    XMAT[s,sp] = (self.L/np.pi**2)*(
                            (np.cos(np.pi*(sp-s))-1)/(sp-s)**2
                            - (np.cos(np.pi*(sp+s+2))-1)/(sp+s+2)**2)
                    XMAT[sp,s] = XMAT[s,sp]
        else:
            raise ValueError("Unknown basis type")
        Tev, Tm = np.linalg.eigh(XMAT)
        return Tm, Tev

class eig:
    def __init__(self,num_basis,fVx,basis_type = "Sine",hbar=1,mass=1,**kwargs):
        fbr
        self.nb = num_basis
        try:
            self.ndvr = kwargs["num_dvr"]
        except:
            self.ndvr = self.nb
        self.fVx = fVx
        self.hbar = hbar
        self.mass = mass
        self.basis_type = basis_type
        self.fbr = fbr(self.nb,basis_type="Sine",**kwargs)
        self.kinetic_factor = self.hbar**2/(2*self.mass)    
        self.Vx = self.fVx(self.fbr.Tev)
        self.KM = self.kinetic_factor * self.fbr.dx2
        self.VM = np.dot(self.fbr.Tm,np.dot(np.diag(self.Vx),self.fbr.Tm.transpose()))
        self.HM = self.KM + self.VM
        self.ev, self.EFC = np.linalg.eigh(self.HM)
        
        
if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    """
    fbr = fbr(200,basis_type="Sine",L=30,space=1)
    fVx = lambda x: (x**2/2 - 0.8)*np.exp(-0.1*x**2)
    Vx = fVx(fbr.Tev)
    hbar = 1
    mass = 1
    kinetic_factor = hbar**2/(2*mass)    
    KM = kinetic_factor * fbr.dx2
    VM = np.dot(fbr.Tm,np.dot(np.diag(Vx),fbr.Tm.transpose()))
    HM = KM + VM
    ev, EFC = np.linalg.eigh(HM)
    plt.plot(fbr.Tev,Vx)
    plt.plot([-3,3],[ev,ev])
    plt.ylim([np.min(Vx),np.max(Vx)+1])
    plt.show()
    """
    fVx = lambda x: (x**2/2 - 0.8)*np.exp(-0.1*x**2)
    eigen = eig(200,fVx,basis_type="Sine",L=30,space=1)
    ev = eigen.ev
    EFC = eigen.EFC
    Vx = eigen.Vx
    x = eigen.fbr.Tev
    plt.plot(x,Vx)
    plt.plot([-3,3],[ev,ev])
    plt.ylim([np.min(Vx),np.max(Vx)+1])
    plt.show()
    
    
