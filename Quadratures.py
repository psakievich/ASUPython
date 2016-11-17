#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 19:23:32 2016

Quadratures for numerical integration

@author: psakievich
"""
import numpy as np
import scipy.special as ss

class GaussLobattoChebyshev():
    def Point(self,n,i):
        return -np.cos(np.pi*i/(n-1.0))
    def Weight(self,n,i):
        n1=int(n)
        i1=int(i)
        if (i1==0 or i1==n1-1):
            return np.pi/(n-1)*0.5*np.sqrt(1.0-self.Point(n,i)**2)
        else:
            return np.pi/(n-1)*np.sqrt(1.0-self.Point(n,i)**2)
    def Points(self,n):
        a=np.empty(int(n))
        for i in range(int(n)):
            a[i]=self.Point(n,i)
        return a
    def Weights(self,n):
        a=np.empty(int(n))
        for i in range(int(n)):
            a[i]=self.Weight(n,i)
        return a
class GaussChebyshev():
    def Point(self,n,i):
        return -np.cos((2.0*float(i)+1.0)*np.pi/(2.0*float(n)))
    def Weight(self,n,i):
        p=self.Point(int(n),i)
        w=np.pi/float(n)*np.sqrt(1.0-p**2)
        return w
    def Points(self,n):
        a=np.empty(int(n))
        for i in range(int(n)):
            a[i]=self.Point(n,i)
        return a
    def Weights(self,n):
        w=np.empty(int(n))
        for i in range(int(n)):
            w[i]=self.Weight(n,i)
        return w
class GaussLegendre():
    def Points(self,n):
        p,w=ss.p_roots(n)
        return p
    def Weights(self,n):
        p,w=ss.p_roots(n)
        return w
    def Point(self,n,i):
        return self.Points(n)[i]
    def Weight(self,n,i):
        return self.Weights(n)[i]
class GaussLobattoLegendre():
    def Points(self,n):
        p=np.empty(n)
        w=np.empty(n)
        p[1:n-1],w[1:n-1]=ss.j_roots(n-2,1.0,1.0)
        p[0]=-1.0
        p[n-1]=1.0
        return p
    def Weights(self,n):
        w=np.empty(n)
        p=self.Points(n)
        for i in range(n):
            w[i]=2.0/(n*(n-1)*ss.eval_legendre(n-1,p[i])**2)
        return w
    def Point(self,n,i):
        return self.Points(n)[i]
    def Weight(self,n,i):
        return self.Weights(n)[i]
    
def Converge(QuadratureClass,start=5,stop=50,inc=5):
    error=[]
    pnts=[]
    for i in range(start,stop,inc):
        pnts.append(i)
        exact=2.0/3.0 #x^2
        test=np.sum(QuadratureClass.Points(i)**2*QuadratureClass.Weights(i))
        error.append(abs(exact-test))
    return pnts,error
    