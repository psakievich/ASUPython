#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 19:23:32 2016

Quadratures for numerical integration

@author: psakievich
"""
import numpy as np

class ChebyshevGaussLobatto():
    def Point(self,n,i):
        return -np.cos(np.pi*i/(n-1))
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
class ChebyshevGauss():
    def Point(self,n,i):
        return -np.cos((2.0*(i+1)-1.0)*np.pi/(2.0*n))
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