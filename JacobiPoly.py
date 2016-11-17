#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 17:16:29 2016

All things Jacobi Polynomials 

ROUTINES TRANSLATED INTO PYTHON FROM speclib.f by Einar Malvin Ronquist
Modifications made based on Karniadakis, George EM. and Sherwin, Spencer,
"Spectral/hp Element Methods for Computational Fluid Dynamics", Second Edition,
pg 586 

@author: psakievich
"""
import numpy as np
import Quadratures as QD

#VERIFIED
def JacobiPoly(x,N,alpha,beta):
    '''Jacobi Polynomial
    Returns the jacobi polynomial and it's derivative of degree N at x
    alpha and beta must be greater than -1
    '''
    #P=0
    P=1.0 #Poly val at x
    dP=0.0# der Poly val at X
    if (N==0):
        return P,dP
    #P=1
    PL=P  #previous poly val at x
    dPL=dP
    apb=alpha+beta
    P=0.5*(alpha-beta+(apb+2.0)*x)
    dP=0.5*(apb+2.0)
    if (N==1):
        return P,dP
    for i in range(2,N+1):
        di=float(i)
        a1=2.0*(di)*(di+apb)*(2.0*di+apb-2.0)
        a2=(2.0*di-1.0+apb)*(alpha**2-beta**2)
        b3=(2.0*di+apb-2.0)
        a3=b3*(b3+1.0)*(b3+2.0)
        a4=2.0*(di+alpha-1.0)*(di+beta-1.0)*(2.0*di+apb)
        PN =((a2+a3*x)*P -a4*PL)/a1
        dPN=((a2+a3*x)*dP-a4*dPL+a3*P)/a1
        PL=P
        dPL=dP
        P=PN
        dP=dPN
    return P, dP
    
def GaussPoints(NP,alpha,beta,eps=1e-12,nMaxIter=10):
    N=NP-1
    XLAST=0.0
    CHG=np.pi*(2.0*N+2.0) #chebyshev seed
    GP=QD.ChebyshevGauss().Points(NP)#np.zeros(NP)
    for j in range(NP):
        if(j==0):
            X=np.cos((2.0*j+1.0)*CHG)
        else:
            X1=np.cos((2.0*j+1.0)*CHG)
            X2=XLAST
            X=0.5*(X1+X2)
        for k in range(nMaxIter):
            P,PD=JacobiPoly(X,NP,alpha,beta)
            RECSUM=0.0
            for i in range(j-1):
                RECSUM=RECSUM+1.0/(X-GP[N-i])
            DELX=-P/(PD-RECSUM*P)
            X=X+DELX
            if(abs(DELX)<eps):
                break
        GP[N-j]=X
        XLAST=X
    return GP
        
        
    
