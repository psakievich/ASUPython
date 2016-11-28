#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 15:52:17 2016

Gather spacially integrated Fourier coefficients

@author: psakievich
"""

import MrImaginaryVtk as MIV
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

def GetIntegratedFouriereCoefficients(PATH,FILEBASE,WAVENUM,TRANGE,SZE=10):
    FTEMP=PATH+FILEBASE+'_{0}_TSTEP_{1}.vts'
    result=np.empty([TRANGE[1]-TRANGE[0],SZE])
    for i in range(TRANGE[0],TRANGE[1]):
        h=MIV.MrVtkVecHandle(FTEMP.format(WAVENUM,i))
        g=h.get()
        result[i,:]=g.integrated_values()
        #result.append(g.integrated_values())
    return result

PATH='analysis/{}/'
FILEBASE='SymS1'
FILERANGE=[0,227]
WAVE=3
wave3=GetIntegratedFouriereCoefficients(PATH.format(WAVE),FILEBASE,WAVE,FILERANGE)


YLIM=(np.min(wave3),np.max(wave3))
XLIM=YLIM

fig=plt.figure()
ax=plt.axes(xlim=XLIM,ylim=YLIM)
ax.grid(True)
velT,=ax.plot([],[],'-o',label=r'$v_{\theta}$')
velR,=ax.plot([],[],'-o',label='$v_r$')
velZ,=ax.plot([],[],'-o',label='$v_z$')
temp,=ax.plot([],[],'-o',label='$T$')
time_text=ax.text(0.05,0.9,'',transform=ax.transAxes)
plt.legend(handles=[velR,velT,velZ,temp])

def init():
    velR.set_data([],[])
    velT.set_data([],[])
    velZ.set_data([],[])
    temp.set_data([],[])
    time_text.set_text('')
    return velR,velT,velZ,temp,time_text
    
def animate(i):
    velR.set_data(wave3.T[1,0:i],wave3.T[0,0:i]) 
    velT.set_data(wave3.T[3,0:i],wave3.T[2,0:i]) 
    velZ.set_data(wave3.T[5,0:i],wave3.T[4,0:i]) 
    temp.set_data(wave3.T[9,0:i],wave3.T[8,0:i]) 
    time_text.set_text('SNAPSHOT={}'.format(i))
    return velR,velT,velZ,temp,time_text
    
ani=animation.FuncAnimation(fig,
                            animate,
                            init_func=init,
                            frames=FILERANGE[1], #generator, iterable or number of frames
                            interval=100, #time delay between drawings (ms)
                            repeat=0) 
plt.show()