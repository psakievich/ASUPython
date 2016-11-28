#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 15:37:52 2016

Process Energy Report
Plot energy report from each step in time series

@author: psakievich
"""
import FileRepair as FR 
from matplotlib import pyplot as plt
from matplotlib import animation

PATH='analysis/'
BASENAME='EnergyReport_'
FILEEXT='.dat'
PLOTMODES=[5,50] #[min,max]
FILERANGE=[0,227] #[min,max]
YLIM=(0,1e-3)

fig=plt.figure()
ax=plt.axes(xlim=(PLOTMODES[0]-1,PLOTMODES[1]+1),ylim=YLIM)
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
    txt=FR.ReadBrokenFile(PATH+BASENAME+'{}'.format(i)+FILEEXT,3)
    result=FR.PackBrokenFileIntoNumpy(txt,8)
    velR.set_data(result.T[0,PLOTMODES[0]:PLOTMODES[1]], \
                  result.T[1,PLOTMODES[0]:PLOTMODES[1]])
    velT.set_data(result.T[0,PLOTMODES[0]:PLOTMODES[1]], \
                  result.T[2,PLOTMODES[0]:PLOTMODES[1]])
    velZ.set_data(result.T[0,PLOTMODES[0]:PLOTMODES[1]], \
                  result.T[3,PLOTMODES[0]:PLOTMODES[1]])
    temp.set_data(result.T[0,PLOTMODES[0]:PLOTMODES[1]], \
                  result.T[5,PLOTMODES[0]:PLOTMODES[1]])
    time_text.set_text('SNAPSHOT={}'.format(i))
    return velR,velT,velZ,temp,time_text
    
ani=animation.FuncAnimation(fig,
                            animate,
                            init_func=init,
                            frames=FILERANGE[1], #generator, iterable or number of frames
                            interval=100, #time delay between drawings (ms)
                            repeat=0) 
plt.show()

