#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 15:48:16 2017

create vtk grid for post processing

@author: psakievich
"""

from vtk import vtkPoints, vtkStructuredGrid,vtkXMLStructuredGridWriter
import Quadratures as Quad
import os,sys
import numpy as np
gl=Quad.GaussLegendre()

nR,nT,nZ=80,512,32
nT+=1

R,T,Z=(gl.Points(nR)+1.0)*0.5*3.15,np.linspace(0,2*np.pi,nT),gl.Points(nZ)*0.5


points=vtkPoints()
for z in range(nZ):
    for t in range(nT):
        for r in range(nR):
            x=R[r]*np.cos(T[t])
            y=R[r]*np.sin(T[t])
            points.InsertNextPoint(x,y,Z[z])

grid=vtkStructuredGrid()
grid.SetDimensions(nR,nT,nZ)
grid.SetPoints(points)

writer=vtkXMLStructuredGridWriter()

writer.SetFileName(os.environ["HOME"]+r"/Desktop/RBC/JFM/Grids/FFTGrid.vts")
writer.SetInputData(grid)
writer.Write()