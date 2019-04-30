#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 21:41:29 2019

@author: psakiev
"""

import unittest
import ParallelMapping as pm
import MrRealVtk as mrv
import MrImaginaryVtk as miv
from vtk import vtkDoubleArray
from vtk import vtkPoints, vtkStructuredGrid
from vtk.util import numpy_support
import Quadratures as Quad
import numpy as np

def CreateStructuredGrid(nSize, writeRank, gBlock):
    gl=Quad.GaussLegendre()
    nT,nR,nZ = nSize
    
    R,T,Z=((gl.Points(nR*gBlock[1])+1.0)*0.5*3.15,
        np.linspace(0,2*np.pi,nT*gBlock[0]),
        gl.Points(nZ*gBlock[2])*0.5)

    points=vtkPoints()
    for z in range(nZ):
        for r in range(nR):
            for t in range(nT):
                tG, rG, zG = pm.Local2Global((t,r,z), writeRank, gBlock, nSize)
                x=R[rG]*np.cos(T[tG])
                y=R[rG]*np.sin(T[tG])
                points.InsertNextPoint(x,y,Z[zG])
    
    grid=vtkStructuredGrid()
    grid.SetDimensions(nR,nT,nZ)
    grid.SetPoints(points)
    
    return grid

def CreateMrReal(nSize, writeRank, gBlock, np_arrays):
    grid = CreateStructuredGrid(nSize, writeRank, gBlock)
    
    vel3d, tem3d, pre3d = (numpy_support.numpy_to_vtk(np_arrays[0], 1, 10),
                           numpy_support.numpy_to_vtk(np_arrays[1], 1, 10),
                           numpy_support.numpy_to_vtk(np_arrays[2], 1, 10))
    
    vel3d.SetName('velocity')
    tem3d.SetName('temperature')
    pre3d.SetName('pressure')
    
    grid.GetPointData().SetVectors(vel3d)
    grid.GetPointData().AddArray(pre3d)
    grid.GetPointData().SetScalars(tem3d)
    return mrv.MrVtkVector(grid)

def CreateMrImaginary(nSize, writeRank, gBlock, np_arrays):
    grid = CreateStructuredGrid(nSize, writeRank, gBlock)
    
    cvel3d, ctem3d, cpre3d = (numpy_support.numpy_to_vtk(np_arrays[0], 1, 10),
                           numpy_support.numpy_to_vtk(np_arrays[1], 1, 10),
                           numpy_support.numpy_to_vtk(np_arrays[2], 1, 10))
    
    cvel3d.SetName('C_velocity')
    ctem3d.SetName('C_temperature')
    cpre3d.SetName('C_pressure')
    
    rvel3d, rtem3d, rpre3d = (numpy_support.numpy_to_vtk(np_arrays[0], 1, 10),
                           numpy_support.numpy_to_vtk(np_arrays[1], 1, 10),
                           numpy_support.numpy_to_vtk(np_arrays[2], 1, 10))
    
    rvel3d.SetName('R_velocity')
    rtem3d.SetName('R_temperature')
    rpre3d.SetName('R_pressure')
    
    grid.GetPointData().SetVectors(cvel3d)
    grid.GetPointData().AddArray(cpre3d)
    grid.GetPointData().SetScalars(ctem3d)
    grid.GetPointData().AddArray(rvel3d)
    grid.GetPointData().AddArray(rpre3d)
    grid.GetPointData().AddArray(rtem3d)
    return miv.MrVtkVector(grid)



class TestComputingStats(unittest.TestCase):
    def setUp(self):
        self.MeanSize = (1, 64, 1024)
        temp = np.arange(0,64*1024)
        vel = np.array([temp,temp*2.0,temp*4.0],dtype=float).T
        self.np_arrays=[vel,temp*0.5,temp*0.25]
    
    def test_ExtractRangeFromMeanShape(self):
        mean = CreateMrImaginary(self.MeanSize, 0, (1,1,1), self.np_arrays)
        v,p,t = pm.ExtractRangeFromMean(mean, (0,0,0), (8,0,1))
        self.assertEqual((8,0,1,3),v.shape)
        self.assertEqual((8,0,1),p.shape)
        self.assertEqual((8,0,1),t.shape)
        
    def test_ExtractRangeFromMeanValues(self):
        mean = CreateMrImaginary(self.MeanSize, 0, (1,1,1), self.np_arrays)
        v,p,t = pm.ExtractRangeFromMean(mean, (0,0,0), (8,1,1))
        vExact = self.np_arrays[0].reshape(64,1,1024,3)
        self.assertTrue(vExact[0:8,0:1,0:1,0].all()==v[:,:,:,0].all())
        

if __name__=='__main__':
    unittest.main()