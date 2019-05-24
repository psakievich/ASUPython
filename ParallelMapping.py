#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 21:14:30 2019

@author: psakiev
"""

import MrRealVtk as mrv
import MrImaginaryVtk as miv
from vtk import vtkPoints, vtkStructuredGrid
from vtk.util import numpy_support
from vtk.numpy_interface import dataset_adapter as dsa
import Quadratures as Quad
import numpy as np
from operator import mul
from mpi4py import MPI

def Local2Global(indLocal, writingRank, gBlock, lBlock):
    kG = indLocal[2] + (writingRank//gBlock[1])*lBlock[2]
    jG = indLocal[1] + (writingRank%gBlock[1]) *lBlock[1]
    iG = indLocal[0]
    return (iG, jG, kG)

def ExtractRangeFromMean(mean, minPoint, maxPoint):
    mean_math = dsa.WrapDataObject(mean.data)
    vel = mean_math.GetPointData().GetArray(0)
    #pre = mean_math.GetPointData().GetArray(1)
    tem = mean_math.GetPointData().GetArray(2)
    shape = [d for d in mean.data.GetDimensions()[::-1]]
    shapeV = [d for d in mean.data.GetDimensions()[::-1]]
    shapeV.append(3)
    vel = vel.reshape(shapeV)
    #pre = pre.reshape(shape)
    tem = tem.reshape(shape)
    v=vel[minPoint[2]:maxPoint[2],
               minPoint[1]:maxPoint[1],
               minPoint[0]:maxPoint[0],
               :]
    #p=pre[minPoint[2]:maxPoint[2],
    #           minPoint[1]:maxPoint[1],
    #           minPoint[0]:maxPoint[0]]
    t=tem[minPoint[2]:maxPoint[2],
               minPoint[1]:maxPoint[1],
               minPoint[0]:maxPoint[0]]
    return v,t
               
                   
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
                tG, rG, zG = Local2Global((t,r,z), writeRank, gBlock, nSize)
                x=R[rG]*np.cos(T[tG])
                y=R[rG]*np.sin(T[tG])
                points.InsertNextPoint(x,y,Z[zG])
    
    grid=vtkStructuredGrid()
    grid.SetDimensions(nT,nR,nZ)
    grid.SetPoints(points)
    
    return grid

def CreateMrReal(nSize, writeRank, gBlock, np_arrays):
    grid = CreateStructuredGrid(nSize, writeRank, gBlock)
    
    vel3d, tem3d = (numpy_support.numpy_to_vtk(np_arrays[0], 1, 10),
                    numpy_support.numpy_to_vtk(np_arrays[1], 1, 10))
    
    vel3d.SetName('velocity')
    tem3d.SetName('temperature')
    #pre3d.SetName('pressure')
    
    grid.GetPointData().SetVectors(vel3d)
    #grid.GetPointData().AddArray(pre3d)
    grid.GetPointData().SetScalars(tem3d)
    return mrv.MrVtkVector(grid)

def CreateMrImaginary(nSize, writeRank, gBlock, np_arrays):
    grid = CreateStructuredGrid(nSize, writeRank, gBlock)
    
    cvel3d, cpre3d, ctem3d = (numpy_support.numpy_to_vtk(np_arrays[0], 1, 10),
                           numpy_support.numpy_to_vtk(np_arrays[1], 1, 10),
                           numpy_support.numpy_to_vtk(np_arrays[2], 1, 10))
    
    cvel3d.SetName('C_velocity')
    ctem3d.SetName('C_temperature')
    cpre3d.SetName('C_pressure')
    
    rvel3d, rpre3d, rtem3d = (numpy_support.numpy_to_vtk(np_arrays[0], 1, 10),
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

def CreateLocalMeanMr(mean, inst, writeRank, gBlock):
    # setup sizes
    InsSize = inst.data.GetDimensions()
    InsTotal = reduce(mul, InsSize)
    corner0 = Local2Global((0,0,0),writeRank, gBlock, InsSize)
    corner1 = Local2Global(InsSize,writeRank, gBlock, InsSize)
    corner1 = (1,corner1[1],corner1[2])
    
    arrays = ExtractRangeFromMean(mean,corner0,corner1)
    r,c =InsTotal, InsTotal/InsSize[0]
    copyArray = np.zeros([r,c])
    for i in range(c):
        copyArray[i*InsSize[0]:(i+1)*InsSize[0],i]=1.0
    
    vmL=np.array([np.matmul(copyArray,arrays[0][:,:,:,0].flatten()),
                       np.matmul(copyArray,arrays[0][:,:,:,1].flatten()),
                       np.matmul(copyArray,arrays[0][:,:,:,2].flatten()),
                       ]).T
    #pmL=np.matmul(copyArray,arrays[1].flatten())
    tmL=np.matmul(copyArray,arrays[1].flatten())

    return CreateMrReal(InsSize,writeRank,gBlock,(vmL,tmL)),corner0,corner1

def AzimuthalAverage(snap, writeRank, gBlock):
    localDim = [d for d in snap.data.GetDimensions()[::-1]]
    localDimV = [d for d in snap.data.GetDimensions()[::-1]]
    localDimV.append(3)
    snap_math = dsa.WrapDataObject(snap.data)
    v,t = (snap_math.PointData[0].reshape(localDimV),
             snap_math.PointData[1].reshape(localDim))
    
    vA,tA = (np.mean(v,axis=2), np.mean(t,axis=2))
    #vA = np.array([vA[:,:,0].flatten(),vA[:,:,1].flatten(),vA[:,:,2].flatten()])
    tA = tA.flatten()
    vA=vA.reshape([len(tA),3])
    newSize = localDim[::-1]
    newSize[0] = 1
    return CreateMrReal(newSize, writeRank, gBlock, (vA, tA))

def AddSubset(totalSet, subSet, corner0, corner1):
    c0=corner0[::-1]
    c1=corner1[::-1]
    tsDim = totalSet.data.GetDimensions()[::-1]
    tsDimV = [d for d in totalSet.data.GetDimensions()[::-1]]
    tsDimV.append(3)
    ts_math = dsa.WrapDataObject(totalSet.data)
    vT,pT,tT = (ts_math.PointData[0].reshape(tsDimV),
                ts_math.PointData[1].reshape(tsDim),
                ts_math.PointData[2].reshape(tsDim))
    
    
    ssDim = subSet.data.GetDimensions()[::-1]
    ssDimV = [d for d in subSet.data.GetDimensions()[::-1]]
    ssDimV.append(3)
    ss_math = dsa.WrapDataObject(subSet.data)
    vS,tS = (ss_math.PointData[0].reshape(ssDimV),
                ss_math.PointData[1].reshape(ssDim))
    
    vT[c0[0]:c1[0],c0[1]:c1[1],c0[2]:c1[2]] += vS
    #pT[c0[0]:c1[0],c0[1]:c1[1],c0[2]:c1[2]] += pS
    tT[c0[0]:c1[0],c0[1]:c1[1],c0[2]:c1[2]] += tS
    
    totalPoints = totalSet.data.GetNumberOfPoints()
    
    pT=pT.flatten()
    tT=tT.flatten()
    vT=vT.reshape([len(pT),3])
    
    for i in range(totalPoints):
        totalSet.data.GetPointData().GetArray(3).SetTuple3(i,vT[i,0],vT[i,1],vT[i,2])
        #totalSet.data.GetPointData().GetArray(4).SetTuple1(i,pT[i])
        totalSet.data.GetPointData().GetArray(5).SetTuple1(i,tT[i])
    
    
def ReduceGridToRank0(grid,comm):
    math_grid=dsa.WrapDataObject(grid.data)
    pd = grid.data.GetPointData()
    if comm.Get_size() == 1:
        return
    if comm.Get_rank() == 0:
        for i in range(grid.data.GetPointData().GetNumberOfArrays()):
            ar = numpy_support.vtk_to_numpy(pd.GetArray(i))
            b = ar.copy()
            comm.Reduce(b,ar, op = MPI.SUM, root=0)
    else:
        for i in range(grid.data.GetPointData().GetNumberOfArrays()):
            ar = numpy_support.vtk_to_numpy(pd.GetArray(i))
            b = ar.copy()
            comm.Reduce(ar,b, op=MPI.SUM, root=0)
        
        
def Cart2Cyl(points, velocity):
    # assume these are vtkArrays
    totalPoints = points.GetNumberOfPoints()
    if(totalPoints != velocity.GetNumberOfTuples()):
        raise Exception("Cart2Cyl::Points and velocity dimensions don't match.")
        
    for i in range(totalPoints):
        x,y,z = points.GetPoint(i)
        u,v,w = velocity.GetTuple3(i)
        
        r = np.sqrt(x**2+y**2)
        theta = np.angle(x+y*1.0j)
        
        vr = u*np.cos(theta)+v*np.sin(theta)
        vtheta = -u*np.sin(theta)+v*np.cos(theta)
        
        points.SetPoint(i,r,theta,z)
        velocity.SetTuple3(i,vr,vtheta,w)
    

def Cyl2Cart(points, velocity):
    # assume these are vtkArrays
    totalPoints = points.GetNumberOfPoints()
    if(totalPoints != velocity.GetNumberOfTuples()):
        raise Exception("Cart2Cyl::Points and velocity dimensions don't match.")
        
    for i in range(totalPoints):
        r,theta,z = points.GetPoint(i)
        vr,vtheta,w = velocity.GetTuple3(i)
        
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        
        u = vr*np.cos(theta)-vtheta*np.sin(theta)
        v = vr*np.sin(theta)+vtheta*np.cos(theta)
        
        points.SetPoint(i,x,y,z)
        velocity.SetTuple3(i,u,v,w)