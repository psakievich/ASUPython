#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 21:14:30 2019

@author: psakiev
"""

from vtk.numpy_interface import dataset_adapter as dsa

def Local2Global(indLocal, writingRank, gBlock, lBlock):
    kG = indLocal[2] + (writingRank//gBlock[2])*lBlock[2]
    jG = indLocal[1] + (writingRank%gBlock[1]) *lBlock[1]
    iG = indLocal[0]
    return (iG, jG, kG)

def ExtractRangeFromMean(mean, minPoint, maxPoint):
    mean_math = dsa.WrapDataObject(mean.data)
    vel = mean_math.GetPointData().GetArray(3)
    pre = mean_math.GetPointData().GetArray(4)
    tem = mean_math.GetPointData().GetArray(5)
    shape = [d for d in mean.data.GetDimensions()]
    shapeV = [d for d in mean.data.GetDimensions()]
    shapeV.append(3)
    vel = vel.reshape(shapeV)
    pre = pre.reshape(shape)
    tem = tem.reshape(shape)
    v=vel[minPoint[0]:maxPoint[0],
               minPoint[1]:maxPoint[1],
               minPoint[2]:maxPoint[2],
               :]
    p=pre[minPoint[0]:maxPoint[0],
               minPoint[1]:maxPoint[1],
               minPoint[2]:maxPoint[2]]
    t=tem[minPoint[0]:maxPoint[0],
               minPoint[1]:maxPoint[1],
               minPoint[2]:maxPoint[2]]
    return v,p,t
               
                   
     