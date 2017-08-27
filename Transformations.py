#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 13:06:10 2017

@author: psakievich
"""

import MrImaginaryVtk as MIV

def TransformVertState(snapshot):
    newSnap=snapshot.DeepCopy()
    npNew=MIV.dsa.WrapDataObject(newSnap.data)
    npSnap=MIV.dsa.WrapDataObject(snapshot.data)
    shape=snapshot.data.GetDimensions()
    npNew.PointData[0][-shape[0]::]=npSnap.PointData[0][0:shape[0]]*(1.0,1.0,-1.0)
    npNew.PointData[3][-shape[0]::]=npSnap.PointData[3][0:shape[0]]*(1.0,1.0,-1.0)
    npNew.PointData[2][-shape[0]::]=-npSnap.PointData[2][0:shape[0]]
    npNew.PointData[5][-shape[0]::]=-npSnap.PointData[5][0:shape[0]]
    npNew.PointData[1][-shape[0]::]=npSnap.PointData[1][0:shape[0]]
    npNew.PointData[4][-shape[0]::]=npSnap.PointData[4][0:shape[0]]
    for i in range(1,shape[2]):
        npNew.PointData[0][-shape[0]*(i+1):-shape[0]*i]= \
                        npSnap.PointData[0][shape[0]*i:shape[0]*(i+1)]*(1.0,1.0,-1.0)
        npNew.PointData[3][-shape[0]*(i+1):-shape[0]*i]= \
                        npSnap.PointData[3][shape[0]*i:shape[0]*(i+1)]*(1.0,1.0,-1.0)
        npNew.PointData[2][-shape[0]*(i+1):-shape[0]*i]= \
                        -npSnap.PointData[2][shape[0]*i:shape[0]*(i+1)]
        npNew.PointData[5][-shape[0]*(i+1):-shape[0]*i]= \
                        -npSnap.PointData[5][shape[0]*i:shape[0]*(i+1)] 
        npNew.PointData[1][-shape[0]*(i+1):-shape[0]*i]= \
                        npSnap.PointData[1][shape[0]*i:shape[0]*(i+1)]
        npNew.PointData[4][-shape[0]*(i+1):-shape[0]*i]= \
                        npSnap.PointData[4][shape[0]*i:shape[0]*(i+1)]                
    return newSnap