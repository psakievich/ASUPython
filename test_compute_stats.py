#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 21:41:29 2019

@author: psakiev
"""

import unittest
import numpy as np
import ParallelMapping as pm
from vtk.numpy_interface import dataset_adapter as dsa
from operator import mul

class TestComputingStats(unittest.TestCase):
    def setUp(self):
        # mean def
        self.MeanSize = (1, 64, 160) #swaped for construction
        self.MeanTotal = reduce(mul, self.MeanSize)
        temp = np.arange(0,self.MeanTotal)
        vel = np.array([temp,temp*2.0,temp*4.0],dtype=float).T
        self.np_mean_arrays=[vel,temp,temp*0.25]
        #ins def
        self.thetaSample=5
        self.realBlocks = (1,8,160)
        self.InsSize = (self.thetaSample/self.realBlocks[0],
                        self.MeanSize[1]/self.realBlocks[1],
                        self.MeanSize[2]/self.realBlocks[2])
        self.InsTotal = reduce(mul, self.InsSize)
        self.writeRank = 0

        temp = np.arange(0,self.InsTotal)
        vel = np.array([temp,temp*2.0,temp*4.0],dtype=float).T
        self.np_ins_arrays=[vel,temp,temp*0.25]
        self.corner0 = pm.Local2Global((0,0,0),self.writeRank, self.realBlocks,
                                       self.InsSize)
        corner1 = pm.Local2Global(self.InsSize,self.writeRank,
                                       self.realBlocks, self.InsSize)
        self.corner1 = (1,corner1[1],corner1[2])
        self.mean = pm.CreateMrImaginary(self.MeanSize, 0, (1,1,1), 
                                      self.np_mean_arrays)
        #test changing fastest index to get arrays to align with inst. snaps
        self.mean.data.SetDimensions(self.MeanSize)  
        self.inst = pm.CreateMrReal(self.InsSize, self.writeRank, self.realBlocks,
                                 self.np_ins_arrays)
    #@unittest.skipIf(self.writeRank !=0)    
    def test_GetCorrectCorners(self):
        if self.writeRank == 0:
            self.assertEqual((0,0,0),self.corner0)
            self.assertEqual((1,self.InsSize[1],self.InsSize[2]),self.corner1)
    
    def test_ExtractRangeFromMeanShape(self):
        v,p,t = pm.ExtractRangeFromMean(self.mean, self.corner0, self.corner1)
        self.assertEqual((1,8,1,3),v.shape)
        self.assertEqual((1,8,1),p.shape)
        self.assertEqual((1,8,1),t.shape)
        
    def test_AreValuesIndexedCorrectly(self):
        if self.writeRank != 0: 
            return
        v,p,t = pm.ExtractRangeFromMean(self.mean,self.corner0,self.MeanSize)
        for z in range(self.MeanSize[2]):
            for r in range(self.MeanSize[1]):
                vNp = self.np_mean_arrays[1][r+z*self.MeanSize[1]]
                vVtk = self.mean.data.GetPointData().GetArray(4).GetTuple(
                        r+z*self.MeanSize[1])
                self.assertEqual(vNp,vVtk,(r,z,vNp,vVtk))
                self.assertEqual(vNp,p[z,r,0],(r,z,vNp,p[z,r,0]))
                
    #@unittest.skip("Skipping extraction")   
    def test_ExtractRangeFromMeanValues(self):
        v,p,t = pm.ExtractRangeFromMean(self.mean, self.corner0, self.corner1)
        vShape = [r for r in self.MeanSize[::-1]]
        vShape.append(3)
        vExact = self.np_mean_arrays[0].reshape(vShape)
        self.assertTrue((vExact[self.corner0[0]:self.corner1[0]
                                ,self.corner0[1]:self.corner1[1]
                                ,self.corner0[2]:self.corner1[2]
                                ,0]==v[:,:,:,0]).all())
        
    #@unittest.skip("Skipping extraction")    
    def test_MatchValuesFromMeanAndInst(self):
        vM,pM,tM = pm.ExtractRangeFromMean(self.mean,self.corner0,self.corner1)
        m_ins = dsa.WrapDataObject(self.inst.data)
        
        pI=m_ins.PointData[1].reshape(self.inst.data.GetDimensions()[::-1])
        self.assertTrue(pI[:,:,0].shape==pM[:,:,0].shape,
                        (pI[:,:,0].shape,pM[:,:,0].shape))
        if self.writeRank != 0:
            return
        self.assertTrue((pI[:,:,0]/self.InsSize[0]==pM[:,:,0]).all(),
                        (pI[:,:,0],pM[:,:,0]))
        
    def test_ComputingFluctuations(self):
        meanL,corner0,corner1 = pm.CreateLocalMeanMr(self.mean,
                                                     self.inst,
                                                     self.writeRank,
                                                     self.realBlocks)
        self.assertEqual(self.corner0,corner0)
        self.assertEqual(self.corner1,corner1)
        arrays = pm.ExtractRangeFromMean(self.mean,self.corner0,self.corner1)
        fluc = self.inst-meanL

        # check against and explict computation of the fluctuation
        for i in range(self.InsSize[2]):
            for j in range(self.InsSize[1]):
                for k in range(self.InsSize[0]):
                    # velocity matches
                    self.assertTrue((
                            fluc.data.GetPointData().GetArray(0).GetTuple3(
                            k+(j+i*self.InsSize[1])*self.InsSize[0])==
                             self.inst.data.GetPointData().GetArray(0).GetTuple3(
                                     k+(j+i*self.InsSize[1])*self.InsSize[0])
                             -arrays[0][i,j,0]).all())
                    # pressure matches
                    self.assertEqual(
                            fluc.data.GetPointData().GetArray(1).GetTuple(
                            k+(j+i*self.InsSize[1])*self.InsSize[0]),
                             self.inst.data.GetPointData().GetArray(1).GetTuple(
                                     k+(j+i*self.InsSize[1])*self.InsSize[0])
                             -arrays[1][i,j,0])
                             
                    #temperature matches
                    self.assertEqual(
                            fluc.data.GetPointData().GetArray(2).GetTuple(
                            k+(j+i*self.InsSize[1])*self.InsSize[0]),
                             self.inst.data.GetPointData().GetArray(2).GetTuple(
                                     k+(j+i*self.InsSize[1])*self.InsSize[0])
                             -arrays[2][i,j,0])
                             
    def test_AzimuthalAverage(self):
        meanL,corner0,corner1 = pm.CreateLocalMeanMr(self.mean,
                                                     self.inst,
                                                     self.writeRank,
                                                     self.realBlocks)
        fluc = self.inst - meanL
        avg = pm.AzimuthalAverage(fluc, self.writeRank, self.realBlocks)
        self.assertEquals(avg.data.GetDimensions(),(1,8,1))
    
    def test_AssembleInSpace(self):
        newGrid = self.mean*0.0
        # test for no overlap
        for i in range(4):
            meanL,corner0,corner1 = pm.CreateLocalMeanMr(self.mean,
                                                     self.inst,
                                                     i,
                                                     self.realBlocks)
            aAvg = pm.AzimuthalAverage(meanL, self.writeRank, self.realBlocks)
            pm.AddSubset(newGrid,aAvg,corner0,corner1)
        ng_math = dsa.WrapDataObject(newGrid.data)
        self.assertFalse(np.sum(ng_math.PointData[1]==0)==8*4)
        

if __name__=='__main__':
    unittest.main()