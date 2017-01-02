#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 15:30:20 2017

Compute Integral Time Scale in Parallel

This computes the integral time scale with each processor computing the 
results for a different wave number.  The details of how this is done and the
mathematical theory can be found in the python documents section.

The flucuations are defined with respect to the combined temporal and azimuthal
 average.
 
 This file should be copied into the directory where you would like the output
 files written.  It should be located on a NFS so all processors can write
 their data for the individual wavenumber's integral time scale

@author: psakievich
"""
#where the custom python libraries are kept
sourcePath=''
#where the files that you will be processing 
#are stored
filePath=''
snapsFileTemplate='{0}/SymS1_{0}_TSTEP_{1}.vts'
avgFile='{0}/S1Wave_{0}_TAVG_{1}-{2}.vts'
outputTemplate='IntegralTSWave_{}_TSTEP_{}-{}'
totalOutTemp='IntegralTTotal_TSTEP_{}-{}'
numFiles=600
index1=0
index2=index1+numFiles-1

# import python modules
import numpy as np
import sys
#add in correct path
sys.path.append(sourcePath)
from mpi4py import MPI
import MrImaginaryVtk as MIV
import vtk
from vtk.numpy_interface import dataset_adapter as dsa

#Set up mpi data
comm=MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()

#Set up fluctuating field
handleAvg=MIV.MrVtkVecHandle(filePath+avgFile.format(rank,index1,index2))
avgField=handleAvg.get()
fluField=[]
#zero out average for all wave numbers except zero
if rank==0:   
    for i in range(index1,index2+1):
        hTemp=MIV.MrVtkVecHandle(filePath+snapsFileTemplate.format(rank,i))
        fluField.append(hTemp.get()+-1.0*avgField)
else:
    for i in range(index1,index2+1):
        hTemp=MIV.MrVtkVecHandle(filePath+snapsFileTemplate.format(rank,i))
        fluField.append(hTemp.get())

#set up variables
nPoints=avgField.data.GetNumberOfPoints()
Rii=np.zeros([numFiles//2,nPoints],dtype=np.float64)
corVel=Rii.copy()

for i in range(numFiles//2):
    #assuming the signal is not periodic so we compute for seperations up to N/2
    #range decreases as seperation time grows
    for j in range(numFiles-i):
        #compute correlations with seperation i for individual fields
        temp=MIV.point_product(fluField[j],fluField[j+i].complex_conjugate())
        #sum velocity components
        mTemp=dsa.WrapDataObject(temp.data)
        Rii[i,:]+=np.sum(mTemp.PointData[3],axis=1)+np.sum(mTemp.PointData[0],axis=1)
    Rii[i,:]=Rii[i,:]/(numFiles-i)
    
if rank >0:
    #double the size to acount for the negative wavenumbers
    Rii*=2.0
    
#Need to do MPI summation here for all wavenumbers
comm.reduce(Rii,corVel,opp=MPI.SUM,root=0)
for i in range(numFiles//2-1,-1,-1):    
    #Convert to correlation coefficient bounded by -1:1
    Rii[i,:]=Rii[i,:]/Rii[0,:]
    corVel[i,:]=corVel[i,:]/corVel[i,:]

#setup output variables for integral times
timeLocal=Rii[0,:]*0.0
if rank==0:
    timeTotal=corVel[0,:]*0.0

for i in range(numFiles//2):
    timeLocal+=Rii[i,:]
if rank==0:
    for i in range(numFiles//2):
        timeTotal+=corVel[i,:]

np.save(outputTemplate.format(rank,index1,index2),allow_pickle=False)
if rank==0:
    np.save(totalOutTemp.format(index1,index2),allow_pickle=False)
    



