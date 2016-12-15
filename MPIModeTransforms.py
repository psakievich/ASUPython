#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 15:05:25 2016

Parallel Processing of iFFT

@author: psakievich
"""

#Dependencies
from mpi4py import MPI
import numpy as np

DEBUG=True
comm=MPI.COMM_WORLD
size=comm.Get_size() #num processors
rank=comm.Get_rank() #process id

#fileTemp='A10SymS1_{}_TSTEP_{}.vts'
#fileDir='../Modes/'
#fileList=[fileTemp.format(i,0) for i in range(4)]
#if DEBUG:
#   if rank == 0:
#      print 'The files listed are:'
#      for f in fileList:
#         print f

def DefineLocalSizes(GlobalSize,Partitions,NumCores):
   '''
   This routine defines the segmentation of an object into chunks that can be 
   sent to other processors
   '''
   nC=int(NumCores)
   if len(GlobalSize)!=len(Partitions):
       print 'ERROR:: Partion size must be declared for each dimension in GlobalSize'
       return None
   if np.prod(Partitions)!=nC:
       print 'ERROR:: Product of partions must match Number of cores'
       return None
   divisor=GlobalSize/Partitions
   remainder=GlobalSize%Partitions
   sizeLocal=np.outer(np.ones(nC,dtype=int),divisor.T)
   #add +1 if remainders
   oddVals=sizeLocal.copy()*0
   
   for i in range(nC):
       #loop over processors
       for j in range(len(Partitions)):
           #loop over dimensions to split
           if len(Partitions)>1:
               if(j==0):
                   index=i/Partitions[j]
               else:
                   index=i%np.prod(Partitions[0:j])
           else:
               index=i%Partitions[0]
             
           if index<remainder[j]:
               oddVals[i,j]+=1
   sizeLocal+=oddVals
   return sizeLocal
   
def DefineExtents(GlobalSize,Partitions,NumCores):
    '''
    Routine to determine the extents of a partioned grid
    
    Extents are left inclusive i.e. range =[). Extents match the local size 
    determined for this partion map by the DefineLocalSizes routine
    '''
    sizeLocal=DefineLocalSizes(GlobalSize,Partitions,NumCores)
    nC=int(NumCores)
    nPrt=len(Partitions)
    extent=np.zeros([nC,nPrt*2],dtype=int)
    cur=np.zeros(nPrt,dtype=int)
    old=cur.copy()
    icur=cur.copy()
    iold=cur.copy()
    for i in range(nC):
        for j in range(nPrt):
            #determine location in the index map
            if nPrt>1:
                if j==0:
                    icur[j]=i/Partitions[j]
                else:
                    icur[j]=i%np.prod(Partitions[0:j])
            else:
                icur[j]=i%Partitions[0]
            #update extent based off changes in index mapp
            if icur[j]==0:
                #if index resets
                iold[j]=0
                cur[j]=0
            elif icur[j]>iold[j]:
                #if index advances
                iold[j]=icur[j]
                cur[j]=old[j]
            extent[i,j*2]=cur[j]
            extent[i,j*2+1]=cur[j]+sizeLocal[i,j]
            old[j]=extent[i,j*2+1]

    return extent
            
   