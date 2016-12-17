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
import vtk
from vtk.numpy_interface import dataset_adapter as dsa

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
   
def DefineExtents(GlobalSize,Partitions,NumCores,sizeLocal):
    '''
    Routine to determine the extents of a partioned grid
    
    Extents are left inclusive i.e. range =[). Extents match the local size 
    determined for this partion map by the DefineLocalSizes routine
    '''
    #sizeLocal=DefineLocalSizes(GlobalSize,Partitions,NumCores)
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

def DivideFileList(fileList,rank,size):
    '''
    Determine which files to read in for this rank
    '''
    numFiles=len(fileList)
    if rank>numFiles:
        return None
    else:
        globalSize=np.array(numFiles,dtype=int)
        readers=max(numFiles,size)
        partitions=np.array(readers,dtype=int)
        return DefineExtents(globalSize,partitions,readers)[rank] 

def RzFlatten(r,z,SIZES,MAJOR='r'):
    '''
    Flaten index from 2d to 1d using the specified major ordering 
    (Default major direction is r)
    '''
    if MAJOR.lower()=='r':
        for i in range(SIZES[1]):
            for j in range(SIZES[0]):
                return j+i*SIZES[0]
    elif MAJOR.lower()=='z':
        for i in range(SIZES[0]):
            for j in range(SIZES[1]):
                return j+i*SIZES[1]
     
    
def Foruier2Real(fileList,aRes,rzDims,rzPars):
    '''
    Main routine for transforming the fourier modes to a real data file
    return vtkStructuredGrid
    
    @param fileList list of files to be read to populate the data set
    @param aRes     azimuthal resolution of the real grid
    @param rzDims   a numpy integer array containing the global r,z resolution
    @param rzPars   a numpy integer array showing how many partitions in each direction
    '''    
    comm=MPI.COMM_WORLD
    size=comm.Get_size()
    rank=comm.Get_rank()
    #split files among processors
    files2Read=DivideFileList(fileList,rank,size)
    #define extents for each processor
    rzLDims=DefineLocalSizes(rzDims,rzPars,size)
    rzExtents=DefineExtents(rzDims,rzPars,size,rzLDims)
    #set buffers
    valueBuffer=np.empty(0,dtype=complex)
    infoBuffer =[]
    '''
    Variables are stored in grids as follows:
    grid[0,:,:,:]-->radial velocity
    grid[1,:,:,:]-->azimuthal velocity
    grid[2,:,:,:]-->vertical velocity
    grid[3,:,:,:]-->pressure
    grid[4,:,:,:]-->temperature
    '''
    fGrid=np.empty(5,[rzLDims[rank][1],rzLDims[rank][0],aRes/2+1],dtype=complex)
    rGrid=np.empty(5,[rzLDims[rank][1],rzLDims[rank][1],aRes])
    #read in data and generate buffers
    if files2Read is not None:
        #read in data for each file in the list
        for i in range(files2Read[0],files2Read[1]):
            reader = vtk.vtkXMLStructuredGridReader()
            reader.SetFileName(fileList[i])
            waveNumber=int(fileList[i].split('_')[1])
            reader.Update()
            waveGrid=reader.GetOutput()
            numpWG=dsa.WrapDataObject(waveGrid)
            #fill in my data
            rRange=range(rzExtents[rank][0],rzExtents[rank][1])
            zRange=range(rzExtents[rank][2],rzExtents[rank][2])
            '''
            Move data from vtk grid to pure numpy arrays for ease of transport
            and post processing.  A little cache thrashing now to eliminate it 
            a lot later.  Not really a way to get around it because the way vtk 
            stores vectors.
            '''
            # account for negative wave number by turning
            # them into positive wave numbers
            if waveNumber<0:
                ccAdjust=-1.0
                waveNumber=abs(waveNumber)
            else:
                ccAdjust=1.0
            #-----------------------------------------------------------------
            '''
            Some functions to make errors less prone
            '''
            def FillFGrid(fGridIndex,realIndex,complexIndex,vecIndex):
                '''
                function for populating fGrid
                '''
                if vecIndex is not None:
                    for ii in range(rzLDims[rank][1]):
                        for jj in range(rzLDims[rank][0]):
                            index=RzFlatten(rRange(jj),zRange(ii),rzDims)
                            fGrid[fGridIndex,ii,jj,waveNumber]= \
                                complex(numpWG.PointData[realIndex][index,vecIndex], \
                                numpWG.PointData[complexIndex][index,vecIndex] \
                                *ccAdjust)
                else:
                    for ii in range(rzLDims[rank][1]):
                        for jj in range(rzLDims[rank][0]):
                            index=RzFlatten(rRange(jj),zRange(ii),rzDims)
                            fGrid[fGridIndex,ii,jj,waveNumber]= \
                                complex(numpWG.PointData[realIndex][index], \
                                numpWG.PointData[complexIndex][index]*ccAdjust)
            def FillBuffers(realIndex,complexIndex,vecIndex,recvId):
                '''
                function to populate the crystal router buffers
                '''
                localBuffer=np.empty(np.prod(rzLDims[recvId]),dtype=complex)
                if vecIndex is not None:
                    for ii in range(rzLDims[recvId][1]):
                        for jj in range(rzLDims[recvId][0]):
                            index=RzFlatten(rRange(jj),zRange(ii),rzDims)
                            iL=RzFlatten(jj,ii,rzLDims[recvId])
                            localBuffer[iL]= \
                                complex(numpWG.PointData[realIndex][index,vecIndex], \
                                numpWG.PointData[complexIndex][index,vecIndex] \
                                *ccAdjust)
                else:
                    for ii in range(rzLDims[recvId][1]):
                        for jj in range(rzLDims[recvId][0]):
                            index=RzFlatten(rRange(jj),zRange(ii),rzDims)
                            iL=RzFlatten(jj,ii,rzLDims[recvId])
                            localBuffer[iL]= \
                                complex(numpWG.PointData[realIndex][index], \
                                numpWG.PointData[complexIndex][index] \
                                *ccAdjust)
                np.append(valueBuffer,localBuffer)
                infoBuffer.append((recvId,waveNumber))
            #------------------------------------------------------------------
            for k in range(size):
                #loop over partitions and fill in buffer or local data depending 
                #on partition id
                if k==rank:
                    FillFGrid(0,3,0,0)#ur
                    FillFGrid(1,3,0,1)#utheta
                    FillFGrid(2,3,0,2)#uz
                    FillFGrid(4,4,1,None)#Pressure
                    FillFGrid(5,5,2,None)#temperature
                else:
                    FillBuffers(3,0,0,k)
                    FillBuffers(3,0,1,k)
                    FillBuffers(3,0,2,k)
                    FillBuffers(4,1,None,k)
                    FillBuffers(5,2,None,k)
            

    #transfer data
    #populate fourier grid from buffers
    #perform ifft
    #return value
    
            
            
   