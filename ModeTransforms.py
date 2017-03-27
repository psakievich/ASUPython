# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 10:21:45 2016

@author: psakievich
"""

from vtk import vtkStructuredGrid, \
vtkXMLStructuredGridReader, \
vtkXMLStructuredGridWriter, \
vtkPoints, \
vtkDoubleArray
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa
'''
Routine for documenting individual modes
Transform Fourier coefficients back to real
space. One period of the mode is documented
over a user specified angle
'''
def FourierToRealDoc(fileName,iFFTsize,outputFile,ang,modeNumber=1):
    reader=vtkXMLStructuredGridReader()
    #load the grid into memory
    reader=vtkXMLStructuredGridReader()
    reader.SetFileName(fileName)
    reader.Update()
    fourierGrid=(dsa.WrapDataObject(reader.GetOutput()))
    #clean up
    del reader
    #print(modeNumber)
    #setup mesh dimensions
    nR=fourierGrid.GetDimensions()[0]
    nTheta=iFFTsize+1
    nZ=fourierGrid.GetDimensions()[2]
    nTotal=nR*nZ*nTheta
    keys=fourierGrid.PointData.keys()
    theta=np.linspace(0,ang,nTheta)
    #Create 3D grid
    grid3d=vtkStructuredGrid()
    grid3d.SetDimensions(nTheta,nR,nZ)
    
    points=vtkPoints()
    points.Allocate(nTotal)
    for k in range(nZ):
        for j in range(nR):
            for i in range(nTheta):
                x=fourierGrid.Points[j+k*nR,0]*np.cos(theta[i])
                y=fourierGrid.Points[j+k*nR,0]*np.sin(theta[i])
                z=fourierGrid.Points[j+k*nR,2]
                #print(x,y,z,theta[i],sGrid.Points[j+k*nR,1],sGrid.Points[j+k*nR,2])
                points.InsertNextPoint(x,y,z)
    #Set Up 3D grid
    grid3d.SetPoints(points)
    vel3d=vtkDoubleArray()
    tem3d=vtkDoubleArray()
    pre3d=vtkDoubleArray()
    
    vel3d.SetName('velocity')
    tem3d.SetName('temperature')
    pre3d.SetName('pressure')
    
    vel3d.SetNumberOfComponents(3)
    tem3d.SetNumberOfComponents(1)
    pre3d.SetNumberOfComponents(1)
    
    vel3d.SetNumberOfTuples(nTotal)
    tem3d.SetNumberOfTuples(nTotal)
    pre3d.SetNumberOfTuples(nTotal)
    #Set up iffts and populate 3d Grid
    for k in range(nZ):
        for i in range(nR):
            temp=np.zeros(iFFTsize/2+1,dtype=complex)
            v=np.array([temp.copy(),temp.copy(),temp.copy()])
            press=temp.copy()
            #set up vectors for ifft

            for ii in range(3):
                v[ii,modeNumber]=complex( \
                    fourierGrid.PointData[keys[3]][i+k*nR,ii], \
                    fourierGrid.PointData[keys[0]][i+k*nR,ii])
            press[modeNumber]=complex( \
                fourierGrid.PointData[keys[4]][i+k*nR], \
                fourierGrid.PointData[keys[1]][i+k*nR])
            temp[modeNumber]=complex( \
                fourierGrid.PointData[keys[5]][i+k*nR], \
                fourierGrid.PointData[keys[2]][i+k*nR])
            #scale by grid size
            #print(i,k,temp)
            v=v*iFFTsize
            press=press*iFFTsize
            temp=temp*iFFTsize
            temp1=np.zeros(nTheta)
            press1=np.zeros(nTheta)
            v1=np.array([np.zeros(nTheta),np.zeros(nTheta),np.zeros(nTheta)])
            #conduct iFFT's
            for ii in range(3):
                v1[ii,0:nTheta-1]=np.fft.irfft(v[ii,:])
                v1[ii,nTheta-1]=v1[ii,0]
            press1[0:nTheta-1]=np.fft.irfft(press)
            temp1[0:nTheta-1]=np.fft.irfft(temp)
            press1[nTheta-1]=press1[0]
            temp1[nTheta-1]=temp1[0]
            #Translate data to cartesiant coordinates for visualization purposes
            for jj in range(nTheta):
                index=jj+i*nTheta+k*nR*nTheta
                vx=v1[0,jj]*np.cos(theta[jj])-v1[1,jj]*np.sin(theta[jj])
                vy=v1[0,jj]*np.sin(theta[jj])+v1[1,jj]*np.cos(theta[jj])
                vel3d.SetTuple3(index,vx,vy,v1[2,jj])
                tem3d.SetTuple1(index,temp1[jj])
                pre3d.SetTuple1(index,press1[jj])
    grid3d.GetPointData().SetVectors(vel3d)
    grid3d.GetPointData().AddArray(pre3d)
    grid3d.GetPointData().SetScalars(tem3d)

    writer=vtkXMLStructuredGridWriter()
    writer.SetFileName(outputFile)
    writer.SetInputData(grid3d)
    writer.Write()
'''
Transform a list of Fourier modes back to real space.
Multiple input files are specifie through fileNames
and output as one file outputFile.

fileNames must be a list i.e. []
iFFTsize must be 2 times larger than the largest 
wave number.
'''
def FourierToReal(fileNames,iFFTsize,outputFile):
    numModes=len(fileNames)
    fourierGrids=[]
    modeNumber=[]
    reader=vtkXMLStructuredGridReader()
    #load each of the grids into memory
    for i in range(numModes):
        reader=vtkXMLStructuredGridReader()
        reader.SetFileName(fileNames[i])
        reader.Update()
        fourierGrids.append(dsa.WrapDataObject(reader.GetOutput()))
        tempVar=fileNames[i].split('_')
        modeNumber.append(int(tempVar[1]))
        #clean up
    del reader, tempVar
    #print(modeNumber)
    #setup mesh dimensions
    nR=fourierGrids[0].GetDimensions()[0]
    nTheta=iFFTsize+1
    nZ=fourierGrids[0].GetDimensions()[2]
    nTotal=nR*nZ*nTheta
    keys=fourierGrids[0].PointData.keys()
    theta=np.linspace(0,2*np.pi,nTheta)
    #Create 3D grid
    grid3d=vtkStructuredGrid()
    grid3d.SetDimensions(nTheta,nR,nZ)
    
    points=vtkPoints()
    points.Allocate(nTotal)
    for k in range(nZ):
        for j in range(nR):
            for i in range(nTheta):
                x=fourierGrids[0].Points[j+k*nR,0]*np.cos(theta[i])
                y=fourierGrids[0].Points[j+k*nR,0]*np.sin(theta[i])
                z=fourierGrids[0].Points[j+k*nR,2]
                #print(x,y,z,theta[i],sGrid.Points[j+k*nR,1],sGrid.Points[j+k*nR,2])
                points.InsertNextPoint(x,y,z)
    #Set Up 3D grid
    grid3d.SetPoints(points)
    vel3d=vtkDoubleArray()
    tem3d=vtkDoubleArray()
    pre3d=vtkDoubleArray()
    
    vel3d.SetName('velocity')
    tem3d.SetName('temperature')
    pre3d.SetName('pressure')
    
    vel3d.SetNumberOfComponents(3)
    tem3d.SetNumberOfComponents(1)
    pre3d.SetNumberOfComponents(1)
    
    vel3d.SetNumberOfTuples(nTotal)
    tem3d.SetNumberOfTuples(nTotal)
    pre3d.SetNumberOfTuples(nTotal)
    #Set up iffts and populate 3d Grid
    for k in range(nZ):
        for i in range(nR):
            temp=np.zeros(iFFTsize/2+1,dtype=complex)
            v=np.array([temp.copy(),temp.copy(),temp.copy()])
            press=temp.copy()
            #set up vectors for ifft
            for kk in range(numModes):
                for ii in range(3):
                    v[ii,modeNumber[kk]]=complex( \
                        fourierGrids[kk].PointData[keys[3]][i+k*nR,ii], \
                        fourierGrids[kk].PointData[keys[0]][i+k*nR,ii])
                press[modeNumber[kk]]=complex( \
                    fourierGrids[kk].PointData[keys[4]][i+k*nR], \
                    fourierGrids[kk].PointData[keys[1]][i+k*nR])
                temp[modeNumber[kk]]=complex( \
                    fourierGrids[kk].PointData[keys[5]][i+k*nR], \
                    fourierGrids[kk].PointData[keys[2]][i+k*nR])
            #scale by grid size
            #print(i,k,temp)
            v=v*iFFTsize
            press=press*iFFTsize
            temp=temp*iFFTsize
            temp1=np.zeros(nTheta)
            press1=np.zeros(nTheta)
            v1=np.array([np.zeros(nTheta),np.zeros(nTheta),np.zeros(nTheta)])
            #conduct iFFT's
            for ii in range(3):
                v1[ii,0:nTheta-1]=np.fft.irfft(v[ii,:])
                v1[ii,nTheta-1]=v1[ii,0]
            press1[0:nTheta-1]=np.fft.irfft(press)
            temp1[0:nTheta-1]=np.fft.irfft(temp)
            press1[nTheta-1]=press1[0]
            temp1[nTheta-1]=temp1[0]
            #Translate data to cartesiant coordinates for visualization purposes
            for jj in range(nTheta):
                index=jj+i*nTheta+k*nR*nTheta
                vx=v1[0,jj]*np.cos(theta[jj])-v1[1,jj]*np.sin(theta[jj])
                vy=v1[0,jj]*np.sin(theta[jj])+v1[1,jj]*np.cos(theta[jj])
                vel3d.SetTuple3(index,vx,vy,v1[2,jj])
                tem3d.SetTuple1(index,temp1[jj])
                pre3d.SetTuple1(index,press1[jj])
    grid3d.GetPointData().SetVectors(vel3d)
    grid3d.GetPointData().AddArray(pre3d)
    grid3d.GetPointData().SetScalars(tem3d)

    writer=vtkXMLStructuredGridWriter()
    writer.SetFileName(outputFile)
    writer.SetInputData(grid3d)
    writer.Write()