# -*- coding: utf-8 -*-
"""
Created on Sat Oct 29 11:16:59 2016

@author: psakievich
"""
import modred as mr
import numpy as np
import Quadratures
#VTK RELATED STUFFS
from vtk import vtkXMLStructuredGridReader,vtkXMLStructuredGridWriter, \
    vtkStructuredGrid
from vtk.numpy_interface import dataset_adapter as dsa
'''
Vector class
'''
class MrVtkVector(mr.Vector):
    #use to define which datasets for inner product
    __MyRealData=[3,5]
    __MyImagData=[0,2]
    def __init__(self,vtkStrGrid):
        self.data=vtkStrGrid
    def __add__(self, other):
        """Return an object that is self+other for all fields 
        """
        new_data=vtkStructuredGrid()
        new_data.DeepCopy(self.data)
        math_me=dsa.WrapDataObject(self.data)
        math_data=dsa.WrapDataObject(new_data)
        math_other=dsa.WrapDataObject(other.data)
        numFlds=len(math_me.PointData.keys())
        for i in range(numFlds):
            math_data.PointData[i][:]= \
                math_me.PointData[i][:]+ \
                math_other.PointData[i][:]

        return MrVtkVector(new_data)
        
    def __mul__(self,scalar):
        """Return an object that is self*scalar for all fields  
        """
        new_data=vtkStructuredGrid()
        new_data.DeepCopy(self.data)
        math_data=dsa.WrapDataObject(new_data) 
        math_me=dsa.WrapDataObject(self.data)
        numFlds=len(math_me.PointData.keys())
        numReal=int(numFlds/2)  
        for i in range(numReal):
            math_data.PointData[i+numReal][:]= \
                math_me.PointData[i+numReal][:]*np.real(scalar)- \
                math_me.PointData[i][:]*np.imag(scalar)
            math_data.PointData[i][:]= \
                math_me.PointData[i][:]*np.real(scalar)+ \
                math_me.PointData[i+numReal][:]*np.imag(scalar)
        return MrVtkVector(new_data)
    
    def inner_product(self,other):
        weighted_me=self.weighted_copy()
        math_me=dsa.WrapDataObject(weighted_me.data)
        math_other=dsa.WrapDataObject(other.data)
        IP=0.0
        for i in range(len(self.__MyImagData)):
            IP=IP+np.vdot(math_me.PointData[self.__MyRealData[i]][:]+ \
                1j*math_me.PointData[self.__MyImagData[i]][:], \
                math_other.PointData[self.__MyRealData[i]][:]+ \
                1j*math_other.PointData[self.__MyImagData[i]][:])
        return IP
        
    def complex_conjugate(self):
        new_data=vtkStructuredGrid()
        new_data.DeepCopy(self.data)
        math_data=dsa.WrapDataObject(new_data)
        math_me=dsa.WrapDataObject(self.data)
        numFlds=len(math_me.PointData.keys())
        for i in range(numFlds/2):
            math_data.PointData[i][:]*=-1.0
        return MrVtkVector(new_data)
        
    def integrated_values(self):
        weighted_me=self.weighted_copy()
        math_me=dsa.WrapDataObject(weighted_me.data)
        numFlds=len(math_me.PointData.keys())
        k=0
        for i in range(numFlds):
           if(len(math_me.PointData[i].shape)>1):
               k=k+math_me.PointData[i].shape[1]
           else:
               k=k+1 
        result=np.empty(k)
        j=0
        for i in range(numFlds):
            if(len(math_me.PointData[i].shape)>1):
                for k in range(math_me.PointData[i].shape[1]):
                    result[j]=np.sum(math_me.PointData[i][:,k])
                    j=j+1
            else:
                result[j]=np.sum(math_me.PointData[i][:])
                j=j+1
        return result
        
    def get_rc_lists(self):
        return (self.__MyRealData,self.__MyImagData)
        
    def weight_matrix(self,QD=Quadratures.GaussLegendre()):
        dims=self.data.GetDimensions()
        bounds=self.data.GetPoints().GetBounds()
        B=np.array([bounds[1]-bounds[0],bounds[3]-bounds[2],bounds[5]-bounds[4]])
        wz=QD.Weights(dims[2])
        wr=QD.Weights(dims[0])
        weights=np.outer(wz,wr) #r is fastest varying in dataset
        weights=np.reshape(weights,dims[0]*dims[2])
        math_me=dsa.WrapDataObject(self.data)
        weights=weights*math_me.Points[:,0] #multiply by R
        weights=weights*0.25*B[0]*B[2] #multiply by jacobian
        return weights

    def weighted_copy(self):
        new_data=vtkStructuredGrid()
        new_data.DeepCopy(self.data)
        math_new=dsa.WrapDataObject(new_data)
        w=self.weight_matrix()
        nFields=len(math_new.PointData.keys())
        for i in range(nFields):
            if(len(math_new.PointData[i].shape)>1):
                for j in range(math_new.PointData[i].shape[1]):
                    math_new.PointData[i][:,j]=math_new.PointData[i][:,j]*w
            else:
                math_new.PointData[i][:]=math_new.PointData[i][:]*w
        return MrVtkVector(new_data)
'''
Vector handle
'''        
class MrVtkVecHandle(mr.VecHandle):
    def __init__(self, vec_path, base_handle=None, scale=None):
        mr.VecHandle.__init__(self,base_handle,scale)
        self.vec_path=vec_path
        
    def _get(self):
        reader=vtkXMLStructuredGridReader()
        reader.SetFileName(self.vec_path)
        reader.Update()
        return(MrVtkVector(reader.GetOutput()))
        
    def _put(self,vec):
        writer=vtkXMLStructuredGridWriter()
        writer.SetInputData(vec.data)
        writer.SetFileName(self.vec_path)
        writer.Write()
class MrVtkVecHandleCreateFluctuation(mr.VecHandle):
    def __init__(self, vec_path_inst, vec_path_mean,base_handle=None, scale=None):
        mr.VecHandle.__init__(self,base_handle,scale)
        self.vec_path=vec_path_inst
        self.vec_path_mean=vec_path_mean
        
    def _get(self):
        reader=vtkXMLStructuredGridReader()
        reader.SetFileName(self.vec_path)
        reader.Update()
        hMean=MrVtkVecHandle(self.vec_path_mean)
        return(MrVtkVector(reader.GetOutput())+ \
               -1.0*hMean.get())
        
    def _put(self,vec):
        writer=vtkXMLStructuredGridWriter()
        writer.SetInputData(vec.data)
        writer.SetFileName(self.vec_path)
        writer.Write()

class MrVtkVecHandleOperateOnFluctuation(mr.VecHandle):
    def __init__(self, vec_path_inst, vec_path_mean,base_handle=None, scale=None):
        mr.VecHandle.__init__(self,base_handle,scale)
        self.vec_path=vec_path_inst
        self.vec_path_mean=vec_path_mean
        
    def _get(self):
        reader=vtkXMLStructuredGridReader()
        reader.SetFileName(self.vec_path)
        reader.Update()
        hMean=MrVtkVecHandle(self.vec_path_mean)
        return(MrVtkVector(reader.GetOutput())+ \
               -1.0*hMean.get())
        
    def _put(self,vec):
        hMean=MrVtkVecHandle(self.vec_path_mean)
        vec+=hMean.get()
        writer=vtkXMLStructuredGridWriter()
        writer.SetInputData(vec.data)
        writer.SetFileName(self.vec_path)
        writer.Write()
'''
 Namespace functions
'''       
def inner_product(v1,v2):
    return v1.inner_product(v2)

def point_product(v1,v2):
    new_data=vtkStructuredGrid()
    new_data.DeepCopy(v1.data)
    math_new=dsa.WrapDataObject(new_data)
    math_v1=dsa.WrapDataObject(v1.data)
    math_v2=dsa.WrapDataObject(v2.data)
    nFields=len(math_new.PointData.keys())
    offset=nFields//2
    for i in range(nFields//2):
        math_new.PointData[i+offset][:]= \
            math_v1.PointData[i+offset][:]*math_v2.PointData[i+offset][:]- \
            math_v1.PointData[i][:]*math_v2.PointData[i][:]
        math_new.PointData[i][:]= \
            math_v1.PointData[i][:]*math_v2.PointData[i+offset][:]+ \
            math_v1.PointData[i+offset][:]*math_v2.PointData[i][:]
    return MrVtkVector(new_data)
    
def point_division(v1,v2):
    new_data=point_product(v1,v2.complex_conjugate())
    divisor=point_product(v2,v2.complex_conjugate())
    math_new=dsa.WrapDataObject(new_data.data)
    math_div=dsa.WrapDataObject(divisor.data)
    nFields=len(math_new.PointData.keys())
    offset=nFields//2
    for i in range(offset):
        math_new.PointData[i+offset][:]= \
            math_new.PointData[i+offset][:]/math_div.PointData[i+offset][:]
        math_new.PointData[i][:]= \
            math_new.PointData[i][:]/math_div.PointData[i+offset][:]
    return new_data