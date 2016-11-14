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
    #use to define which dataset to perform operations on
    __MyRealData=[3,5]
    __MyImagData=[0,2]
    def __init__(self,vtkStrGrid):
        self.data=vtkStrGrid
    def __add__(self, other):
        """Return an object that is self+other for fields identified by 
        __MyRealData and __MyImagData"""
        new_data=vtkStructuredGrid()
        new_data.DeepCopy(self.data)
        math_me=dsa.WrapDataObject(self.data)
        math_data=dsa.WrapDataObject(new_data)
        math_other=dsa.WrapDataObject(other.data)
        for i in range(len(self.__MyRealData)):
            math_data.PointData[self.__MyRealData[i]][:]= \
                math_me.PointData[self.__MyRealData[i]][:]+ \
                math_other.PointData[self.__MyRealData[i]][:]
            math_data.PointData[self.__MyImagData[i]][:]= \
                math_me.PointData[self.__MyImagData[i]][:]+ \
                math_other.PointData[self.__MyImagData[i]][:]
        return MrVtkVector(new_data)
        
    def __mul__(self,scalar):
        """Return an object that is self*scalar for fields identified by 
        __MyRealData and __MyImagData"""
        new_data=vtkStructuredGrid()
        new_data.DeepCopy(self.data)
        math_data=dsa.WrapDataObject(new_data)
        math_me=dsa.WrapDataObject(self.data)
        for i in range(len(self.__MyRealData)):
            math_data.PointData[self.__MyRealData[i]][:]= \
                math_me.PointData[self.__MyRealData[i]][:]*np.real(scalar)- \
                math_me.PointData[self.__MyImagData[i]][:]*np.imag(scalar)
            math_data.PointData[self.__MyImagData[i]][:]= \
                math_me.PointData[self.__MyImagData[i]][:]*np.real(scalar)+ \
                math_me.PointData[self.__MyRealData[i]][:]*np.imag(scalar)
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
        for i in range(len(self.__MyRealData)):
            math_data.PointData[self.__MyImagData[i]][:]= \
                math_me.PointData[self.__MyImagData[i]][:]*-1.0
        return MrVtkVector(new_data)
        
    def get_rc_lists(self):
        return (self.__MyRealData,self.__MyImagData)
        
    def weight_matrix(self):
        dims=self.data.GetDimensions()
        total=self.data.GetNumberOfPoints()
        bounds=self.data.GetPoints().GetBounds()
        B=np.array([bounds[1]-bounds[0],bounds[3]-bounds[2],bounds[5]-bounds[4]])
        weights=np.empty(total)
        QD=Quadratures.ChebyshevGauss()
        for k in range(dims(2)):
            for j in range(dims(1)):
                for i in range(dims(0)):
                    ii=i+j*dims(0)+k*(dims(0)*dims(1))
                    weights[ii]=QD.Weight(dims(0),i)*QD.Weight(dims(2),k)
        return weights*0.25*B[0]*B[2]

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
    RL,CL=v1.get_rc_lists()
    for i in range(len(RL)):
        math_new.PointData[RL[i]][:]= \
            math_v1.PointData[RL[i]][:]*math_v2.PointData[RL[i]][:]- \
            math_v1.PointData[RL[i]][:]*math_v2.PointData[CL[i]][:]
        math_new.PointData[CL[i]][:]= \
            math_v1.PointData[CL[i]][:]*math_v2.PointData[RL[i]][:]+ \
            math_v1.PointData[RL[i]][:]*math_v2.PointData[CL[i]][:]
    return MrVtkVector(new_data)
    
def point_division(v1,v2):
    new_data=point_product(v1,v2.complex_conjugate())
    divisor=point_product(v2,v2.complex_conjugate())
    math_new=dsa.WrapDataObject(new_data.data)
    math_div=dsa.WrapDataObject(divisor.data)
    RL,CL=v1.get_rc_lists()
    for i in range(len(RL)):
        math_new.PointData[RL[i]][:]= \
            math_new.PointData[RL[i]][:]/math_div.PointData[RL[i]][:]
        math_new.PointData[CL[i]][:]= \
            math_new.PointData[CL[i]][:]/math_div.PointData[RL[i]][:]
    return new_data