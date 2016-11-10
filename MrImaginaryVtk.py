# -*- coding: utf-8 -*-
"""
Created on Sat Oct 29 11:16:59 2016

@author: psakievich
"""
import modred as mr
import numpy as np
#VTK RELATED STUFFS
from vtk import vtkXMLStructuredGridReader,vtkXMLStructuredGridWriter, \
    vtkStructuredGrid
from vtk.numpy_interface import dataset_adapter as dsa


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
        math_me=dsa.WrapDataObject(self.data)
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