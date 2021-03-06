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
from vtk.util import numpy_support as ns
'''
Vector class
This class operates on the flow field variables 
as a single, flattened vector.  The vector 
interfaces with the VTK structured grid.  
Scalar*Vector, Vector+Vector and 
(Vector,Vector).  The actual variables that 
are used in the in the inner product are defined
by the variables __MyRealData and __MyImagData.
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
    def DeepCopy(self):
        temp=vtkStructuredGrid()
        temp.DeepCopy(self.data)
        return MrVtkVector(temp)   
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
          iR=self.__MyRealData[i]
          iI=self.__MyImagData[i]
          tempR=ns.vtk_to_numpy(weighted_me.data.GetPointData().GetArray(iR))
          tempI=ns.vtk_to_numpy(weighted_me.data.GetPointData().GetArray(iI))*1j
          me=tempR-tempI
          #print type(me),type(me[0])
        #  tempR=ns.vtk_to_numpy(other.data.GetPointData().GetArray(iR))
        #  tempI=ns.vtk_to_numpy(other.data.GetPointData().GetArray(iI))*1j
        #  ot=tempR+tempI
        #  IP+=np.sum(me*ot)
          IP=IP+np.vdot(math_me.PointData[iR][:]+ \
               1j*math_me.PointData[iI][:], \
               math_other.PointData[iR][:]+ \
               1j*math_other.PointData[iI][:])

         # if(len(math_me.PointData[iI].shape)>1):
         #   for j in range(math_me.PointData[iI].shape[1]):
         #      IP+=np.sum((math_me.PointData[iR][:,j]-1j*math_me.PointData[iI][:,j])*\
         #                 (math_other.PointData[iR][:,j]+1j*math_other.PointData[iI][:,j]).T)
         # else:
         #   IP+=np.sum((math_me.PointData[iR][:]-1j*math_me.PointData[iI][:])*\
         #              (math_other.PointData[iR][:]+1j*math_other.PointData[iI][:]).T)
        #print(type(IP),IP)
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
        
    def weight_matrix(self,QD=Quadratures.GaussLegendre(),rBool=True,zBool=True):
        '''
        Weighting matrix for the numerical integration. Different 
        Quadratures can be specified
        '''
        dims=self.data.GetDimensions()
        bounds=self.data.GetPoints().GetBounds()
        #B=np.array([bounds[1]-bounds[0],bounds[3]-bounds[2],bounds[5]-bounds[4]])
        #B=[3.15,0.0,1.0]
        if zBool:
          wz=QD.Weights(dims[2])
        else:
          wz=np.ones(dims[2])
        if rBool:
          wr=QD.Weights(dims[0])
        else:
          wr=np.ones(dims[0])
        pr=QD.Points(dims[0])
        pz=QD.Points(dims[2])
        B=np.array([(bounds[1]-bounds[0])/(pr[dims[0]-1]-pr[0])*2.0, \
                    (bounds[3]-bounds[2]), \
                    (bounds[5]-bounds[4])/(pz[dims[2]-1]-pz[0])*2.0])
        weights=np.outer(wz,wr) #r is fastest varying in dataset
        weights=np.reshape(weights,dims[0]*dims[2])
        math_me=dsa.WrapDataObject(self.data)
        jacobian=1.0
        if rBool:
          weights*=math_me.Points[:,0] #multiply by R
          jacobian*=0.5*B[0]
        if zBool:
          jacobian*=0.5*B[2]
        weights=weights*jacobian #multiply by jacobian
        return weights
    
    def weighted_copy(self, rBool=True, zBool=True):
        new_data=vtkStructuredGrid()
        new_data.DeepCopy(self.data)
        math_new=dsa.WrapDataObject(new_data)
        w=self.weight_matrix(rBool=rBool,zBool=zBool)
        nFields=len(math_new.PointData.keys())
        for i in range(nFields):
            if(len(math_new.PointData[i].shape)>1):
                for j in range(math_new.PointData[i].shape[1]):
                    math_new.PointData[i][:,j]=math_new.PointData[i][:,j]*w
            else:
                math_new.PointData[i][:]=math_new.PointData[i][:]*w
        return MrVtkVector(new_data)
    
    def power(self,power):
        new_data=vtkStructuredGrid()
        new_data.DeepCopy(self.data)
        math_me=dsa.WrapDataObject(self.data)
        math_new=dsa.WrapDataObject(new_data)
        numFlds=len(math_me.PointData.keys())
        for i in range(numFlds):
            math_new.PointData[i][:]=math_me.PointData[i][:]**power
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