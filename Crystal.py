#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 21:00:40 2016
Crystal Router
@author: psakievich
"""
from mpi4py import MPI
import numpy as np

class CrystalRouter():
    '''
    All-to-all message passing interface
    values must be numpy arrays
    '''
    def __init__(self,MpiComm,dataType=np.float64):
        self.comm=MpiComm
        self.rank=self.comm.Get_rank()
        self.size=self.comm.Get_size()
        self.sizeHyper=int(np.ceil(np.log2(self.size)))
        self.dataType=dataType
        self.iValBuff=np.empty(0,dtype=self.dataType)
        self.oValBuff=np.empty(0,dtype=self.dataType)
        self.nValBuff=np.empty(0,dtype=self.dataType)
        self.iInfoBuff=[]
        self.nInfoBuff=[]
        self.oInfoBuff=[]
    def AddMessage(self,Value,Info):
        '''
        Add message to the out buffer
        Value is a 1-d numpy array
        Info is a tuple of data typically (d_rank,message size)
        '''
        np.append(self.oValBuff,Value)
        self.oInfoBuff.append(Info)
    def RetrieveMessages(self):
        '''
        Return message information meant for the local rank
        '''
        return self.iValBuff, self.iInfoBuff
    def ResetBuffers(self,dataType=np.float64):
        '''
        Blank out all buffers and possibly reset numpy data type
        '''
        self.__init__(dataType)
    def GetHyperCubePartners(self,rank,size,level):
        '''
        Get Communication partner for n-level of crystal router communication
        n is defined in range 1:log_2(sizeLocal) where sizeLocal is 2^n
        '''
        # round up to next largest hypercube size if not perfect hypercube
        sizeLocal=2**np.ceil(np.log2(size))
        
        def SetOffset(size,level):
            return size//2**level
        offset=SetOffset(sizeLocal,level)
        shift=rank//SetOffset(sizeLocal,level-1)*SetOffset(sizeLocal,level-1)
        if rank<offset+shift:
            partner = int(rank+offset)
        else:
            partner = int(rank-offset)
        if partner < size:
            return partner
        else:
            return None
            
        
    def SendRecv(self,size,rank,dim):
        partner=self.GetHyperCubePartners(rank,size,dim)
        state=MPI.Status
        if partner is not None:
            if rank<partner:
                #info send and receive
                self.comm.send(self.oInfoBuff,dest=partner,tag=rank)
                tempInfo=self.comm.recv(source=partner,tag=partner)
                #val send and receives
                self.comm.Send(self.oValBuff,dest=partner,tag=rank)
                #use probe to get size of numpy array to be received
                self.comm.Probe(source=partner,tag=partner,status=state)
                tempVal=np.empty(state.Get_elements(self.dataType),dtype=self.dataType)
                self.comm.Recv(tempVal,source=partner,tag=partner)
                
            else:
                #info receive and send
                tempInfo=self.comm.recv(source=partner,tag=partner)
                self.comm.send(self.oInfoBuff,dest=partner,tag=rank)
                #val receive and send
                self.comm.Probe(source=partner,tag=partner,status=state)
                tempVal=np.empty(state.Get_elements(self.dataType),dtype=self.dataType)
                self.comm.Recv(tempVal,source=partner,tag=partner)
                self.comm.Send(self.oValBuff,dest=partner,tag=rank)
                
            
