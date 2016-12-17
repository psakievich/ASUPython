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
    def __init__(self,dataType=np.float64):
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
        Info is a tuple of data typically (s_rank,d_rank,message size)
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
    def SendRecv(self,size,rank,comm,dim):
        offset=size/2**dim
        status=MPI.Status
        if rank<offset:
            comm.send(self.oInfoBuff,dest=rank+offset,tag=rank)
            temp=comm.recv(source=rank+offset,tag=rank+offset)
            comm.probe(status)
            
