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
    Main routines are AddMessage and Comminicate
    '''
    def __init__(self,MpiComm,dataType=np.float64):
        self.comm=MpiComm
        self.rank=self.comm.Get_rank()
        self.size=self.comm.Get_size()
        self.sizeHyper=int(np.ceil(np.log2(self.size)))
        self.level=1
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
        if self.IsNext(self.rank,Info[0],self.size,1):
            np.append(self.nValBuff,Value)
            self.nInfoBuff.append(Info)
        else:
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
            
    def IsNext(self,rank,partner,size,level):
        # round up to next largest hypercube size if not perfect hypercube
        sizeLocal=2**np.ceil(np.log2(size))
        
        def SetOffset(size,level):
            return size//2**level
        offset=SetOffset(sizeLocal,level)
        shift=rank//SetOffset(sizeLocal,level-1)*SetOffset(sizeLocal,level-1)
        if rank<offset+shift:
            if partner>=offset+shift:
                return True
            else:
                return False
        else:
            if partner<offset+shift:
                return True
            else:
                return False
    def SortBuffers(self,size,rank,dim):
        tempVal=np.empty(0,dtype=self.dataType)
        tempInfo=[]
        while len(self.oInfoBuff)>0:
            partner=self.oInfoBuff[0][0]
            length=self.oInfoBuff[0][1]
            if self.IsNext(rank,partner,size,dim):
                # if message will be transfered next
                np.append(self.nValBuff,self.oValBuff[0:length])
                self.nInfoBuff.append(self.oInfoBuff[0])
            else:
                # if needs to stay in out buffer move to temporary buffer
                np.append(tempVal,self.oInfoBuff[0:length])
                tempInfo.append(self.oInfoBuff[0])
            #clear the message to continue marching through the buffer    
            del self.oInfoBuff[0]
            np.delete(self.oValBuff,range(length))
        #assign left overs back to buffers
        self.oInfoBuff=tempInfo
        self.oValBuff=tempVal
        
    def SendRecv(self,size,rank,dim):
        partner=self.GetHyperCubePartners(rank,size,dim)
        state=MPI.Status
        #communicate
        if partner is not None:
            if rank<partner:
                #info send and receive
                self.comm.send(self.nInfoBuff,dest=partner,tag=rank)
                tempInfo=self.comm.recv(source=partner,tag=partner)
                #val send and receives
                self.comm.Send(self.nValBuff,dest=partner,tag=rank)
                #use probe to get size of numpy array to be received
                self.comm.Probe(source=partner,tag=partner,status=state)
                tempVal=np.empty(state.Get_elements(self.dataType),dtype=self.dataType)
                self.comm.Recv(tempVal,source=partner,tag=partner)
                
            else:
                #info receive and send
                tempInfo=self.comm.recv(source=partner,tag=partner)
                self.comm.send(self.nInfoBuff,dest=partner,tag=rank)
                #val receive and send
                self.comm.Probe(source=partner,tag=partner,status=state)
                tempVal=np.empty(state.Get_elements(self.dataType),dtype=self.dataType)
                self.comm.Recv(tempVal,source=partner,tag=partner)
                self.comm.Send(self.nValBuff,dest=partner,tag=rank)
            #clear nextBuffers
            del self.nInfoBuff[::]
            np.delete(self.nValBuff,range(len(self.nValBuff)))
            #fill buffers for next round of communication
            while len(tempInfo)>0:
                endPnt=tempInfo[0][1]
                nPartner=tempInfo[0][0]
                if nPartner==rank:
                    # if message destination is this rank append to inBuffers
                    np.append(self.iValBuff,tempVal[0:endPnt])
                    self.iInfoBuff.append(tempInfo[0])
                elif self.IsNext(rank,nPartner,size,dim+1):
                    # if message will be transfered in the next communication
                    np.append(self.nValBuff,tempVal[0:endPnt])
                    self.nInfoBuff.append(tempInfo[0])
                else:
                    # if message needs to be stored during next communication
                    np.append(self.oValBuff,tempVal[0:endPnt])
                    self.oInfoBuff.append(tempInfo[0])
                #delete transfered data from temp buffers
                del tempInfo[0]
                np.delete(tempVal,range(endPnt))
    
    def Communicate(self):
        '''
        March over the length of the hypercube and perform the send recev
        operations
        '''
        for i in range(self.sizeHyper):
            self.SortBuffers(self.size,self.rank,i)
            self.SendRecv(self.size,self.rank,i)
                    
                    
                
                
            
