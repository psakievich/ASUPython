#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 15:26:42 2016

File repair and plotting from the fracturing that 
Stampede Fortran modules create when writing ASCII
data 

@author: psakievich
"""
import numpy as np

def ReadBrokenFile(fName,catEveryXLines):
    dataSet=[]
    with open(fName) as f:
        i=0
        temp=''
        for line in f:
            line=line.strip()
            temp=temp+' '+line
            i=i+1
            if(i%catEveryXLines==0):
                dataSet.append(temp)
                temp=''
    return dataSet

def PackBrokenFileIntoNumpy(dataSet,lineLen):
    len1=len(dataSet)
    arr=np.empty([len1,lineLen])
    for i in range(len1):
        arr[i,:]=np.array(dataSet[i].split(),dtype=float)
    return arr