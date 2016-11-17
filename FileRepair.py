#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 15:26:42 2016

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

def PackBrokenFileIntoNumpy(dataSet):
    len1=len(dataSet)
    len2=len(dataSet[0])
    arr=np.empty([len1,len2])
    for i in range(len1):
        arr[i,:]=np.array(dataSet[i].strip(),dtype=float)
    return arr