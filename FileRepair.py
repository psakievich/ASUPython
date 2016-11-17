#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 15:26:42 2016

@author: psakievich
"""
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