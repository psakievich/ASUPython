#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 16:10:07 2016

TEST POD RUN

@author: psakievich
"""

import modred as mr
import MrImaginaryVtk as MIV

PATH='analysis/{}/'
FILEBASE='SymS1'
FILESTD='_{}_TSTEP_{}.vts'

NSNAPS=10
NINC=10
WAVE=3

inHandles=[MIV.MrVtkVecHandle(PATH.format(WAVE)+FILEBASE+FILESTD.format(WAVE,i)) \
           for i in range(0,NSNAPS*NINC,NINC)]

outHandles=[MIV.MrVtkVecHandle(PATH.format(WAVE)+'PODMODE_{}.vts'.format(i)) \
            for i in range(10)]

myPOD=mr.PODHandles(MIV.inner_product)
myPOD.compute_decomp(inHandles)
myPOD.compute_modes(range(NSNAPS),outHandles)


