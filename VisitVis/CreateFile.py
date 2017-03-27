import os,sys
path=sys.argv[1]
fIn=sys.argv[2]
fOut=sys.argv[3]
sys.path.append(path)
import ModeTransforms as MT
import numpy as np
MT.FourierToRealDoc(fIn,32,fOut,np.pi*0.5)