import os,sys
sys.path.append(os.environ["HOME"]+"/Research/ASUPython/")
import modred as mr
import MrImaginaryVtk as MIV
import numpy as np

myMode=int(sys.argv[1])
numModes=int(sys.argv[2])

modes=[MIV.MrVtkVecHandle("{0}/POD_{0}_Mode_{1}.vts".format(myMode,i)) for i in range(numModes)]
vec_space = mr.VectorSpaceHandles(inner_product=MIV.inner_product)
IP_mat = vec_space.compute_symmetric_inner_product_mat(modes)
if not np.allclose(IP_mat, np.eye(len(modes))):
    print('Warning: modes are not orthonormal', IP_mat)
