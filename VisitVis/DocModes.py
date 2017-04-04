import sys,os,subprocess
from mpi4py import MPI
sys.path.append(os.environ["HOME"]+"/Research/ASUPython/")
import ModeTransforms as MT
import numpy as np

comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()
script=os.environ["HOME"]+"/Research/ASUPython/VisitVis/ClientSave.py"

Modes=[0,1,2,3,4,5,6,7,8,9,10,20]
nModes=len(Modes)
nLoop=nModes//size
nMod=nModes%size
path=os.environ["HOME"]+"/Desktop/JFM/{}/"
myFiles=[]
myTotal=[]
for i in Modes:
  myPath=path.format(i)
  for file in os.listdir(myPath):
    if file.endswith(".vts") and int(file.strip(".vts").split('_')[3])<11:
      myFiles.append(file)
      myTotal.append(myPath+file)
def Process(i):
   MT.FourierToRealDoc(myTotal[i],32,"R"+myFiles[i],np.pi*0.5)
   subprocess.call(["visit","-cli","-nowin","-s",script,"R"+myFiles[i]])
nModes=len(myTotal)
nLoop=nModes//size
nMod=nModes%size

for i in range(nLoop):
  myMode=rank+i*size
  Process(myMode) 

if (rank<nMod):
  myMode=rank+size*nLoop
  Process(myMode)
