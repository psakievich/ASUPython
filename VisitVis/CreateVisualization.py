import os,sys,subprocess
sys.path.append(os.environ["HOME"]+"/Documents/workspace/ASUPython/")
import ModeTransforms as MT
import numpy as np

fName=sys.argv[1]
MT.FourierToRealDoc(fName,32,"temp.vts",np.pi*0.5)
visit="/Applications/VisIt.app/Contents/Resources/bin/visit"
subprocess.call([visit,"-cli","-s","ClientPlotting.py","temp.vts"])
os.remove("temp.vts")