import sys,os,time
sys.path.append(os.environ["HOME"]+"/Research/ASUPython/VisitVis/")
import DocumentPOD as DP
fName=sys.argv[1]
AddWindow()
if fName is not None:
  obj=DP.POD_Plot(fName)
  obj.OpenDB()
  obj.GeneratePlot1()
  obj.GeneratePlot2()
  obj.Save()
  obj.ClearPlots()
  obj.CloseDB()
  sys.exit()
