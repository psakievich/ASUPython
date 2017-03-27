import sys,os,time
sys.path.append(os.environ["HOME"]+"/Documents/workspace/ASUPython/VisitVis/")
import DocumentPOD as DP
print "\n",sys.argv
fName=sys.argv[1]
AddWindow()
if fName is not None:
  obj=DP.POD_Plot(fName)
  obj.OpenDB()
  obj.GeneratePlot1()
  obj.GeneratePlot2()
  time.sleep(10)
 # obj.ClearPlots()
 # obj.CloseDB()
 # sys.exit()