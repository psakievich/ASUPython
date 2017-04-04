import sys,os
sys.path.append(os.environ["HOME"]+"/Research/ASUPython/VisitVis/")
import LatexDoc as LD

f=LD.TeXFile("MyDictionary.tex")
f.AddPackage("graphicx")
f.AddPackage("auto-pst-pdf")
f.AddPackage("nopageno")
f.BeginDocument()
PicTemp="DocK_{}_m_{}"
CapTemp="k={} m={}:{}"
for i in range(11):
   f.AddFigure(PicTemp.format(1,i))
for i in range(2):
   f.AddFigure(PicTemp.format(2,i))
for i in range(9):
   f.AddFigure(PicTemp.format(3,i))
for i in range(11):
   f.AddFigure(PicTemp.format(4,i))
f.Write("\clearpage")
for i in range(11):
   f.AddFigure(PicTemp.format(5,i))
for i in range(11):
   f.AddFigure(PicTemp.format(6,i))
for i in range(11):
   f.AddFigure(PicTemp.format(7,i))
for i in range(11):
   f.AddFigure(PicTemp.format(8,i))
f.Write("\clearpage")
for i in range(11):
   f.AddFigure(PicTemp.format(9,i))
for i in range(11):
   f.AddFigure(PicTemp.format(10,i))
for i in range(11):
   f.AddFigure(PicTemp.format(20,i))


f.EndDocument()
