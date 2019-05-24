'''
Batch process the plots from POD modes
'''

from visit import *

dbTemp='{0}/POD_{0}_Mode_*.vts database'
picTemp='TotalEnergyPOD_{0}_Mode'

class Mode:
  def __init__(self,Database,PicFamily):
    self.db=Database
    self.picTemp=PicFamily
  def OpenDB(self):
    OpenDatabase(self.db)
  def CloseDB(self):
    CloseDatabase(self.db)
  def GenerateVar(self):
    DefineScalarExpression("TotalEnergy", \
     "R_velocity[0]^2+C_velocity[0]^2+"
     "R_velocity[1]^2+C_velocity[1]^2+"
     "R_velocity[2]^2+C_velocity[2]^2+"
     "R_temperature[0]^2+C_temperature[0]^2")
  def SetView(self):
    v=GetView3D()
    v.viewNormal=(0,-1,0)
    v.viewUp=(0,0,1)
    v.imageZoom=0.75
    SetView3D(v)
  def SetAnnotations(self):
    a=AnnotationAttributes()
    a.axes3D.visible=0
    a.axes3D.triadFlag=0
    a.axes3D.bboxFlag=0
    a.userInfoFlag=0
    a.databaseInfoFlag=0   
    SetAnnotationAttributes(a)
  def GeneratePlots(self,nWaves=None):
    AddPlot("Pseudocolor","TotalEnergy")
    DrawPlots()
    p=PseudocolorAttributes()
    p.colorTableName="hot_desaturated"
    p.minFlag=True
    p.maxFlag=True
    SetPlotOptions(p)
    plotName=GetPlotList().GetPlots(0).plotName
    legend=GetAnnotationObject(plotName)
    legend.managePosition=0
    legend.position=(0.01,0.67)
    legend.fontBold=1
    legend.fontHeight=0.02
    legend.yScale=1.25
    s=GetSaveWindowAttributes()
    s.fileName=self.picTemp
    s.format=s.POSTSCRIPT #post script
    s.outputDirectory="RZEnergy"
    s.outputToCurrentDirectory=0
    SetSaveWindowAttributes(s)
    if nWaves is not None:
       loopVals=min(TimeSliderGetNStates(),nWaves)
    else:
       loopVals=TimeSliderGetNStates()
   
    for state in range(loopVals):
      SetTimeSliderState(state)
      SaveWindow()
    DeleteAllPlots()
  def Initialize(self):
    self.OpenDB()
    self.GenerateVar()
    self.SetView()
    self.SetAnnotations()

def ProcessBatch(myModes,myPics,waveLim=None):
   for i in range(len(myModes)):
      m=Mode(myModes[i],myPics[i])
      m.OpenDB()
      m.Initialize()
      m.GeneratePlots(nWaves=waveLim)
      m.CloseDB()
    
  
