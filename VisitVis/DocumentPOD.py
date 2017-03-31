'''
Batch process POD modes 
'''
from visit import *
class POD_Plot:
   def __init__(self,fName):
      self.fileName=fName
      self.label=None
   def OpenDB(self):
      OpenDatabase(self.fileName)
   def SetWindowProperties(self):
      a=AnnotationAttributes()
      a.axes3D.visible=0
      a.axes3D.bboxFlag=0
      a.userInfoFlag=0
      a.databaseInfoFlag=0
      SetAnnotationAttributes(a)
   def GetLegend(self):
      plotName=GetPlotList().GetPlots(0).plotName
      legend=GetAnnotationObject(plotName)
      legend.fontBold=1
      legend.fontHeight=0.02
      legend.yScale=1.25
      legend.xScale=1.55
      return legend 
   def GeneratePlot2(self):
      SetActiveWindow(2)
      self.OpenDB()
      AddPlot("Pseudocolor","temperature")
      self.SetWindowProperties()
      p=PseudocolorAttributes()
      p.colorTableName="hot_desaturated"
      SetPlotOptions(p)
      v=GetView3D()
      v.viewNormal=(0,-0.9,0.5)
      v.centerOfRotation = (1.57491, 1.57491, 0)
      legend=self.GetLegend()
      legend.orientation=legend.HorizontalBottom
      legend.managePosition=0
      legend.position=(0.3,0.25)
      SetView3D(v)
      DrawPlots()
   def GeneratePlot1(self):
      SetActiveWindow(1)
      self.OpenDB()
      #stream line plots
      AddPlot('Streamline','velocity')
      self.SetWindowProperties()
      sa=StreamlineAttributes()
      sa.sourceType=sa.SpecifiedBox
      sa.useWholeBox=0
      sa.boxExtents=(0,3.15,0,3.15,-0.3,0.3)
      sa.sampleDensity0=15
      sa.sampleDensity1=15
      sa.sampleDensity2=3
      sa.coloringMethod=sa.ColorByVariable
      sa.coloringVariable='temperature'
      sa.colorTableName='hot_desaturated'
      sa.displayMethod=sa.Tubes
      sa.geomDisplayQuality=sa.Super
      sa.showSeeds=0
      sa.integrationDirection=sa.Both
      sa.issueTerminationWarnings=0
      sa.issueStiffnessWarnings=0
      sa.issueCriticalPointsWarnings=0
      SetPlotOptions(sa)
      #view
      v=GetView3D()
      v.viewNormal=(0.762321, 0.466435, 0.44867)
      v.viewUp=(-0.341935, -0.250559, 0.905704)
      SetView3D(v)
      #setup legend
      legend=self.GetLegend()
      legend.orientation=legend.HorizontalBottom
      legend.managePosition=0
      legend.position=(0.3,0.25)
      #setup label
      self.label=CreateAnnotationObject("Text2D")
      self.label.position=(0.35,0.8)
      self.label.height=0.04
      self.label.fontBold=1
      temp=self.fileName.strip(".vts").split('_')
      self.label.text=r'k={} m={}'.format(temp[1],int(temp[3])+1)
      DrawPlots()
   def ClearPlots(self):
      #self.label.Delete()
      SetActiveWindow(1)
      DeleteAllPlots()
      SetActiveWindow(2)
      DeleteAllPlots()
   def Save(self):
      s=SaveWindowAttributes()
      temp=self.fileName.strip(".vts").split('_')
      s.fileName='DocK_{}_m_{}'.format(temp[1],temp[3])
      s.family=0
      s.format=s.POSTSCRIPT
      s.saveTiled=1
      SetSaveWindowAttributes(s)
      SaveWindow()
   def CloseDB(self):
      SetActiveWindow(1)
      CloseDatabase(self.fileName)
      SetActiveWindow(2)
      CloseDatabase(self.fileName)
