'''
We will generate a class for creating latex files
that document a large number of pictures
'''

class TeXFile:
  def __init__(self,texFileName):
     self.myFile=open(texFileName,'w')
     w=self.myFile
     w.write("\documentclass[12pt,letterpaper]{article}\n")
  def Write(self,text):
     self.myFile.write(text+'\n')
  def AddPackage(self,packageName):
     self.Write(r'\usepackage{{{}}}'.format(packageName))
  def StartDocument(self):
     self.Write(r'\begin{{document}}')
  def AddFigure(self,imageName,caption=None,sizing=r'width=0.9\textwidth'):
     self.Write(r'\begin{{figure}}')
     self.Write(r'\centering')
     self.Write(r'[{}]'.format(sizing)+r'{{{}}}'.format(imageName))
     if caption is not None:
        self.Write(r'\caption{{{}}}'.format(caption))
     self.Write(r'\end{{figure}}')
  def EndDocument(self)
     self.Write(r'\end{{document}}')

