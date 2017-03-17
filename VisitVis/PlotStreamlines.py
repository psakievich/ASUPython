AddPlot('Streamline','velocity')
sa=StreamlineAttributes()
sa.sourceType=4      #plane source
sa.sampleDensity0=10 #x dir
sa.sampleDensity1=10 #y dir
sa.coloringMethod=6  #color by variable
sa.colorTableName='hot_desaturated'
sa.integrationDirection=2 #forward and backward'
sa.coloringVariable='temperature'
sa.displayMethod=1   #tubes
sa.showSeeds=0
sa.geomDisplayQuality=2
sa.sampleDistance0=3.15
sa.sampleDistance1=3.15
sa.issueTerminationWarnings=0
sa.issueStiffnessWarnings=0
sa.issueCriticalPointsWarnings=0
SetPlotOptions(sa)
DrawPlots()
