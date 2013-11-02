from visit import *
DeleteAllPlots()

wat = GetSaveWindowAttributes()
wat.fileName = 'scatter-bv-'
wat.width = 1024
wat.height = 768
SetSaveWindowAttributes(wat)

AddPlot("Scatter", "particle/r", 1, 0)
SetActivePlots(0)
ScatterAtts = ScatterAttributes()
ScatterAtts.var1 = "particle/B"
ScatterAtts.var1Role = ScatterAtts.Coordinate0  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var1MinFlag = 0
ScatterAtts.var1MaxFlag = 0
ScatterAtts.var1Min = 0
ScatterAtts.var1Max = 1
ScatterAtts.var1Scaling = ScatterAtts.Linear  # Linear, Log, Skew
ScatterAtts.var1SkewFactor = 1
ScatterAtts.var2Role = ScatterAtts.Coordinate1  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var2 = "particle/V"
ScatterAtts.var2MinFlag = 0
ScatterAtts.var2MaxFlag = 0
ScatterAtts.var2Min = 0
ScatterAtts.var2Max = 1
ScatterAtts.var2Scaling = ScatterAtts.Linear  # Linear, Log, Skew
ScatterAtts.var2SkewFactor = 1
ScatterAtts.var3Role = ScatterAtts.None  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var3 = "default"
ScatterAtts.var3MinFlag = 0
ScatterAtts.var3MaxFlag = 0
ScatterAtts.var3Min = 0
ScatterAtts.var3Max = 1
ScatterAtts.var3Scaling = ScatterAtts.Linear  # Linear, Log, Skew
ScatterAtts.var3SkewFactor = 1
ScatterAtts.var4Role = ScatterAtts.Color  # Coordinate0, Coordinate1, Coordinate2, Color, None
ScatterAtts.var4 = "particle/V"
ScatterAtts.var4MinFlag = 0
ScatterAtts.var4MaxFlag = 0
ScatterAtts.var4Min = 0
ScatterAtts.var4Max = 1
ScatterAtts.var4Scaling = ScatterAtts.Linear  # Linear, Log, Skew
ScatterAtts.var4SkewFactor = 1
ScatterAtts.pointSize = 0.05
ScatterAtts.pointSizePixels = 5
ScatterAtts.pointType = ScatterAtts.Sphere  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
ScatterAtts.scaleCube = 1
ScatterAtts.colorType = ScatterAtts.ColorByColorTable  # ColorByForegroundColor, ColorBySingleColor, ColorByColorTable
ScatterAtts.singleColor = (255, 0, 0, 255)
ScatterAtts.colorTableName = "Default"
ScatterAtts.invertColorTable = 0
ScatterAtts.legendFlag = 1
SetPlotOptions(ScatterAtts)

aat = AnnotationAttributes()
aat.userInfoFlag = 0
aat.databaseInfoFlag = 0
aat.timeInfoFlag = 0
SetAnnotationAttributes(aat)

plotName = GetPlotList().GetPlots(0).plotName
legend = GetAnnotationObject(plotName)
legend.orientation = legend.VerticalLeft  # VerticalRight, VerticalLeft, HorizontalTop, Horiz
legend.fontBold = 1

banner = CreateAnnotationObject("Text2D")
banner.text = "scatter b vs. v t=$time"
banner.position = (0.2, 0.93)
banner.fontBold = 1
banner.height = 0.02

DrawPlots()
SaveWindow()

banner.Delete()
DeleteAllPlots()
ResetView()
