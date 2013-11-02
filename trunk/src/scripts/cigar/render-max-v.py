from visit import *
DeleteAllPlots()
ResetView()

wat = GetSaveWindowAttributes()
wat.fileName = 'binning-maxv'
wat.width = 1024
wat.height = 768
SetSaveWindowAttributes(wat)

AddPlot("Pseudocolor", "particle/V", 1, 1)
SetActivePlots(0)

AddOperator("DataBinning", 1)
DataBinningAtts = DataBinningAttributes()
DataBinningAtts.numDimensions = DataBinningAtts.Two  # One, Two, Three
DataBinningAtts.dim1BinBasedOn = DataBinningAtts.X  # X, Y, Z, Variable
DataBinningAtts.dim1Var = "default"
DataBinningAtts.dim1SpecifyRange = 0
DataBinningAtts.dim1MinRange = 0
DataBinningAtts.dim1MaxRange = 1
DataBinningAtts.dim1NumBins = 50
DataBinningAtts.dim2BinBasedOn = DataBinningAtts.Y  # X, Y, Z, Variable
DataBinningAtts.dim2Var = "default"
DataBinningAtts.dim2SpecifyRange = 0
DataBinningAtts.dim2MinRange = 0
DataBinningAtts.dim2MaxRange = 1
DataBinningAtts.dim2NumBins = 50
DataBinningAtts.dim3BinBasedOn = DataBinningAtts.Variable  # X, Y, Z, Variable
DataBinningAtts.dim3Var = "default"
DataBinningAtts.dim3SpecifyRange = 0
DataBinningAtts.dim3MinRange = 0
DataBinningAtts.dim3MaxRange = 1
DataBinningAtts.dim3NumBins = 50
DataBinningAtts.outOfBoundsBehavior = DataBinningAtts.Clamp  # Clamp, Discard
DataBinningAtts.reductionOperator = DataBinningAtts.Maximum  # Average, Minimum, Maximum, StandardDeviation, Variance, Sum, Count, RMS, PDF
DataBinningAtts.varForReduction = "particle/V"
DataBinningAtts.emptyVal = float('NaN')
DataBinningAtts.outputType = DataBinningAtts.OutputOnBins  # OutputOnBins, OutputOnInputMesh
DataBinningAtts.removeEmptyValFromCurve = 1
SetOperatorOptions(DataBinningAtts, 1)

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
banner.text = "max(V(x,y))  t=$time"
banner.position = (0.2, 0.93)
banner.fontBold = 1
banner.height = 0.02

DrawPlots()
SaveWindow()

banner.Delete()
DeleteAllPlots()
ResetView()
