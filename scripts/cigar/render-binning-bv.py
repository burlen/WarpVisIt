from visit import *
DeleteAllPlots()
ResetView()

wat = GetSaveWindowAttributes()
wat.fileName = 'binning-bv'
wat.width = 1024
wat.height = 768
SetSaveWindowAttributes(wat)

AddPlot("Pseudocolor", "particle/vz", 1, 1)
AddOperator("DataBinning", 1)
SetActivePlots(0)
DataBinningAtts = DataBinningAttributes()
DataBinningAtts.numDimensions = DataBinningAtts.Two  # One, Two, Three
DataBinningAtts.dim1BinBasedOn = DataBinningAtts.Variable  # X, Y, Z, Variable
DataBinningAtts.dim1Var = "particle/B"
DataBinningAtts.dim1SpecifyRange = 0
DataBinningAtts.dim1MinRange = 0
DataBinningAtts.dim1MaxRange = 1
DataBinningAtts.dim1NumBins = 50
DataBinningAtts.dim2BinBasedOn = DataBinningAtts.Variable  # X, Y, Z, Variable
DataBinningAtts.dim2Var = "particle/V"
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
DataBinningAtts.reductionOperator = DataBinningAtts.Count  # Average, Minimum, Maximum, StandardDeviation, Variance, Sum, Count, RMS, PDF
DataBinningAtts.varForReduction = "particle/vz"
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
banner.text = "number of v vs. b  t=$time"
banner.position = (0.2, 0.93)
banner.fontBold = 1
banner.height = 0.02

DrawPlots()
SaveWindow()

banner.Delete()
DeleteAllPlots()
ResetView()
