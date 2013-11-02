from visit import *
DeleteAllPlots()

wat = GetSaveWindowAttributes()
wat.fileName = 'cigar-v-of-x-y-'
SetSaveWindowAttributes(wat)

AddPlot("Pseudocolor", "particle/V", 1, 1)
DrawPlots()
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

DrawPlots()
SaveWindow()
