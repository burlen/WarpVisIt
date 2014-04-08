from visit import *

#----------------------------------------------------------------------------
def meshPlot(winId,meshName):
    """VisIt Du gehst mir auf die Eier."""
    SetActiveWindow(winId)
    AddPlot('Mesh',meshName)
    plotId = GetNumPlots()-1
    clearAnnotations(winId, plotId)
    DrawPlots()
    return plotId

#----------------------------------------------------------------------------
def meshValid(winId,meshName):
    """VisIt Du gehst mir auf die Eier."""
    SetActiveWindow(winId)
    plotId = meshPlot(winId,meshName)
    SetActivePlots(plotId)
    SetQueryFloatFormat("%g")
    SetQueryOutputToValue()
    numNodes = Query("NumNodes")
    sys.stderr.write('%s has %g nodes\n'%(meshName,numNodes))
    SetActivePlots(plotId)
    DeleteActivePlots()
    return numNodes>0

#----------------------------------------------------------------------------
def deletePlots(winId):
    SetActiveWindow(winId)
    DeleteAllPlots()
    return

#----------------------------------------------------------------------------
def clearWindow(winId):
    """Hari suggested this, doesn't work"""
    SetActiveWindow(winId)
    InvertBackgroundColor()
    InvertBackgroundColor()
    return

#----------------------------------------------------------------------------
def plotEmpty(winId, plotId):
    """check for empty plot by query num nodes"""

    SetActiveWindow(winId)
    SetActivePlots(plotId)

    SetQueryFloatFormat("%g")
    SetQueryOutputToValue()
    numNodes = Query("NumNodes")

    sys.stderr.write('NumNodes=%g'%(numNodes))

    return numNodes<1

#----------------------------------------------------------------------------
def setView(winId, plotId):

    SetActiveWindow(winId)
    SetActivePlots(plotId)

    ResetView()

    View3DAtts = GetView3D()
    View3DAtts.viewNormal = (0.82, 0.53, 0.25)
    View3DAtts.viewUp = (-0.19, -0.16, 0.97)
    SetView3D(View3DAtts)

    return

#----------------------------------------------------------------------------
def plotPColorElec0Uz(winId):

    SetActiveWindow(winId)

    AddPlot("Pseudocolor", "Electron-0/uz", 1, 0)
    plotId = GetNumPlots() - 1
    SetActivePlots(plotId)

    AddOperator("Project", 0)
    ProjectAtts = ProjectAttributes()
    ProjectAtts.projectionType = ProjectAtts.XZCartesian  # ZYCartesian, XZCartesian, XYCartesian, XRCylindrical, YRCylindrical, ZRCylindrical
    ProjectAtts.vectorTransformMethod = ProjectAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
    SetOperatorOptions(ProjectAtts, 1)

    AddOperator("Elevate", 0)
    ElevateAtts = ElevateAttributes()
    ElevateAtts.useXYLimits = 1
    ElevateAtts.limitsMode = ElevateAtts.OriginalData  # OriginalData, CurrentPlot
    ElevateAtts.scaling = ElevateAtts.Linear  # Linear, Log, Skew
    ElevateAtts.skewFactor = 1
    ElevateAtts.minFlag = 1
    ElevateAtts.min = -1e+09
    ElevateAtts.maxFlag = 1
    ElevateAtts.max = 1e+09
    ElevateAtts.zeroFlag = 0
    ElevateAtts.variable = "Electron-0/uy"
    SetOperatorOptions(ElevateAtts, 1)

    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
    PseudocolorAtts.skewFactor = 1
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = -7e+08
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = 7e+08
    PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
    PseudocolorAtts.colorTableName = "hot_and_cold"
    PseudocolorAtts.invertColorTable = 0
    PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
    PseudocolorAtts.opacityVariable = ""
    PseudocolorAtts.opacity = 1
    PseudocolorAtts.opacityVarMin = 0
    PseudocolorAtts.opacityVarMax = 1
    PseudocolorAtts.opacityVarMinFlag = 0
    PseudocolorAtts.opacityVarMaxFlag = 0
    PseudocolorAtts.pointSize = 0.05
    PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    PseudocolorAtts.pointSizeVarEnabled = 0
    PseudocolorAtts.pointSizeVar = "default"
    PseudocolorAtts.pointSizePixels = 2
    PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
    PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    PseudocolorAtts.lineWidth = 0
    PseudocolorAtts.tubeDisplayDensity = 10
    PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
    PseudocolorAtts.tubeRadiusAbsolute = 0.125
    PseudocolorAtts.tubeRadiusBBox = 0.005
    PseudocolorAtts.varyTubeRadius = 0
    PseudocolorAtts.varyTubeRadiusVariable = ""
    PseudocolorAtts.varyTubeRadiusFactor = 10
    PseudocolorAtts.endPointType = PseudocolorAtts.None  # None, Tails, Heads, Both
    PseudocolorAtts.endPointStyle = PseudocolorAtts.Spheres  # Spheres, Cones
    PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
    PseudocolorAtts.endPointRadiusAbsolute = 1
    PseudocolorAtts.endPointRadiusBBox = 0.005
    PseudocolorAtts.endPointRatio = 2
    PseudocolorAtts.renderSurfaces = 1
    PseudocolorAtts.renderWireframe = 0
    PseudocolorAtts.renderPoints = 0
    PseudocolorAtts.smoothingLevel = 0
    PseudocolorAtts.legendFlag = 1
    PseudocolorAtts.lightingFlag = 1
    SetPlotOptions(PseudocolorAtts)

    setAnnotations(winId, plotId, xAxisName='X', yAxisName='Z', zAxisName='Uy', showDB=1)

    DrawPlots()
    setView(winId, plotId)

    return plotId

#----------------------------------------------------------------------------
def plotPColorElec0Uy(winId):

    SetActiveWindow(winId)

    AddPlot("Pseudocolor", "Electron-0/uy", 1, 0)
    plotId = GetNumPlots() - 1
    SetActivePlots(plotId)

    AddOperator("Project", 0)
    ProjectAtts = ProjectAttributes()
    ProjectAtts.projectionType = ProjectAtts.XZCartesian  # ZYCartesian, XZCartesian, XYCartesian, XRCylindrical, YRCylindrical, ZRCylindrical
    ProjectAtts.vectorTransformMethod = ProjectAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
    SetOperatorOptions(ProjectAtts, 1)

    AddOperator("Elevate", 0)
    ElevateAtts = ElevateAttributes()
    ElevateAtts.useXYLimits = 1
    ElevateAtts.limitsMode = ElevateAtts.OriginalData  # OriginalData, CurrentPlot
    ElevateAtts.scaling = ElevateAtts.Linear  # Linear, Log, Skew
    ElevateAtts.skewFactor = 1
    ElevateAtts.minFlag = 1
    ElevateAtts.min = -1e+09
    ElevateAtts.maxFlag = 1
    ElevateAtts.max = 1e+09
    ElevateAtts.zeroFlag = 0
    ElevateAtts.variable = "Electron-0/uz"
    SetOperatorOptions(ElevateAtts, 1)

    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
    PseudocolorAtts.skewFactor = 1
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = -7e+08
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = 7e+08
    PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
    PseudocolorAtts.colorTableName = "hot_and_cold"
    PseudocolorAtts.invertColorTable = 0
    PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
    PseudocolorAtts.opacityVariable = ""
    PseudocolorAtts.opacity = 1
    PseudocolorAtts.opacityVarMin = 0
    PseudocolorAtts.opacityVarMax = 1
    PseudocolorAtts.opacityVarMinFlag = 0
    PseudocolorAtts.opacityVarMaxFlag = 0
    PseudocolorAtts.pointSize = 0.05
    PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    PseudocolorAtts.pointSizeVarEnabled = 0
    PseudocolorAtts.pointSizeVar = "default"
    PseudocolorAtts.pointSizePixels = 2
    PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
    PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    PseudocolorAtts.lineWidth = 0
    PseudocolorAtts.tubeDisplayDensity = 10
    PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
    PseudocolorAtts.tubeRadiusAbsolute = 0.125
    PseudocolorAtts.tubeRadiusBBox = 0.005
    PseudocolorAtts.varyTubeRadius = 0
    PseudocolorAtts.varyTubeRadiusVariable = ""
    PseudocolorAtts.varyTubeRadiusFactor = 10
    PseudocolorAtts.endPointType = PseudocolorAtts.None  # None, Tails, Heads, Both
    PseudocolorAtts.endPointStyle = PseudocolorAtts.Spheres  # Spheres, Cones
    PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
    PseudocolorAtts.endPointRadiusAbsolute = 1
    PseudocolorAtts.endPointRadiusBBox = 0.005
    PseudocolorAtts.endPointRatio = 2
    PseudocolorAtts.renderSurfaces = 1
    PseudocolorAtts.renderWireframe = 0
    PseudocolorAtts.renderPoints = 0
    PseudocolorAtts.smoothingLevel = 0
    PseudocolorAtts.legendFlag = 1
    PseudocolorAtts.lightingFlag = 1
    SetPlotOptions(PseudocolorAtts)

    setAnnotations(winId, plotId, xAxisName='X', yAxisName='Z', zAxisName='Uz', showDB=1)

    DrawPlots()
    setView(winId, plotId)

    return plotId

#----------------------------------------------------------------------------
def plotPColorElec1Uz(winId):

    SetActiveWindow(winId)

    AddPlot("Pseudocolor", "Electron-1/uz", 1, 0)
    plotId = GetNumPlots() - 1
    SetActivePlots(plotId)

    AddOperator("Project", 0)
    ProjectAtts = ProjectAttributes()
    ProjectAtts.projectionType = ProjectAtts.XZCartesian  # ZYCartesian, XZCartesian, XYCartesian, XRCylindrical, YRCylindrical, ZRCylindrical
    ProjectAtts.vectorTransformMethod = ProjectAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
    SetOperatorOptions(ProjectAtts, 0)

    AddOperator("Elevate", 0)
    ElevateAtts = ElevateAttributes()
    ElevateAtts.useXYLimits = 1
    ElevateAtts.limitsMode = ElevateAtts.OriginalData  # OriginalData, CurrentPlot
    ElevateAtts.scaling = ElevateAtts.Linear  # Linear, Log, Skew
    ElevateAtts.skewFactor = 1
    ElevateAtts.minFlag = 0
    ElevateAtts.min = -5e+08
    ElevateAtts.maxFlag = 0
    ElevateAtts.max = 5e+08
    ElevateAtts.zeroFlag = 0
    ElevateAtts.variable = "Electron-1/ux"
    SetOperatorOptions(ElevateAtts, 0)

    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
    PseudocolorAtts.skewFactor = 1
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
    PseudocolorAtts.minFlag = 0
    PseudocolorAtts.min = -5e+08
    PseudocolorAtts.maxFlag = 0
    PseudocolorAtts.max = 5e+08
    PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
    PseudocolorAtts.colorTableName = "Spectral"
    PseudocolorAtts.invertColorTable = 1
    PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
    PseudocolorAtts.opacityVariable = ""
    PseudocolorAtts.opacity = 1
    PseudocolorAtts.opacityVarMin = 0
    PseudocolorAtts.opacityVarMax = 1
    PseudocolorAtts.opacityVarMinFlag = 0
    PseudocolorAtts.opacityVarMaxFlag = 0
    PseudocolorAtts.pointSize = 0.05
    PseudocolorAtts.pointType = PseudocolorAtts.Sphere  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    PseudocolorAtts.pointSizeVarEnabled = 0
    PseudocolorAtts.pointSizeVar = "default"
    PseudocolorAtts.pointSizePixels = 10
    PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
    PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    PseudocolorAtts.lineWidth = 0
    PseudocolorAtts.tubeDisplayDensity = 10
    PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
    PseudocolorAtts.tubeRadiusAbsolute = 0.125
    PseudocolorAtts.tubeRadiusBBox = 0.005
    PseudocolorAtts.varyTubeRadius = 0
    PseudocolorAtts.varyTubeRadiusVariable = ""
    PseudocolorAtts.varyTubeRadiusFactor = 10
    PseudocolorAtts.endPointType = PseudocolorAtts.None  # None, Tails, Heads, Both
    PseudocolorAtts.endPointStyle = PseudocolorAtts.Spheres  # Spheres, Cones
    PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
    PseudocolorAtts.endPointRadiusAbsolute = 1
    PseudocolorAtts.endPointRadiusBBox = 0.005
    PseudocolorAtts.endPointRatio = 2
    PseudocolorAtts.renderSurfaces = 1
    PseudocolorAtts.renderWireframe = 0
    PseudocolorAtts.renderPoints = 0
    PseudocolorAtts.smoothingLevel = 0
    PseudocolorAtts.legendFlag = 1
    PseudocolorAtts.lightingFlag = 1
    SetPlotOptions(PseudocolorAtts)

    setAnnotations(winId, plotId, xAxisName='X', yAxisName='Z', zAxisName='Ux', showDB=1)

    DrawPlots()
    setView(winId, plotId)

    return plotId

#----------------------------------------------------------------------------
def plotBinElec0Weights(winId, nbx, nbz):

    SetActiveWindow(winId)

    # VisIt complains if we use a variable that doesn't exist
    # plot the pcolor on a variable that does, fix it later.
    #AddPlot("Pseudocolor", "operators/DataBinning/2D/Electron-0", 1, 1)
    AddPlot("Pseudocolor", "Electron-0/weights", 1, 0)
    plotId = GetNumPlots() - 1
    SetActivePlots(plotId)

    DrawPlots()
    setView(winId, plotId)

    AddOperator("DataBinning",0)
    DataBinningAtts = DataBinningAttributes()
    DataBinningAtts.numDimensions = DataBinningAtts.Two  # One, Two, Three
    DataBinningAtts.dim1BinBasedOn = DataBinningAtts.X  # X, Y, Z, Variable
    DataBinningAtts.dim1Var = "default"
    DataBinningAtts.dim1SpecifyRange = 0
    DataBinningAtts.dim1MinRange = 0
    DataBinningAtts.dim1MaxRange = 1
    DataBinningAtts.dim1NumBins = nbx
    DataBinningAtts.dim2BinBasedOn = DataBinningAtts.Z  # X, Y, Z, Variable
    DataBinningAtts.dim2Var = "default"
    DataBinningAtts.dim2SpecifyRange = 0
    DataBinningAtts.dim2MinRange = 0
    DataBinningAtts.dim2MaxRange = 1
    DataBinningAtts.dim2NumBins = nbz
    DataBinningAtts.dim3BinBasedOn = DataBinningAtts.Variable  # X, Y, Z, Variable
    DataBinningAtts.dim3Var = "default"
    DataBinningAtts.dim3SpecifyRange = 0
    DataBinningAtts.dim3MinRange = 0
    DataBinningAtts.dim3MaxRange = 1
    DataBinningAtts.dim3NumBins = 50
    DataBinningAtts.outOfBoundsBehavior = DataBinningAtts.Clamp  # Clamp, Discard
    DataBinningAtts.reductionOperator = DataBinningAtts.Sum  # Average, Minimum, Maximum, StandardDeviation, Variance, Sum, Count, RMS, PDF
    DataBinningAtts.varForReduction = "Electron-0/weights"
    DataBinningAtts.emptyVal = 0
    #DataBinningAtts.outputType = DataBinningAtts.OutputOnInputMesh  # OutputOnBins, OutputOnInputMesh
    DataBinningAtts.outputType = DataBinningAtts.OutputOnBins  # OutputOnBins, OutputOnInputMesh
    DataBinningAtts.removeEmptyValFromCurve = 1
    SetOperatorOptions(DataBinningAtts, 0)

    #AddOperator("Project", 0)
    #ProjectAtts = ProjectAttributes()
    #ProjectAtts.projectionType = ProjectAtts.XZCartesian  # ZYCartesian, XZCartesian, XYCartesian, XRCylindrical, YRCylindrical, ZRCylindrical
    #ProjectAtts.vectorTransformMethod = ProjectAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
    #SetOperatorOptions(ProjectAtts, 0)

    AddOperator("Elevate", 0)
    ElevateAtts = ElevateAttributes()
    ElevateAtts.useXYLimits = 0
    ElevateAtts.limitsMode = ElevateAtts.OriginalData  # OriginalData, CurrentPlot
    ElevateAtts.scaling = ElevateAtts.Linear  # Linear, Log, Skew
    ElevateAtts.skewFactor = 1
    ElevateAtts.minFlag = 0
    ElevateAtts.min = 0
    ElevateAtts.maxFlag = 0
    ElevateAtts.max = 1
    ElevateAtts.zeroFlag = 1
    ElevateAtts.variable = "default"
    SetOperatorOptions(ElevateAtts, 0)

    # fix pcolor variable
    ChangeActivePlotsVar("operators/DataBinning/2D/Electron-0")

    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
    PseudocolorAtts.skewFactor = 1
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = 0
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = 6e11
    #PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
    PseudocolorAtts.centering = PseudocolorAtts.Nodal  # Natural, Nodal, Zonal
    PseudocolorAtts.colorTableName = "Spectral"
    PseudocolorAtts.invertColorTable = 1
    PseudocolorAtts.opacityType = PseudocolorAtts.Constant  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
    PseudocolorAtts.opacityVariable = ""
    #PseudocolorAtts.opacity = 0.505882
    PseudocolorAtts.opacity = 1.0
    PseudocolorAtts.opacityVarMin = 0
    PseudocolorAtts.opacityVarMax = 1
    PseudocolorAtts.opacityVarMinFlag = 0
    PseudocolorAtts.opacityVarMaxFlag = 0
    PseudocolorAtts.pointSize = 0.05
    PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    PseudocolorAtts.pointSizeVarEnabled = 0
    PseudocolorAtts.pointSizeVar = "default"
    PseudocolorAtts.pointSizePixels = 2
    PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
    PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    PseudocolorAtts.lineWidth = 0
    PseudocolorAtts.tubeDisplayDensity = 10
    PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
    PseudocolorAtts.tubeRadiusAbsolute = 0.125
    PseudocolorAtts.tubeRadiusBBox = 0.005
    PseudocolorAtts.varyTubeRadius = 0
    PseudocolorAtts.varyTubeRadiusVariable = ""
    PseudocolorAtts.varyTubeRadiusFactor = 10
    PseudocolorAtts.endPointType = PseudocolorAtts.None  # None, Tails, Heads, Both
    PseudocolorAtts.endPointStyle = PseudocolorAtts.Spheres  # Spheres, Cones
    PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
    PseudocolorAtts.endPointRadiusAbsolute = 1
    PseudocolorAtts.endPointRadiusBBox = 0.005
    PseudocolorAtts.endPointRatio = 2
    PseudocolorAtts.renderSurfaces = 1
    PseudocolorAtts.renderWireframe = 0
    PseudocolorAtts.renderPoints = 0
    PseudocolorAtts.smoothingLevel = 0
    PseudocolorAtts.legendFlag = 1
    #PseudocolorAtts.lightingFlag = 1
    PseudocolorAtts.lightingFlag = 0
    SetPlotOptions(PseudocolorAtts)

    setAnnotations(winId, plotId, xAxisName='X', yAxisName='Z', zAxisName='', showDB=1)

    DrawPlots()
    setView(winId, plotId)

    return plotId

#----------------------------------------------------------------------------
def plotPColorElec1XZUz(winId):

    SetActiveWindow(winId)

    AddPlot("Pseudocolor", "Electron-1/uz", 1, 0)
    plotId = GetNumPlots() - 1
    SetActivePlots(plotId)

    AddOperator("Project", 0)
    ProjectAtts = ProjectAttributes()
    ProjectAtts.projectionType = ProjectAtts.XZCartesian  # ZYCartesian, XZCartesian, XYCartesian, XRCylindrical, YRCylindrical, ZRCylindrical
    ProjectAtts.vectorTransformMethod = ProjectAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
    SetOperatorOptions(ProjectAtts, 0)

    DefineScalarExpression("Electron-1-uz", "<Electron-1/uz>*1e-14")

    AddOperator("Elevate", 0)
    ElevateAtts = ElevateAttributes()
    ElevateAtts.useXYLimits = 0
    ElevateAtts.limitsMode = ElevateAtts.OriginalData  # OriginalData, CurrentPlot
    ElevateAtts.scaling = ElevateAtts.Linear  # Linear, Log, Skew
    ElevateAtts.skewFactor = 1
    ElevateAtts.minFlag = 0
    ElevateAtts.min = 0
    ElevateAtts.maxFlag = 0
    ElevateAtts.max = 7e+08
    ElevateAtts.zeroFlag = 0
    ElevateAtts.variable = "Electron-1-uz"
    SetOperatorOptions(ElevateAtts, 0)

    #AddOperator("Transform", 0)
    #TransformAtts = TransformAttributes()
    #TransformAtts.doRotate = 0
    #TransformAtts.rotateOrigin = (0, 0, 0)
    #TransformAtts.rotateAxis = (0, 0, 1)
    #TransformAtts.rotateAmount = 0
    #TransformAtts.rotateType = TransformAtts.Deg  # Deg, Rad
    #TransformAtts.doScale = 1
    #TransformAtts.scaleOrigin = (0, 0, 0)
    #TransformAtts.scaleX = 1
    #TransformAtts.scaleY = 1
    #TransformAtts.scaleZ = 10
    #TransformAtts.doTranslate = 0
    #TransformAtts.translateX = 0
    #TransformAtts.translateY = 0
    #TransformAtts.translateZ = 0
    #TransformAtts.transformType = TransformAtts.Similarity  # Similarity, Coordinate, Linear
    #TransformAtts.inputCoordSys = TransformAtts.Cartesian  # Cartesian, Cylindrical, Spherical
    #TransformAtts.outputCoordSys = TransformAtts.Spherical  # Cartesian, Cylindrical, Spherical
    #TransformAtts.m00 = 1
    #TransformAtts.m01 = 0
    #TransformAtts.m02 = 0
    #TransformAtts.m03 = 0
    #TransformAtts.m10 = 0
    #TransformAtts.m11 = 1
    #TransformAtts.m12 = 0
    #TransformAtts.m13 = 0
    #TransformAtts.m20 = 0
    #TransformAtts.m21 = 0
    #TransformAtts.m22 = 1
    #TransformAtts.m23 = 0
    #TransformAtts.m30 = 0
    #TransformAtts.m31 = 0
    #TransformAtts.m32 = 0
    #TransformAtts.m33 = 1
    #TransformAtts.invertLinearTransform = 0
    #TransformAtts.vectorTransformMethod = TransformAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
    #TransformAtts.transformVectors = 1
    #SetOperatorOptions(TransformAtts, 0)

    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
    PseudocolorAtts.skewFactor = 1
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
    PseudocolorAtts.minFlag = 0
    PseudocolorAtts.min = 0
    PseudocolorAtts.maxFlag = 0
    PseudocolorAtts.max = 1
    PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
    PseudocolorAtts.colorTableName = "hot_desaturated"
    PseudocolorAtts.invertColorTable = 0
    PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
    PseudocolorAtts.opacityVariable = ""
    PseudocolorAtts.opacity = 1
    PseudocolorAtts.opacityVarMin = 0
    PseudocolorAtts.opacityVarMax = 1
    PseudocolorAtts.opacityVarMinFlag = 0
    PseudocolorAtts.opacityVarMaxFlag = 0
    PseudocolorAtts.pointSize = 0.05
    PseudocolorAtts.pointType = PseudocolorAtts.Sphere  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    PseudocolorAtts.pointSizeVarEnabled = 0
    PseudocolorAtts.pointSizeVar = "default"
    PseudocolorAtts.pointSizePixels = 5
    PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
    PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    PseudocolorAtts.lineWidth = 0
    PseudocolorAtts.tubeDisplayDensity = 10
    PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
    PseudocolorAtts.tubeRadiusAbsolute = 0.125
    PseudocolorAtts.tubeRadiusBBox = 0.005
    PseudocolorAtts.varyTubeRadius = 0
    PseudocolorAtts.varyTubeRadiusVariable = ""
    PseudocolorAtts.varyTubeRadiusFactor = 10
    PseudocolorAtts.endPointType = PseudocolorAtts.None  # None, Tails, Heads, Both
    PseudocolorAtts.endPointStyle = PseudocolorAtts.Spheres  # Spheres, Cones
    PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
    PseudocolorAtts.endPointRadiusAbsolute = 1
    PseudocolorAtts.endPointRadiusBBox = 0.005
    PseudocolorAtts.endPointRatio = 2
    PseudocolorAtts.renderSurfaces = 1
    PseudocolorAtts.renderWireframe = 0
    PseudocolorAtts.renderPoints = 0
    PseudocolorAtts.smoothingLevel = 0
    PseudocolorAtts.legendFlag = 1
    PseudocolorAtts.lightingFlag = 1
    SetPlotOptions(PseudocolorAtts)

    setAnnotations(winId, plotId, xAxisName='X', yAxisName='Z', zAxisName='Uz', showDB=1)

    DrawPlots()
    setView(winId, plotId)

    return plotId

#----------------------------------------------------------------------------
def clearAnnotations(winId, plotId):
    SetActiveWindow(winId)
    SetActivePlots(plotId)
    AnnotationAtts = AnnotationAttributes()
    AnnotationAtts.axes2D.visible = 0
    AnnotationAtts.axes2D.tickLocation = AnnotationAtts.axes2D.Outside  # Inside, Outside, Both
    AnnotationAtts.axes2D.tickAxes = AnnotationAtts.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
    AnnotationAtts.axes2D.xAxis.title.font.font = AnnotationAtts.axes2D.xAxis.title.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes2D.xAxis.label.font.font = AnnotationAtts.axes2D.xAxis.label.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes2D.yAxis.title.font.font = AnnotationAtts.axes2D.yAxis.title.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes2D.yAxis.label.font.font = AnnotationAtts.axes2D.yAxis.label.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes3D.visible = 0
    AnnotationAtts.axes3D.tickLocation = AnnotationAtts.axes3D.Inside  # Inside, Outside, Both
    AnnotationAtts.axes3D.axesType = AnnotationAtts.axes3D.ClosestTriad  # ClosestTriad, FurthestTriad, OutsideEdges, StaticTriad, StaticEdges
    AnnotationAtts.axes3D.triadFlag = 0
    AnnotationAtts.axes3D.bboxFlag = 0
    AnnotationAtts.axes3D.xAxis.title.font.font = AnnotationAtts.axes3D.xAxis.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.xAxis.label.font.font = AnnotationAtts.axes3D.xAxis.label.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.yAxis.title.font.font = AnnotationAtts.axes3D.yAxis.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.yAxis.label.font.font = AnnotationAtts.axes3D.yAxis.label.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.zAxis.title.font.font = AnnotationAtts.axes3D.zAxis.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.zAxis.label.font.font = AnnotationAtts.axes3D.zAxis.label.font.Arial  # Arial, Courier, Times
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.userInfoFont.font = AnnotationAtts.userInfoFont.Arial  # Arial, Courier, Times
    AnnotationAtts.databaseInfoFlag = 0
    AnnotationAtts.databaseInfoFont.font = AnnotationAtts.databaseInfoFont.Arial  # Arial, Courier, Times
    AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.File  # File, Directory, Full, Smart, SmartDirectory
    AnnotationAtts.legendInfoFlag = 0
    AnnotationAtts.gradientBackgroundStyle = AnnotationAtts.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
    AnnotationAtts.backgroundMode = AnnotationAtts.Solid  # Solid, Gradient, Image, ImageSphere
    AnnotationAtts.axesArray.axes.title.font.font = AnnotationAtts.axesArray.axes.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axesArray.axes.label.font.font = AnnotationAtts.axesArray.axes.label.font.Arial  # Arial, Courier, Times
    SetAnnotationAttributes(AnnotationAtts)
    return

#----------------------------------------------------------------------------
def setAnnotations(winId, plotId, xAxisName='X', yAxisName='Z', zAxisName='Uz', showDB=0):
    SetActiveWindow(winId)
    SetActivePlots(plotId)
    AnnotationAtts = AnnotationAttributes()
    AnnotationAtts.axes2D.tickLocation = AnnotationAtts.axes2D.Outside  # Inside, Outside, Both
    AnnotationAtts.axes2D.tickAxes = AnnotationAtts.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
    AnnotationAtts.axes2D.xAxis.title.font.font = AnnotationAtts.axes2D.xAxis.title.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes2D.xAxis.title.title = xAxisName
    AnnotationAtts.axes2D.xAxis.label.font.font = AnnotationAtts.axes2D.xAxis.label.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes2D.yAxis.title.font.font = AnnotationAtts.axes2D.yAxis.title.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes2D.yAxis.title.title = yAxisName
    AnnotationAtts.axes2D.yAxis.label.font.font = AnnotationAtts.axes2D.yAxis.label.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes3D.tickLocation = AnnotationAtts.axes3D.Inside  # Inside, Outside, Both
    AnnotationAtts.axes3D.axesType = AnnotationAtts.axes3D.OutsideEdges  # ClosestTriad, FurthestTriad, OutsideEdges, StaticTriad, StaticEdges
    AnnotationAtts.axes3D.bboxFlag = 0
    AnnotationAtts.axes3D.xAxis.title.font.font = AnnotationAtts.axes3D.xAxis.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.xAxis.title.font.scale = 1.2
    AnnotationAtts.axes3D.xAxis.title.font.bold = 1
    AnnotationAtts.axes3D.xAxis.title.userTitle = 1
    AnnotationAtts.axes3D.xAxis.title.title = xAxisName
    AnnotationAtts.axes3D.xAxis.label.visible = 0
    AnnotationAtts.axes3D.xAxis.label.font.font = AnnotationAtts.axes3D.xAxis.label.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.yAxis.title.font.font = AnnotationAtts.axes3D.yAxis.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.yAxis.title.font.scale = 1.2
    AnnotationAtts.axes3D.yAxis.title.font.bold = 1
    AnnotationAtts.axes3D.yAxis.title.userTitle = 1
    AnnotationAtts.axes3D.yAxis.title.title = yAxisName
    AnnotationAtts.axes3D.yAxis.label.visible = 0
    AnnotationAtts.axes3D.yAxis.label.font.font = AnnotationAtts.axes3D.yAxis.label.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.zAxis.title.font.font = AnnotationAtts.axes3D.zAxis.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.zAxis.title.font.scale = 1.2
    AnnotationAtts.axes3D.zAxis.title.font.bold = 1
    AnnotationAtts.axes3D.zAxis.title.userTitle = 1
    AnnotationAtts.axes3D.zAxis.title.title = zAxisName
    AnnotationAtts.axes3D.zAxis.label.visible = 0
    AnnotationAtts.axes3D.zAxis.label.font.font = AnnotationAtts.axes3D.zAxis.label.font.Arial  # Arial, Courier, Times
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.userInfoFont.font = AnnotationAtts.userInfoFont.Arial  # Arial, Courier, Times
    AnnotationAtts.databaseInfoFlag = showDB
    AnnotationAtts.databaseInfoFont.font = AnnotationAtts.databaseInfoFont.Arial  # Arial, Courier, Times
    AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.File  # File, Directory, Full, Smart, SmartDirectory
    AnnotationAtts.gradientBackgroundStyle = AnnotationAtts.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
    AnnotationAtts.backgroundMode = AnnotationAtts.Solid  # Solid, Gradient, Image, ImageSphere
    AnnotationAtts.axesArray.axes.title.font.font = AnnotationAtts.axesArray.axes.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axesArray.axes.label.font.font = AnnotationAtts.axesArray.axes.label.font.Arial  # Arial, Courier, Times
    SetAnnotationAttributes(AnnotationAtts)
    return

#----------------------------------------------------------------------------
def saveWindows(imw,imh):
    nrow=2
    ncol=2

    imw/=ncol
    imh/=nrow

    row = lambda i: (int(i)/ncol)
    col = lambda i: (int(i)%ncol)

    SaveWindowAtts = SaveWindowAttributes()

    i=0
    SaveWindowAtts.subWindowAtts.win3.position = (col(i)*imw, row(i)*imh)
    i+=1
    SaveWindowAtts.subWindowAtts.win4.position = (col(i)*imw, row(i)*imh)
    i+=1
    SaveWindowAtts.subWindowAtts.win1.position = (col(i)*imw, row(i)*imh)
    i+=1
    SaveWindowAtts.subWindowAtts.win2.position = (col(i)*imw, row(i)*imh)
    i+=1

    SaveWindowAtts.subWindowAtts.win1.size = (imw,imh)
    SaveWindowAtts.subWindowAtts.win2.size = (imw,imh)
    SaveWindowAtts.subWindowAtts.win3.size = (imw,imh)
    SaveWindowAtts.subWindowAtts.win4.size = (imw,imh)

    SaveWindowAtts.outputToCurrentDirectory = 1
    SaveWindowAtts.outputDirectory = "./"
    SaveWindowAtts.fileName = "lpa-4-view"
    SaveWindowAtts.family = 1
    SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
    SaveWindowAtts.width = imw*nrow
    SaveWindowAtts.height = imh*ncol

    SaveWindowAtts.screenCapture = 0
    SaveWindowAtts.saveTiled = 0
    SaveWindowAtts.quality = 80
    SaveWindowAtts.progressive = 0
    SaveWindowAtts.binary = 0
    SaveWindowAtts.stereo = 0
    SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
    SaveWindowAtts.forceMerge = 0
    SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint  # NoConstraint, EqualWidthHeight, ScreenProportions
    SaveWindowAtts.advancedMultiWindowSave = 1

    SetSaveWindowAttributes(SaveWindowAtts)

    SaveWindow()
    return


# do the vis!
SetWindowLayout(1)
elec0 = meshValid(1,'Electron-0')
elec1 = meshValid(1,'Electron-1')

SetWindowLayout(4)
i=1
while i<=4:
    deletePlots(i)
    i+=1

if elec0:
    plotPColorElec0Uz(1)
    plotPColorElec0Uy(2)
    plotBinElec0Weights(4,200,400)
else:
    meshPlot(1,'Electron-0')
    meshPlot(2,'Electron-0')
    meshPlot(4,'Electron-0')

if elec1:
    plotPColorElec1Uz(3)
    plotPColorElec1XZUz(4)
else:
    meshPlot(3,'Electron-1')

saveWindows(1920,1080)

SetQueryOutputToValue()
memUse = Query("Memory Usage")

try:
    memUse = sorted(memUse)
    n = len(memUse)
    if n%2 == 0:
        medMem = (memUse[n/2] + memUse[n/2-1])/2.0
    else:
        medMem = memUse[n/2+1]
    sys.stderr.write('VisItMemUse = [%g %g %g]\n'%(min(memUse), medMem, max(memUse)))
except:
    sys.stderr.write('VisItMemUse = %g\n'%(memUse))

i=1
while i<=4:
    deletePlots(i)
    i+=1
