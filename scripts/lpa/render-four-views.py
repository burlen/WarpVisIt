from visit import *
SetWindowLayout(4)


def hide_if_empty():
    SetQueryFloatFormat("%g")
    SetQueryOutputToValue()
    numNodes = Query("NumNodes")
    if numNodes < 10:
        HideActivePlots()

def set_view():
    ResetView()
    View3DAtts = GetView3D()
    View3DAtts.viewNormal = (0.82, 0.53, 0.25)
    View3DAtts.viewUp = (-0.19, -0.16, 0.97)
    SetView3D(View3DAtts)


def setup_plot1():
    AddPlot("Pseudocolor", "Electron-0/uz", 1, 1)
    SetActivePlots(0)
    SetActivePlots(0)
    AddOperator("Project", 1)
    ProjectAtts = ProjectAttributes()
    ProjectAtts.projectionType = ProjectAtts.XZCartesian  # ZYCartesian, XZCartesian, XYCartesian, XRCylindrical, YRCylindrical, ZRCylindrical
    ProjectAtts.vectorTransformMethod = ProjectAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
    SetOperatorOptions(ProjectAtts, 1)
    AddOperator("Elevate", 1)
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
    PseudocolorAtts.min = -1e+09
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = 1e+09
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
    DrawPlots()


def setup_plot2():

    AddPlot("Pseudocolor", "Electron-0/uy", 1, 1)
    SetActivePlots(0)
    SetActivePlots(0)
    AddOperator("Project", 1)
    ProjectAtts = ProjectAttributes()
    ProjectAtts.projectionType = ProjectAtts.XZCartesian  # ZYCartesian, XZCartesian, XYCartesian, XRCylindrical, YRCylindrical, ZRCylindrical
    ProjectAtts.vectorTransformMethod = ProjectAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
    SetOperatorOptions(ProjectAtts, 1)
    AddOperator("Elevate", 1)
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
    PseudocolorAtts.min = -1e+09
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = 1e+09
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
    DrawPlots()


def setup_plot3():

    AddPlot("Pseudocolor", "Electron-1/uz", 1, 1)
    SetActivePlots(0)
    SetActivePlots(0)
    AddOperator("Project", 1)
    ProjectAtts = ProjectAttributes()
    ProjectAtts.projectionType = ProjectAtts.XZCartesian  # ZYCartesian, XZCartesian, XYCartesian, XRCylindrical, YRCylindrical, ZRCylindrical
    ProjectAtts.vectorTransformMethod = ProjectAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
    SetOperatorOptions(ProjectAtts, 1)
    AddOperator("Elevate", 1)
    ElevateAtts = ElevateAttributes()
    ElevateAtts.useXYLimits = 1
    ElevateAtts.limitsMode = ElevateAtts.OriginalData  # OriginalData, CurrentPlot
    ElevateAtts.scaling = ElevateAtts.Linear  # Linear, Log, Skew
    ElevateAtts.skewFactor = 1
    ElevateAtts.minFlag = 1
    ElevateAtts.min = -5e+08
    ElevateAtts.maxFlag = 1
    ElevateAtts.max = 5e+08
    ElevateAtts.zeroFlag = 0
    ElevateAtts.variable = "Electron-0/ux"
    SetOperatorOptions(ElevateAtts, 1)
    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
    PseudocolorAtts.skewFactor = 1
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
    PseudocolorAtts.minFlag = 0
    PseudocolorAtts.min = -5e+08
    PseudocolorAtts.maxFlag = 0
    PseudocolorAtts.max = 5e+08
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
    DrawPlots()


def setup_plot4():
    #########################################
    # elec density                          #
    #########################################
    AddPlot("Pseudocolor", "operators/DataBinning/2D/Electron-0", 1, 1)
    DataBinningAtts = DataBinningAttributes()
    DataBinningAtts.numDimensions = DataBinningAtts.Two  # One, Two, Three
    DataBinningAtts.dim1BinBasedOn = DataBinningAtts.X  # X, Y, Z, Variable
    DataBinningAtts.dim1Var = "default"
    DataBinningAtts.dim1SpecifyRange = 0
    DataBinningAtts.dim1MinRange = 0
    DataBinningAtts.dim1MaxRange = 1
    DataBinningAtts.dim1NumBins = 200
    DataBinningAtts.dim2BinBasedOn = DataBinningAtts.Z  # X, Y, Z, Variable
    DataBinningAtts.dim2Var = "default"
    DataBinningAtts.dim2SpecifyRange = 0
    DataBinningAtts.dim2MinRange = 0
    DataBinningAtts.dim2MaxRange = 1
    DataBinningAtts.dim2NumBins = 400
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
    DataBinningAtts.outputType = DataBinningAtts.OutputOnInputMesh  # OutputOnBins, OutputOnInputMesh
    DataBinningAtts.removeEmptyValFromCurve = 1
    SetOperatorOptions(DataBinningAtts, 1)
    AddOperator("Project", 1)
    ProjectAtts = ProjectAttributes()
    ProjectAtts.projectionType = ProjectAtts.XZCartesian  # ZYCartesian, XZCartesian, XYCartesian, XRCylindrical, YRCylindrical, ZRCylindrical
    ProjectAtts.vectorTransformMethod = ProjectAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
    SetOperatorOptions(ProjectAtts, 1)
    AddOperator("Elevate", 1)
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
    SetOperatorOptions(ElevateAtts, 1)
    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
    PseudocolorAtts.skewFactor = 1
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
    PseudocolorAtts.minFlag = 0
    PseudocolorAtts.min = 0
    PseudocolorAtts.maxFlag = 0
    PseudocolorAtts.max = 1
    PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
    PseudocolorAtts.colorTableName = "Spectral"
    PseudocolorAtts.invertColorTable = 1
    PseudocolorAtts.opacityType = PseudocolorAtts.Constant  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
    PseudocolorAtts.opacityVariable = ""
    PseudocolorAtts.opacity = 0.505882
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
    DrawPlots()
    ##############################################
    #   ielec uz particles                       #
    ##############################################
    AddPlot("Pseudocolor", "Electron-1/uz", 1, 0)
    AddOperator("Project", 0)
    SetActivePlots(1)
    SetActivePlots(1)
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
    ElevateAtts.min = 0
    ElevateAtts.maxFlag = 1
    ElevateAtts.max = 7e+08
    ElevateAtts.zeroFlag = 0
    ElevateAtts.variable = "default"
    SetOperatorOptions(ElevateAtts, 0)
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
    DrawPlots()


def set_annotations(xAxisName='X', yAxisName='Z', zAxisName='Uz', showDB=0):

    # Logging for SetAnnotationObjectOptions is not implemented yet.
    AnnotationAtts = AnnotationAttributes()
    AnnotationAtts.axes2D.visible = 1
    AnnotationAtts.axes2D.autoSetTicks = 1
    AnnotationAtts.axes2D.autoSetScaling = 1
    AnnotationAtts.axes2D.lineWidth = 0
    AnnotationAtts.axes2D.tickLocation = AnnotationAtts.axes2D.Outside  # Inside, Outside, Both
    AnnotationAtts.axes2D.tickAxes = AnnotationAtts.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
    AnnotationAtts.axes2D.xAxis.title.visible = 1
    AnnotationAtts.axes2D.xAxis.title.font.font = AnnotationAtts.axes2D.xAxis.title.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes2D.xAxis.title.font.scale = 1
    AnnotationAtts.axes2D.xAxis.title.font.useForegroundColor = 1
    AnnotationAtts.axes2D.xAxis.title.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes2D.xAxis.title.font.bold = 1
    AnnotationAtts.axes2D.xAxis.title.font.italic = 1
    AnnotationAtts.axes2D.xAxis.title.userTitle = 0
    AnnotationAtts.axes2D.xAxis.title.userUnits = 0
    AnnotationAtts.axes2D.xAxis.title.title = xAxisName
    AnnotationAtts.axes2D.xAxis.title.units = ""
    AnnotationAtts.axes2D.xAxis.label.visible = 1
    AnnotationAtts.axes2D.xAxis.label.font.font = AnnotationAtts.axes2D.xAxis.label.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes2D.xAxis.label.font.scale = 1
    AnnotationAtts.axes2D.xAxis.label.font.useForegroundColor = 1
    AnnotationAtts.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes2D.xAxis.label.font.bold = 1
    AnnotationAtts.axes2D.xAxis.label.font.italic = 1
    AnnotationAtts.axes2D.xAxis.label.scaling = 0
    AnnotationAtts.axes2D.xAxis.tickMarks.visible = 1
    AnnotationAtts.axes2D.xAxis.tickMarks.majorMinimum = 0
    AnnotationAtts.axes2D.xAxis.tickMarks.majorMaximum = 1
    AnnotationAtts.axes2D.xAxis.tickMarks.minorSpacing = 0.02
    AnnotationAtts.axes2D.xAxis.tickMarks.majorSpacing = 0.2
    AnnotationAtts.axes2D.xAxis.grid = 0
    AnnotationAtts.axes2D.yAxis.title.visible = 1
    AnnotationAtts.axes2D.yAxis.title.font.font = AnnotationAtts.axes2D.yAxis.title.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes2D.yAxis.title.font.scale = 1
    AnnotationAtts.axes2D.yAxis.title.font.useForegroundColor = 1
    AnnotationAtts.axes2D.yAxis.title.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes2D.yAxis.title.font.bold = 1
    AnnotationAtts.axes2D.yAxis.title.font.italic = 1
    AnnotationAtts.axes2D.yAxis.title.userTitle = 0
    AnnotationAtts.axes2D.yAxis.title.userUnits = 0
    AnnotationAtts.axes2D.yAxis.title.title = yAxisName
    AnnotationAtts.axes2D.yAxis.title.units = ""
    AnnotationAtts.axes2D.yAxis.label.visible = 1
    AnnotationAtts.axes2D.yAxis.label.font.font = AnnotationAtts.axes2D.yAxis.label.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes2D.yAxis.label.font.scale = 1
    AnnotationAtts.axes2D.yAxis.label.font.useForegroundColor = 1
    AnnotationAtts.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes2D.yAxis.label.font.bold = 1
    AnnotationAtts.axes2D.yAxis.label.font.italic = 1
    AnnotationAtts.axes2D.yAxis.label.scaling = 0
    AnnotationAtts.axes2D.yAxis.tickMarks.visible = 1
    AnnotationAtts.axes2D.yAxis.tickMarks.majorMinimum = 0
    AnnotationAtts.axes2D.yAxis.tickMarks.majorMaximum = 1
    AnnotationAtts.axes2D.yAxis.tickMarks.minorSpacing = 0.02
    AnnotationAtts.axes2D.yAxis.tickMarks.majorSpacing = 0.2
    AnnotationAtts.axes2D.yAxis.grid = 0
    AnnotationAtts.axes3D.visible = 1
    AnnotationAtts.axes3D.autoSetTicks = 1
    AnnotationAtts.axes3D.autoSetScaling = 1
    AnnotationAtts.axes3D.lineWidth = 0
    AnnotationAtts.axes3D.tickLocation = AnnotationAtts.axes3D.Inside  # Inside, Outside, Both
    AnnotationAtts.axes3D.axesType = AnnotationAtts.axes3D.OutsideEdges  # ClosestTriad, FurthestTriad, OutsideEdges, StaticTriad, StaticEdges
    AnnotationAtts.axes3D.triadFlag = 1
    AnnotationAtts.axes3D.bboxFlag = 0
    AnnotationAtts.axes3D.xAxis.title.visible = 1
    AnnotationAtts.axes3D.xAxis.title.font.font = AnnotationAtts.axes3D.xAxis.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.xAxis.title.font.scale = 1.2
    AnnotationAtts.axes3D.xAxis.title.font.useForegroundColor = 1
    AnnotationAtts.axes3D.xAxis.title.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes3D.xAxis.title.font.bold = 1
    AnnotationAtts.axes3D.xAxis.title.font.italic = 0
    AnnotationAtts.axes3D.xAxis.title.userTitle = 1
    AnnotationAtts.axes3D.xAxis.title.userUnits = 0
    AnnotationAtts.axes3D.xAxis.title.title = xAxisName
    AnnotationAtts.axes3D.xAxis.title.units = ""
    AnnotationAtts.axes3D.xAxis.label.visible = 0
    AnnotationAtts.axes3D.xAxis.label.font.font = AnnotationAtts.axes3D.xAxis.label.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.xAxis.label.font.scale = 1
    AnnotationAtts.axes3D.xAxis.label.font.useForegroundColor = 1
    AnnotationAtts.axes3D.xAxis.label.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes3D.xAxis.label.font.bold = 0
    AnnotationAtts.axes3D.xAxis.label.font.italic = 0
    AnnotationAtts.axes3D.xAxis.label.scaling = 0
    AnnotationAtts.axes3D.xAxis.tickMarks.visible = 1
    AnnotationAtts.axes3D.xAxis.tickMarks.majorMinimum = 0
    AnnotationAtts.axes3D.xAxis.tickMarks.majorMaximum = 1
    AnnotationAtts.axes3D.xAxis.tickMarks.minorSpacing = 0.02
    AnnotationAtts.axes3D.xAxis.tickMarks.majorSpacing = 0.2
    AnnotationAtts.axes3D.xAxis.grid = 0
    AnnotationAtts.axes3D.yAxis.title.visible = 1
    AnnotationAtts.axes3D.yAxis.title.font.font = AnnotationAtts.axes3D.yAxis.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.yAxis.title.font.scale = 1.2
    AnnotationAtts.axes3D.yAxis.title.font.useForegroundColor = 1
    AnnotationAtts.axes3D.yAxis.title.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes3D.yAxis.title.font.bold = 1
    AnnotationAtts.axes3D.yAxis.title.font.italic = 0
    AnnotationAtts.axes3D.yAxis.title.userTitle = 1
    AnnotationAtts.axes3D.yAxis.title.userUnits = 0
    AnnotationAtts.axes3D.yAxis.title.title = yAxisName
    AnnotationAtts.axes3D.yAxis.title.units = ""
    AnnotationAtts.axes3D.yAxis.label.visible = 0
    AnnotationAtts.axes3D.yAxis.label.font.font = AnnotationAtts.axes3D.yAxis.label.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.yAxis.label.font.scale = 1
    AnnotationAtts.axes3D.yAxis.label.font.useForegroundColor = 1
    AnnotationAtts.axes3D.yAxis.label.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes3D.yAxis.label.font.bold = 0
    AnnotationAtts.axes3D.yAxis.label.font.italic = 0
    AnnotationAtts.axes3D.yAxis.label.scaling = 0
    AnnotationAtts.axes3D.yAxis.tickMarks.visible = 1
    AnnotationAtts.axes3D.yAxis.tickMarks.majorMinimum = 0
    AnnotationAtts.axes3D.yAxis.tickMarks.majorMaximum = 1
    AnnotationAtts.axes3D.yAxis.tickMarks.minorSpacing = 0.02
    AnnotationAtts.axes3D.yAxis.tickMarks.majorSpacing = 0.2
    AnnotationAtts.axes3D.yAxis.grid = 0
    AnnotationAtts.axes3D.zAxis.title.visible = 1
    AnnotationAtts.axes3D.zAxis.title.font.font = AnnotationAtts.axes3D.zAxis.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.zAxis.title.font.scale = 1.2
    AnnotationAtts.axes3D.zAxis.title.font.useForegroundColor = 1
    AnnotationAtts.axes3D.zAxis.title.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes3D.zAxis.title.font.bold = 1
    AnnotationAtts.axes3D.zAxis.title.font.italic = 0
    AnnotationAtts.axes3D.zAxis.title.userTitle = 1
    AnnotationAtts.axes3D.zAxis.title.userUnits = 0
    AnnotationAtts.axes3D.zAxis.title.title = zAxisName
    AnnotationAtts.axes3D.zAxis.title.units = ""
    AnnotationAtts.axes3D.zAxis.label.visible = 0
    AnnotationAtts.axes3D.zAxis.label.font.font = AnnotationAtts.axes3D.zAxis.label.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.zAxis.label.font.scale = 1
    AnnotationAtts.axes3D.zAxis.label.font.useForegroundColor = 1
    AnnotationAtts.axes3D.zAxis.label.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes3D.zAxis.label.font.bold = 0
    AnnotationAtts.axes3D.zAxis.label.font.italic = 0
    AnnotationAtts.axes3D.zAxis.label.scaling = 0
    AnnotationAtts.axes3D.zAxis.tickMarks.visible = 1
    AnnotationAtts.axes3D.zAxis.tickMarks.majorMinimum = 0
    AnnotationAtts.axes3D.zAxis.tickMarks.majorMaximum = 1
    AnnotationAtts.axes3D.zAxis.tickMarks.minorSpacing = 0.02
    AnnotationAtts.axes3D.zAxis.tickMarks.majorSpacing = 0.2
    AnnotationAtts.axes3D.zAxis.grid = 0
    AnnotationAtts.axes3D.setBBoxLocation = 0
    AnnotationAtts.axes3D.bboxLocation = (0, 1, 0, 1, 0, 1)
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.userInfoFont.font = AnnotationAtts.userInfoFont.Arial  # Arial, Courier, Times
    AnnotationAtts.userInfoFont.scale = 1
    AnnotationAtts.userInfoFont.useForegroundColor = 1
    AnnotationAtts.userInfoFont.color = (0, 0, 0, 255)
    AnnotationAtts.userInfoFont.bold = 0
    AnnotationAtts.userInfoFont.italic = 0
    AnnotationAtts.databaseInfoFlag = showDB
    AnnotationAtts.timeInfoFlag = 1
    AnnotationAtts.databaseInfoFont.font = AnnotationAtts.databaseInfoFont.Arial  # Arial, Courier, Times
    AnnotationAtts.databaseInfoFont.scale = 1
    AnnotationAtts.databaseInfoFont.useForegroundColor = 1
    AnnotationAtts.databaseInfoFont.color = (0, 0, 0, 255)
    AnnotationAtts.databaseInfoFont.bold = 0
    AnnotationAtts.databaseInfoFont.italic = 0
    AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.File  # File, Directory, Full, Smart, SmartDirectory
    AnnotationAtts.databaseInfoTimeScale = 1
    AnnotationAtts.databaseInfoTimeOffset = 0
    AnnotationAtts.legendInfoFlag = 1
    AnnotationAtts.backgroundColor = (255, 255, 255, 255)
    AnnotationAtts.foregroundColor = (0, 0, 0, 255)
    AnnotationAtts.gradientBackgroundStyle = AnnotationAtts.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
    AnnotationAtts.gradientColor1 = (0, 0, 255, 255)
    AnnotationAtts.gradientColor2 = (0, 0, 0, 255)
    AnnotationAtts.backgroundMode = AnnotationAtts.Solid  # Solid, Gradient, Image, ImageSphere
    AnnotationAtts.backgroundImage = ""
    AnnotationAtts.imageRepeatX = 1
    AnnotationAtts.imageRepeatY = 1
    AnnotationAtts.axesArray.visible = 1
    AnnotationAtts.axesArray.ticksVisible = 1
    AnnotationAtts.axesArray.autoSetTicks = 1
    AnnotationAtts.axesArray.autoSetScaling = 1
    AnnotationAtts.axesArray.lineWidth = 0
    AnnotationAtts.axesArray.axes.title.visible = 1
    AnnotationAtts.axesArray.axes.title.font.font = AnnotationAtts.axesArray.axes.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axesArray.axes.title.font.scale = 1
    AnnotationAtts.axesArray.axes.title.font.useForegroundColor = 1
    AnnotationAtts.axesArray.axes.title.font.color = (0, 0, 0, 255)
    AnnotationAtts.axesArray.axes.title.font.bold = 0
    AnnotationAtts.axesArray.axes.title.font.italic = 0
    AnnotationAtts.axesArray.axes.title.userTitle = 0
    AnnotationAtts.axesArray.axes.title.userUnits = 0
    AnnotationAtts.axesArray.axes.title.title = ""
    AnnotationAtts.axesArray.axes.title.units = ""
    AnnotationAtts.axesArray.axes.label.visible = 1
    AnnotationAtts.axesArray.axes.label.font.font = AnnotationAtts.axesArray.axes.label.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axesArray.axes.label.font.scale = 1
    AnnotationAtts.axesArray.axes.label.font.useForegroundColor = 1
    AnnotationAtts.axesArray.axes.label.font.color = (0, 0, 0, 255)
    AnnotationAtts.axesArray.axes.label.font.bold = 0
    AnnotationAtts.axesArray.axes.label.font.italic = 0
    AnnotationAtts.axesArray.axes.label.scaling = 0
    AnnotationAtts.axesArray.axes.tickMarks.visible = 1
    AnnotationAtts.axesArray.axes.tickMarks.majorMinimum = 0
    AnnotationAtts.axesArray.axes.tickMarks.majorMaximum = 1
    AnnotationAtts.axesArray.axes.tickMarks.minorSpacing = 0.02
    AnnotationAtts.axesArray.axes.tickMarks.majorSpacing = 0.2
    AnnotationAtts.axesArray.axes.grid = 0
    SetAnnotationAttributes(AnnotationAtts)


def save_window_matrix():

    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.outputToCurrentDirectory = 0
    SaveWindowAtts.outputDirectory = "/Users/oruebel/Devel"
    SaveWindowAtts.fileName = "lpa_testimage_"
    SaveWindowAtts.family = 1
    SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
    SaveWindowAtts.width = 3200
    SaveWindowAtts.height = 2200
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
    SaveWindowAtts.subWindowAtts.win1.position = (0, 1100)
    SaveWindowAtts.subWindowAtts.win1.size = (1600, 1100)
    SaveWindowAtts.subWindowAtts.win1.layer = 0
    SaveWindowAtts.subWindowAtts.win1.transparency = 0
    SaveWindowAtts.subWindowAtts.win1.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win2.position = (1600, 1100)
    SaveWindowAtts.subWindowAtts.win2.size = (1600, 1100)
    SaveWindowAtts.subWindowAtts.win2.layer = 0
    SaveWindowAtts.subWindowAtts.win2.transparency = 0
    SaveWindowAtts.subWindowAtts.win2.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win3.position = (0, 0)
    SaveWindowAtts.subWindowAtts.win3.size = (1600, 1100)
    SaveWindowAtts.subWindowAtts.win3.layer = 0
    SaveWindowAtts.subWindowAtts.win3.transparency = 0
    SaveWindowAtts.subWindowAtts.win3.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win4.position = (1600, 0)
    SaveWindowAtts.subWindowAtts.win4.size = (1600, 1100)
    SaveWindowAtts.subWindowAtts.win4.layer = 0
    SaveWindowAtts.subWindowAtts.win4.transparency = 0
    SaveWindowAtts.subWindowAtts.win4.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win5.position = (0, 0)
    SaveWindowAtts.subWindowAtts.win5.size = (128, 128)
    SaveWindowAtts.subWindowAtts.win5.layer = 0
    SaveWindowAtts.subWindowAtts.win5.transparency = 0
    SaveWindowAtts.subWindowAtts.win5.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win6.position = (0, 0)
    SaveWindowAtts.subWindowAtts.win6.size = (128, 128)
    SaveWindowAtts.subWindowAtts.win6.layer = 0
    SaveWindowAtts.subWindowAtts.win6.transparency = 0
    SaveWindowAtts.subWindowAtts.win6.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win7.position = (0, 0)
    SaveWindowAtts.subWindowAtts.win7.size = (128, 128)
    SaveWindowAtts.subWindowAtts.win7.layer = 0
    SaveWindowAtts.subWindowAtts.win7.transparency = 0
    SaveWindowAtts.subWindowAtts.win7.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win8.position = (0, 0)
    SaveWindowAtts.subWindowAtts.win8.size = (128, 128)
    SaveWindowAtts.subWindowAtts.win8.layer = 0
    SaveWindowAtts.subWindowAtts.win8.transparency = 0
    SaveWindowAtts.subWindowAtts.win8.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win9.position = (0, 0)
    SaveWindowAtts.subWindowAtts.win9.size = (128, 128)
    SaveWindowAtts.subWindowAtts.win9.layer = 0
    SaveWindowAtts.subWindowAtts.win9.transparency = 0
    SaveWindowAtts.subWindowAtts.win9.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win10.position = (0, 0)
    SaveWindowAtts.subWindowAtts.win10.size = (128, 128)
    SaveWindowAtts.subWindowAtts.win10.layer = 0
    SaveWindowAtts.subWindowAtts.win10.transparency = 0
    SaveWindowAtts.subWindowAtts.win10.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win11.position = (0, 0)
    SaveWindowAtts.subWindowAtts.win11.size = (128, 128)
    SaveWindowAtts.subWindowAtts.win11.layer = 0
    SaveWindowAtts.subWindowAtts.win11.transparency = 0
    SaveWindowAtts.subWindowAtts.win11.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win12.position = (0, 0)
    SaveWindowAtts.subWindowAtts.win12.size = (128, 128)
    SaveWindowAtts.subWindowAtts.win12.layer = 0
    SaveWindowAtts.subWindowAtts.win12.transparency = 0
    SaveWindowAtts.subWindowAtts.win12.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win13.position = (0, 0)
    SaveWindowAtts.subWindowAtts.win13.size = (128, 128)
    SaveWindowAtts.subWindowAtts.win13.layer = 0
    SaveWindowAtts.subWindowAtts.win13.transparency = 0
    SaveWindowAtts.subWindowAtts.win13.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win14.position = (0, 0)
    SaveWindowAtts.subWindowAtts.win14.size = (128, 128)
    SaveWindowAtts.subWindowAtts.win14.layer = 0
    SaveWindowAtts.subWindowAtts.win14.transparency = 0
    SaveWindowAtts.subWindowAtts.win14.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win15.position = (0, 0)
    SaveWindowAtts.subWindowAtts.win15.size = (128, 128)
    SaveWindowAtts.subWindowAtts.win15.layer = 0
    SaveWindowAtts.subWindowAtts.win15.transparency = 0
    SaveWindowAtts.subWindowAtts.win15.omitWindow = 0
    SaveWindowAtts.subWindowAtts.win16.position = (0, 0)
    SaveWindowAtts.subWindowAtts.win16.size = (128, 128)
    SaveWindowAtts.subWindowAtts.win16.layer = 0
    SaveWindowAtts.subWindowAtts.win16.transparency = 0
    SaveWindowAtts.subWindowAtts.win16.omitWindow = 0
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()


SetActiveWindow(1)
SetActiveWindow(1)
DeleteAllPlots()
setup_plot1()
set_annotations(xAxisName='X', yAxisName='Z', zAxisName='Uy', showDB=1)
set_view()
hide_if_empty()

SetActiveWindow(2)
SetActiveWindow(2)
DeleteAllPlots()
setup_plot2()
set_annotations(xAxisName='X', yAxisName='Z', zAxisName='Uz', showDB=0)
set_view()
hide_if_empty()

SetActiveWindow(3)
SetActiveWindow(3)
DeleteAllPlots()
setup_plot3()
set_annotations(xAxisName='X', yAxisName='Z', zAxisName='Ux', showDB=0)
set_view()
hide_if_empty()

SetActiveWindow(4)
SetActiveWindow(4)
DeleteAllPlots()
setup_plot4()
set_annotations(xAxisName='X', yAxisName='Z', zAxisName='Uz', showDB=0)
set_view()
hide_if_empty()

save_window_matrix()
  

SetActiveWindow(1)
DeleteAllPlots()
SetActiveWindow(2)
DeleteAllPlots()
SetActiveWindow(3)
DeleteAllPlots()
SetActiveWindow(4)
DeleteAllPlots()





