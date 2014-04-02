from visit import *

def hide_if_empty():
    sys.stderr.write('=============hide_if_empty\n')
    np = GetNumPlots()
    i = 0
    while(i<np):
        SetActivePlots(i)
        SetQueryFloatFormat("%g")
        SetQueryOutputToValue()
        numNodes = Query("NumNodes")
        if numNodes < 10:
            #HideActivePlots()
            #DeleteActivePlots()
            DeleteAllPlots()
            return True
        i = i + 1
    return False

def set_view():
    sys.stderr.write('=============set_view\n')
    ResetView()
    View3DAtts = GetView3D()
    View3DAtts.viewNormal = (0.82, 0.53, 0.25)
    View3DAtts.viewUp = (-0.19, -0.16, 0.97)
    SetView3D(View3DAtts)


def setup_plot1():
    sys.stderr.write('=============setup_plot1\n')
    AddPlot("Pseudocolor", "Electron-0/uz", 1, 0)
    SetActivePlots(0)
    SetActivePlots(0)
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
    DrawPlots()

def setup_plot2():
    sys.stderr.write('=============setup_plot2\n')
    AddPlot("Pseudocolor", "Electron-0/uy", 1, 0)
    SetActivePlots(0)
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
    DrawPlots()

def setup_plot3():
    sys.stderr.write('=============setup_plot3\n')
    AddPlot("Pseudocolor", "Electron-1/uz", 1, 0)
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
    ElevateAtts.minFlag = 0
    ElevateAtts.min = -5e+08
    ElevateAtts.maxFlag = 0
    ElevateAtts.max = 5e+08
    ElevateAtts.zeroFlag = 0
    ElevateAtts.variable = "Electron-1/ux"
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
    DrawPlots()

def setup_plot4():
    sys.stderr.write('=============setup_plot4\n')
    #########################################
    # elec density                          #
    #########################################
    #AddPlot("Pseudocolor", "operators/DataBinning/2D/Electron-0", 1, 1)
    AddPlot("Pseudocolor", "Electron-0/weights", 1, 0)
    #
    AddOperator("DataBinning",0)
    SetActivePlots(0)
    #
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
    #DataBinningAtts.outputType = DataBinningAtts.OutputOnInputMesh  # OutputOnBins, OutputOnInputMesh
    DataBinningAtts.outputType = DataBinningAtts.OutputOnBins  # OutputOnBins, OutputOnInputMesh
    DataBinningAtts.removeEmptyValFromCurve = 1
    SetOperatorOptions(DataBinningAtts, 1)
    #AddOperator("Project", 1)
    #ProjectAtts = ProjectAttributes()
    #ProjectAtts.projectionType = ProjectAtts.XZCartesian  # ZYCartesian, XZCartesian, XYCartesian, XRCylindrical, YRCylindrical, ZRCylindrical
    #ProjectAtts.vectorTransformMethod = ProjectAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
    #SetOperatorOptions(ProjectAtts, 1)
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
    #
    ChangeActivePlotsVar("operators/DataBinning/2D/Electron-0")
    #
    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
    PseudocolorAtts.skewFactor = 1
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
    PseudocolorAtts.minFlag = 0
    PseudocolorAtts.min = 0
    PseudocolorAtts.maxFlag = 0
    PseudocolorAtts.max = 1
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
    ElevateAtts.maxFlag = 0
    ElevateAtts.max = 7e+08
    ElevateAtts.zeroFlag = 0
    ElevateAtts.variable = "default"
    SetOperatorOptions(ElevateAtts, 0)
    AddOperator("Transform", 0)
    TransformAtts = TransformAttributes()
    TransformAtts.doRotate = 0
    TransformAtts.rotateOrigin = (0, 0, 0)
    TransformAtts.rotateAxis = (0, 0, 1)
    TransformAtts.rotateAmount = 0
    TransformAtts.rotateType = TransformAtts.Deg  # Deg, Rad
    TransformAtts.doScale = 1
    TransformAtts.scaleOrigin = (0, 0, 0)
    TransformAtts.scaleX = 1
    TransformAtts.scaleY = 1
    TransformAtts.scaleZ = 10
    TransformAtts.doTranslate = 0
    TransformAtts.translateX = 0
    TransformAtts.translateY = 0
    TransformAtts.translateZ = 0
    TransformAtts.transformType = TransformAtts.Similarity  # Similarity, Coordinate, Linear
    TransformAtts.inputCoordSys = TransformAtts.Cartesian  # Cartesian, Cylindrical, Spherical
    TransformAtts.outputCoordSys = TransformAtts.Spherical  # Cartesian, Cylindrical, Spherical
    TransformAtts.m00 = 1
    TransformAtts.m01 = 0
    TransformAtts.m02 = 0
    TransformAtts.m03 = 0
    TransformAtts.m10 = 0
    TransformAtts.m11 = 1
    TransformAtts.m12 = 0
    TransformAtts.m13 = 0
    TransformAtts.m20 = 0
    TransformAtts.m21 = 0
    TransformAtts.m22 = 1
    TransformAtts.m23 = 0
    TransformAtts.m30 = 0
    TransformAtts.m31 = 0
    TransformAtts.m32 = 0
    TransformAtts.m33 = 1
    TransformAtts.invertLinearTransform = 0
    TransformAtts.vectorTransformMethod = TransformAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
    TransformAtts.transformVectors = 1
    SetOperatorOptions(TransformAtts, 0)
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
    sys.stderr.write('=============set_annotations\n')
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

def save_window_matrix(omitWin1, omitWin2, omitWin3, omitWin4):
    sys.stderr.write('=============omit_window_matrix\n')
    if (omitWin1 and omitWin2 and omitWin3 and omitWin4):
        return

    row = lambda i: (int(i)%2)
    col = lambda i: (int(i)/2)

    x = 0
    y = 0
    i = int(0)
    dx = 1600
    dy = 1100

    SaveWindowAtts = SaveWindowAttributes()

    if (not omitWin3):
        SaveWindowAtts.subWindowAtts.win3.position = (row(i)*dx, col(i)*dy)
        SaveWindowAtts.subWindowAtts.win3.size = (dx,dy)
        SaveWindowAtts.subWindowAtts.win3.layer = 0
        SaveWindowAtts.subWindowAtts.win3.transparency = 0
        SaveWindowAtts.subWindowAtts.win3.omitWindow = 0
        i+=1
    if (not omitWin4):
        SaveWindowAtts.subWindowAtts.win4.position = (row(i)*dx, col(i)*dy)
        SaveWindowAtts.subWindowAtts.win4.size = (dx,dy)
        SaveWindowAtts.subWindowAtts.win4.layer = 0
        SaveWindowAtts.subWindowAtts.win4.transparency = 0
        SaveWindowAtts.subWindowAtts.win4.omitWindow = 0
        i+=1
    if (not omitWin1):
        SaveWindowAtts.subWindowAtts.win1.position = (row(i)*dx, col(i)*dy)
        SaveWindowAtts.subWindowAtts.win1.size = (dx,dy)
        SaveWindowAtts.subWindowAtts.win1.layer = 0
        SaveWindowAtts.subWindowAtts.win1.transparency = 0
        SaveWindowAtts.subWindowAtts.win1.omitWindow = 0
        i+=1
    if (not omitWin2):
        SaveWindowAtts.subWindowAtts.win2.position = (row(i)*dx, col(i)*dy)
        SaveWindowAtts.subWindowAtts.win2.size = (dx,dy)
        SaveWindowAtts.subWindowAtts.win2.layer = 0
        SaveWindowAtts.subWindowAtts.win2.transparency = 0
        SaveWindowAtts.subWindowAtts.win2.omitWindow = 0
        i+=1

    imw = 1600
    if (i>0):
        imw = 3200

    imh = 1100
    if (i>2):
        imh = 2200

    SaveWindowAtts.outputToCurrentDirectory = 1
    SaveWindowAtts.outputDirectory = "./"
    SaveWindowAtts.fileName = "lpa-4-view"
    SaveWindowAtts.family = 1
    SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
    SaveWindowAtts.width = imw
    SaveWindowAtts.height = imh

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

def delete_plots():
    SetActiveWindow(1)
    SetActiveWindow(1)
    DeleteAllPlots()
    SetActiveWindow(2)
    SetActiveWindow(2)
    DeleteAllPlots()
    SetActiveWindow(3)
    SetActiveWindow(3)
    DeleteAllPlots()
    SetActiveWindow(4)
    SetActiveWindow(4)
    DeleteAllPlots()


#SendSimulationCommand('localhost', simFile, 'pause')
SetWindowLayout(4)
delete_plots()
#SetWindowLayout(1)
#SetWindowLayout(4)

SetActiveWindow(1)
SetActiveWindow(1)
setup_plot1()
omitWin1 = hide_if_empty()
if (not omitWin1):
    set_annotations(xAxisName='X', yAxisName='Z', zAxisName='Uy', showDB=1)
    set_view()

SetActiveWindow(2)
SetActiveWindow(2)
setup_plot2()
omitWin2 = hide_if_empty()
if (not omitWin2):
    set_annotations(xAxisName='X', yAxisName='Z', zAxisName='Uz', showDB=1)
    set_view()

SetActiveWindow(3)
SetActiveWindow(3)
setup_plot3()
omitWin3 = hide_if_empty()
if (not omitWin3):
    set_annotations(xAxisName='X', yAxisName='Z', zAxisName='Ux', showDB=1)
    set_view()

SetActiveWindow(4)
SetActiveWindow(4)
setup_plot4()
omitWin4 = hide_if_empty()
if (not omitWin4):
    set_annotations(xAxisName='X', yAxisName='Z', zAxisName='Uz', showDB=1)
    set_view()

sys.stderr.flush()

save_window_matrix(omitWin1, omitWin2, omitWin3, omitWin4)
delete_plots()
