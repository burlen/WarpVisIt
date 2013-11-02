from visit import *
DeleteAllPlots()
ResetView()

wat = GetSaveWindowAttributes()
wat.fileName = 'volume-phi-'
wat.width = 1024
wat.height = 768
SetSaveWindowAttributes(wat)

AddPlot("Volume", "grid/phi", 1, 0)
SetActivePlots(0)

AddOperator("Reflect", 0)
ReflectAtts = ReflectAttributes()
ReflectAtts.octant = ReflectAtts.PXPYPZ  # PXPYPZ, NXPYPZ, PXNYPZ, NXNYPZ, PXPYNZ, NXPYNZ, PXNYNZ, NXNYNZ
ReflectAtts.useXBoundary = 1
ReflectAtts.specifiedX = 0
ReflectAtts.useYBoundary = 1
ReflectAtts.specifiedY = 0
ReflectAtts.useZBoundary = 1
ReflectAtts.specifiedZ = 0
ReflectAtts.reflections = (1, 1, 1, 1, 0, 0, 0, 0)
SetOperatorOptions(ReflectAtts, 0)

VolumeAtts = VolumeAttributes()
VolumeAtts.legendFlag = 1
VolumeAtts.lightingFlag = 1
VolumeAtts.colorControlPoints.GetControlPoints(0).colors = (0, 0, 255, 255)
VolumeAtts.colorControlPoints.GetControlPoints(0).position = 0
VolumeAtts.colorControlPoints.GetControlPoints(1).colors = (0, 255, 255, 255)
VolumeAtts.colorControlPoints.GetControlPoints(1).position = 0.25
VolumeAtts.colorControlPoints.GetControlPoints(2).colors = (0, 255, 0, 255)
VolumeAtts.colorControlPoints.GetControlPoints(2).position = 0.5
VolumeAtts.colorControlPoints.GetControlPoints(3).colors = (255, 255, 0, 255)
VolumeAtts.colorControlPoints.GetControlPoints(3).position = 0.75
VolumeAtts.colorControlPoints.GetControlPoints(4).colors = (255, 0, 0, 255)
VolumeAtts.colorControlPoints.GetControlPoints(4).position = 1
VolumeAtts.colorControlPoints.smoothing = VolumeAtts.colorControlPoints.Linear  # None, Linear, CubicSpline
VolumeAtts.colorControlPoints.equalSpacingFlag = 0
VolumeAtts.colorControlPoints.discreteFlag = 0
VolumeAtts.colorControlPoints.externalFlag = 0
VolumeAtts.opacityAttenuation = 1
VolumeAtts.opacityMode = VolumeAtts.FreeformMode  # FreeformMode, GaussianMode, ColorTableMode
#controlPoints does not contain any GaussianControlPoint objects.
VolumeAtts.resampleFlag = 1
VolumeAtts.resampleTarget = 50000
VolumeAtts.opacityVariable = "default"
VolumeAtts.compactVariable = "default"
VolumeAtts.freeformOpacity = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255)
VolumeAtts.useColorVarMin = 0
VolumeAtts.colorVarMin = 0
VolumeAtts.useColorVarMax = 0
VolumeAtts.colorVarMax = 0
VolumeAtts.useOpacityVarMin = 0
VolumeAtts.opacityVarMin = 0
VolumeAtts.useOpacityVarMax = 0
VolumeAtts.opacityVarMax = 0
VolumeAtts.smoothData = 0
VolumeAtts.samplesPerRay = 500
VolumeAtts.rendererType = VolumeAtts.RayCasting  # Splatting, Texture3D, RayCasting, RayCastingIntegration, SLIVR, RayCastingSLIVR, Tuvok
VolumeAtts.gradientType = VolumeAtts.SobelOperator  # CenteredDifferences, SobelOperator
VolumeAtts.num3DSlices = 200
VolumeAtts.scaling = VolumeAtts.Linear  # Linear, Log, Skew
VolumeAtts.skewFactor = 1
VolumeAtts.limitsMode = VolumeAtts.OriginalData  # OriginalData, CurrentPlot
VolumeAtts.sampling = VolumeAtts.Rasterization  # KernelBased, Rasterization, Trilinear
VolumeAtts.rendererSamples = 3
#transferFunction2DWidgets does not contain any TransferFunctionWidget objects.
VolumeAtts.transferFunctionDim = 1
VolumeAtts.lowGradientLightingReduction = VolumeAtts.Lower  # Off, Lowest, Lower, Low, Medium, High, Higher, Highest
VolumeAtts.lowGradientLightingClampFlag = 0
VolumeAtts.lowGradientLightingClampValue = 1
VolumeAtts.materialProperties = (0.4, 0.75, 0, 15)
SetPlotOptions(VolumeAtts)

aat = AnnotationAttributes()
aat.userInfoFlag = 0
aat.databaseInfoFlag = 0
aat.timeInfoFlag = 0
aat.axes3D.visible = 0
SetAnnotationAttributes(aat)

plotName = GetPlotList().GetPlots(0).plotName
legend = GetAnnotationObject(plotName)
legend.orientation = legend.VerticalLeft  # VerticalRight, VerticalLeft, HorizontalTop, Horiz
legend.fontBold = 1

banner = CreateAnnotationObject("Text2D")
banner.text = "volume phi  t=$time"
banner.position = (0.2, 0.93)
banner.fontBold = 1
banner.height = 0.02

View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (-0.677904, 0.539731, 0.499136)
View3DAtts.focus = (0.017225, 0, 0)
View3DAtts.viewUp = (0.220775, 0.797084, -0.562064)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 0.72103
View3DAtts.nearPlane = -1.44206
View3DAtts.farPlane = 1.44206
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1.0
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (0.017225, 0, 0)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
SetView3D(View3DAtts)

DrawPlots()
SaveWindow()

banner.Delete()
DeleteAllPlots()
ResetView()
