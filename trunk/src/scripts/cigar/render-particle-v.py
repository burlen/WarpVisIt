from visit import *
DeleteAllPlots()
ResetView()

wat = GetSaveWindowAttributes()
wat.fileName = 'particle-v-'
wat.width = 1024
wat.height = 768
SetSaveWindowAttributes(wat)

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

RenderingAtts = GetRenderingAttributes()
RenderingAtts.antialiasing = 1
RenderingAtts.scalableActivationMode = RenderingAtts.Always
RenderingAtts.doShadowing = 1
RenderingAtts.shadowStrength = 0.5
SetRenderingAttributes(RenderingAtts)

light0 = LightAttributes()
light0.enabledFlag = 1
light0.type = light0.Camera
light0.direction = (0.2, -0.2, -0.8)
light0.color = (255, 255, 255, 255)
light0.brightness = 1
SetLight(0, light0)

AddPlot("Pseudocolor", "particle/V")
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.legendFlag = 1
PseudocolorAtts.lightingFlag = 1
PseudocolorAtts.pointSize = 0.005
PseudocolorAtts.pointType = PseudocolorAtts.Icosahedron
PseudocolorAtts.opacity = 1.0
PseudocolorAtts.colorTableName = "hot_desaturated"
PseudocolorAtts.pointSizePixels = 2
PseudocolorAtts.opacityType = PseudocolorAtts.Constant
SetPlotOptions(PseudocolorAtts)

aat = AnnotationAttributes()
aat.userInfoFlag = 0
aat.databaseInfoFlag = 0
aat.timeInfoFlag = 0
aat.axes3D.visible = 1
aat.axes3D.bboxFlag = 1
aat.userInfoFlag = 0
SetAnnotationAttributes(aat)

plotName = GetPlotList().GetPlots(0).plotName
legend = GetAnnotationObject(plotName)
legend.orientation = legend.VerticalLeft  # VerticalRight, VerticalLeft, HorizontalTop, Horiz
legend.fontBold = 1

banner = CreateAnnotationObject("Text2D")
banner.text = "particle V  t=$time"
banner.position = (0.2, 0.93)
banner.fontBold = 1
banner.height = 0.02

DrawPlots()
SaveWindow()

banner.Delete()
DeleteAllPlots()
ResetView()
