from visit import *
DeleteAllPlots()

wat = GetSaveWindowAttributes()
wat.fileName = 'cigar-part-and-field-'
SetSaveWindowAttributes(wat)

AnnotationAtts = GetAnnotationAttributes()
AnnotationAtts.axes3D.visible = 0
AnnotationAtts.axes3D.bboxFlag = 0
AnnotationAtts.userInfoFlag = 0
SetAnnotationAttributes(AnnotationAtts)

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

AddPlot("Pseudocolor", "particle/vx")
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

AddPlot("Pseudocolor", "grid/rho")
#AddOperator("Clip", 0)
#ClipAtts = ClipAttributes()
#ClipAtts.funcType = ClipAtts.Plane
#ClipAtts.plane1Status = 1
#ClipAtts.plane2Status = 0
#ClipAtts.plane3Status = 0
#ClipAtts.plane1Origin = (-0.15, 0, 0)
#ClipAtts.plane1Normal = (1, 0, 0)
#ClipAtts.planeInverse = 1
#ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1
#ClipAtts.center = (0, 0, 0)
#ClipAtts.radius = 1
#ClipAtts.sphereInverse = 0
#SetOperatorOptions(ClipAtts, 0)

view = GetView3D()
view.viewNormal = (-0.54, 0.32, 0.78)
view.focus = (0.0098, 0.0098, 0)
view.viewUp = (0.04248, 0.933108, -0.35707)
view.nearPlane = -1.4
view.farPlane = 1.4
SetView3D(view)

DrawPlots()
SaveWindow()
