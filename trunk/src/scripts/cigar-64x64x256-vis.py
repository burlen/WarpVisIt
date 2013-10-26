import time
from visit import *

def Render():
    """
    Example of script to control VisIt rendering.
    """
    visit.DeleteAllPlots()

    wat = visit.GetSaveWindowAttributes()
    wat.fileName = 'cigar'
    visit.SetSaveWindowAttributes(wat)

    AnnotationAtts = visit.GetAnnotationAttributes()
    AnnotationAtts.axes3D.visible = 0
    AnnotationAtts.axes3D.bboxFlag = 0
    AnnotationAtts.userInfoFlag = 0
    visit.SetAnnotationAttributes(AnnotationAtts)

    RenderingAtts = visit.GetRenderingAttributes()
    RenderingAtts.antialiasing = 1
    RenderingAtts.scalableActivationMode = RenderingAtts.Always
    RenderingAtts.doShadowing = 1
    RenderingAtts.shadowStrength = 0.5
    visit.SetRenderingAttributes(RenderingAtts)

    light0 = visit.LightAttributes()
    light0.enabledFlag = 1
    light0.type = light0.Camera
    light0.direction = ( 0.2, -0.2, -0.8)
    light0.color = (255, 255, 255, 255)
    light0.brightness = 1
    visit.SetLight(0, light0)

    visit.AddPlot("Pseudocolor", "particle/vx")
    PseudocolorAtts = visit.PseudocolorAttributes()
    PseudocolorAtts.legendFlag = 1
    PseudocolorAtts.lightingFlag = 1
    PseudocolorAtts.pointSize = 0.005
    PseudocolorAtts.pointType = PseudocolorAtts.Icosahedron
    PseudocolorAtts.opacity = 1.0
    PseudocolorAtts.colorTableName = "hot_desaturated"
    PseudocolorAtts.pointSizePixels = 2
    PseudocolorAtts.opacityType = PseudocolorAtts.Constant
    visit.SetPlotOptions(PseudocolorAtts)

    visit.AddPlot("Pseudocolor", "grid/rho")
    #visit.AddOperator("Clip", 0)
    #ClipAtts = visit.ClipAttributes()
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
    #visit.SetOperatorOptions(ClipAtts, 0)

    view = visit.GetView3D()
    view.viewNormal = (-0.54, 0.32, 0.78)
    view.focus = (0.0098, 0.0098, 0)
    view.viewUp = (0.04248, 0.933108, -0.35707)
    view.nearPlane = -1.4
    view.farPlane = 1.4
    visit.SetView3D(view)

    visit.DrawPlots()

    visit.SaveWindow()
