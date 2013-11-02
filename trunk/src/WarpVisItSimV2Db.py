import numpy as np
import warp
import parallel
import simV2
from WarpVisItUtil import pError
from WarpVisItUtil import pDebug

# VisIt has a nice feature that if you use
# / in a variable name you can construct nested
# menus and organize variables by mesh. However
# this breaks the data binning operator.
varsep = '/'

#-----------------------------------------------------------------------------
def getTime():
    """Return current simulation time"""
    t = warp.top.time
    pDebug('time = %f'%(t))
    return t

#-----------------------------------------------------------------------------
def getIteration():
    """Return current simulation iteration"""
    it = warp.top.it
    pDebug('iteration = %i'%(it))
    return it

#-----------------------------------------------------------------------------
def getMetaData(userData):
    """
    Callback function used to provide visit with metadata
    """
    valid = lambda x : x != simV2.VISIT_INVALID_HANDLE
    pDebug('getMetaData')

    # VisIt assumes that each process has some data
    # so you must have at least as many domains as
    # processes. you just return invalid handle in
    # getMesh/getVar on those processes that do not
    nDomains = parallel.number_of_PE()

    simmd = simV2.VisIt_SimulationMetaData_alloc()
    if not valid(simmd):
        pError('VisIt_SimulationMetaData_alloc failed')
        return None

    # time
    simV2.VisIt_SimulationMetaData_setCycleTime(simmd,
          warp.top.it,
          warp.top.time)

    # particle mesh
    meshmd = simV2.VisIt_MeshMetaData_alloc()
    if not valid(meshmd):
        pError('VisIt_MeshMetaData_alloc failed')
        return None
    simV2.VisIt_MeshMetaData_setName(meshmd, 'particle')
    simV2.VisIt_MeshMetaData_setMeshType(meshmd, simV2.VISIT_MESHTYPE_POINT)
    simV2.VisIt_MeshMetaData_setTopologicalDimension(meshmd, 0)
    simV2.VisIt_MeshMetaData_setSpatialDimension(meshmd, 3)
    simV2.VisIt_MeshMetaData_setNumDomains(meshmd, nDomains)
    simV2.VisIt_MeshMetaData_setDomainTitle(meshmd, 'Domains')
    simV2.VisIt_MeshMetaData_setDomainPieceName(meshmd, 'domain')
    simV2.VisIt_MeshMetaData_setNumGroups(meshmd, 0)
    simV2.VisIt_MeshMetaData_setXUnits(meshmd, 'm')
    simV2.VisIt_MeshMetaData_setYUnits(meshmd, 'm')
    simV2.VisIt_MeshMetaData_setZUnits(meshmd, 'm')
    simV2.VisIt_MeshMetaData_setXLabel(meshmd, 'x')
    simV2.VisIt_MeshMetaData_setYLabel(meshmd, 'y')
    simV2.VisIt_MeshMetaData_setZLabel(meshmd, 'z')
    simV2.VisIt_SimulationMetaData_addMesh(simmd, meshmd)

    # field mesh
    meshmd = simV2.VisIt_MeshMetaData_alloc()
    if not valid(meshmd):
        pError('VisIt_MeshMetaData_alloc failed')
        return None
    simV2.VisIt_MeshMetaData_setName(meshmd, 'grid')
    simV2.VisIt_MeshMetaData_setMeshType(meshmd, simV2.VISIT_MESHTYPE_RECTILINEAR)
    simV2.VisIt_MeshMetaData_setTopologicalDimension(meshmd, 3)
    simV2.VisIt_MeshMetaData_setSpatialDimension(meshmd, 3)
    simV2.VisIt_MeshMetaData_setNumDomains(meshmd, nDomains)
    simV2.VisIt_MeshMetaData_setDomainTitle(meshmd, 'Domains')
    simV2.VisIt_MeshMetaData_setDomainPieceName(meshmd, 'domain')
    simV2.VisIt_MeshMetaData_setNumGroups(meshmd, 0)
    simV2.VisIt_MeshMetaData_setXUnits(meshmd, 'm')
    simV2.VisIt_MeshMetaData_setYUnits(meshmd, 'm')
    simV2.VisIt_MeshMetaData_setZUnits(meshmd, 'm')
    simV2.VisIt_MeshMetaData_setXLabel(meshmd, 'x')
    simV2.VisIt_MeshMetaData_setYLabel(meshmd, 'y')
    simV2.VisIt_MeshMetaData_setZLabel(meshmd, 'z')
    simV2.VisIt_SimulationMetaData_addMesh(simmd, meshmd)

    # CSG mesh
    meshmd = simV2.VisIt_MeshMetaData_alloc()
    if not valid(meshmd):
        pError('VisIt_MeshMetaData_alloc failed')
        return None
    simV2.VisIt_MeshMetaData_setName(meshmd, 'csg')
    simV2.VisIt_MeshMetaData_setMeshType(meshmd, simV2.VISIT_MESHTYPE_CSG)
    simV2.VisIt_MeshMetaData_setTopologicalDimension(meshmd, 3)
    simV2.VisIt_MeshMetaData_setSpatialDimension(meshmd, 3)
    simV2.VisIt_MeshMetaData_setNumDomains(meshmd, nDomains)
    simV2.VisIt_MeshMetaData_setDomainTitle(meshmd, 'Domains')
    simV2.VisIt_MeshMetaData_setDomainPieceName(meshmd, 'domain')
    simV2.VisIt_MeshMetaData_setNumGroups(meshmd, 0)
    simV2.VisIt_MeshMetaData_setXUnits(meshmd, 'm')
    simV2.VisIt_MeshMetaData_setYUnits(meshmd, 'm')
    simV2.VisIt_MeshMetaData_setZUnits(meshmd, 'm')
    simV2.VisIt_MeshMetaData_setXLabel(meshmd, 'x')
    simV2.VisIt_MeshMetaData_setYLabel(meshmd, 'y')
    simV2.VisIt_MeshMetaData_setZLabel(meshmd, 'z')
    simV2.VisIt_SimulationMetaData_addMesh(simmd, meshmd)

    # particle data arrays
    for var in ['pid','r','rp','theta','w','vdrifts',
        'vx','vy','vz','vr','ux','uy','uz','etheta',
        'er','ex','ey','ez','bx','by','bz','xp','yp',
        'zp','tp','rank']:
        vmd = simV2.VisIt_VariableMetaData_alloc()
        if not valid(vmd):
            pError('VisIt_VariableMetaData_alloc failed')
            return None
        varname = 'particle%s%s'%(varsep,var)
        simV2.VisIt_VariableMetaData_setName(vmd, varname)
        simV2.VisIt_VariableMetaData_setMeshName(vmd, 'particle')
        simV2.VisIt_VariableMetaData_setType(vmd, simV2.VISIT_VARTYPE_SCALAR)
        simV2.VisIt_VariableMetaData_setCentering(vmd, simV2.VISIT_VARCENTERING_NODE)
        simV2.VisIt_SimulationMetaData_addVariable(simmd, vmd)

    # field data arrays
    for var in ['ax','ay','az','A',
        'bx','by','bz','B','J','jx','jy',
        'jz','phi','rho','rank']:
        vmd = simV2.VisIt_VariableMetaData_alloc()
        if not valid(vmd):
            pError('VisIt_VariableMetaData_alloc failed')
            return None
        varname = 'grid%s%s'%(varsep,var)
        simV2.VisIt_VariableMetaData_setName(vmd, varname)
        simV2.VisIt_VariableMetaData_setMeshName(vmd, 'grid')
        simV2.VisIt_VariableMetaData_setType(vmd, simV2.VISIT_VARTYPE_SCALAR)
        simV2.VisIt_VariableMetaData_setCentering(vmd, simV2.VISIT_VARCENTERING_NODE)
        simV2.VisIt_SimulationMetaData_addVariable(simmd, vmd)

    # expressions
    expname = lambda x : x[0]
    expstr = lambda x : x[1]
    exptype = lambda x : x[2]
    for exp in [('{0}/v','{{<{0}/vx>, <{0}/vy>, <{0}/vz>}}',simV2.VISIT_VARTYPE_VECTOR),
          ('{0}/V','sqrt(<{0}/vx>*<{0}/vx>+<{0}/vy>*<{0}/vy>+<{0}/vz>*<{0}/vz>)',simV2.VISIT_VARTYPE_SCALAR),
          ('{0}/B','sqrt(<{0}/bx>*<{0}/bx>+<{0}/by>*<{0}/by>+<{0}/bz>*<{0}/bz>)',simV2.VISIT_VARTYPE_SCALAR),
          ('{0}/E','sqrt(<{0}/ex>*<{0}/ex>+<{0}/ey>*<{0}/ey>+<{0}/ez>*<{0}/ez>)',simV2.VISIT_VARTYPE_SCALAR) ]:
#          ('const','coord(csg)[0]*0',simV2.VISIT_VARTYPE_SCALAR)]:
#          ('particle proc id','procid(particles)',simV2.VISIT_VARTYPE_SCALAR),
#          ('field proc id','procid(fields)',simV2.VISIT_VARTYPE_SCALAR)]:

        expmd = simV2.VisIt_ExpressionMetaData_alloc()
        if not valid(expmd):
            pError('VisIt_ExpressionMetaData_alloc failed')
            return None
        simV2.VisIt_ExpressionMetaData_setName(expmd, expname(exp).format('particle'))
        simV2.VisIt_ExpressionMetaData_setDefinition(expmd, expstr(exp).format('particle'))
        simV2.VisIt_ExpressionMetaData_setType(expmd, exptype(exp))
        simV2.VisIt_SimulationMetaData_addExpression(simmd, expmd)

    # commands
    for cmd in ['halt', 'step', 'run', 'run_without_update', 'end']:
        cmdmd = simV2.VisIt_CommandMetaData_alloc()
        if not valid(cmdmd):
            pError('VisIt_CommandMetaData_alloc failed')
            return None
        simV2.VisIt_CommandMetaData_setName(cmdmd, cmd)
        simV2.VisIt_SimulationMetaData_addGenericCommand(simmd, cmdmd)

    return simmd

#-----------------------------------------------------------------------------
def getMesh(domain, name, userData):
    """
    Callback function used to send mesh data (e.g., particles) to VisIt.
    """
    # VisIt assumes that each process has some data
    # so you must have at least as many domains as
    # processes. you just return invalid handle in
    # getMesh/getVar on those processes that do not
    valid = lambda x : x != simV2.VISIT_INVALID_HANDLE

    pDebug('getMesh %i %s'%(domain, name))

    # particle mesh
    if name == 'particle':
        x = warp.getx(gather=0)
        if not x.size:
            return simV2.VISIT_INVALID_HANDLE
        pDebug(str(x))
        n = x.size
        xCopy = x.tolist() # FIXME -- makes a copy
        xvd = simV2.VisIt_VariableData_alloc()
        if not valid(xvd):
            pError('VisIt_VariableData_alloc failed')
            return None
        simV2.VisIt_VariableData_setDataF(xvd, simV2.VISIT_OWNER_SIM, 1, n, xCopy)

        y = warp.gety(gather=0)
        pDebug(str(y))
        yCopy = y.tolist() # FIXME
        yvd = simV2.VisIt_VariableData_alloc()
        if not valid(yvd):
            pError('VisIt_VariableData_alloc failed')
            return None
        simV2.VisIt_VariableData_setDataF(yvd, simV2.VISIT_OWNER_SIM, 1, n, yCopy)

        z = warp.getz(gather=0)
        pDebug(str(z))
        zCopy = z.tolist() # FIXME
        zvd = simV2.VisIt_VariableData_alloc()
        if not valid(zvd):
            pError('VisIt_VariableData_alloc failed')
            return None
        simV2.VisIt_VariableData_setDataF(zvd, simV2.VISIT_OWNER_SIM, 1, n, zCopy)

        mesh = simV2.VisIt_PointMesh_alloc()
        if not valid(mesh):
            pError('VisIt_PointMesh_alloc failed')
            return None

        simV2.VisIt_PointMesh_setCoordsXYZ(mesh, xvd, yvd, zvd)
        return mesh

    # uniform mesh
    if name == 'grid':
        size = getGridSize()
        coords = getGridCoordinates()

        xvd = []
        i = 0
        while i < 3:
            vd = simV2.VisIt_VariableData_alloc()
            if not valid(vd):
                pError('VisIt_VariableData_alloc failed')
                return None
            simV2.VisIt_VariableData_setDataF(vd,
                simV2.VISIT_OWNER_SIM, 1, int(size[i]),
                coords[i].tolist()) # FIXME -- again copy
            xvd.append(vd)
            i += 1

        mesh = simV2.VisIt_RectilinearMesh_alloc()
        if not valid(mesh):
            pError('VisIt_RectilinearMesh_alloc failed')
            return None

        simV2.VisIt_RectilinearMesh_setCoordsXYZ(mesh, xvd[0], xvd[1], xvd[2])
        return mesh

    # CSG mesh
    if name == 'csg':
        # FIXME -- parameters must be user setable
        piperad_outer = 3.445e-2 + 0.3
        piperad_inner = piperad_outer * 0.8
        pipelength = 2.0 #4.

        zmin = warp.top.zbeam - pipelength/2.
        zmax = warp.top.zbeam + pipelength/2.
        ext_min = [-piperad_outer , -piperad_outer , zmin]
        ext_max = [ piperad_outer ,  piperad_outer , zmax]
        bound_types = [simV2.VISIT_CSG_CYLINDER_PNLR, simV2.VISIT_CSG_CYLINDER_PNLR]
        bound_coeffs = [0., 0., zmin, 0., 0., 1., pipelength, piperad_outer,
                        0., 0., zmax, 0., 0., 1., pipelength, piperad_inner]
        region_operators = [simV2.VISIT_CSG_INNER,
              simV2.VISIT_CSG_OUTER,
              simV2.VISIT_CSG_INTERSECT]
        leftids =  [ 0,  1, 0]
        rightids = [-1, -1, 1]
        zonelist = [2]

        mesh = simV2.VisIt_CSGMesh_alloc()
        if not valid(mesh):
            pError('VisIt_CSGMesh_alloc failed')
            return None

        cbt = simV2.VisIt_VariableData_alloc()
        if not valid(cbt):
            pError('VisIt_VariableData_alloc failed')
            return None
        simV2.VisIt_VariableData_setDataI(cbt, simV2.VISIT_OWNER_COPY, 1, 2, bound_types)
        simV2.VisIt_CSGMesh_setBoundaryTypes(mesh, cbt)

        cbc = simV2.VisIt_VariableData_alloc()
        if not valid(cbc):
            pError('VisIt_VariableData_alloc failed')
            return None
        simV2.VisIt_VariableData_setDataF(cbc, simV2.VISIT_OWNER_COPY, 1, 16, bound_coeffs)
        simV2.VisIt_CSGMesh_setBoundaryCoeffs(mesh, cbc)

        simV2.VisIt_CSGMesh_setExtents(mesh, ext_min, ext_max)

        cro = simV2.VisIt_VariableData_alloc()
        if not valid(cro):
            pError('VisIt_VariableData_alloc failed')
            return None
        simV2.VisIt_VariableData_setDataI(cro, simV2.VISIT_OWNER_COPY, 1, 3, region_operators)

        cli = simV2.VisIt_VariableData_alloc()
        if not valid(cli):
            pError('VisIt_VariableData_alloc failed')
            return None
        simV2.VisIt_VariableData_setDataI(cli, simV2.VISIT_OWNER_COPY, 1, 3, leftids)

        cri = simV2.VisIt_VariableData_alloc()
        if not valid(cri):
            pError('VisIt_VariableData_alloc failed')
            return None
        simV2.VisIt_VariableData_setDataI(cri, simV2.VISIT_OWNER_COPY, 1, 3, rightids)

        simV2.VisIt_CSGMesh_setRegions(mesh, cro , cli, cri)

        czl = simV2.VisIt_VariableData_alloc()
        if not valid(czl):
            pError('VisIt_VariableData_alloc failed')
            return None
        simV2.VisIt_VariableData_setDataI(czl, simV2.VISIT_OWNER_COPY, 1, 1, zonelist)
        simV2.VisIt_CSGMesh_setZonelist(mesh, czl)
        return mesh
    # invalid
    pError('Unrecognized mesh name %s'%(name))
    return None

#-----------------------------------------------------------------------------
def getVar(domain, varid, userData):
    """
    Callback function used to send variable data (e.g., vx) to VisIt.
    """
    valid = lambda x : x != simV2.VISIT_INVALID_HANDLE
    # VisIt assumes that each process has some data
    # so you must have at least as many domains as
    # processes. you just return invalid handle in
    # getMesh/getVar on those processes that do not
    pDebug('getVar %i %s'%(domain, varid))

    tok = varid.split(varsep)
    meshname = tok[0]
    varname = tok[1]

    try:

        # particle data
        if meshname == 'particle':
            # pid
            if varname == 'pid':
                return passParticleData(warp.getpid(gather=0))
            # w
            if varname == 'w':
                return passParticleData(warp.getw(gather=0))
            # r
            if varname == 'r':
                return passParticleData(warp.getr(gather=0))
            # rp
            if varname == 'rp':
                return passParticleData(warp.getrp(gather=0))
            # theta
            if varname == 'theta':
                return passParticleData(warp.gettheta(gather=0))
            # vdrifts
            if varname == 'vdrifts':
                return passParticleData(warp.getvdrifts(gather=0))
            # vx
            if varname == 'vx':
                return passParticleData(warp.getvx(gather=0))
            # vy
            if varname == 'vy':
                return passParticleData(warp.getvy(gather=0))
            # vz
            if varname == 'vz':
                return passParticleData(warp.getvz(gather=0))
            # vr
            if varname == 'vr':
                return passParticleData(warp.getvr(gather=0))
            # ux
            if varname == 'ux':
                return passParticleData(warp.getux(gather=0))
            # uy
            if varname == 'uy':
                return passParticleData(warp.getuy(gather=0))
            # uz
            if varname == 'uz':
                return passParticleData(warp.getuz(gather=0))
            # etheta
            if varname == 'etheta':
                return passParticleData(warp.getetheta(gather=0))
            # er
            if varname == 'er':
                return passParticleData(warp.geter(gather=0))
            # ex
            if varname == 'ex':
                return passParticleData(warp.getex(gather=0))
            # ey
            if varname == 'ey':
                return passParticleData(warp.getey(gather=0))
            # ez
            if varname == 'ez':
                return passParticleData(warp.getez(gather=0))
            # bx
            if varname == 'bx':
                return passParticleData(warp.getbx(gather=0))
            # by
            if varname == 'by':
                return passParticleData(warp.getby(gather=0))
            # bz
            if varname == 'bz':
                return passParticleData(warp.getbz(gather=0))
            # xp
            if varname == 'xp':
                return passParticleData(warp.getxp(gather=0))
            # yp
            if varname == 'yp':
                return passParticleData(warp.getyp(gather=0))
            # zp
            if varname == 'zp':
                return passParticleData(warp.getzp(gather=0))
            # tp
            if varname == 'tp':
                return passParticleData(warp.gettp(gather=0))
            # rank
            if varname == 'rank':
                n = warp.getx(gather=0).size # FIXME -- must be defined somehwere?
                if not n:
                    return simV2.VISIT_INVALID_HANDLE
                rank = [float(parallel.get_rank())] * n # FIXME -- list ?
                vd = simV2.VisIt_VariableData_alloc()
                if not valid(vd):
                    pError('VisIt_VariableData_alloc failed')
                    return None
                simV2.VisIt_VariableData_setDataD(vd, simV2.VISIT_OWNER_SIM, 1, n, rank)
                return vd


        # grided data
        elif meshname == 'grid':
            # ax
            if varname == 'ax':
                return passGridData(warp.geta(comp='x',local=1))
            # ay
            if varname == 'ay':
                return passGridData(warp.geta(comp='y',local=1))
            # az
            if varname == 'az':
                return passGridData(warp.geta(comp='z',local=1))
            # A
            if varname == 'A':
                return passGridData(warp.geta(comp='A',local=1))
            # b
            if varname == 'B':
                return passGridData(warp.getb(comp='B',local=1))
            # bx
            if varname == 'bx':
                return passGridData(warp.getb(comp='x',local=1))
            # by
            if varname == 'by':
                return passGridData(warp.getb(comp='y',local=1))
            # bz
            if varname == 'bz':
                return passGridData(warp.getb(comp='z',local=1))
            # j
            if varname == 'J':
                return passGridData(warp.getj(comp='J',local=1))
            # jx
            if varname == 'jx':
                return passGridData(warp.getj(comp='x',local=1))
            # jy
            if varname == 'jy':
                return passGridData(warp.getj(comp='y',local=1))
            # jz
            if varname == 'jz':
                return passGridData(warp.getj(comp='z',local=1))
            # phi
            if varname == 'phi':
                return passGridData(warp.getphi(local=1))
            # rho
            if varname == 'rho':
                return passGridData(warp.getrho(local=1))
            # rank
            if varname == 'rank':
                rho = warp.getrho(local=1) # FIXME -- hack to verify
                rhoCopy = np.reshape(rho, rho.size, order='F').tolist() # FIXME
                rhoCopy = [float(parallel.get_rank())] * len(rhoCopy)
                vd = simV2.VisIt_VariableData_alloc()
                if not valid(vd):
                    pError('VisIt_VariableData_alloc failed')
                    return None
                simV2.VisIt_VariableData_setDataD(vd, simV2.VISIT_OWNER_COPY, 1, rho.size, rhoCopy)
                return vd

    # the simulation crashed when we asked for
    # that variable
    except:
        pError('Warp failed to produce %s'%(varname))
        return simV2.VISIT_INVALID_HANDLE

    # invalid
    pError('Unrecognized variable requested %s'%(varid))
    return simV2.VISIT_INVALID_HANDLE

#-----------------------------------------------------------------------------
def getDomains(name, userData):
    """
    Callback function used to get the list of domains handled by this process
    """
    # VisIt assumes that each process has some data
    # so you must have at least as many domains as
    # processes. you just return invalid handle in
    # getMesh/getVar on those processes that do not
    valid = lambda x : x != simV2.VISIT_INVALID_HANDLE
    pDebug('getDomains %s'%(name))

    rank = parallel.get_rank()
    numPE = parallel.number_of_PE()

    vd = simV2.VisIt_VariableData_alloc()
    if not valid(vd):
        pError('VisIt_VariableData_alloc failed')
        return None
    simV2.VisIt_VariableData_setDataI(vd, simV2.VISIT_OWNER_VISIT, 1, 1, [rank])

    doms = simV2.VisIt_DomainList_alloc()
    if not valid(doms):
        pError('VisIt_DomainList_alloc failed')
        return None
    simV2.VisIt_DomainList_setDomains(doms, numPE, vd)

    return doms

#-----------------------------------------------------------------------------
def passGridData(data):
    """Helper for passing data to VisIt"""
    pDebug(data.shape)
    valid = lambda x : x != simV2.VISIT_INVALID_HANDLE
    dataCopy = np.reshape(data, data.size, order='F').tolist() # FIXME -- copy*2
    vd = simV2.VisIt_VariableData_alloc()
    if not valid(vd):
        pError('VisIt_VariableData_alloc failed')
        return None
    simV2.VisIt_VariableData_setDataD(vd, simV2.VISIT_OWNER_COPY, 1, data.size, dataCopy)
    return vd

#-----------------------------------------------------------------------------
def passParticleData(data):
    """Helper for passing data to VisIt"""
    valid = lambda x : x != simV2.VISIT_INVALID_HANDLE
    if not data.size:
        return simV2.VISIT_INVALID_HANDLE
    dataCopy = data.tolist() # FIXME -- this makes a copy
    vd = simV2.VisIt_VariableData_alloc()
    if not valid(vd):
        pError('VisIt_VariableData_alloc failed')
        return None
    simV2.VisIt_VariableData_setDataD(vd, simV2.VISIT_OWNER_SIM, 1, data.size, dataCopy)
    return vd

#-----------------------------------------------------------------------------
def getDecomp():
    """Get the domain decomposition from Warp"""
    solver = (warp.getregisteredsolver() or warp.w3d)
    if solver == warp.w3d:
        decomp = warp.top.fsdecomp
    else:
        decomp = solver.fsdecomp
    return decomp

#-----------------------------------------------------------------------------
def getGridSize(cells=False):
    """
    Get this rank's grid size in units of cells
    (or units of points if cells is False).
    """
    decomp = getDecomp()

    i = 0 if cells else 1

    size = [[],[],[]]
    size[0] = decomp.nx[decomp.ixproc]+i
    size[1] = decomp.ny[decomp.iyproc]+i
    size[2] = decomp.nz[decomp.izproc]+i

    pDebug(size)
    return size

#-----------------------------------------------------------------------------
def getGridExtent(cells=False):
    """
    Get this rank's index space grid extent in units of cells
    (or units of points if cells is False).
    """
    decomp = getDecomp()

    i = 1 if cells else 0

    ext = [[],[],[]]
    ext[0] = [decomp.ix[decomp.ixproc], decomp.ix[decomp.ixproc]+decomp.nx[decomp.ixproc]-i]
    ext[1] = [decomp.iy[decomp.iyproc], decomp.iy[decomp.iyproc]+decomp.ny[decomp.iyproc]-i]
    ext[2] = [decomp.iz[decomp.izproc], decomp.iz[decomp.izproc]+decomp.nz[decomp.izproc]-i]

    return ext

#-----------------------------------------------------------------------------
def getGridIndices(cells=False):
    """
    Get this rank's grid indices in units of cells
    (or units of points if cells is False).
    """
    decomp = getDecomp()

    i = 0 if cells else 1

    ids = []
    ids.append(np.arange(decomp.ix[decomp.ixproc], decomp.ix[decomp.ixproc]+decomp.nx[decomp.ixproc]+i, dtype=np.float32))
    ids.append(np.arange(decomp.iy[decomp.iyproc], decomp.iy[decomp.iyproc]+decomp.ny[decomp.iyproc]+i, dtype=np.float32))
    ids.append(np.arange(decomp.iz[decomp.izproc], decomp.iz[decomp.izproc]+decomp.nz[decomp.izproc]+i, dtype=np.float32))

    return ids

#-----------------------------------------------------------------------------
def getGridSpacing():
    """Get the grid cell spacing"""
    dx = [[],[],[]]
    dx[0] = warp.w3d.dx
    dx[1] = warp.w3d.dy
    dx[2] = warp.w3d.dz
    return dx

#-----------------------------------------------------------------------------
def getGridOrigin():
    """Get the lower left corner of the grid"""
    x0 = [[],[],[]]
    x0[0] = warp.w3d.xmmin
    x0[1] = warp.w3d.ymmin
    x0[2] = warp.w3d.zmmin
    return x0

#-----------------------------------------------------------------------------
def getGridCoordinates(cells=False):
    """
    Get the coordinates as if the grid were a stretched mesh.
    This is a fundemental flaw in VisIt - it can't handle uniform
    meshes.
    """
    x0 = getGridOrigin()
    dx = getGridSpacing()
    coords = getGridIndices()
    i = 0
    while i<3:
        coords[i] *= dx[i]
        coords[i] += x0[i]
        i += 1
    return coords
