# MagnetSweepDemo.py  (REAL QH geometry)
# ---------------------------------------------------------------------------------------------------
# Fusion 360 SCRIPT: distributed radial-build + magnet-sweep cutaway using the REAL QH (Landreman-Paul)
# plasma boundary (wout_qh.nc) and the REAL 20 modular coils (coils_qh), embedded in qh_real_geom.py.
#   * plasma lofted all the way around (real QH boundary),
#   * first wall on field periods 1 & 3 (opposite), blanket (breeder/shield/vessel) on period 1,
#   * finished REAL modular coils on periods 2 & 4 (circular conductor swept along the real filament),
#   * one REAL coil on period 1 left as centroid + guide rail + cross-section for a live Create->Sweep.
# Everything is scaled so the plasma major radius = target_R0 (dialog). Fusion internal unit = cm.
# ---------------------------------------------------------------------------------------------------
import adsk.core, adsk.fusion, traceback, math, os, sys

_here = os.path.dirname(os.path.abspath(__file__))
if _here not in sys.path:
    sys.path.insert(0, _here)
import qh_real_geom as G          # NFP, R0, XM, XN, RMNC, ZMNS, COILS (metres)

_app = None
_ui = None
_handlers = []
CMD_ID = 'ParaStellMagnetSweepRealCmd'
S = 1.0            # metres->cm scale, set on execute
COIL_R = 0.5       # conductor radius (cm), set on execute
COIL_COLOR = (205, 140, 45)
LAYER_COLORS = {'Plasma': (235, 70, 70), 'First_Wall': (110, 110, 115), 'Breeder': (70, 135, 205),
                'Shield': (160, 55, 55), 'Vacuum_Vessel': (60, 60, 68)}


# ---------- generic helpers ----------
def _P(t):
    return adsk.core.Point3D.create(t[0], t[1], t[2])

def outward_normals(uv):
    n = len(uv); out = []
    for j in range(n):
        jm, jp = (j-1) % n, (j+1) % n
        tu = uv[jp][0]-uv[jm][0]; tv = uv[jp][1]-uv[jm][1]
        nu, nv = tv, -tu
        if nu*uv[j][0] + nv*uv[j][1] < 0.0:
            nu, nv = -nu, -nv
        L = math.hypot(nu, nv) or 1.0
        out.append((nu/L, nv/L))
    return out

def meridional_plane(comp, phi):
    planes = comp.constructionPlanes
    pin = planes.createInput()
    pin.setByAngle(comp.zConstructionAxis, adsk.core.ValueInput.createByReal(phi), comp.xZConstructionPlane)
    pl = planes.add(pin)
    try:
        pl.isLightBulbOn = False
    except Exception:
        pass
    return pl

def closed_profile_on_plane(comp, plane, world_pts, name, visible):
    sk = comp.sketches.add(plane)
    sk.isComputeDeferred = True
    coll = adsk.core.ObjectCollection.create()
    for w in world_pts:
        coll.add(sk.modelToSketchSpace(_P(w)))
    spline = sk.sketchCurves.sketchFittedSplines.add(coll)
    try:
        spline.isClosed = True
    except Exception:
        pass
    sk.isComputeDeferred = False
    sk.name = name
    try:
        sk.isVisible = visible
    except Exception:
        pass
    return (sk.profiles.item(0) if sk.profiles.count > 0 else None), sk

def spline_3d(comp, world_pts, closed, name, visible):
    sk = comp.sketches.add(comp.xYConstructionPlane)
    sk.isComputeDeferred = True
    coll = adsk.core.ObjectCollection.create()
    for w in world_pts:
        coll.add(_P(w))
    spline = sk.sketchCurves.sketchFittedSplines.add(coll)
    try:
        spline.isClosed = closed
    except Exception:
        pass
    sk.isComputeDeferred = False
    sk.name = name
    try:
        sk.isVisible = visible
    except Exception:
        pass
    return spline, sk

def loft_profiles(comp, profiles, closed, name):
    profiles = [p for p in profiles if p is not None]
    if len(profiles) < 2:
        raise RuntimeError('loft "%s": <2 profiles' % name)
    lofts = comp.features.loftFeatures
    li = lofts.createInput(adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    for pr in profiles:
        li.loftSections.add(pr)
    li.isSolid = True
    try:
        li.isClosed = closed
    except Exception:
        pass
    body = lofts.add(li).bodies.item(0)
    body.name = name
    return body

def _cut(comp, target, tool, keep):
    try:
        combines = comp.features.combineFeatures
        tools = adsk.core.ObjectCollection.create(); tools.add(tool)
        ci = combines.createInput(target, tools)
        ci.operation = adsk.fusion.FeatureOperations.CutFeatureOperation
        ci.isKeepToolBodies = keep
        combines.add(ci); return True
    except Exception:
        return False

def try_color(design, body, rgb):
    try:
        lib = _app.materialLibraries.itemByName('Fusion 360 Appearance Library')
        base = None
        for nm in ('Paint - Enamel Glossy (Generic)', 'Plastic - Matte (Generic)', 'Steel - Satin'):
            base = lib.appearances.itemByName(nm)
            if base:
                break
        if not base:
            return
        appr = design.appearances.addByCopy(base, body.name + '_appr')
        for prop in appr.appearanceProperties:
            if isinstance(prop, adsk.core.ColorProperty):
                prop.value = adsk.core.Color.create(rgb[0], rgb[1], rgb[2], 255); break
        body.appearance = appr
    except Exception:
        pass


# ---------- REAL QH geometry ----------
def bnd_RZ(theta, phi):
    """Real plasma boundary (metres) from the embedded VMEC Fourier coefficients."""
    R = 0.0; Z = 0.0
    for m, n, rc, zs in zip(G.XM, G.XN, G.RMNC, G.ZMNS):
        ang = m*theta - n*phi
        R += rc*math.cos(ang); Z += zs*math.sin(ang)
    return R, Z

def real_station(phi, ntheta):
    """Boundary cross-section at phi: local (u,v) in cm about the section centroid, plus centroid world."""
    thetas = [2*math.pi*j/ntheta for j in range(ntheta)]
    RZ = [bnd_RZ(th, phi) for th in thetas]
    Rc = sum(r for r, _ in RZ)/ntheta
    Zc = sum(z for _, z in RZ)/ntheta
    uv = [((R-Rc)*S, (Z-Zc)*S) for (R, Z) in RZ]
    center = (Rc*S*math.cos(phi), Rc*S*math.sin(phi), Zc*S)
    return uv, outward_normals(uv), center

def offset_world(uv, normals, center, phi, off):
    cphi, sphi = math.cos(phi), math.sin(phi)
    out = []
    for j in range(len(uv)):
        u = uv[j][0] + off*normals[j][0]
        v = uv[j][1] + off*normals[j][1]
        out.append((center[0] + u*cphi, center[1] + u*sphi, center[2] + v))
    return out

def offset_profile(comp, phi, off, ntheta, name, visible):
    uv, normals, center = real_station(phi, ntheta)
    plane = meridional_plane(comp, phi)
    return closed_profile_on_plane(comp, plane, offset_world(uv, normals, center, phi, off), name, visible)


# ---------- stage builders ----------
def cum_offsets(p):
    o = p['t_sol']; o1 = o + p['t_fw']; o2 = o1 + p['t_breeder']; o3 = o2 + p['t_shield']; o4 = o3 + p['t_vv']
    return dict(sol=o, fw=o1, breeder=o2, shield=o3, vv=o4)

def build_plasma_full(comp, p, design):
    n = p['Nphi_period']*G.NFP
    phis = [2*math.pi*i/n for i in range(n)]
    profs, sks = [], []
    for phi in phis:
        pr, sk = offset_profile(comp, phi, 0.0, p['Ntheta'], 'Plasma_%03d' % (int(round(math.degrees(phi))) % 360), False)
        profs.append(pr); sks.append(sk)
    try:
        body = loft_profiles(comp, profs, True, 'Plasma')
    except Exception:
        body = loft_profiles(comp, profs, False, 'Plasma')
    try_color(design, body, LAYER_COLORS['Plasma'])
    for sk in sks:
        try:
            sk.isVisible = False
        except Exception:
            pass
    return body

def build_shell_sector(comp, p, design, outer_off, inner_off, phi0, phi1, name, color):
    nst = p['Nphi_period'] + 1
    phis = [phi0 + (phi1-phi0)*i/(nst-1) for i in range(nst)]
    outer, inner, sks = [], [], []
    for phi in phis:
        deg = int(round(math.degrees(phi))) % 360
        po, so = offset_profile(comp, phi, outer_off, p['Ntheta'], '%s_out_%03d' % (name, deg), False)
        pi, si = offset_profile(comp, phi, inner_off, p['Ntheta'], '%s_in_%03d' % (name, deg), False)
        outer.append(po); inner.append(pi); sks += [so, si]
    ob = loft_profiles(comp, outer, False, name + '_outer_solid')
    ib = loft_profiles(comp, inner, False, name + '_inner_solid')
    _cut(comp, ob, ib, keep=False)
    ob.name = name
    try_color(design, ob, color)
    for sk in sks:
        try:
            sk.isVisible = False
        except Exception:
            pass
    return ob

def _coil_centroid3d(pts):
    n = len(pts)
    return (sum(p[0] for p in pts)/n, sum(p[1] for p in pts)/n, sum(p[2] for p in pts)/n)

def build_real_coil(comp, design, pts_m, name):
    """Sweep a circular conductor along a REAL coil filament (non-planar)."""
    pts = [(x*S, y*S, z*S) for (x, y, z) in pts_m]
    cspline, csk = spline_3d(comp, pts, True, name + '_centroid', False)
    planes = comp.constructionPlanes
    pin = planes.createInput()
    pin.setByDistanceOnPath(cspline, adsk.core.ValueInput.createByReal(0.0))   # plane _|_ filament at start
    secplane = planes.add(pin)
    try:
        secplane.isLightBulbOn = False
    except Exception:
        pass
    sk = comp.sketches.add(secplane)
    sk.sketchCurves.sketchCircles.addByCenterRadius(sk.modelToSketchSpace(_P(pts[0])), COIL_R)
    sk.name = name + '_section'
    prof = sk.profiles.item(0)
    path = comp.features.createPath(cspline, True)
    sweeps = comp.features.sweepFeatures
    si = sweeps.createInput(prof, path, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    body = sweeps.add(si).bodies.item(0)
    body.name = name
    try_color(design, body, COIL_COLOR)
    try:
        sk.isVisible = False
    except Exception:
        pass
    return body

def build_demo_coil(comp, design, pts_m):
    """A REAL coil in period 1 left un-swept: visible centroid + guide rail + cross-section for live sweep."""
    pts = [(x*S, y*S, z*S) for (x, y, z) in pts_m]
    _, csk = spline_3d(comp, pts, True, 'DemoCoil_Centroid', True)
    cen = _coil_centroid3d(pts)
    gl = []
    for pnt in pts:
        d = (pnt[0]-cen[0], pnt[1]-cen[1], pnt[2]-cen[2])
        L = math.sqrt(d[0]**2 + d[1]**2 + d[2]**2) or 1.0
        gl.append((pnt[0] + COIL_R*d[0]/L, pnt[1] + COIL_R*d[1]/L, pnt[2] + COIL_R*d[2]/L))
    gspline, gsk = spline_3d(comp, gl, True, 'DemoCoil_GuideRail', True)
    # rectangular cross-section on the plane perpendicular to the filament at its start
    planes = comp.constructionPlanes
    pin = planes.createInput()
    pin.setByDistanceOnPath(csk.sketchCurves.sketchFittedSplines.item(0),
                            adsk.core.ValueInput.createByReal(0.0))
    secplane = planes.add(pin)
    try:
        secplane.isLightBulbOn = False
    except Exception:
        pass
    sk = comp.sketches.add(secplane)
    c2 = sk.modelToSketchSpace(_P(pts[0]))
    corner = adsk.core.Point3D.create(c2.x + COIL_R, c2.y + 1.5*COIL_R, 0.0)
    sk.sketchCurves.sketchLines.addCenterPointRectangle(c2, corner)
    sk.name = 'DemoCoil_Section'
    sk.isVisible = True
    return True

def coil_phi_deg(pts_m):
    n = len(pts_m)
    x = sum(p[0] for p in pts_m)/n; y = sum(p[1] for p in pts_m)/n
    return math.degrees(math.atan2(y, x)) % 360


# ---------- command ----------
def _defaults():
    return dict(target_R0_mm=600.0, Ntheta=64, Nphi_period=8,
                t_sol_mm=8.0, t_fw_mm=5.0, t_breeder_mm=22.0, t_shield_mm=16.0, t_vv_mm=8.0,
                coil_r_mm=18.0)


class CommandExecuteHandler(adsk.core.CommandEventHandler):
    def __init__(self):
        super().__init__()

    def notify(self, args):
        global S, COIL_R
        try:
            inp = args.command.commandInputs
            design = adsk.fusion.Design.cast(_app.activeProduct)
            if not design:
                _ui.messageBox('Open a Fusion DESIGN document first, then run this script.')
                return
            design.designType = adsk.fusion.DesignTypes.ParametricDesignType
            gv = lambda i: inp.itemById(i).value
            target_R0 = gv('target_R0')                    # cm
            S = target_R0 / G.R0                           # metres -> cm scale
            COIL_R = gv('coil_r')
            p = dict(Ntheta=inp.itemById('Ntheta').value, Nphi_period=inp.itemById('Nphi_period').value,
                     t_sol=gv('t_sol'), t_fw=gv('t_fw'), t_breeder=gv('t_breeder'),
                     t_shield=gv('t_shield'), t_vv=gv('t_vv'))
            root = design.rootComponent
            try:
                root.name = 'ParaStell_MagnetSweep_RealQH'
            except Exception:
                pass
            o = cum_offsets(p)
            per = 2*math.pi/G.NFP
            def rng(k):
                return ((k-1)*per, k*per)
            report = []

            try:
                build_plasma_full(root, p, design); report.append('Plasma (real QH, full torus): OK')
            except Exception as e:
                report.append('Plasma: FAILED (%s)' % e)
            for k in (1, 3):
                try:
                    a, b = rng(k)
                    build_shell_sector(root, p, design, o['fw'], o['sol'], a, b, 'FirstWall_P%d' % k, LAYER_COLORS['First_Wall'])
                    report.append('First wall period %d: OK' % k)
                except Exception as e:
                    report.append('First wall period %d: FAILED (%s)' % (k, e))
            a1, b1 = rng(1)
            for nm, oo, io, col in (('Breeder_P1', o['breeder'], o['fw'], LAYER_COLORS['Breeder']),
                                    ('Shield_P1', o['shield'], o['breeder'], LAYER_COLORS['Shield']),
                                    ('VacVessel_P1', o['vv'], o['shield'], LAYER_COLORS['Vacuum_Vessel'])):
                try:
                    build_shell_sector(root, p, design, oo, io, a1, b1, nm, col)
                    report.append('%s: OK' % nm)
                except Exception as e:
                    report.append('%s: FAILED (%s)' % (nm, e))

            # sort the 20 real coils into field periods; sweep those in periods 2 & 4; demo one in period 1
            per_deg = 360.0/G.NFP
            demo_idx, made = None, 0
            best = 1e9
            for i, pts in enumerate(G.COILS):
                phi = coil_phi_deg(pts)
                k = int(phi // per_deg) + 1
                if k in (2, 4):
                    try:
                        build_real_coil(root, design, pts, 'Coil_P%d_%02d' % (k, i)); made += 1
                    except Exception:
                        pass
                elif k == 1:
                    dd = abs(phi - 0.5*per_deg)
                    if dd < best:
                        best = dd; demo_idx = i
            report.append('Real coils swept on periods 2 & 4: %d' % made)
            if demo_idx is not None:
                try:
                    build_demo_coil(root, design, G.COILS[demo_idx])
                    report.append('Demo magnet (real coil, period 1): centroid + guide rail + section')
                except Exception as e:
                    report.append('Demo magnet: FAILED (%s)' % e)

            try:
                _app.activeViewport.fit()
            except Exception:
                pass
            _ui.messageBox('Real-QH magnet sweep demo generated.\n\n' + '\n'.join(report) +
                           '\n\nLive sweep (period 1): Create -> Sweep, profile = DemoCoil_Section,\n'
                           'path = DemoCoil_Centroid, optional Guide Rail = DemoCoil_GuideRail.')
        except:
            if _ui:
                _ui.messageBox('Execute failed:\n{}'.format(traceback.format_exc()))


class CommandCreatedHandler(adsk.core.CommandCreatedEventHandler):
    def __init__(self):
        super().__init__()

    def notify(self, args):
        try:
            cmd = args.command
            cmd.setDialogInitialSize(420, 520)
            inp = cmd.commandInputs
            d = _defaults()
            def valg(i, name, mm, unit='mm'):
                inp.addValueInput(i, name, unit, adsk.core.ValueInput.createByString('%g %s' % (mm, unit)))
            valg('target_R0', 'Plasma major radius R0', d['target_R0_mm'])
            inp.addIntegerSpinnerCommandInput('Ntheta', 'Poloidal points', 24, 160, 1, d['Ntheta'])
            inp.addIntegerSpinnerCommandInput('Nphi_period', 'Stations per field period', 4, 24, 1, d['Nphi_period'])
            valg('t_sol', 'SOL / vacuum', d['t_sol_mm'])
            valg('t_fw', 'First wall', d['t_fw_mm'])
            valg('t_breeder', 'Breeder', d['t_breeder_mm'])
            valg('t_shield', 'Shield', d['t_shield_mm'])
            valg('t_vv', 'Vacuum vessel', d['t_vv_mm'])
            valg('coil_r', 'Conductor radius', d['coil_r_mm'])
            onExec = CommandExecuteHandler(); cmd.execute.add(onExec); _handlers.append(onExec)
            onDestroy = CommandDestroyedHandler(); cmd.destroy.add(onDestroy); _handlers.append(onDestroy)
        except:
            if _ui:
                _ui.messageBox('Command setup failed:\n{}'.format(traceback.format_exc()))


class CommandDestroyedHandler(adsk.core.CommandEventHandler):
    def __init__(self):
        super().__init__()

    def notify(self, args):
        adsk.terminate()


def run(context):
    global _app, _ui
    try:
        _app = adsk.core.Application.get()
        _ui = _app.userInterface
        design = adsk.fusion.Design.cast(_app.activeProduct)
        if not design:
            _ui.messageBox('Open a Fusion DESIGN document first, then run this script.')
            return
        old = _ui.commandDefinitions.itemById(CMD_ID)
        if old:
            old.deleteMe()
        cmdDef = _ui.commandDefinitions.addButtonDefinition(
            CMD_ID, 'Magnet Sweep Demo (Real QH)',
            'Real QH plasma + real modular coils: cutaway with a live-sweep magnet on period 1.')
        onCreated = CommandCreatedHandler()
        cmdDef.commandCreated.add(onCreated)
        _handlers.append(onCreated)
        cmdDef.execute()
        adsk.autoTerminate(False)
    except:
        if _ui:
            _ui.messageBox('Script start failed:\n{}'.format(traceback.format_exc()))
