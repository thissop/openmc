# ParaStellDemo.py
# ---------------------------------------------------------------------------------------------------
# Autodesk Fusion 360 SCRIPT: step-by-step ParaStell teaching model.
#   * First 3 field periods (0-270 deg): each layer lofted into a colored solid; sketches hidden.
#   * 4th field period (270-360 deg): cross-section sketches (non-uniform breeder<->shield trade)
#     left VISIBLE and un-lofted so you can sweep/loft them live during the talk.
# (Moreno, Bader & Wilson, Front. Nucl. Eng. 2024).
#
# What it builds (one click), all inside a parent component "ParaStell_Geometry_Demo":
#   A. Reference_Curves      -- toroidal magnetic-axis-like centerline + fixed-theta guide rails.
#   B. Simple_Swept_Plasma_Demo  -- ONE bean profile swept along the centerline (constant section).
#   C. Plasma_Loft          -- many CHANGING poloidal sections lofted (a genuinely 3-D plasma).
#   D/E. Uniform_Radial_Build and Nonuniform_Radial_Build -- nested lofted layers (plasma, SOL,
#        first wall, breeder, shield, vacuum vessel); the non-uniform build locally trades breeder
#        thickness for shield thickness at a fixed total (fixed outer envelope).
#   F. Coils                -- 8 swept modular-coil rings (2 per field period for nfp=4).
#
# TEACHING MAPPING (Fusion feature  <->  the paper's step):
#   * changing poloidal profile sketches   <->  VMEC-like equilibrium: profile shape depends on phi.
#   * offset points along the local normal  <->  ParaStell radial build:  o_l = sum_i t_i, moved along n_hat.
#   * LOFT through the profiles             <->  ParaStell's CadQuery loft across poloidal profiles.
#   * SWEEP one profile along the axis path <->  a coil-style sweep / the "single profile" contrast.
#
# ROBUSTNESS NOTES (why the geometry is built the way it is):
#   * Fusion internal length unit is CENTIMETERS. The dialog takes mm/deg and stores internal cm/rad,
#     so all math below is already in cm/rad (input.value returns cm for lengths, rad for angles).
#   * Poloidal planes are meridional planes (constant phi, containing Z), made with the documented
#     ConstructionPlaneInput.setByAngle(zAxis, phi, xZplane). This keeps plane creation reliable while
#     still letting the cross-section SHAPE change with phi (the whole point of the loft-vs-sweep demo).
#   * Profile points are computed in world coordinates, then converted with sketch.modelToSketchSpace()
#     before being added as a closed SketchFittedSpline -> a clean planar profile for lofting.
#   * Each major stage is wrapped in try/except so a failure in one stage still leaves the others.
# ---------------------------------------------------------------------------------------------------

import adsk.core, adsk.fusion, traceback, math

_app = None
_ui = None
_handlers = []          # keep handler objects alive (classic Fusion gotcha)
CMD_ID = 'ParaStellStepByStepCmd'
PANEL_ID = 'SolidScriptsAddinsPanel'     # the "ADD-INS" panel in the SOLID tab


# ============================== small vector helpers (plain tuples, cm) =============================
def _add(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def _mul(a, s): return (a[0]*s, a[1]*s, a[2]*s)
def _P(t): return adsk.core.Point3D.create(t[0], t[1], t[2])

def _ang_diff(a, b):
    """Smallest signed angular difference a-b, wrapped into [-pi, pi]."""
    d = (a - b + math.pi) % (2.0*math.pi) - math.pi
    return d


# ============================== the analytic stellarator-like plasma ================================
# axis (magnetic-axis-like) curve: a circle of radius R0 with a small vertical excursion.
def axis_point(phi, R0, A_axis, Nfp):
    return (R0*math.cos(phi), R0*math.sin(phi), A_axis*math.sin(Nfp*phi))

def e_radial(phi):
    """Unit "horizontal" direction of the meridional plane at phi (points outward in major radius)."""
    return (math.cos(phi), math.sin(phi), 0.0)

E_VERT = (0.0, 0.0, 1.0)     # "vertical" direction of every meridional plane

def r_theta(theta, phi, a0, e1, e2, Nfp):
    """Minor radius vs poloidal angle theta, deliberately phi-dependent so the SHAPE changes toroidally
    (this is what a single sweep cannot reproduce)."""
    return a0 * (1.0 + e1*math.cos(Nfp*phi)*math.cos(theta)
                     + e2*math.sin(Nfp*phi)*math.sin(2.0*theta))

def local_uv(theta, phi, p):
    """Local (u,v) coordinates of a plasma-boundary point in the meridional plane at phi.
    u is along e_radial, v is along E_VERT. The whole cross-section is rotated by alpha(phi)."""
    r = r_theta(theta, phi, p['a0'], p['e1'], p['e2'], p['Nfp'])
    u = r*math.cos(theta)
    v = p['elong']*r*math.sin(theta)
    a = p['alpha0']*math.sin(p['Nfp']*phi)         # cross-section rotation vs phi
    ca, sa = math.cos(a), math.sin(a)
    return (u*ca - v*sa, u*sa + v*ca)

def uv_to_world(u, v, phi, axisPt):
    """Map local (u,v) in the meridional plane at phi to a 3-D world point (cm)."""
    er = e_radial(phi)
    return (axisPt[0] + u*er[0], axisPt[1] + u*er[1], axisPt[2] + v)   # v is along +Z

def outward_normals(uv):
    """Unit outward normals of a closed 2-D profile given as an ordered list of (u,v) points.
    Numerical tangent  t_j ~ uv[j+1]-uv[j-1];  normal = rotate t by 90 deg; flip to point away from
    the (approximately centered) origin. This is the 2-D version of ParaStell's first-wall normal."""
    n = len(uv)
    out = []
    for j in range(n):
        jm, jp = (j-1) % n, (j+1) % n
        tu = uv[jp][0] - uv[jm][0]
        tv = uv[jp][1] - uv[jm][1]
        nu, nv = tv, -tu                                   # rotate tangent by -90 deg
        if nu*uv[j][0] + nv*uv[j][1] < 0.0:                # ensure it points outward
            nu, nv = -nu, -nv
        L = math.hypot(nu, nv) or 1.0
        out.append((nu/L, nv/L))
    return out


# ============================== radial-build thicknesses (the trade) ================================
def local_delta(phi, theta, p):
    """Smooth localized breeder->shield trade amount delta(phi,theta) >= 0 (a 2-D Gaussian window)."""
    dphi = _ang_diff(phi, p['phi0'])
    dth = _ang_diff(theta, p['theta0'])
    g = math.exp(-0.5*((dphi/p['sig_phi'])**2 + (dth/p['sig_theta'])**2))
    return p['delta_max']*g

def layer_thicknesses(phi, theta, p, nonuniform):
    """Return (t_sol, t_fw, t_breeder, t_shield, t_vv) in cm.
    Non-uniform case: t_breeder = nom - delta, t_shield = nom + delta, so breeder+shield is CONSTANT
    -> the outer envelope is identical to the uniform build (the key teaching point)."""
    t_sol, t_fw, t_vv = p['t_sol'], p['t_fw'], p['t_vv']
    if nonuniform:
        d = local_delta(phi, theta, p)
        d = min(d, p['t_breeder'])                          # clamp so breeder never goes negative
        return (t_sol, t_fw, p['t_breeder'] - d, p['t_shield'] + d, t_vv)
    return (t_sol, t_fw, p['t_breeder'], p['t_shield'], t_vv)


# ============================== Fusion geometry helpers =============================================
def meridional_plane(comp, phi):
    """Construction plane = XZ plane rotated about the Z axis by phi (a constant-phi meridional plane).
    Uses the documented setByAngle(linearEntity=Z axis, angle, planarEntity=XZ plane)."""
    planes = comp.constructionPlanes
    pin = planes.createInput()
    pin.setByAngle(comp.zConstructionAxis,
                   adsk.core.ValueInput.createByReal(phi),
                   comp.xZConstructionPlane)
    return planes.add(pin)

def closed_profile_on_plane(comp, plane, world_pts):
    """Create a sketch on `plane`, add a CLOSED fitted spline through the given world points, and return
    (profile, sketch). world_pts are 3-D tuples (cm) that already lie in the plane."""
    sk = comp.sketches.add(plane)
    sk.isComputeDeferred = True                            # faster + avoids intermediate solves
    coll = adsk.core.ObjectCollection.create()
    for w in world_pts:
        coll.add(sk.modelToSketchSpace(_P(w)))             # world -> sketch-local (in-plane, z~0)
    spline = sk.sketchCurves.sketchFittedSplines.add(coll)
    try:
        spline.isClosed = True
    except Exception:
        pass
    sk.isComputeDeferred = False
    prof = sk.profiles.item(0) if sk.profiles.count > 0 else None
    return prof, sk

def spline_3d(comp, world_pts, closed):
    """Create a genuinely 3-D fitted spline (on an XY sketch) through non-coplanar world points."""
    sk = comp.sketches.add(comp.xYConstructionPlane)
    sk.isComputeDeferred = True
    coll = adsk.core.ObjectCollection.create()
    for w in world_pts:
        coll.add(_P(w))                                    # 3-D points -> 3-D spline in the sketch
    spline = sk.sketchCurves.sketchFittedSplines.add(coll)
    try:
        spline.isClosed = closed
    except Exception:
        pass
    sk.isComputeDeferred = False
    return spline, sk

def loft_profiles(comp, profiles, closed, name):
    """LOFT (solid) through an ordered list of profiles. Set isClosed to wrap a full torus."""
    profiles = [p for p in profiles if p is not None]
    if len(profiles) < 2:
        raise RuntimeError('loft "%s" needs >=2 valid profiles (got %d)' % (name, len(profiles)))
    lofts = comp.features.loftFeatures
    li = lofts.createInput(adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    for prof in profiles:
        li.loftSections.add(prof)
    li.isSolid = True
    try:
        li.isClosed = closed
    except Exception:
        pass
    feat = lofts.add(li)
    body = feat.bodies.item(0)
    body.name = name
    return body

def try_cut(comp, target_body, tool_body):
    """Boolean-subtract tool_body from target_body, KEEPING tool_body (so it can be reused as the next
    target). Returns True on success. Non-fatal on failure."""
    try:
        combines = comp.features.combineFeatures
        tools = adsk.core.ObjectCollection.create()
        tools.add(tool_body)
        ci = combines.createInput(target_body, tools)
        ci.operation = adsk.fusion.FeatureOperations.CutFeatureOperation
        ci.isKeepToolBodies = True
        combines.add(ci)
        return True
    except Exception:
        return False

def try_color(design, body, rgb, opacity=255):
    """Best-effort appearance/color. Fully guarded -- geometry must never fail because of appearance."""
    try:
        lib = _app.materialLibraries.itemByName('Fusion 360 Appearance Library')
        base = None
        for nm in ('Paint - Enamel Glossy (Generic)', 'Plastic - Matte (Generic)', 'Steel - Satin'):
            base = lib.appearances.itemByName(nm)
            if base:
                break
        if base is None:
            return
        appr = design.appearances.addByCopy(base, body.name + '_appr')
        for prop in appr.appearanceProperties:
            if isinstance(prop, adsk.core.ColorProperty):
                prop.value = adsk.core.Color.create(rgb[0], rgb[1], rgb[2], opacity)
                break
        body.appearance = appr
    except Exception:
        pass


# ============================== stage builders =====================================================
LAYER_NAMES = ['SOL', 'First_Wall', 'Breeder', 'Shield', 'Vacuum_Vessel']
LAYER_COLORS = {'Plasma': (235, 70, 70), 'SOL': (190, 190, 190), 'First_Wall': (95, 95, 95),
                'Breeder': (70, 135, 205), 'Shield': (160, 55, 55), 'Vacuum_Vessel': (55, 55, 62)}
PHI_SEAM = 2.0*math.pi*3.0/4.0      # 270 deg: end of the 3 pre-built field periods (start of the demo)


def _station_profiles(comp, phi, p, prefix, visible):
    """Create the plasma + every layer closed-profile sketch at ONE toroidal station.
    Returns (plasma_profile, {layer: profile}, [sketches]). The construction plane is always hidden;
    the sketches' visibility is set by `visible` (hidden for the pre-built blanket, shown for the demo)."""
    deg = int(round(math.degrees(phi))) % 360
    axisPt = axis_point(phi, p['R0'], p['A_axis'], p['Nfp'])
    plane = meridional_plane(comp, phi)
    plane.name = '%s_Plane_%03d' % (prefix, deg)
    try:
        plane.isLightBulbOn = False
    except Exception:
        pass
    thetas = [2*math.pi*j/p['Ntheta'] for j in range(p['Ntheta'])]
    plasma_uv = [local_uv(th, phi, p) for th in thetas]
    normals = outward_normals(plasma_uv)
    sks = []
    pworld = [uv_to_world(u, v, phi, axisPt) for (u, v) in plasma_uv]
    pprof, sk = closed_profile_on_plane(comp, plane, pworld)
    sk.name = '%s_Plasma_%03d' % (prefix, deg)
    try:
        sk.isVisible = visible
    except Exception:
        pass
    sks.append(sk)
    lprofs = {}
    cum = [0.0]*p['Ntheta']
    for li, nm in enumerate(LAYER_NAMES):
        for j, th in enumerate(thetas):
            cum[j] += layer_thicknesses(phi, th, p, True)[li]     # non-uniform everywhere (localized)
        lw = [uv_to_world(plasma_uv[j][0] + cum[j]*normals[j][0],
                          plasma_uv[j][1] + cum[j]*normals[j][1], phi, axisPt)
              for j in range(p['Ntheta'])]
        prof2, sk2 = closed_profile_on_plane(comp, plane, lw)
        sk2.name = '%s_%s_%03d' % (prefix, nm, deg)
        try:
            sk2.isVisible = visible
        except Exception:
            pass
        sks.append(sk2)
        lprofs[nm] = prof2
    return pprof, lprofs, sks


def build_prebuilt_blanket(comp, p, design):
    """First 3 field periods (0..270 deg): loft each layer into a solid, color it per layer, and HIDE
    every construction sketch/plane. This is the finished blanket the audience sees from the start."""
    n = max(3, p['Nphi_built'])
    phis = [PHI_SEAM*i/(n-1) for i in range(n)]      # 0 .. 270 deg inclusive
    plasma_profiles = []
    layer_profiles = {nm: [] for nm in LAYER_NAMES}
    all_sketches = []
    for phi in phis:
        pprof, lprofs, sks = _station_profiles(comp, phi, p, 'Blanket', visible=False)
        plasma_profiles.append(pprof)
        for nm in LAYER_NAMES:
            layer_profiles[nm].append(lprofs[nm])
        all_sketches += sks
    # OPEN loft (0..270 deg, not closed back around the torus)
    solids = {}
    solids['Plasma'] = loft_profiles(comp, plasma_profiles, False, 'Blanket_Plasma')
    for nm in LAYER_NAMES:
        solids[nm] = loft_profiles(comp, layer_profiles[nm], False, 'Blanket_' + nm + '_solid')
    # shells via ordered cut cascade (outer -> inner; keep tool bodies). Non-fatal.
    order = ['Vacuum_Vessel', 'Shield', 'Breeder', 'First_Wall', 'SOL']
    inner = {'Vacuum_Vessel': 'Shield', 'Shield': 'Breeder', 'Breeder': 'First_Wall',
             'First_Wall': 'SOL', 'SOL': 'Plasma'}
    shelled = True
    for nm in order:
        if not try_cut(comp, solids[nm], solids[inner[nm]]):
            shelled = False
    solids['Plasma'].name = 'Blanket_Plasma'
    for nm in LAYER_NAMES:
        solids[nm].name = 'Blanket_' + nm
    for nm, body in solids.items():
        try_color(design, body, LAYER_COLORS.get(nm, (128, 128, 128)),
                  opacity=150 if not shelled else 255)
    for sk in all_sketches:                          # hide the construction sketches
        try:
            sk.isVisible = False
        except Exception:
            pass
    return shelled


def build_demo_sector(comp, p, design):
    """4th field period (270..360 deg): create the cross-section sketches (plasma + every layer) with
    the NON-UNIFORM breeder->shield trade, left VISIBLE and UN-lofted so they can be walked through
    live. Also drop a visible axis arc to use as the plasma-sweep path."""
    m = max(2, p['Ndemo'])
    phis = [PHI_SEAM + (2*math.pi - PHI_SEAM)*i/(m-1) for i in range(m)]    # 270 .. 360 deg inclusive
    for phi in phis:
        _station_profiles(comp, phi, p, 'Demo', visible=True)
    arc = [axis_point(PHI_SEAM + (2*math.pi - PHI_SEAM)*i/60.0, p['R0'], p['A_axis'], p['Nfp'])
           for i in range(61)]
    _, sk = spline_3d(comp, arc, False)
    sk.name = 'Demo_Axis_Sweep_Path'
    try:
        sk.isVisible = True
    except Exception:
        pass
    return m


# ============================== command (dialog + execute) =========================================
def _defaults():
    # entered in mm/deg; stored internally as cm/rad by Fusion.
    return dict(
        R0_mm=100.0, a0_mm=22.0, Nfp=4, Ntheta=48, Nphi_built=19, Ndemo=5,
        elong=1.35, e1=0.28, e2=0.16, alpha0_deg=18.0, A_axis_mm=8.0,
        t_sol_mm=4.0, t_fw_mm=3.0, t_breeder_mm=12.0, t_shield_mm=10.0, t_vv_mm=4.0,
        delta_max_mm=6.0, phi0_deg=315.0, theta0_deg=0.0, sig_phi_deg=15.0, sig_theta_deg=40.0,
    )


class CommandExecuteHandler(adsk.core.CommandEventHandler):
    def __init__(self):
        super().__init__()

    def notify(self, args):
        try:
            inputs = args.command.commandInputs
            design = adsk.fusion.Design.cast(_app.activeProduct)
            if not design:
                _ui.messageBox('Open a Fusion DESIGN document first, then run this script.')
                return
            design.designType = adsk.fusion.DesignTypes.ParametricDesignType
            gv = lambda i: inputs.itemById(i).value          # cm for length, rad for angle
            p = dict(
                R0=gv('R0'), a0=gv('a0'), A_axis=gv('A_axis'), Nfp=inputs.itemById('Nfp').value,
                elong=inputs.itemById('elong').value, e1=inputs.itemById('e1').value,
                e2=inputs.itemById('e2').value, alpha0=gv('alpha0'),
                Ntheta=inputs.itemById('Ntheta').value,
                Nphi_built=inputs.itemById('Nphi_built').value, Ndemo=inputs.itemById('Ndemo').value,
                t_sol=gv('t_sol'), t_fw=gv('t_fw'), t_breeder=gv('t_breeder'),
                t_shield=gv('t_shield'), t_vv=gv('t_vv'),
                delta_max=gv('delta_max'), phi0=gv('phi0'), theta0=gv('theta0'),
                sig_phi=gv('sig_phi'), sig_theta=gv('sig_theta'),
            )
            root = design.rootComponent
            try:
                root.name = 'ParaStell_StepByStep'
            except Exception:
                pass
            report = []
            try:
                shelled = build_prebuilt_blanket(root, p, design)
                report.append('Pre-built blanket (3 periods): OK%s' % ('' if shelled else ' (solids, no shell)'))
            except Exception as e:
                report.append('Pre-built blanket: FAILED (%s)' % e)
            try:
                m = build_demo_sector(root, p, design)
                report.append('Demo sector (4th period): %d visible stations' % m)
            except Exception as e:
                report.append('Demo sector: FAILED (%s)' % e)
            try:
                _app.activeViewport.fit()
            except Exception:
                pass
            msg = ('Step-by-step model generated.\n\n' + '\n'.join(report) +
                   '\n\nHidden: all Blanket_* construction sketches & planes.\n'
                   'Visible: Demo_* cross-sections + Demo_Axis_Sweep_Path.\n\n'
                   'Walk-through order:\n'
                   '  1) Sweep Demo_Plasma_270 along Demo_Axis_Sweep_Path  -> constant section.\n'
                   '  2) Loft the Demo_Plasma_* sketches  -> the shape CHANGES with phi.\n'
                   '  3) Loft each Demo_<layer>_* set  -> watch breeder thin & shield thicken.')
            _ui.messageBox(msg)
        except:
            if _ui:
                _ui.messageBox('Execute failed:\n{}'.format(traceback.format_exc()))


class CommandCreatedHandler(adsk.core.CommandCreatedEventHandler):
    def __init__(self):
        super().__init__()

    def notify(self, args):
        try:
            cmd = args.command
            cmd.setDialogInitialSize(420, 640)
            inp = cmd.commandInputs
            d = _defaults()

            def valg(children, i, name, mm, unit='mm'):
                children.addValueInput(i, name, unit, adsk.core.ValueInput.createByString('%g %s' % (mm, unit)))

            g = inp.addGroupCommandInput('g_shape', 'Plasma shape'); gs = g.children
            valg(gs, 'R0', 'Major radius R0', d['R0_mm'])
            valg(gs, 'a0', 'Plasma minor radius a0', d['a0_mm'])
            valg(gs, 'A_axis', 'Axis vertical excursion', d['A_axis_mm'])
            gs.addIntegerSpinnerCommandInput('Nfp', 'Field periods Nfp', 1, 8, 1, d['Nfp'])
            gs.addValueInput('elong', 'Elongation', '', adsk.core.ValueInput.createByReal(d['elong']))
            gs.addValueInput('e1', 'Shaping e1', '', adsk.core.ValueInput.createByReal(d['e1']))
            gs.addValueInput('e2', 'Shaping e2', '', adsk.core.ValueInput.createByReal(d['e2']))
            valg(gs, 'alpha0', 'Section rotation amp.', d['alpha0_deg'], 'deg')

            g2 = inp.addGroupCommandInput('g_res', 'Resolution'); r = g2.children
            r.addIntegerSpinnerCommandInput('Ntheta', 'Poloidal points N_theta', 12, 128, 1, d['Ntheta'])
            r.addIntegerSpinnerCommandInput('Nphi_built', 'Built stations (0-270 deg)', 4, 96, 1, d['Nphi_built'])
            r.addIntegerSpinnerCommandInput('Ndemo', 'Demo stations (270-360 deg)', 2, 24, 1, d['Ndemo'])

            g3 = inp.addGroupCommandInput('g_build', 'Radial build (constant offsets)'); b = g3.children
            valg(b, 't_sol', 'SOL / vacuum', d['t_sol_mm'])
            valg(b, 't_fw', 'First wall', d['t_fw_mm'])
            valg(b, 't_breeder', 'Breeder (nominal)', d['t_breeder_mm'])
            valg(b, 't_shield', 'Shield (nominal)', d['t_shield_mm'])
            valg(b, 't_vv', 'Vacuum vessel', d['t_vv_mm'])

            g4 = inp.addGroupCommandInput('g_trade', 'Breeder->shield trade (in the demo sector)'); t = g4.children
            valg(t, 'delta_max', 'Max trade delta', d['delta_max_mm'])
            valg(t, 'phi0', 'Region center phi', d['phi0_deg'], 'deg')
            valg(t, 'theta0', 'Region center theta', d['theta0_deg'], 'deg')
            valg(t, 'sig_phi', 'Region width sigma_phi', d['sig_phi_deg'], 'deg')
            valg(t, 'sig_theta', 'Region width sigma_theta', d['sig_theta_deg'], 'deg')

            onExec = CommandExecuteHandler()
            cmd.execute.add(onExec)
            _handlers.append(onExec)
            onDestroy = CommandDestroyedHandler()
            cmd.destroy.add(onDestroy)
            _handlers.append(onDestroy)
        except:
            if _ui:
                _ui.messageBox('Command setup failed:\n{}'.format(traceback.format_exc()))


# ============================== script entry point =================================================
class CommandDestroyedHandler(adsk.core.CommandEventHandler):
    def __init__(self):
        super().__init__()

    def notify(self, args):
        adsk.terminate()        # dialog closed (OK or Cancel) -> end the script


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
            CMD_ID, 'ParaStell Step-By-Step',
            'Pre-built 3-period blanket (colored, sketches hidden) + a live-build 4th-period demo sector.')
        onCreated = CommandCreatedHandler()
        cmdDef.commandCreated.add(onCreated)
        _handlers.append(onCreated)
        cmdDef.execute()
        adsk.autoTerminate(False)
    except:
        if _ui:
            _ui.messageBox('Script start failed:\n{}'.format(traceback.format_exc()))
