# ParaStellDemo.py
# ---------------------------------------------------------------------------------------------------
# Autodesk Fusion 360 add-in: an educational CAD demonstration of the geometry ideas behind ParaStell
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
CMD_ID = 'ParaStellGeometryDemoCmd'
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
def build_reference_curves(comp, p):
    """A. toroidal centerline + four fixed-theta guide rails (teaching curves; not required by the loft)."""
    Nphi_dense = 120
    axis_pts = [axis_point(2*math.pi*i/Nphi_dense, p['R0'], p['A_axis'], p['Nfp'])
                for i in range(Nphi_dense)]
    _, sk_axis = spline_3d(comp, axis_pts, True)
    sk_axis.name = 'Axis_Centerline'
    for theta in (0.0, math.pi/2, math.pi, 3*math.pi/2):
        rail = []
        for i in range(Nphi_dense):
            phi = 2*math.pi*i/Nphi_dense
            u, v = local_uv(theta, phi, p)
            rail.append(uv_to_world(u, v, phi, axis_point(phi, p['R0'], p['A_axis'], p['Nfp'])))
        _, sk = spline_3d(comp, rail, True)
        sk.name = 'GuideRail_theta_%03d' % int(round(math.degrees(theta)))
    return sk_axis

def build_sweep_demo(comp, p, axis_sketch):
    """B. ONE bean profile (at phi=0) swept along the centerline path -> constant cross-section body.
    This is intentionally the *wrong* method; it contrasts with the loft."""
    phi = 0.0
    axisPt = axis_point(phi, p['R0'], p['A_axis'], p['Nfp'])
    thetas = [2*math.pi*j/p['Ntheta'] for j in range(p['Ntheta'])]
    pts = [uv_to_world(*local_uv(th, phi, p), phi=phi, axisPt=axisPt) for th in thetas]
    plane = meridional_plane(comp, phi)
    plane.name = 'Plane_phi_000_sweepdemo'
    prof, sk = closed_profile_on_plane(comp, plane, pts)
    sk.name = 'SweepDemo_Profile_phi0'
    if prof is None:
        raise RuntimeError('sweep-demo profile is not closed')
    # path = the toroidal centerline spline
    path_curve = axis_sketch.sketchCurves.sketchFittedSplines.item(0)
    path = comp.features.createPath(path_curve, True)
    sweeps = comp.features.sweepFeatures
    si = sweeps.createInput(prof, path, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    try:
        si.profileScaling = adsk.fusion.SweepProfileScalingOptions.SweepProfileNoScalingOption
    except Exception:
        pass
    feat = sweeps.add(si)
    body = feat.bodies.item(0)
    body.name = 'Simple_Swept_Plasma_Demo'
    return body

def build_plasma_and_layers(parent_occ_comp, p, nonuniform, design):
    """C/D/E. Build the plasma loft and the nested layer lofts inside a dedicated sub-component."""
    tag = 'Nonuniform' if nonuniform else 'Uniform'
    comp = parent_occ_comp          # Part docs allow ONE component -> build in root, distinguish by name prefix

    Nphi, Ntheta = p['Nphi'], p['Ntheta']
    close = p['close_full_torus']
    phis = [2*math.pi*i/Nphi for i in range(Nphi)]
    layer_names = ['SOL', 'First_Wall', 'Breeder', 'Shield', 'Vacuum_Vessel']

    plasma_profiles = []
    layer_profiles = {nm: [] for nm in layer_names}

    for phi in phis:
        axisPt = axis_point(phi, p['R0'], p['A_axis'], p['Nfp'])
        plane = meridional_plane(comp, phi)
        plane.name = '%s_Plane_phi_%03d' % (tag, int(round(math.degrees(phi))))
        thetas = [2*math.pi*j/Ntheta for j in range(Ntheta)]
        plasma_uv = [local_uv(th, phi, p) for th in thetas]
        normals = outward_normals(plasma_uv)

        # plasma profile (offset 0)
        pworld = [uv_to_world(u, v, phi, axisPt) for (u, v) in plasma_uv]
        prof, sk = closed_profile_on_plane(comp, plane, pworld)
        sk.name = '%s_Plasma_phi_%03d' % (tag, int(round(math.degrees(phi))))
        plasma_profiles.append(prof)

        # cumulative outer boundary for each layer, offset along the local normal
        cum = [0.0]*Ntheta
        for li, nm in enumerate(layer_names):
            for j, th in enumerate(thetas):
                t = layer_thicknesses(phi, th, p, nonuniform)[li]
                cum[j] += t
            lw = []
            for j in range(Ntheta):
                u = plasma_uv[j][0] + cum[j]*normals[j][0]
                v = plasma_uv[j][1] + cum[j]*normals[j][1]
                lw.append(uv_to_world(u, v, phi, axisPt))
            prof2, sk2 = closed_profile_on_plane(comp, plane, lw)
            sk2.name = '%s_%s_phi_%03d' % (tag, nm, int(round(math.degrees(phi))))
            layer_profiles[nm].append(prof2)

    # loft each boundary into a solid
    solids = {}
    solids['Plasma'] = loft_profiles(comp, plasma_profiles, close, tag + '_Plasma_solid')
    for nm in layer_names:
        solids[nm] = loft_profiles(comp, layer_profiles[nm], close, tag + '_' + nm + '_solid')

    # shells via ordered Boolean cascade (outer -> inner; each solid is target exactly once)
    order = ['Vacuum_Vessel', 'Shield', 'Breeder', 'First_Wall', 'SOL']
    inner = {'Vacuum_Vessel': 'Shield', 'Shield': 'Breeder', 'Breeder': 'First_Wall',
             'First_Wall': 'SOL', 'SOL': 'Plasma'}
    shelled = True
    for nm in order:
        if not try_cut(comp, solids[nm], solids[inner[nm]]):
            shelled = False
    # tidy names + colors
    colors = {'Plasma': (255, 80, 80), 'SOL': (170, 170, 170), 'First_Wall': (90, 90, 90),
              'Breeder': (70, 130, 200), 'Shield': (150, 60, 60), 'Vacuum_Vessel': (60, 60, 60)}
    solids['Plasma'].name = tag + '_Plasma'
    for nm in layer_names:
        solids[nm].name = tag + '_' + nm
    for nm, body in solids.items():
        try_color(design, body, colors.get(nm, (128, 128, 128)),
                  opacity=140 if not shelled else 255)     # translucent if we couldn't shell
    return comp, shelled

def build_coils(comp, p, design):
    """F. 8 swept modular-coil rings (2 per field period, nfp=4). Each = a circular conductor swept
    along a planar 'picture-frame' coil centerline that encircles the plasma cross-section."""
    ncoils = 2*p['Nfp']
    # generous coil minor radius so the ring clears the whole radial build
    build = p['t_sol']+p['t_fw']+p['t_breeder']+p['t_shield']+p['t_vv']
    env = p['a0']*(1.0+abs(p['e1'])+abs(p['e2']))*p['elong'] + build
    coil_a = env + p['coil_gap']            # horizontal (major-radius) semi-size
    coil_b = env*1.05 + p['coil_gap']       # vertical semi-size
    made = 0
    for k in range(ncoils):
        try:
            phi = 2*math.pi*k/ncoils
            axisPt = axis_point(phi, p['R0'], p['A_axis'], p['Nfp'])
            er = e_radial(phi)
            # coil centerline: an ellipse in the meridional plane at phi, encircling the cross-section
            M = 80
            cl = []
            for m in range(M):
                t = 2*math.pi*m/M
                u = coil_a*math.cos(t)
                v = coil_b*math.sin(t)
                cl.append((axisPt[0]+u*er[0], axisPt[1]+u*er[1], axisPt[2]+v))
            cspline, csk = spline_3d(comp, cl, True)
            csk.name = 'Coil_%02d_centerline' % k
            # conductor cross-section: a small circle in the horizontal plane at the t=0 point,
            # where the coil tangent is vertical (so the circle is perpendicular to the path there).
            start = (axisPt[0]+coil_a*er[0], axisPt[1]+coil_a*er[1], axisPt[2])
            planes = comp.constructionPlanes
            pin = planes.createInput()
            pin.setByOffset(comp.xYConstructionPlane, adsk.core.ValueInput.createByReal(start[2]))
            hplane = planes.add(pin)
            hplane.name = 'Coil_%02d_section_plane' % k
            sk = comp.sketches.add(hplane)
            center2d = sk.modelToSketchSpace(_P(start))
            sk.sketchCurves.sketchCircles.addByCenterRadius(center2d, p['coil_r'])
            sk.name = 'Coil_%02d_section' % k
            prof = sk.profiles.item(0)
            path = comp.features.createPath(cspline, True)
            sweeps = comp.features.sweepFeatures
            si = sweeps.createInput(prof, path, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
            feat = sweeps.add(si)
            body = feat.bodies.item(0)
            body.name = 'Coil_%02d' % k
            try_color(design, body, (210, 160, 40))
            made += 1
        except Exception:
            pass                                            # skip a coil that fails; keep the rest
    return made


# ============================== the command (dialog + execute) =====================================
def _defaults():
    # dialog values are entered in mm/deg; stored internally as cm/rad by Fusion.
    return dict(
        R0_mm=100.0, a0_mm=22.0, Nfp=4, Nphi=12, Ntheta=48,
        elong=1.35, e1=0.28, e2=0.16, alpha0_deg=18.0, A_axis_mm=8.0,
        t_sol_mm=4.0, t_fw_mm=3.0, t_breeder_mm=12.0, t_shield_mm=10.0, t_vv_mm=4.0,
        delta_max_mm=6.0, phi0_deg=60.0, theta0_deg=0.0, sig_phi_deg=25.0, sig_theta_deg=40.0,
        coil_gap_mm=10.0, coil_r_mm=5.0,
    )


class CommandExecuteHandler(adsk.core.CommandEventHandler):
    def __init__(self):
        super().__init__()

    def notify(self, args):
        try:
            inputs = args.command.commandInputs
            design = adsk.fusion.Design.cast(_app.activeProduct)
            if not design:
                _ui.messageBox('No active Fusion design. Open a Design document and try again.')
                return
            design.designType = adsk.fusion.DesignTypes.ParametricDesignType

            gv = lambda i: inputs.itemById(i).value      # value: cm for length, rad for angle
            iv = lambda i: inputs.itemById(i).value      # int spinners
            bv = lambda i: inputs.itemById(i).value      # bool checkbox

            # assemble the parameter dict in internal units (cm, rad)
            p = dict(
                R0=gv('R0'), a0=gv('a0'), Nfp=iv('Nfp'), Nphi=iv('Nphi'), Ntheta=iv('Ntheta'),
                elong=inputs.itemById('elong').value, e1=inputs.itemById('e1').value,
                e2=inputs.itemById('e2').value, alpha0=gv('alpha0'), A_axis=gv('A_axis'),
                t_sol=gv('t_sol'), t_fw=gv('t_fw'), t_breeder=gv('t_breeder'),
                t_shield=gv('t_shield'), t_vv=gv('t_vv'),
                delta_max=gv('delta_max'), phi0=gv('phi0'), theta0=gv('theta0'),
                sig_phi=gv('sig_phi'), sig_theta=gv('sig_theta'),
                coil_gap=gv('coil_gap'), coil_r=gv('coil_r'),
                close_full_torus=bv('close_torus'),
            )
            do_sweep, do_uniform = bv('do_sweep'), bv('do_uniform')
            do_nonuniform, do_coils = bv('do_nonuniform'), bv('do_coils')

            root = design.rootComponent
            parent = root                       # Part docs allow only one component -> build in root
            try:
                root.name = 'ParaStell_Geometry_Demo'
            except Exception:
                pass

            report = []
            # A. reference curves
            try:
                axis_sketch = build_reference_curves(parent, p)
                report.append('Reference_Curves: OK')
            except Exception as e:
                axis_sketch = None
                report.append('Reference_Curves: FAILED (%s)' % e)
            # B. sweep demo
            if do_sweep and axis_sketch is not None:
                try:
                    build_sweep_demo(parent, p, axis_sketch)
                    report.append('Simple_Swept_Plasma_Demo: OK')
                except Exception as e:
                    report.append('Simple_Swept_Plasma_Demo: FAILED (%s)' % e)
            # C/D. uniform radial build (includes the plasma loft)
            if do_uniform:
                try:
                    _, shelled = build_plasma_and_layers(parent, p, False, design)
                    report.append('Uniform_Radial_Build: OK%s' % ('' if shelled else ' (solids, no shell)'))
                except Exception as e:
                    report.append('Uniform_Radial_Build: FAILED (%s)' % e)
            # E. non-uniform radial build (breeder-for-shield trade)
            if do_nonuniform:
                try:
                    _, shelled = build_plasma_and_layers(parent, p, True, design)
                    report.append('Nonuniform_Radial_Build: OK%s' % ('' if shelled else ' (solids, no shell)'))
                except Exception as e:
                    report.append('Nonuniform_Radial_Build: FAILED (%s)' % e)
            # F. coils
            if do_coils:
                try:
                    n = build_coils(parent, p, design)
                    report.append('Coils: %d of %d built' % (n, 2*p['Nfp']))
                except Exception as e:
                    report.append('Coils: FAILED (%s)' % e)

            try:
                _app.activeViewport.fit()
            except Exception:
                pass
            _ui.messageBox('ParaStell demo generated.\n\n' + '\n'.join(report))
        except:
            if _ui:
                _ui.messageBox('Execute failed:\n{}'.format(traceback.format_exc()))


class CommandCreatedHandler(adsk.core.CommandCreatedEventHandler):
    def __init__(self):
        super().__init__()

    def notify(self, args):
        try:
            cmd = args.command
            cmd.setDialogInitialSize(420, 720)
            inp = cmd.commandInputs
            d = _defaults()

            def val(i, name, mm, unit='mm'):
                inp.addValueInput(i, name, unit, adsk.core.ValueInput.createByString('%g %s' % (mm, unit)))
            def plain(i, name, v):
                inp.addValueInput(i, name, '', adsk.core.ValueInput.createByReal(v))   # unitless
            def spin(i, name, lo, hi, val_):
                inp.addIntegerSpinnerCommandInput(i, name, lo, hi, 1, val_)
            def chk(i, name, v):
                inp.addBoolValueInput(i, name, True, '', v)

            g = inp.addGroupCommandInput('g_shape', 'Plasma shape'); gs = g.children
            def valg(children, i, name, mm, unit='mm'):
                children.addValueInput(i, name, unit, adsk.core.ValueInput.createByString('%g %s' % (mm, unit)))
            valg(gs, 'R0', 'Major radius R0', d['R0_mm'])
            valg(gs, 'a0', 'Plasma minor radius a0', d['a0_mm'])
            valg(gs, 'A_axis', 'Axis vertical excursion', d['A_axis_mm'])
            gs.addIntegerSpinnerCommandInput('Nfp', 'Field periods Nfp', 1, 8, 1, d['Nfp'])
            gs.addValueInput('elong', 'Elongation', '', adsk.core.ValueInput.createByReal(d['elong']))
            gs.addValueInput('e1', 'Shaping e1', '', adsk.core.ValueInput.createByReal(d['e1']))
            gs.addValueInput('e2', 'Shaping e2', '', adsk.core.ValueInput.createByReal(d['e2']))
            valg(gs, 'alpha0', 'Section rotation amp.', d['alpha0_deg'], 'deg')

            g2 = inp.addGroupCommandInput('g_res', 'Resolution'); r = g2.children
            r.addIntegerSpinnerCommandInput('Nphi', 'Toroidal sections N_phi', 4, 96, 1, d['Nphi'])
            r.addIntegerSpinnerCommandInput('Ntheta', 'Poloidal points N_theta', 12, 128, 1, d['Ntheta'])
            r.addBoolValueInput('close_torus', 'Close full torus (else open)', True, '', True)

            g3 = inp.addGroupCommandInput('g_build', 'Radial build (constant offsets)'); b = g3.children
            valg(b, 't_sol', 'SOL / vacuum', d['t_sol_mm'])
            valg(b, 't_fw', 'First wall', d['t_fw_mm'])
            valg(b, 't_breeder', 'Breeder (nominal)', d['t_breeder_mm'])
            valg(b, 't_shield', 'Shield (nominal)', d['t_shield_mm'])
            valg(b, 't_vv', 'Vacuum vessel', d['t_vv_mm'])

            g4 = inp.addGroupCommandInput('g_trade', 'Breeder->shield trade (non-uniform)'); t = g4.children
            valg(t, 'delta_max', 'Max trade delta', d['delta_max_mm'])
            valg(t, 'phi0', 'Region center phi', d['phi0_deg'], 'deg')
            valg(t, 'theta0', 'Region center theta', d['theta0_deg'], 'deg')
            valg(t, 'sig_phi', 'Region width sigma_phi', d['sig_phi_deg'], 'deg')
            valg(t, 'sig_theta', 'Region width sigma_theta', d['sig_theta_deg'], 'deg')

            g5 = inp.addGroupCommandInput('g_coils', 'Coils'); c = g5.children
            valg(c, 'coil_gap', 'Coil standoff gap', d['coil_gap_mm'])
            valg(c, 'coil_r', 'Conductor radius', d['coil_r_mm'])

            g6 = inp.addGroupCommandInput('g_stages', 'Stages to build'); s = g6.children
            s.addBoolValueInput('do_sweep', 'Simple sweep demo', True, '', True)
            s.addBoolValueInput('do_uniform', 'Uniform radial build', True, '', True)
            s.addBoolValueInput('do_nonuniform', 'Non-uniform radial build', True, '', True)
            s.addBoolValueInput('do_coils', 'Coils (swept)', True, '', True)

            onExec = CommandExecuteHandler()
            cmd.execute.add(onExec)
            _handlers.append(onExec)
        except:
            if _ui:
                _ui.messageBox('Command setup failed:\n{}'.format(traceback.format_exc()))


# ============================== add-in lifecycle ===================================================
def run(context):
    global _app, _ui
    try:
        _app = adsk.core.Application.get()
        _ui = _app.userInterface
        cmdDef = _ui.commandDefinitions.itemById(CMD_ID)
        if not cmdDef:
            cmdDef = _ui.commandDefinitions.addButtonDefinition(
                CMD_ID, 'ParaStell Geometry Demo',
                'Generate the stellarator geometry teaching demo (sweep vs loft, radial build, coils).')
        onCreated = CommandCreatedHandler()
        cmdDef.commandCreated.add(onCreated)
        _handlers.append(onCreated)
        panel = _ui.allToolbarPanels.itemById(PANEL_ID)
        if panel and not panel.controls.itemById(CMD_ID):
            panel.controls.addCommand(cmdDef)
        # Uncomment to auto-open the dialog when the add-in starts:
        # cmdDef.execute()
    except:
        if _ui:
            _ui.messageBox('Add-in start failed:\n{}'.format(traceback.format_exc()))


def stop(context):
    try:
        panel = _ui.allToolbarPanels.itemById(PANEL_ID)
        if panel:
            ctrl = panel.controls.itemById(CMD_ID)
            if ctrl:
                ctrl.deleteMe()
        cmdDef = _ui.commandDefinitions.itemById(CMD_ID)
        if cmdDef:
            cmdDef.deleteMe()
        _handlers.clear()
    except:
        if _ui:
            _ui.messageBox('Add-in stop failed:\n{}'.format(traceback.format_exc()))
