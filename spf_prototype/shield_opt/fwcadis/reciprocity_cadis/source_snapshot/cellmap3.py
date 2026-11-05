import openmc, numpy as np, openmc.lib
def build_materials():
    w=openmc.Material();w.add_element("W",1.0);w.set_density("g/cm3",19.3)
    rafm=openmc.Material()
    for el,wo in [("Fe",.885),("Cr",.09),("W",.011),("Mn",.006),("V",.002),("C",.006)]: rafm.add_element(el,wo,"wo")
    rafm.set_density("g/cm3",7.9)
    fw=openmc.Material.mix_materials([w,rafm],[0.2/3.2,3.0/3.2],"vo");fw.name="first_wall"
    mult=openmc.Material(name="multiplier");mult.add_element("Be",1.0);mult.set_density("g/cm3",1.85)
    br=openmc.Material(name="breeder");br.add_element("Li",2.0,enrichment=60,enrichment_target="Li6",enrichment_type="ao");br.add_element("Be",1.0);br.add_element("F",4.0);br.set_density("g/cm3",1.94)
    bw=openmc.Material(name="back_wall")
    for el,wo in [("Fe",.885),("Cr",.09),("W",.011),("Mn",.006),("V",.002),("C",.006)]: bw.add_element(el,wo,"wo")
    bw.set_density("g/cm3",7.9)
    sh=openmc.Material(name="shield");sh.add_element("W",1.0);sh.add_element("C",1.0);sh.set_density("g/cm3",15.6)
    gap=openmc.Material(name="gap");gap.add_nuclide("H1",1.0);gap.set_density("g/cm3",1e-8)
    vv=openmc.Material(name="vac_vessel")
    for el,wo in [("Fe",.62),("Cr",.16),("Ni",.11),("Mo",.02),("Mn",.015),("B",.03),("O",.03),("H",.005)]: vv.add_element(el,wo,"wo")
    vv.set_density("g/cm3",7.9)
    ts=openmc.Material(name="thermal_shield");ts.add_element("Fe",.89,"wo");ts.add_element("Cr",.09,"wo");ts.add_element("Mn",.02,"wo");ts.set_density("g/cm3",7.9)
    mg=openmc.Material(name="magnets");mg.add_element("Cu",1.0);mg.set_density("g/cm3",8.5)
    vac=openmc.Material(name="Vacuum");vac.add_nuclide("H1",1.0);vac.set_density("g/cm3",1e-10)
    return openmc.Materials([fw,mult,br,bw,sh,gap,vv,ts,mg,vac])
mats=build_materials()
dag=openmc.DAGMCUniverse("dagmc_qh_op_uniform.h5m", auto_geom_ids=True)
bb=dag.bounding_box; ext=float(np.max(np.abs(np.concatenate([np.array(bb.lower_left),np.array(bb.upper_right)]))))
sph=openmc.Sphere(r=ext*1.5+10, boundary_type="vacuum")
geom=openmc.Geometry([openmc.Cell(region=-sph, fill=dag)])
s=openmc.Settings(); s.run_mode="fixed source"; s.batches=1; s.particles=1
s.source=[openmc.IndependentSource(space=openmc.stats.Point((1367,0,0)))]
m=openmc.Model(geometry=geom, materials=mats, settings=s)
m.export_to_model_xml(); openmc.lib.init(["-s","1"])
# radial ray at phi=45 deg, z=0, from R=1150 outward to 1650
prev=None
for R in np.arange(1150.0, 1660.0, 2.0):
    x,y = R*np.cos(np.pi/4), R*np.sin(np.pi/4)
    try:
        cell,inst = openmc.lib.find_cell((x,y,0.0))
        cid = cell.id
        matname = cell.fill.name if (cell.fill is not None and hasattr(cell.fill,"name")) else "VOID/None"
    except Exception:
        cid, matname = None, "OUTSIDE"
    key=(cid,matname)
    if key!=prev:
        print("R=%.0f  cell=%s  mat=%s" % (R, cid, matname)); prev=key
openmc.lib.finalize()
