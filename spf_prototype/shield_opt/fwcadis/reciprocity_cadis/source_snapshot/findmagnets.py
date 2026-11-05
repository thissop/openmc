import openmc, numpy as np, openmc.lib
exec(open("build_materials_snippet.py").read())
# read first 5 coils (period 1) from coils_qh
def read_makegrid(path):
    lines=[l.rstrip("\n") for l in open(path)]
    start=next(i for i,l in enumerate(lines) if l.strip()=="mirror NIL")+1
    end=next(i for i,l in enumerate(lines) if l.strip()=="end")
    coils=[]; cur=[]
    for l in lines[start:end]:
        t=l.split()
        if len(t)<4: continue
        cur.append([float(x) for x in t[:4]])
        if len(t)>4: coils.append(np.array(cur)); cur=[]
    return coils
coils=read_makegrid("/burg-archive/home/tjk2147/pstl_test/coils_qh")[:5]   # period-1 coils
mats=build_materials()
dag=openmc.DAGMCUniverse("dagmc_qh_op_uniform.h5m", auto_geom_ids=True)
bb=dag.bounding_box; ext=float(np.max(np.abs(np.concatenate([np.array(bb.lower_left),np.array(bb.upper_right)]))))
sph=openmc.Sphere(r=ext*1.5+10, boundary_type="vacuum")
geom=openmc.Geometry([openmc.Cell(region=-sph, fill=dag)])
s=openmc.Settings(); s.run_mode="fixed source"; s.batches=1; s.particles=1
s.source=[openmc.IndependentSource(space=openmc.stats.Point((1367,50,0)))]
m=openmc.Model(geometry=geom, materials=mats, settings=s)
m.export_to_model_xml(); openmc.lib.init(["-s","1"])
mag=set()
for k,c in enumerate(coils):
    for idx in (0, len(c)//4, len(c)//2, 3*len(c)//4):    # a few points around each filament
        x,y,z = c[idx,0]*100, c[idx,1]*100, c[idx,2]*100   # m -> cm
        try:
            cell,inst=openmc.lib.find_cell((x,y,z))
            nm=cell.fill.name if (cell.fill is not None and hasattr(cell.fill,"name")) else "VOID"
            if nm=="magnets": mag.add(cell.id)
        except Exception as e:
            pass
print("MAGNET_CELLS", sorted(mag))
openmc.lib.finalize()
