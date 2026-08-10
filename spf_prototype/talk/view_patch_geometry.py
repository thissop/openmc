"""Interactive 3D view of the QH stellarator wall patches + coils.

Run:  python3 view_patch_geometry.py
A window opens that you can drag to rotate, scroll to zoom. Close it when done.
Nothing is saved. To save an angle you like, note the elevation/azimuth printed on close,
or use the toolbar's save button.
"""
import matplotlib.pyplot as plt   # uses your default interactive backend (a draggable window)
from gen_patch_fig import draw

fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111, projection="3d")
draw(ax)
ax.view_init(elev=34, azim=-58)

def _on_close(_):
    print(f"final view  elev={ax.elev:.0f}  azim={ax.azim:.0f}")
fig.canvas.mpl_connect("close_event", _on_close)

print("drag to rotate, scroll to zoom, close window to exit")
plt.show()
