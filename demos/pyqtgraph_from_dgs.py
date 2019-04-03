import sys
import os
from io import BytesIO

import numpy as np

from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl


import dg_file as dgf
import dg_metadata as dgm
import dg_image

from spatialnde.coordframes import coordframe
from spatialnde.ndeobj import ndepart,ndeassembly
from spatialnde.exporters.x3d import X3DSerialization
from spatialnde.cadpart.appearance import simple_material


# unitsperintensity is the contrast setting
# offset is the brightness offset (as in dg_scope)
# colormap is the colormap ("hot", "gray", "colorful")
# frameno is the frame index (0-based; note that dg_scope frame number 
# displayed is one based, so subtract one from that to get this parameter)

dgsname=sys.argv[1]

proj_channel_name = sys.argv[2]

frameno = int(sys.argv[3])

colormap = sys.argv[4]

unitsperintensity = float(sys.argv[5])

offset=float(sys.argv[6])


# Usage: python pyqtgraph_from_dgs.py <dgsfilename> <projectionchannel> <frameno> <colormap> <unitsperintensity> <offset>
# e.g. 
# python /tmp/pyqtgraph_from_dgs.py 0012-flash_collect_data_spatial_orient__flashdata_dgsfile.dgs ProjAccumDiffStack_T17_F7611_6F 67 hot 10 4.9

# This one's only a single shot, but it does illustrate compressed loading
# and raw data (may require Python 3 and updates to other packages)
# python /tmp/pyqtgraph_from_dgs.py 0004-flash_collect_data_spatialT17-F7611-4H_orient_front_convex_flashdata_dgsfile.dgs.bz2 ProjDiffStack_T17_F7611_ideal 67 hot 10 4.9

# Post-analysis and data fusion:
# python /tmp/pyqtgraph_from_dgs.py 0005-flash_collect_data_spatialT17-F7611-4H_orient_back_flashdata_dgsfile_greensinversion_tik_7.5e-11_fused.dgs Projgreensinversion_DiffStack_T17_F7611_ideal 4 hot 5000 1800

(junkmd,wfmdict)=dgf.loadsnapshot(dgsname,memmapok=True)


# Trace through cross referenes to underlying geometry channel
geomchan = proj_channel_name
texchanprefix=""
while len(dgm.GetMetaDatumWIStr(wfmdict[geomchan],"X3DGeomRef","")) > 0:
    texchanprefix+=dgm.GetMetaDatumWIStr(wfmdict[geomchan],"TexChanPrefix","")
    geomchan = dgm.GetMetaDatumWIStr(wfmdict[geomchan],"X3DGeomRef","")
    pass


# Read geometry from metadata into an ndepart object                        
X3DGeom = dgm.GetMetaDatumWIStr(wfmdict[geomchan],"X3DGeom","")
X3DGeom_fh = BytesIO(X3DGeom.encode('utf-8'))

objframe=coordframe()
X3DObj = ndepart.fromx3d(objframe,None,X3DGeom_fh,tol=1e-6)
X3DGeom_fh.close()


# Assume a triangular mesh, extract vertexes and texture vertexes
assert(X3DObj.implpart.surfaces[0].vertexidx.shape[0] % 4 == 0)
n_tris = X3DObj.implpart.surfaces[0].vertexidx.shape[0] // 4
indices_reshape=X3DObj.implpart.surfaces[0].vertexidx.reshape(n_tris,4)
assert((indices_reshape[:,3]==-1).all())  # terminating -1s on face indices

verts=X3DObj.implpart.surfaces[0].vertices[indices_reshape[:,:3].reshape(n_tris*3),:].reshape(n_tris,3,3)

texindices_reshape=X3DObj.implpart.surfaces[0].intrinsicparameterization.texcoordidx.reshape(n_tris,4)
assert((texindices_reshape[:,3]==-1).all())  # terminating -1s on face indices
texverts=X3DObj.implpart.surfaces[0].intrinsicparameterization.texcoord[texindices_reshape[:,:3].reshape(n_tris*3),:].reshape(n_tris,3,2)

normals = X3DObj.implpart.surfaces[0].facetnormals


# Extract texture image from given URL
appearanceurl = X3DObj.implpart.surfaces[0].appearance.texture_url

assert(appearanceurl[0]=='#') # appearance URL should be a fragment indicating a channel
appearance_chan = texchanprefix + appearanceurl[1:]

# Extract texture image, convert to RGBA Image per a particular colormap, offset and intensity scaling

# (Suggest widgets to control frameno, colormap, offset, and unitsperintensity)
# and rerun these next few lines when they are updated. 


if len(wfmdict[appearance_chan].data.shape) > 2:
    frameindex = (frameno,)
    pass
else:
    frameindex = ()  # empty tuple
    pass
PILimg = dg_image.toimage(wfmdict[appearance_chan],wfmdict,frameindex,unitsperintensity,offset,colormap=colormap).convert("RGBA")
appearance_chan_data=np.frombuffer(PILimg.tobytes("raw","RGBA",0,-1),dtype=np.uint8).reshape(PILimg.size[1],PILimg.size[0],4)



app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.show()
w.setWindowTitle('pyqtgraph example: GLTexturedMeshItem')
w.setCameraPosition(distance=40)

g = gl.GLGridItem()
g.scale(2,2,1)
w.addItem(g)


## Mesh item will automatically compute face normals.
m1 = gl.GLTexturedMeshItem(vertexes=verts, normals=normals, texvertexes = texverts, texture = appearance_chan_data, shader="textureMap", faces=None,color=[1.0,1.0,1.0,1.0])

# NOTE: Can update texture image  by calling m1.setTexture() again later (recommended from main event/rendering thread only!) 

#m1 = gl.GLMeshItem(vertexes=verts,normals=normals,faces=None,color=[0.5,0.5,0.5,1.0],shader="shaded")
# m1.translate(5, 5, 0)
#m1.setGLOptions('additive')
w.addItem(m1)


## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
        pass
    pass

