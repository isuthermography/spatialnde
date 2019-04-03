import sys
import os
from io import BytesIO

from spatialnde.coordframes import coordframe
from spatialnde.ndeobj import ndepart,ndeassembly
from spatialnde.exporters.x3d import X3DSerialization
from spatialnde.cadpart.appearance import simple_material

# Example of how to read in an X3D file, apply an appearance, 
# and write it out


x3dname=sys.argv[1]
(x3dpath,x3dfile)=os.path.split(x3dname)
x3dbasename=os.path.splitext(x3dfile)[0]


use_bytesio=True
if use_bytesio:
    fh=open(x3dname,"r")
    x3d_bytes = fh.read()
    x3d_bytesio = BytesIO(x3d_bytes)
    pass

    
objframe=coordframe() 

if use_bytesio:
    obj = ndepart.fromx3d(objframe,None,x3d_bytesio,tol=1e-6)
    pass
else:
    obj = ndepart.fromx3d(objframe,None,x3dname,tol=1e-6)
    pass

# surface data in obj.implpart.surfaces[0] 
# which is a spatialnde.cadpart.polygonalsurface.polygonalsurface object

# parameterization (texture mapping), if present in the x3d file, in
# obj.implpart.surfaces[0].intrinsicparameterization which is a
# spatialnde.cadpart.polygonalsurface_texcoordparameterization.polygonalsurface_texcoordparameterization
# object
# 

# vertex array in obj.implpart.surfaces[0].vertices

# vertex ids per polygon in 
# obj.implpart.surfaces[0].vertexidx_indices and .vertexidx
#
# i.e. obj.implpart.surfaces[0].vertexidx[obj.implpart.surfaces[0].vertexidx_indices[27,:]] gives the 
# first vertex indices involved in mesh element #27. Subsequent vertex
# indices follow, terminated by -1

# To write out X3D: 


x3dwriter=X3DSerialization.tofileorbuffer("/tmp/objout.x3d",x3dnamespace=None)

# You can add additional tags to the X3D by 
# creating the XML tags then appending them to x3dwriter.scene
Appear = simple_material.from_color((.6,.01,.01),DefName="coloring")

obj.assign_appearances(Appear)

obj.X3DWrite(x3dwriter,objframe)
x3dwriter.finish()
