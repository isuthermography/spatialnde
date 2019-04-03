#! /usr/bin/env python
#
#
# Usage: python dgs_update_geometry.py <dgsfile> <geomchannel> <replacement_x3d>

# Makes a backup of dgsfile (.dgs.oldgeometry) and will fail if backup file
# already exists
#
# IMPORTANT: DO NOT REPLACE THE GEOMETRY WITH A DIFFERENT UV PARAMETERIZATION
# (i.e. different texture mapping) as the (u,v) coordinates are how
# stored data and/or landmark coords are mapped to the object. 

import sys
import os
import copy
import dg_file as dgf
import dg_metadata as dgm
from spatialnde.dataguzzler import dg_3d
from spatialnde.coordframes import coordframe
from spatialnde.ndeobj import ndepart,ndeassembly

from spatialnde.exporters.vrml import VRMLSerialization
from spatialnde.exporters.x3d import X3DSerialization

try:
    from cStringIO import StringIO  # python 2.x
    pass
except ImportError:
    from io import StringIO # python 3.x
    pass


dgsfilename=sys.argv[1]
geomchannel=sys.argv[2]

replacement_x3d=sys.argv[3]

backup_dgs = dgsfilename+".oldgeometry"

if os.path.exists(backup_dgs):
    sys.stderr.write("Backup .dgs file %s already exists. Exiting.\n" % (backup_dgs))
    sys.exit(1)
    pass


os.rename(dgsfilename,backup_dgs)

(metadata,wfmdict)=dgf.loadsnapshot(backup_dgs)
geomchan=wfmdict[geomchannel]
assert("X3DGeom" in geomchan.MetaData) # Specify the REAL geometry channel, not all the references to it!!!


# Load replacement geometry

objframe=coordframe() 
newobj = ndepart.fromx3d(objframe,None,replacement_x3d,tol=1e-6)

# need to transfer texure urls, which might be different
(oldobj,texchanprefix) = dg_3d.ndepart_from_dataguzzler_wfm(geomchan,wfmdict,objframe)

# assume new and old surfaces line up
assert(len(newobj.implpart.surfaces)==len(oldobj.implpart.surfaces))

# copy appearances (texture URL's)
for cnt in range(len(newobj.implpart.surfaces)):
    oldappearance = copy.copy(oldobj.implpart.surfaces[cnt].appearance)
    newobj.implpart.surfaces[0].assign_appearance(oldappearance)
    pass


# Update X3D and VRML metadata
VRMLBuf=StringIO()
vrmlwriter = VRMLSerialization.tofileorbuffer(VRMLBuf)

newobj.VRMLWrite(vrmlwriter,objframe,UVparameterization=None)
vrmlwriter.finish()

X3DBuf=StringIO()
x3dnamespace=None # just use default
x3dwriter=X3DSerialization.tofileorbuffer(X3DBuf,x3dnamespace=x3dnamespace)
newobj.X3DWrite(x3dwriter,objframe,UVparameterization=None)
x3dwriter.finish()

dgm.AddMetaDatumWI(geomchan,dgm.MetaDatum("VRML97Geom",VRMLBuf.getvalue()))
dgm.AddMetaDatumWI(geomchan,dgm.MetaDatum("X3DGeom",X3DBuf.getvalue()))
    
assert(not(os.path.exists(dgsfilename)))

dgf.savesnapshot(dgsfilename,wfmdict,Metadata=metadata)
