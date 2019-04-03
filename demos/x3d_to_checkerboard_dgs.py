# Load in a .x3d file, store to a .dgs file with a checkerboard pattern

import sys
import os
import numpy as np
import dataguzzler as dg
import dg_file as dgf
import dg_metadata as dgm

try:
    from cStringIO import StringIO  # python 2.x
    pass
except ImportError:
    from io import StringIO # python 3.x
    pass


from spatialnde.coordframes import coordframe
from spatialnde.ndeobj import ndepart

from spatialnde.exporters.vrml import VRMLAppearance


imageshape=(1024,1024)
cbrows=20
cbcols=20

x3dname=sys.argv[1]
(x3dpath,x3dfile)=os.path.split(x3dname)
x3dbasename=os.path.splitext(x3dfile)[0]


objframe=coordframe() 

obj = ndepart.fromx3d(objframe,None,x3dname,tol=1e-6)
VRMLBuf=StringIO()

VRMLAppear=VRMLAppearance.Simple(DefName="ImageMap")
obj.VRMLFile(VRMLBuf,objframe,UVparameterization=None,appearance=VRMLAppear)
VRMLAppear.reset() # needed if we want to write a new file

wfmdict={}
wfmdict[x3dbasename]=dg.wfminfo()
wfmdict[x3dbasename].Name=x3dbasename
wfmdict[x3dbasename].dimlen=np.array((),dtype='i8')
dgm.AddMetaDatumWI(wfmdict[x3dbasename],dgm.MetaDatum("VRML97Geom",VRMLBuf.getvalue()))

x3dparamname=x3dbasename+"_parameterization"
wfmdict[x3dparamname]=dg.wfminfo()
wfmdict[x3dparamname].Name=x3dparamname
wfmdict[x3dparamname].ndim=2
wfmdict[x3dparamname].dimlen=np.array(imageshape,dtype='i8')
wfmdict[x3dparamname].n=np.prod(wfmdict[x3dparamname].dimlen)
wfmdict[x3dparamname].data=np.zeros(wfmdict[x3dparamname].dimlen,dtype='f')

xpos=np.arange(imageshape[0],dtype='d')
ypos=np.arange(imageshape[1],dtype='d')

xchecker = (xpos//(imageshape[0]*1.0/(cbcols)) % 2).astype(np.bool)
ychecker = (ypos//(imageshape[1]*1.0/(cbrows)) % 2).astype(np.bool)


wfmdict[x3dparamname].data[:,:]=xchecker.reshape(imageshape[0],1) ^ ychecker.reshape(1,imageshape[1])  # XOR operator

dgm.AddMetaDatumWI(wfmdict[x3dbasename],dgm.MetaDatum("TextureChan_ImageMap",x3dparamname+":0"))

dgf.savesnapshot(os.path.join(x3dpath,x3dbasename+".dgs"),wfmdict)

                   



