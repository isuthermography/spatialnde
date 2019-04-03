# Load in a .x3d file into dataguzzler.
# For use with dataguzzler_define_object_landmarks.confm4
import sys
import os
import ast
import numpy as np
import dataguzzler as dg
#import dg_file as dgf
import dg_comm as dgc
import dg_metadata as dgm

try:
    from cStringIO import StringIO  # python 2.x
    pass
except ImportError:
    from io import StringIO # python 3.x
    pass


from spatialnde.coordframes import coordframe
from spatialnde.ndeobj import ndepart,ndeassembly

from spatialnde.cadpart.appearance import texture_url

from spatialnde.dataguzzler.dg_3d import create_x3d_vrml_channel



if len(sys.argv) < 2:
    sys.stderr.write("""
USAGE: python dataguzzler_define_landmarks_load_x3d.py <x3dfilename> <imageshape>  [scalefactorx] [scalefactory]

    imageshape should be a tuple e.g. (1024,1024).
    scalefactorx and scalefactory can be provided 
    in a second run once you have determined the scaling
    from your UV parameterization. See WORKFLOW.txt
""")
    sys.exit(0)
    
    

x3dname=sys.argv[1]
(x3dpath,x3dfile)=os.path.split(x3dname)
x3dbasename=os.path.splitext(x3dfile)[0]

#imageshape=(1024,1024)
imageshape=ast.literal_eval(sys.argv[2])
assert(isinstance(imageshape,tuple))
assert(len(imageshape)==2)
assert(isinstance(imageshape[0],int))
assert(isinstance(imageshape[1],int))


cbrows=20
cbcols=20


# strip trailing _curvature from basename, if present
if x3dbasename.lower().endswith("_curvature"):
    x3dbasename=x3dbasename[:-10]
    pass


# strip trailing _uv from basename, if present
if x3dbasename.endswith("_uv") or x3dbasename.endswith("_UV"):
    x3dbasename=x3dbasename[:-3]
    pass

# ***!!! NOTE: We should really clear out any
# Characters that can't be in dataguzzler identifiers (waveform names)
# from x3dbasename ***!!!

x3dbasename=x3dbasename.replace(" ","_").replace("-","_")



objframe=coordframe() 

texturename = "%s_tex" % (x3dbasename)
texturerawname = "%s_tex_raw" % (x3dbasename)
texturesrcname = "%s_tex_src" % (x3dbasename)
Appear = texture_url.from_url("#"+texturename,DefName="TextureMap")

obj = ndepart.fromx3d(objframe,None,x3dname,tol=1e-6)
obj.assign_appearances(Appear)

assembly = ndeassembly(parts=[obj])

#VRMLBuf=StringIO()
#vrmlwriter=VRMLSerialization.tofileorbuffer(VRMLBuf)
#obj.VRMLWrite(vrmlwriter,objframe,UVparameterization=None)
#vrmlwriter.finish()

#X3DBuf=StringIO()
#x3dwriter=X3DSerialization.tofileorbuffer(X3DBuf)
#obj.X3DWrite(x3dwriter,objframe,UVparameterization=None)
#x3dwriter.finish()


wfmdict={}

# Create a channel holding the X3D+VRML object representation
x3dchan=create_x3d_vrml_channel(wfmdict,x3dbasename,assembly,objframe)
# Program object channel so that keyboard commands in scope display 
# adjust texture channel
dgm.AddMetaDatumWI(x3dchan,dgm.MetaDatum("Scope_KeyboardControlChan",texturename))

#wfmdict["x3dfile"]=dg.wfminfo()
#wfmdict["x3dfile"].Name="x3dfile"
#wfmdict["x3dfile"].dimlen=np.array((1,),dtype='i8')
#wfmdict["x3dfile"].ndim=1

#wfmdict["x3dfile"].wfmrevision=0
#dgm.AddMetaDatumWI(wfmdict["x3dfile"],dgm.MetaDatum("VRML97Geom",VRMLBuf.getvalue()))
#wfmdict["x3dfile"].data=np.array((0,),dtype='f')

#x3dparamname="Param"
#x3dparamrawname="ParamRaw"
wfmdict[texturerawname]=dg.wfminfo()
wfmdict[texturerawname].Name=texturerawname
wfmdict[texturerawname].ndim=2
wfmdict[texturerawname].dimlen=np.array(imageshape,dtype='i8')
wfmdict[texturerawname].wfmrevision=0
wfmdict[texturerawname].n=np.prod(wfmdict[texturerawname].dimlen)
wfmdict[texturerawname].data=np.zeros(wfmdict[texturerawname].dimlen,dtype='f')
dgm.AddMetaDatumWI(wfmdict[texturerawname],dgm.CreateMetaDatumStr("Coord1","U Position"))
dgm.AddMetaDatumWI(wfmdict[texturerawname],dgm.CreateMetaDatumStr("Coord2","V Position"))

if len(sys.argv) > 3:
    ScaleFactorX=float(sys.argv[3])
    ScaleFactorY=float(sys.argv[4])
    dgm.AddMetaDatumWI(wfmdict[texturerawname],dgm.CreateMetaDatumStr("Units1","meters"))
    dgm.AddMetaDatumWI(wfmdict[texturerawname],dgm.CreateMetaDatumDbl("Step1",ScaleFactorX/imageshape[0]))
    dgm.AddMetaDatumWI(wfmdict[texturerawname],dgm.CreateMetaDatumDbl("IniVal1",-0.5*ScaleFactorX))
    dgm.AddMetaDatumWI(wfmdict[texturerawname],dgm.CreateMetaDatumStr("Units2","meters"))
    dgm.AddMetaDatumWI(wfmdict[texturerawname],dgm.CreateMetaDatumDbl("Step2",ScaleFactorY/imageshape[1]))
    dgm.AddMetaDatumWI(wfmdict[texturerawname],dgm.CreateMetaDatumDbl("IniVal2",-0.5*ScaleFactorY))
    pass
else:
    dgm.AddMetaDatumWI(wfmdict[texturerawname],dgm.CreateMetaDatumStr("Units1","pixels"))
    dgm.AddMetaDatumWI(wfmdict[texturerawname],dgm.CreateMetaDatumDbl("Step1",1.0/imageshape[0]))
    dgm.AddMetaDatumWI(wfmdict[texturerawname],dgm.CreateMetaDatumDbl("IniVal1",-0.5/imageshape[0]))
    dgm.AddMetaDatumWI(wfmdict[texturerawname],dgm.CreateMetaDatumStr("Units2","pixels"))
    dgm.AddMetaDatumWI(wfmdict[texturerawname],dgm.CreateMetaDatumDbl("Step2",1.0/imageshape[1]))
    dgm.AddMetaDatumWI(wfmdict[texturerawname],dgm.CreateMetaDatumDbl("IniVal2",-0.5/imageshape[1]))
    pass


xpos=np.arange(imageshape[0],dtype='d')
ypos=np.arange(imageshape[1],dtype='d')

xchecker = (xpos//(imageshape[0]*1.0/(cbcols)) % 2).astype(np.bool)
ychecker = (ypos//(imageshape[1]*1.0/(cbrows)) % 2).astype(np.bool)


wfmdict[texturerawname].data[:,:]=xchecker.reshape(imageshape[0],1) ^ ychecker.reshape(1,imageshape[1])  # XOR operator

#dgm.AddMetaDatumWI(wfmdict[x3dbasename],dgm.MetaDatum("TextureChan_ImageMap",x3dparamname+":-1"))

#dgf.savesnapshot(os.path.join(x3dpath,x3dbasename+".dgs"),wfmdict)

dgch=dgc.client();

for wfmname in wfmdict:
    dgc.uploadwfm(dgch,wfmdict[wfmname])
    pass

dgc.command(dgch,"WFM:COPY %s %s" % (texturerawname,texturesrcname))
dgc.close(dgch)




