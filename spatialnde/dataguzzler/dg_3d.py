import copy
import json



#try:
#    from cStringIO import StringIO  # python 2.x
#    pass
#except ImportError:
#    from io import StringIO # python 3.x
#    pass
from io import BytesIO

import numpy as np
import dg_file as dgf
import dg_metadata as dgm
import dataguzzler as dg

from ..exporters.vrml import VRMLSerialization
from ..exporters.x3d import X3DSerialization
from ..ndeobj import ndepart


def create_x3d_vrml_channel(wfmdict,channame,assembly,coordframe,x3dnamespace=None):
    wfm=dg.wfminfo()
    wfm.Name=channame
    wfm.dimlen=np.array((1,),dtype='i8')
    wfm.ndim=1
    wfm.data=np.array((0,),dtype='f')
    wfm.wfmrevision=0

    VRMLBuf=BytesIO()
    vrmlwriter = VRMLSerialization.tofileorbuffer(VRMLBuf)

    for part in assembly.parts:
        part.VRMLWrite(vrmlwriter,coordframe,UVparameterization=None)
        pass
    vrmlwriter.finish()

    X3DBuf=BytesIO()
    x3dwriter=X3DSerialization.tofileorbuffer(X3DBuf,x3dnamespace=x3dnamespace)
    for part in assembly.parts:
        part.X3DWrite(x3dwriter,coordframe,UVparameterization=None)
        pass
    x3dwriter.finish()

    dgm.AddMetaDatumWI(wfm,dgm.MetaDatum("VRML97Geom",VRMLBuf.getvalue().decode('utf-8')))
    dgm.AddMetaDatumWI(wfm,dgm.MetaDatum("X3DGeom",X3DBuf.getvalue().decode('utf-8')))
    
    TextureNames=set([])
    for part in assembly.parts:
        for surface in part.implpart.surfaces:
            if hasattr(surface.appearance,"texture_url") and surface.appearance.texture_url.startswith("#"):
                TextureNames.add(surface.appearance.texture_url[1:])
                pass
            pass
        pass
    
    # Add metadata listing the texture names
    dgm.AddMetaDatumWI(wfm,dgm.MetaDatum("TextureChans","|".join(TextureNames)))
    
    if wfmdict is not None:
        wfmdict[channame]=wfm
        pass

    return wfm


def ndepartparams_from_landmarked3d(l3d_wfmdict,wfmname_list,TexChanPrefix=""):

    

    landmarkdict2d={}
    landmarkdict3d={}

    UVScalingDict={}
    
    for wfmname in wfmname_list:
        if not(wfmname.startswith(TexChanPrefix)):
            raise ValueError("Error processing texture channel: Texture channel name does not start with specified prefix (%s)" % (wfmname,TexChanPrefix))
        unprefixedname=wfmname[len(TexChanPrefix):]
        
        wfm=l3d_wfmdict[wfmname]

        IniVal1=dgm.GetMetaDatumWIDbl(wfm,"IniVal1",0.0)
        Step1=dgm.GetMetaDatumWIDbl(wfm,"Step1",1.0)
        IniVal2=dgm.GetMetaDatumWIDbl(wfm,"IniVal2",0.0)
        Step2=dgm.GetMetaDatumWIDbl(wfm,"Step2",1.0)

        # lowerleft defined at corner of pixel, not pixel center
        lowerleft_meaningfulunits=(IniVal1-Step1/2.0,IniVal2-Step2/2.0)
        meaningfulunits_per_texcoord = (Step1*wfm.dimlen[0],Step2*wfm.dimlen[1])

        # index must be a string because params must be JSON-compatible
        texurl = "#"+unprefixedname
        UVScalingDict[texurl]=(lowerleft_meaningfulunits,meaningfulunits_per_texcoord)
        
        for mdname in wfm.MetaData:
            if (mdname.startswith("FIDUCIAL_") or mdname.startswith("LANDMARK_")) and mdname.endswith("_X"):
                fiducialname=mdname[9:-2]
                fiducialmdname=mdname[:-2]
                
                landmarkx=dgm.GetMetaDatumWIDbl(wfm,"%s_X" % (fiducialmdname),np.NaN)
                landmarky=dgm.GetMetaDatumWIDbl(wfm,"%s_Y" % (fiducialmdname),np.NaN)
                
                
                landmarkdict2d[fiducialname]=(texurl,landmarkx,
                                              landmarky)
                pass
            pass
        pass

    # Was UV_ScalingParamsBySurfaceNum
    implpartparams = { "UV_ScalingParamsByTexURL": UVScalingDict }
    landmarks_params = (landmarkdict2d,landmarkdict3d)
    ndepartparams = (landmarks_params,implpartparams)
    return ndepartparams
    

def blank_uv_from_landmarked3d(l3d_wfmdict,wfmname):
    # Replace specified waveform with a blank, ready for stuff to be added 
    oldwfm=l3d_wfmdict[wfmname]

    l3d_wfmdict[wfmname]=dg.wfminfo()
    l3d_wfmdict[wfmname].Name=wfmname
    l3d_wfmdict[wfmname].ndim=2
    l3d_wfmdict[wfmname].dimlen=np.array(oldwfm.dimlen,dtype='i8')
    l3d_wfmdict[wfmname].wfmrevision=0
    l3d_wfmdict[wfmname].n=np.prod(l3d_wfmdict[wfmname].dimlen)
    l3d_wfmdict[wfmname].data=np.zeros(l3d_wfmdict[wfmname].dimlen,dtype='f',order='F')
    l3d_wfmdict[wfmname].MetaData = copy.deepcopy(oldwfm.MetaData)
    return l3d_wfmdict[wfmname]

def ndepart_from_dataguzzler_wfm(wfm_3d,wfmdict,objframe):
    # Really, we shouldn't need wvm_uv but as it is
    # we need to extract ndepartparams from it prior to load (!)

    geomstr=""
    TextureChans=dgm.GetMetaDatumWIStr(wfm_3d,"TextureChans","").split("|")  # TextureChans are already prefixed...


    TexChanPrefix=""
    refwfm_3d = wfm_3d
    while len(geomstr)==0:
        geomstr=dgm.GetMetaDatumWIStr(refwfm_3d,"X3DGeom","")
        TexChanPrefix+=dgm.GetMetaDatumWIStr(refwfm_3d,"TexChanPrefix","")

        if len(geomstr)==0:
            geomrefstr=dgm.GetMetaDatumWIStr(refwfm_3d,"X3DGeomRef","")
            if len(geomrefstr)==0:
                break;
            if geomrefstr not in wfmdict:
                raise ValueError("X3DGeomRef metadata refers to nonexisting waveform %s" % (geomrefstr))
            refwfm_3d=wfmdict[geomrefstr]
            
            pass
        pass

    ndepartparams=ndepartparams_from_landmarked3d(wfmdict,TextureChans,TexChanPrefix)
    
    x3dbuf=BytesIO(geomstr.encode('utf-8'))

    obj = ndepart.fromx3d(objframe,ndepartparams,x3dbuf)

    return (obj,TexChanPrefix)
