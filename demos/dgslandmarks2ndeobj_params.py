#! /usr/bin/env python

# Usage: dgslandmarks2ndeobj_params.py <dgsfile> <texchan>
# to load landmarks and scalefactor from a .dgsfile
#
# <texchan> e.g. CR03_SPAR_01H_tex

import sys
import json

import dg_file as dgf
import dg_eval

dgsname=sys.argv[1]
texchan=sys.argv[2]

(md,wfmdict)=dgf.loadsnapshot(dgsname)

(ndim,DimLen,IniVal,Step,bases) = dg_eval.geom(wfmdict[texchan])
scalefactor_x=Step[0]*DimLen[0]
scalefactor_y=Step[1]*DimLen[1]
                                               
LandmarkDict={}

# We don't distinguish here between landmarks and fiducials
# for the time being

# Also assume only a single surface!!!
surfaceid = 0

for mdname in wfmdict[texchan].MetaData.keys():
    if ((mdname.startswith("LANDMARK_")
         or mdname.startswith("FIDUCIAL_")) and 
        mdname.endswith("_X")):
        LandmarkName=mdname[9:-2]
        LandmarkCoords=(surfaceid,wfmdict[texchan].MetaData[mdname].Value,wfmdict[texchan].MetaData[mdname[:-1]+"Y"].Value)
        LandmarkDict[LandmarkName]=LandmarkCoords
        pass
                                               
# Now dump JSON
# ... suitable for ndepart.fromx3d()

uvparam_params=None # so far...
landmarks_params=(LandmarkDict,{}) # 2nd element is landmarkdict3d, which we do not yet generate
# implpartparams are the parameters of the intrinsic parameterization,
# (lowerleft_meaningfulunits, meaningfulunits_per_texcoord)
# Remember that lowerleft in terms of IniVal is the coordinate of the first point,
# whereas lowerleft in terms of texture is the lower left of the first point

implpartparams = ( (IniVal[0]-Step[0]/2.0,IniVal[1]-Step[1]/2.0), (scalefactor_x,scalefactor_y))
                                               
ndepartparams=(uvparam_params,landmarks_params,implpartparams)
                                               
ndepart_paramstr = json.dumps(ndepartparams)
print(ndepart_paramstr)
