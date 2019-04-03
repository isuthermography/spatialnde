#! /usr/bin/env python
#
#
# Usage: python landmarks_from_dgsfile.py <dgsfile> <specimenid_with_dashes_shifted_to_underscores>

import sys
import dg_file as dgf
import dg_metadata as dgm
from spatialnde.dataguzzler import dg_3d

dgsfilename=sys.argv[1]
specimenid=sys.argv[2]

(metadata,wfmdict)=dgf.loadsnapshot(dgsfilename)

# load from specimenname_tex channel
wfmname="%s_tex" % (specimenid)

(landmarks_params,implpartparams) = dg_3d.ndepartparams_from_landmarked3d(wfmdict,[ wfmname ])
(landmarkdict2d,landmarkdict3d) = landmarks_params
print("Resolution: (%d,%d)" % (wfmdict[wfmname].dimlen[0],wfmdict[wfmname].dimlen[1]))
print("ScaleFactors: %.8g %.8g" % (dgm.GetMetaDatumWIDbl(wfmdict[wfmname],"Step1",1.0)*wfmdict[wfmname].dimlen[0],
                                   dgm.GetMetaDatumWIDbl(wfmdict[wfmname],"Step2",1.0)*wfmdict[wfmname].dimlen[1]))
# landmarkdict3d not currently used
print("IniVal: %.8g %.8g" % (dgm.GetMetaDatumWIDbl(wfmdict[wfmname],"IniVal1",0.0),
                             dgm.GetMetaDatumWIDbl(wfmdict[wfmname],"IniVal2",0.0)))
print("Step: %.8g %.8g" % (dgm.GetMetaDatumWIDbl(wfmdict[wfmname],"Step1",1.0),
                           dgm.GetMetaDatumWIDbl(wfmdict[wfmname],"Step2",1.0)))

for landmarkname in landmarkdict2d:
    (surfacenum,landmarku,landmarkv) = landmarkdict2d[landmarkname]
    print("LANDMARK %s (%.12g,%.12g)" % (landmarkname,landmarku,landmarkv))
    pass


