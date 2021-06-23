

import sys, os
import numpy as np
import posixpath
import cv2

try:
    # py2.x
    from urllib import pathname2url
    from urllib import url2pathname
    from urllib import quote
    from urllib import unquote
    pass
except ImportError:
    # py3.x
    from urllib.request import pathname2url
    from urllib.request import url2pathname
    from urllib.parse import quote
    from urllib.parse import unquote
    pass


from spatialnde import intrinsic_calibration
from spatialnde import dg_undistort

import dg_file as dgf
from limatix.dc_value import hrefvalue
from limatix.dc_value import numericunitsvalue

def run(_xmldoc, _tag, dc_dgsfile_href, cic_href,grow_factor_numericunits=numericunitsvalue(1.2)):

    # dc_dgsfile is dgsfile element under measurement element

    dc_dgsfile_path = dc_dgsfile_href.getpath()

    # cic file is the camera intrinsic calibration
    #cic_file = os.path.join("/dataawareness/Camera_Calibration_files","flir_a35_63201313_calib_20150614_20150517sdh4.cic")
    cic_file=cic_href.getpath()
    
    # load in calibration
    calib=intrinsic_calibration.intrinsic_calibration.fromfile(cic_file)
    
    # create undistortion processor
    processor=intrinsic_calibration.undistortion_processor(calib,calib.sizex,calib.sizey,int(calib.sizex*grow_factor_numericunits.value()),int(calib.sizey*grow_factor_numericunits.value()),cv2.INTER_LINEAR,1.0)

    # load in .dgs file
    
    (dgs_metadata,dgs_wfmdict)=dgf.loadsnapshot(dc_dgsfile_path,memmapok=True)
    
    new_wfmdict=dg_undistort.undistort_wfmdict(processor,dgs_wfmdict)
    

    undistorteddgs="%s_undistorted.dgs" % (posixpath.splitext(dc_dgsfile_href.get_bare_unquoted_filename())[0])
    undistortedhref=hrefvalue(quoted(undistorteddgs),contexthref=dc_dgsfile_href)
    dgf_write=dgf.creat(undistortedhref.getpath())
    dgf.startchunk(dgf_write,"SNAPSHOT");

    # provide identical metadata chunk
    dgf.writemetadata(dgf_write,dgs_metadata);
    
    # Write channels in same order as original waveform (dgs_wfmdict is an ordereddict)
    for wfmname in dgs_wfmdict:
        dgf.writenamedwfm(dgf_write,new_wfmdict[wfmname])
        pass
    dgf.endchunk(dgf_write) # SNAPSHOT
    dgf.close(dgf_write)
    
    # ipython interactive execution only works properly if the results
    # are returned at the very bottom of the function

    return {"dc:dgs_undistorted": undistortedhref}

