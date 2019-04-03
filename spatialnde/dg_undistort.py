# dg_undistort.py: Remove distortion from dataguzzler waveforms
#
# ... Should probably be moved into dataguzzler subdirectory

import sys
import copy
import collections
import dataguzzler as dg
import numpy as np


def undistort_wfm_given_emptystruct(processor,wfm,emptystruct):
    # emptystruct should have data,dimlen, and ndim set

    emptystruct.Name=wfm.Name
    emptystruct.MetaData=copy.deepcopy(wfm.MetaData)
    emptystruct.wfmrevision=wfm.wfmrevision
    emptystruct.IsProcExpr=wfm.IsProcExpr
    emptystruct.IsProcRGBA=wfm.IsProcRGBA
    emptystruct.NeedData=wfm.NeedData
    emptystruct.NeedMetaData=wfm.NeedMetaData
    emptystruct.HaveData=wfm.HaveData
    emptystruct.HaveMetaData=wfm.HaveMetaData

    # input shape must match processor
    assert(wfm.data.shape[0]==processor.insizex)
    assert(wfm.data.shape[1]==processor.insizey)

    if wfm.ndim==2:
        toundistort=wfm.data[:,::-1].T
        undistorted_image=processor.undistort_numpy_image(toundistort)
        emptystruct.data[::]=undistorted_image[::-1,:].T
        pass
    else:
        # wfm.ndim==3
        for framecnt in range(wfm.data.shape[2]):
            # Correct for dataguzzler image origins in lower left corner
            toundistort=wfm.data[:,::-1,framecnt].T
            undistorted_image=processor.undistort_numpy_image(toundistort)
            emptystruct.data[:,:,framecnt]=undistorted_image[::-1,:].T
            pass
        pass
    
    return emptystruct

# undistort a single waveform, given a processor
def undistort_wfm(processor,wfm):
    newwfm=dg.wfminfo()
    
    assert(wfm.ndim >= 2 and wfm.ndim <= 3)

    newwfm.ndim=wfm.ndim
    if wfm.ndim==2:
        newwfm.dimlen=np.array((processor.outsizex,processor.outsizey),dtype=np.uint64)
        newwfm.data=np.empty((processor.outsizex,processor.outsizey),dtype=wfm.data.dtype)
        pass
    else:
        # wfm.ndim==3

        newwfm.dimlen=np.array((processor.outsizex,processor.outsizey,wfm.data.shape[2]),dtype=np.uint64)
        newwfm.data=np.empty((processor.outsizex,processor.outsizey,wfm.data.shape[2]),dtype=wfm.data.dtype)
        pass
    return undistort_wfm_given_emptystruct(processor,wfm,newwfm)



def undistort_wfmdict(processor,wfmdict):
    # undistort all channels in a waveform dict that are
    # (a) 2 or 3-dimensional
    # (b) first two axis lengths match processor.insizex and processor.insizey
    # Return new dictionary with all channels:
    #   * unprocessed channels will be references to input dictionary entries
    #   * processed channels will be new waveforms
    # If wfmdict is ordered, order will be preserved. 
    newwfmdict=collections.OrderedDict()
    
    for wfmname in wfmdict:
        if ((wfmdict[wfmname].ndim==2 or wfmdict[wfmname].ndim==3) and
            wfmdict[wfmname].data.shape[0]==processor.insizex and
            wfmdict[wfmname].data.shape[1]==processor.insizey):
            newwfmdict[wfmname]=undistort_wfm(processor,wfmdict[wfmname])
            pass
        else:
            newwfmdict[wfmname]=wfmdict[wfmname]
            pass
        pass
    return newwfmdict

