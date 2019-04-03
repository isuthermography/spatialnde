

import sys, os
import copy
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


import dataguzzler as dg

from spatialnde import intrinsic_calibration
from spatialnde import dg_undistort
from spatialnde.calibration_image import uint_to_float
from spatialnde.calibration_image import float_to_uint

from matplotlib import pyplot as plt

import dg_file as dgf
import dg_metadata as dgm
import dg_eval

from limatix.dc_value import hrefvalue
from limatix.dc_value import numericunitsvalue
from limatix.dc_value import arrayvalue


def numeric_input(prompt,default=None):
    defaultstr=""
    if default is not None:
        defaultstr=" [%s]" % (default)
        pass
    while True:
        try: 
            inp=raw_input("%s:%s " % (prompt,defaultstr)).strip()
            if len(inp)==0 and default is not None:
                inp=default
                pass
            value=numericunitsvalue(inp)
            pass
        except KeyboardInterrupt:
            raise
        except Exception as e: 
            print(e)
            print("Try again... (or press Ctrl-C to break)")
            pass
        else:
            return value
        pass
    pass



class PointGetter(object):
    xdata=None
    ydata=None

    cid=None  # Mysterious return value from mpl_connect. It seems to be an integer (?)

    def __init__(self,fig):
        self.cid = fig.canvas.mpl_connect('button_press_event',self.got_button_press)
        self.xdata=None
        self.ydata=None
        self.button=None
        self.markerplot=None
        pass

    def got_button_press(self,event):
        self.xdata=event.xdata
        self.ydata=event.ydata
        self.button=event.button
        if self.markerplot is not None:
            self.markerplot.remove()
        self.markerplot=plt.plot(self.xdata,self.ydata,'x',markersize=15,markeredgewidth=3)
        pass
    pass


def selection_input(optiondict,prompt,default=None):
    defaultstr=""
    if default is not None:
        defaultstr=" [%s]" % (default)
        pass

    while True:
        try: 
            print(prompt)
            for key in optiondict: 
                print("  %s: %s" % (key,optiondict[key]))
                pass
            inp=raw_input("   -->%s " % (defaultstr)).strip()
            if len(inp)==0 and default is not None:
                inp=default
                pass
            if inp not in optiondict: 
                raise ValueError("Invalid entry: %s" % (inp))
            pass
        except KeyboardInterrupt:
            raise
        except Exception as e: 
            print(e)
            print("Try again... (or press Ctrl-C to break)")
            pass
        else:
            return (inp,optiondict[inp])
        pass
    pass
        
def position_input(imagedata,prompt):
    fig=plt.figure()
    plt.imshow(imagedata)
    plt.title(prompt)
    getter=PointGetter(fig)
    plt.show()
    plt.close(fig)
    return (getter.xdata,getter.ydata)


# perspective single waveform, given a transformation matrix
def perspective_wfm(transformation_matrix,pixelsizex,pixelsizey,new_x_dim,new_y_dim,wfm,wfmgeom):
    (ndim,DimLen,IniVal,Step,bases)=wfmgeom

    assert(wfm.ndim >= 2 and wfm.ndim <= 3)
    
    newwfm=dg.wfminfo()
    newwfm.Name=wfm.Name
    newwfm.MetaData=copy.deepcopy(wfm.MetaData)
    newwfm.wfmrevision=wfm.wfmrevision
    newwfm.ndim=wfm.ndim
    newwfm.IsProcExpr=wfm.IsProcExpr
    newwfm.IsProcRGBA=wfm.IsProcRGBA
    newwfm.NeedData=wfm.NeedData
    newwfm.NeedMetaData=wfm.NeedMetaData
    newwfm.HaveData=wfm.HaveData
    newwfm.HaveMetaData=wfm.HaveMetaData

    # input shape must match processor
    #assert(wfm.data.shape[0]==processor.insizex)
    #assert(wfm.data.shape[1]==processor.insizey)
    NewStep1=pixelsizex.value('m')
    NewStep2=pixelsizey.value('m')
    # Origin is treated as the (minx, miny) image corner. IniVal represents pixel center
    NewIniVal1=pixelsizex.value('m')/2.0
    NewIniVal2=pixelsizey.value('m')/2.0
    dgm.AddMetaDatumWI(newwfm,dgm.MetaDatum("IniVal1",NewIniVal1))
    dgm.AddMetaDatumWI(newwfm,dgm.MetaDatum("IniVal2",NewIniVal2))
    dgm.AddMetaDatumWI(newwfm,dgm.MetaDatum("Step1",NewStep1))
    dgm.AddMetaDatumWI(newwfm,dgm.MetaDatum("Step2",NewStep2))
    dgm.AddMetaDatumWI(newwfm,dgm.MetaDatum("Units1","meters"))
    dgm.AddMetaDatumWI(newwfm,dgm.MetaDatum("Units2","meters"))
    

    if wfm.ndim==2:
        toperspective=wfm.data[:,::-1].T
        #immin=np.nanmin(toperspective)
        #immax=np.nanmax(toperspective)

        #toperspective_uint=float_to_uint(toperspective,immin=immin,immax=immax,bits=16)
        #orthographic_image=cv2.warpPerspective(toperspective_uint,transformation_matrix, (new_x_dim,new_y_dim))
        orthographic_image=cv2.warpPerspective(toperspective,transformation_matrix, (new_x_dim,new_y_dim))
        #plt.imshow(orthographic_image)
        #plt.colorbar()
        #plt.show()

        #newwfm.data=uint_to_float(orthographic_image[::-1,:].T,immin,immax,bits=16,dtype=wfm.data.dtype)
        newwfm.data=orthographic_image[::-1,:].T
        pass
    else:
        # wfm.ndim==3
        newwfm.data=np.zeros((new_x_dim,new_y_dim,wfm.data.shape[2]),dtype=wfm.data.dtype)
        for framecnt in range(wfm.data.shape[2]):
            # Correct for dataguzzler image origins in lower left corner
            toperspective=wfm.data[:,::-1,framecnt].T
            orthographic_image=cv2.warpPerspective(toperspective,transformation_matrix, (new_x_dim,new_y_dim))
            newwfm.data[:,:,framecnt]=orthographic_image[::-1,:].T
            pass
        pass
    newwfm.dimlen=np.array(newwfm.data.shape,dtype=np.uint64)
    
    return newwfm


def perspective_wfmdict(transformation_matrix,pixelsizex,pixelsizey,new_x_dim,new_y_dim,wfmdict,wfmgeom):
    # perspective all channels in a waveform dict that are
    # (a) 2 or 3-dimensional
    # (b) Match geometry specified in wfmgeom (result of dg_eval.geom())
    # Return new dictionary with all channels:
    #   * unprocessed channels will be references to input dictionary entries
    #   * processed channels will be new waveforms
    newwfmdict={}
    (ndim,DimLen,IniVal,Step,bases)=wfmgeom

    for wfmname in wfmdict:
        if ((wfmdict[wfmname].ndim==2 or wfmdict[wfmname].ndim==3) and
            (wfmdict[wfmname].data.shape[0]==DimLen[0]) and
            (wfmdict[wfmname].data.shape[1]==DimLen[1]) and
            (wfmdict[wfmname].MetaData["IniVal1"].Value==IniVal[0]) and
            (wfmdict[wfmname].MetaData["Step1"].Value==Step[0]) and
            (wfmdict[wfmname].MetaData["IniVal2"].Value==IniVal[1]) and
            (wfmdict[wfmname].MetaData["Step2"].Value==Step[1])):
            
            newwfmdict[wfmname]=perspective_wfm(transformation_matrix,pixelsizex,pixelsizey,new_x_dim,new_y_dim,wfmdict[wfmname],wfmgeom)
            pass
        else:
            newwfmdict[wfmname]=wfmdict[wfmname]
            pass
        pass
    return newwfmdict

#!!! TODO: Encode IniVal1, Step1, IniVal2, Step2 with physical coordinates

def run(_xmldoc, _tag, channel_name_str, frame_number_int, pixelsizex_numericunits, pixelsizey_numericunits, dc_dgs_undistorted_href):

    # dc_dgs_undistorted is dgs_undistorted element under measurement element

    dc_dgs_undistorted_path = dc_dgs_undistorted_href.getpath()


    # load in .dgs file
    dgfh=dgf.open(dc_dgs_undistorted_path)
    chunk = dgf.nextchunk(dgfh)
    chunk_metadata, wave_forms, wfmdict = dgf.procSNAPSHOT(dgfh,memmapok=True)
    
    print("")

    xlen=numeric_input("Enter specimen x-axis length and units --> ")
    
    # print("Got xlen: %f mm" % (xlen.value("mm")))

    ylen=numeric_input("Enter specimen y-axis length and units --> ")
    # print("Got ylen: %f mm" % (ylen.value("mm")))

    #face_viewed=selection_input({"f":"front","b":"back"},"Enter face viewed",default="f")

    #zpos=numeric_input("Enter z position viewed and units --> ",default=numericunitsvalue("0 m"))

    gotimage=wfmdict[channel_name_str].data[:,::-1,int(frame_number_int)].T
    #origin_pixels=position_input(gotimage,"Select origin")
    pos1_pixels=position_input(gotimage,"select point at minimum x, minimum y")
    print(" minimum x, minimum y = (%f,%f)" % (pos1_pixels[0],pos1_pixels[1]))

    pos2_pixels=position_input(gotimage,"select point at minimum x, maximum y")
    print(" minimum x, maximum y = (%f,%f)" % (pos2_pixels[0],pos2_pixels[1]))
    
    pos3_pixels=position_input(gotimage,"select point at maximum x, maximum y")
    print(" maximum x, maximum y = (%f,%f)" % (pos3_pixels[0],pos3_pixels[1]))

    pos4_pixels=position_input(gotimage,"select point at maximum x, minimum y")
    print(" maximum x, minimum y = (%f,%f)" % (pos4_pixels[0],pos4_pixels[1]))
    

    new_x_dim=int(xlen.value('m')/pixelsizex_numericunits.value('m'))
    new_y_dim=int(ylen.value('m')/pixelsizey_numericunits.value('m'))

    pts1=np.array((pos1_pixels,pos2_pixels,pos3_pixels,pos4_pixels),dtype='d')
    pts2=np.array(( (0.0,0.0), (0.0,new_y_dim-1), (new_x_dim-1,new_y_dim-1), (new_x_dim-1,0.0)))

    # Reverse sense of pts2 so our transform doesn't give the image an extra flip, considering minimum x is to the lower left
    pts2[:,1]=new_y_dim-1-pts2[:,1]
    
    #vecplusy=pts1[1,:]-pts1[0,:]
    #vecplusx=pts1[3,:]-pts1[0,:]
    #origindist=np.sqrt(np.sum((pts1-origin_pixels)**2.0,1))
    #originidx=np.argmin(origindist)

    #print(origin_pixels)
    #print(origindist)
    #print(originidx)
    #originstrings=("minxminy","minxmaxy","maxxmaxy","maxxminy")

    #if (origindist[originidx] > 30):
    #    print("WARNING... CLICKED ORIGIN DOES NOT SEEM TO LINE UP")
    #    print("WITH A CORNER. USING NEAREST CORNER ANYWAY")
    #    pass

    #pts2=pts2-originshift[originidx]


    transformation_matrix = cv2.getPerspectiveTransform(np.float32(pts1),np.float32(pts2))
    
    
    
    #dst = cv2.warpPerspective(im_undistorted,transformation_matrix,(xlen.value('m')/pixelsizex.value('m'),ylen.value('m')/pixelsizey.value('m')))


    new_wfmdict=perspective_wfmdict(transformation_matrix,pixelsizex_numericunits,pixelsizey_numericunits,new_x_dim,new_y_dim,wfmdict, dg_eval.geom(wfmdict[channel_name_str],raw=True))

    orthographicdgs="%s_orthographic.dgs" % (posixpath.splitext(dc_dgs_undistorted_href.get_bare_unquoted_filename())[0])
    orthographichref=hrefvalue(quote(orthographicdgs),contexthref=dc_dgs_undistorted_href)
    
    dgf_write=dgf.creat(orthographichref.getpath())
    dgf.startchunk(dgf_write,"SNAPSHOT");

    # provide empty metadata chunk
    EmptyMetadata={};
    dgf.writemetadata(dgf_write,EmptyMetadata);
    
    # Write channels in same order as original waveform
    for origwfm in wave_forms:
        dgf.writenamedwfm(dgf_write,new_wfmdict[origwfm.Name])
        pass
    dgf.endchunk(dgf_write) # SNAPSHOT
    dgf.close(dgf_write)
    
    dgf.close(dgfh)
    # ipython interactive execution only works properly if the results
    # are returned at the very bottom of the function

    return {"dc:dgs_orthographic": orthographichref, 
            "dc:dgs_orthographic_xlen": xlen,
            "dc:dgs_orthographic_ylen": ylen,
            "dc:dgs_orthographic_cornerpts_imagepixels_origin_incr_y_then_x": arrayvalue(pts1), 
            "dc:dgs_orthographic_pixelsizex": pixelsizex_numericunits,
            "dc:dgs_orthographic_pixelsizey": pixelsizey_numericunits,}

