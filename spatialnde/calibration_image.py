import numpy as np
import cv2
import sys

import os
import os.path

# import dg_file as dgf   # Now imported when used
from .timestamp import roundtosecond,now

from limatix import xmldoc
from limatix import dc_value

#sys.path.insert(0, "spatialnde/")
#from calibration_image_patterns import AsymCircPattern
#from calibration_image_patterns import AsymCircPattern2
#from calibration_image_patterns import ChessBoardPattern
#from calibration_image_patterns import SymCircPattern

from .calibration_image_patterns import AsymCircPattern
from .calibration_image_patterns import AsymCircPattern2
from .calibration_image_patterns import ChessBoardPattern
from .calibration_image_patterns import SymCircPattern


class calibration_image_params(dc_value.value):
    # Represents input parameters to a distortion correction operation
    input_href=None
    channel_name=None  # applies if inputfile is of a type that has channels
    frame_number=None # applies if inputfile is of a type that has frames
                      # NOTE: for .dgs files frame_number is zero based
                      #       as in direct waveform access, not one-based
                      #       as seen in the scope
    offset=None  # subtract this offset from raw image data prior to analysis. It represents what will ultimately be the middle gray (signal value=0.5)
    units_per_intensity=None  # divide out this scale factor from offset-subtracted raw image data prior to analysis
    invert=None  # Should we invert image prior to analysis?
    cvflags=None     # cvflags are flags to findCirclesGrid or findChessboardCorners
    # i.e. cv2.CALIB_CB_CLUSTERING (for symmetric or asymmetric circles)
    # cv2.CALIB_CB_ADAPTIVE_THRESHOLD (for chessboard)
    # cv2.CALIB_CB_NORMALIZE_IMAGE (for chessboard)
    # cv2.CALIB_CB_FILTER_QUADS (for chessboard)
    # cv2.CALIB_CB_FAST_CHECK (for chessboard, unlikely)

    pattern_type=None
    pattern_num_x=None
    pattern_num_y=None
    pattern_spacing=None
    

    nsmap={

        "snci": "http://thermal.cnde.iastate.edu/spatialnde/calibrationimage",
        "xlink": "http://www.w3.org/1999/xlink",
        "dcv": "http://limatix.org/dcvalue",
        }

    
    def __init__(self,input_href,channel_name,frame_number,offset,units_per_intensity,invert,cvflags,pattern_type,pattern_num_x,pattern_num_y,pattern_spacing):
        assert(isinstance(input_href,dc_value.hrefvalue))
        self.input_href=input_href
        

        self.channel_name=channel_name
        self.frame_number=frame_number
        self.offset=offset
        self.units_per_intensity=units_per_intensity
        self.invert=invert
        self.cvflags=cvflags
        self.pattern_type=pattern_type
        self.pattern_num_x=pattern_num_x
        self.pattern_num_y=pattern_num_y
        assert(isinstance(pattern_spacing,dc_value.numericunitsvalue))
        self.pattern_spacing=pattern_spacing

        # object is considered final once created. 
        self.final=True
        pass

    def xmlrepr(self,xmldocu,element,defunits=None):
        inputel=xmldocu.addelement(element,"snci:input")
        self.input_href.xmlrepr(xmldocu,inputel)

        if self.channel_name is not None:
            xmldocu.addsimpleelement(element,"snci:channel_name",(self.channel_name,))
            pass
        
        if self.frame_number is not None:
            xmldocu.addsimpleelement(element,"snci:frame_number",(str(self.frame_number),))
            pass
        
        
        xmldocu.addsimpleelement(element,"snci:offset",(repr(self.offset),))
        xmldocu.addsimpleelement(element,"snci:units_per_intensity",(repr(self.units_per_intensity),))
        xmldocu.addsimpleelement(element,"snci:invert",(str(self.invert),))
        if self.cvflags is not None and self.cvflags != 0:
            cvflags_el=xmldocu.addelement(element,"snci:cvflags")
            cvflags_strs=[]
            workflags=self.cvflags
            if self.pattern_type=="chessboard":
                if workflags & cv2.CALIB_CB_ADAPTIVE_THRESH:
                    cvflags_strs.append("CALIB_CB_ADAPTIVE_THRESH")
                    workflags &= ~cv2.CALIB_CB_ADAPTIVE_THRESH
                    pass
                
                if workflags & cv2.CALIB_CB_NORMALIZE_IMAGE:
                    cvflags_strs.append("CALIB_CB_NORMALIZE_IMAGE")
                    workflags &= ~cv2.CALIB_CB_NORMALIZE_IMAGE
                    pass

                if workflags & cv2.CALIB_CB_FILTER_QUADS:
                    cvflags_strs.append("CALIB_CB_FILTER_QUADS")
                    workflags &= ~cv2.CALIB_CB_FILTER_QUADS
                    pass

                if workflags & cv2.CALIB_CB_FAST_CHECK:
                    cvflags_strs.append("CALIB_CB_FAST_CHECK")
                    workflags &= ~cv2.CALIB_CB_FAST_CHECK
                    pass
                pass
            else: 
                # circles
                if workflags & cv2.CALIB_CB_CLUSTERING:
                    cvflags_strs.append("CALIB_CB_CLUSTERING")
                    workflags &= ~cv2.CALIB_CB_CLUSTERING
                    pass
                pass
            if workflags != 0:
                raise ValueError("Unknown CV flags: 0x%x" % (workflags))

            xmldocu.settext(cvflags_el,"|".join(cvflags_strs))
            pass
        xmldocu.addsimpleelement(element,"snci:pattern_type",(self.pattern_type,))
        xmldocu.addsimpleelement(element,"snci:pattern_num_x",(str(self.pattern_num_x),))
        xmldocu.addsimpleelement(element,"snci:pattern_num_y",(str(self.pattern_num_y),))

        pattern_spacingel=xmldocu.addelement(element,"snci:pattern_spacing")
        self.pattern_spacing.xmlrepr(xmldocu,pattern_spacingel,defunits=str(self.pattern_spacing.units()))
        
        pass

    @classmethod
    def fromxml(cls,xmldocu,element,defunits=None,xml_attribute=None):
        input_href=None
        channel_name=None
        frame_number=None
        offset=None
        units_per_intensity=None
        invert=None
        cvflags=None
        pattern_type=None
        pattern_spacing=None

        for child in xmldocu.children(element):
            if xmldocu.tag_is(child,"snci:input"):
                input_href=dc_value.hrefvalue.fromxml(xmldocu,child)
                pass
            elif xmldocu.tag_is(child,"snci:channel_name"):
                channel_name=xmldocu.gettext(child)
                pass
            elif xmldocu.tag_is(child,"snci:frame_number"):
                frame_number=int(xmldocu.gettext(child))
                pass
            elif xmldocu.tag_is(child,"snci:offset"):
                offset=float(xmldocu.gettext(child))
                pass
            elif xmldocu.tag_is(child,"snci:units_per_intensity"):
                units_per_intensity=float(xmldocu.gettext(child))
                pass
            elif xmldocu.tag_is(child,"snci:cvflags"):
                cvflags=0
                for flagname in xmldocu.gettext(child).split('|'):
                    cvflags |= getattr(cv2,flagname)
                    pass
                pass
            elif xmldocu.tag_is(child,"snci:invert"):
                invert=xmldocu.gettext(child).lower() in ("true","yes")
                pass
            elif xmldocu.tag_is(child,"snci:pattern_type"):
                pattern_type=xmldocu.gettext(child)
                pass
            elif xmldocu.tag_is(child,"snci:pattern_num_x"):
                pattern_num_x=int(xmldocu.gettext(child))
                pass
            elif xmldocu.tag_is(child,"snci:pattern_num_y"):
                pattern_num_y=int(xmldocu.gettext(child))
                pass
            elif xmldocu.tag_is(child,"snci:pattern_spacing"):
                pattern_spacing=dc_value.numericunitsvalue.fromxml(xmldocu,child)

                pass
            else:
                raise ValueError("Unknown distortion correction parameter: %s" % (xmldocu.gettag(child)))
            pass

        return cls(input_href,channel_name,frame_number,offset,units_per_intensity,invert,cvflags,pattern_type,pattern_num_x,pattern_num_y,pattern_spacing)
        pass
            
    def copy_and_update(self,**kwargs):
        paramdict={
            "input_href": self.input_href,
            "channel_name": self.channel_name,
            "frame_number": self.frame_number,
            "offset": self.offset,
            "units_per_intensity": self.units_per_intensity,
            "invert": self.invert,
            "cvflags": self.cvflags,
            "pattern_type": self.pattern_type,
            "pattern_num_x": self.pattern_num_x,
            "pattern_num_y": self.pattern_num_y,
            "pattern_spacing": self.pattern_spacing
            }
        
        for arg in kwargs.keys():
            assert(arg in paramdict)
            pass
        paramdict.update(kwargs)
        return calibration_image_params(**paramdict)
    
    def __str__(self):
        tempdoc=xmldoc.xmldoc.newdoc("snci:calibration_image_params",nsmap=self.nsmap)
        self.xmlrepr(tempdoc,tempdoc.getroot())
        return tempdoc.tostring(pretty_print=True)

    def __repr__(self):
        return self.__str__()
    
    pass

class calibration_image(dc_value.value):
    # Represents input parameters and results of identifying
    # points on a calibration image
    params=None   # class calibration_image_params
    object_points=None # Numpy array of object points
    image_points=None # Numpy array of image points
    timestamp=None # when objectpoints/image points were generated (ISO-8601 from spatialnde.timestamp)

    numpy_image=None  # not saved; only available when we have just done the point identification
    uint_image=None  # not saved; only available when we have just done the point identification


    nsmap={

        "snci": "http://thermal.cnde.iastate.edu/spatialnde/calibrationimage",
        "snic": "hhtp://thermal.cnde.iastate.edu/spatialnde/intrinsiccalibration",
        "xlink": "http://www.w3.org/1999/xlink",
        "dcv": "http://limatix.org/dcvalue",

        }

    def __init__(self,params,object_points=None,image_points=None,timestamp=None):
        # if object_points and image_points are None
        # will run opencv object/image point measurement based on the 
        # input params

        assert(isinstance(params,calibration_image_params))
        if object_points is None and image_points is None:
            (self.numpy_image,self.uint_image)=load_and_scale_calib_image(params)
            try:
                (object_points,image_points)=evaluate_calibration_image(params,self.uint_image,invert=params.invert,cvflags=params.cvflags)
                pass
            except ValueError:
                # if we couldn't evaluate, object_points and image_points will be None
                pass

            # save timestamp
            timestamp=roundtosecond(now()).isoformat()
            
            pass
        elif object_points is None and image_points is not None:
            (self.numpy_image,self.uint_image)=load_and_scale_calib_image(params)
            try:
                (object_points,image_points)=evaluate_calibration_image(params,self.uint_image,invert=params.invert,cvflags=params.cvflags,image_points=image_points)
                pass
            except ValueError:
                # if we couldn't evaluate, object_points and image_points will be None
                pass

            # save timestamp
            timestamp=roundtosecond(now()).isoformat()
            
            pass
        
        self.params=params
        self.object_points=object_points
        self.image_points=image_points
        self.timestamp=timestamp
        pass

    def xmlrepr(self,xmldocu,element,defunits=None):

        # Distortion correction contains params
        paramsel=xmldocu.addelement(element,"snci:calibration_image_params")
        self.params.xmlrepr(xmldocu,paramsel)
        
        # ... followed by object points array
        if self.object_points is not None:
            objptel=xmldocu.addelement(element,"snci:object_points")
            dc_value.arrayvalue(self.object_points).xmlrepr(xmldocu,objptel)
            pass

        # ... followed by image points array
        if self.image_points is not None:
            camptsel=xmldocu.addelement(element,"snci:image_points")
            dc_value.arrayvalue(self.image_points).xmlrepr(xmldocu,camptsel)
            pass
        
        if self.timestamp is not None:
            xmldocu.addsimpleelement(element,"snci:image_points_timestamp",(self.timestamp,))
            pass
        

        pass

    @classmethod
    def fromxml(cls,xmldocu,element,defunits=None,xml_attribute=None):
        params=None
        object_points=None
        image_points=None
        timestamp=None

        for child in xmldocu.children(element):
            if xmldocu.tag_is(child,"snci:calibration_image_params"):
                params=calibration_image_params.fromxml(xmldocu,child)
                pass
            elif xmldocu.tag_is(child,"snci:object_points"):
                object_points_value=dc_value.arrayvalue.fromxml(xmldocu,child)
                object_points=object_points_value.value().astype(dtype=np.float32) # Line 3168 in calib3d/src/calibration.cpp, objectPoints must be CV_32F 
                pass
            elif xmldocu.tag_is(child,"snci:image_points"):
                image_points_value=dc_value.arrayvalue.fromxml(xmldocu,child)
                image_points=image_points_value.value().astype(dtype=np.float32) # Line 3192 in calib3d/src/calibration.cpp, imagePoints must be CV_32F  
                pass
            elif xmldocu.tag_is(child,"snci:image_points_timestamp"):
                timestamp=xmldocu.gettext(child)
                pass
            pass
            
        return cls(params,object_points,image_points,timestamp) 


    def __str__(self):
        tempdoc=xmldoc.xmldoc.newdoc("snci:calibration_image",nsmap=self.nsmap)
        self.xmlrepr(tempdoc,tempdoc.getroot())
        return tempdoc.tostring(pretty_print=True)

    def __repr__(self):
        return self.__str__()
    


    pass

def extract_image(filename, channel_name=None,frame_number=None):

    """Extract an image from the provided file (.dgs, .png, .jpg, etc.) and channel/frame specification. Return single precision floating point numpy array."""

    inputfilepath=filename
    inputfilename=os.path.split(inputfilepath)[1]
    inputfileext=os.path.splitext(inputfilename)[1]
    
    if inputfileext.lower()==".dgs":
        import dg_file as dgf
        
        dgfh=dgf.open(inputfilepath)
        chunk = dgf.nextchunk(dgfh)
        chunk_metadata, wave_forms, wave_form_dictionary = dgf.procSNAPSHOT(dgfh,memmapok=True)
        wfm = wave_form_dictionary[channel_name]
        if len(wfm.data.shape) > 2:
            retval=wfm.data[:,::-1,frame_number].T  # reverse y order and transpose for consistency with non-dataguzzler image storage conventions
            pass
        else: 
            retval=wfm.data[:,::-1].T # reverse y order and transpose for consistency with non-dataguzzler image storage conventions
            pass
        dgf.close(dgfh)
        pass
    elif inputfileext.lower()==".png":
        import scipy.misc
        retval=np.array(scipy.misc.imread(inputfilepath, flatten=1),dtype='f') 
        pass
    else:
        raise ValueError("Unknown file type %s" % (inputfileext.lower()))
    return retval


def float_to_uint(im, immin=None, immax=None, bits=8):
    
    """Convert an image in floating point to an 8 or 16 bit integer image"""

    assert(im.dtype == np.float64 or im.dtype == np.float32)

    if immin is None:
        immin = np.nanmin(im)
        pass

    if immax is None:
        immax = np.nanmax(im)
        pass

    if bits == 8:
        upper_bound = 255.0
        lower_bound = 1.0
        im_transformed=np.uint8((im-immin)/(immax-immin)*(upper_bound-lower_bound) + lower_bound) ## 8-bit, range from 1 to 256
        im_transformed[im < immin]=lower_bound
        im_transformed[im > immax]=upper_bound
        im_transformed[np.isnan(im)] = 0

    elif bits == 16:
        upper_bound = 65535.0
        lower_bound = 10.0
        im_transformed=np.uint16((im-immin)/(immax-immin)*(upper_bound-lower_bound) + lower_bound) ## 16-bit, range from 10 to 65536
        im_transformed[im < immin]=lower_bound
        im_transformed[im > immax]=upper_bound
        im_transformed[np.isnan(im)] = 0
    
    return im_transformed


def uint_to_float(im, immin, immax, bits=8,dtype=np.dtype('float32')):

    """Convert a 8 or 16 bit image into a floating point image. Invalid pixels are stored as nan"""

    if bits == 8:
        assert(im.dtype == np.uint8)
        upper_bound = 255.0
        lower_bound = 1.0
    elif bits == 16:
        assert(im.dtype == np.uint16)
        upper_bound = 65535.0
        lower_bound = 10.0
    
    im_transformed = (np.array(im,dtype=dtype)-lower_bound)/(upper_bound-lower_bound)*(immax-immin)+immin

    return im_transformed


def normalize_float(im, units_per_intensity, offset):
    normalized_float = (im-offset)/units_per_intensity + 0.5

    return normalized_float


def load_and_scale_calib_image(params):
    numpy_image = extract_image(params.input_href.getpath(),params.channel_name,params.frame_number)

    scaled_image = normalize_float(numpy_image, params.units_per_intensity, params.offset)
    uint_image = float_to_uint(scaled_image, immin=0.0, immax=1.0, bits=8) # must use 8-bit because opencv blob detector doesn't support anything else
    

    return (numpy_image,uint_image)

def evaluate_calibration_image(params,uint_image,invert=False,image_points=None,cvflags=0):
    # cvflags are flags to findCirclesGrid or findChessboardCorners
    # i.e. cv2.CALIB_CB_CLUSTERING (for symmetric or asymmetric circles)
    # cv2.CALIB_CB_ADAPTIVE_THRESHOLD (for chessboard)
    # cv2.CALIB_CB_NORMALIZE_IMAGE (for chessboard)
    # cv2.CALIB_CB_FILTER_QUADS (for chessboard)
    # cv2.CALIB_CB_FAST_CHECK (for chessboard, unlikely)

    assert(uint_image.dtype == np.uint8)

    if cvflags is None: 
        cvflags=0
        pass

    if params.pattern_type=="asymmetric_circles":
        object_points=AsymCircPattern((params.pattern_num_x,params.pattern_num_y),params.pattern_spacing.value())
        findcirclesflags=cv2.CALIB_CB_ASYMMETRIC_GRID
        pass
    elif params.pattern_type=="asymmetric_circles_offset":
        object_points=AsymCircPattern2((params.pattern_num_x,params.pattern_num_y),params.pattern_spacing.value())
        findcirclesflags=cv2.CALIB_CB_ASYMMETRIC_GRID
        pass
    elif params.pattern_type=="chessboard":
        object_points=ChessBoardPattern((params.pattern_num_x,params.pattern_num_y),params.pattern_spacing.value())
        pass
    elif params.pattern_type=="symmetric_circles":
        object_points=SymCircPattern((params.pattern_num_x,params.pattern_num_y),params.pattern_spacing.value())
        findcirclesflags=cv2.CALIB_CB_SYMMETRIC_GRID
        pass
    else :
        raise ValueError("Unknown pattern type")
    

    #cv2.imshow('numpy_image',numpy_image)
    #cv2.imshow('uint_image',uint_image)
    #cv2.waitKey(0);
    #import pylab as pl
    #pl.imshow(uint_image,vmin=0,vmax=255)
    #pl.show()

    if invert:
        # invert is useful for chessboard because per documentation
        # you need a white border around the chessboard. So if it is
        # black, just invert it
        uint_image=255-uint_image
        pass

    if image_points is None:
        if params.pattern_type=="chessboard":
            (found, image_points) = cv2.findChessboardCorners(uint_image, (params.pattern_num_x,params.pattern_num_y),flags=cvflags)
            pass
        else:
            (found, image_points) = cv2.findCirclesGrid(uint_image, (params.pattern_num_x,params.pattern_num_y), flags=findcirclesflags|cvflags)
            pass
    
        #import pylab as pl
        #print(image_points)
        ##image_points_np=np.array(image_points,dtype='d')
        ##print(image_points_np)
        ##print(image_points_np.shape)
        #pl.imshow(uint_image,vmin=0,vmax=255)
        #pl.plot(image_points[:,0,0],image_points[:,0,1],'o',markersize=10)
        #pl.show()

        if not found:
            raise ValueError("All points not found for %d by %d grid" % (params.pattern_num_x,params.pattern_num_y))
        pass
    
    return (object_points, image_points)
