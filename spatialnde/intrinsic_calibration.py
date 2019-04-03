import numpy as np
import cv2
import sys
import copy

import os

from .timestamp import roundtosecond,now

from limatix import xmldoc
from limatix import dc_value

#sys.path.insert(0, "spatialnde/")
#from calibration_image import calibration_image
#from calibration_image import float_to_uint
#from calibration_image import uint_to_float

from .calibration_image import calibration_image
from .calibration_image import float_to_uint
from .calibration_image import uint_to_float

class intrinsic_calibration(dc_value.value):
    # Represents input parameters and results of a camera intrinsic calibration operation
    calibration_images=None   # tuple of class calibration_image
    sizex=None  # horizontal size of image from camera, in pixels
    sizey=None  # vertical size of image from camera, in pixels
    cameramatrix=None  # Numpy
    distcoeffs=None  # Numpy
    _reprojection_error = None # float ... use .reprojection_error property instead 
    timestamp=None  # ISO-8601 from spatialnde.timestamp

    nsmap={
        "snic": "http://thermal.cnde.iastate.edu/spatialnde/intrinsiccalibration",
        "snci": "http://thermal.cnde.iastate.edu/spatialnde/calibrationimage",
        "xlink": "http://www.w3.org/1999/xlink",
        "dcv": "http://limatix.org/dcvalue",
        
        }

    def __init__(self,calibration_images,sizex,sizey,cameramatrix=None,distcoeffs=None,reprojection_error=None,timestamp=None):
        # if cameramatrix and distcoeffs are None
        # will run opencv camera intrinsic calibration measurement based on the 
        # input params

        for calibimage in calibration_images:
            assert(isinstance(calibimage,calibration_image))
            pass
        
        self.sizex=sizex
        self.sizey=sizey

        if cameramatrix is None and distcoeffs is None:
            (cameramatrix,distcoeffs)=evaluate_intrinsic_calibration(calibration_images,sizex,sizey)

            # save timestamp  ... via spatialnde.timestamp
            timestamp=roundtosecond(now()).isoformat()
            
            pass
        
        self.calibration_images=tuple(calibration_images)
        self.cameramatrix=cameramatrix
        self.distcoeffs=distcoeffs

        #self.reprojection_error=evaluate_reprojection_error(self.calibration_images, self.cameramatrix, self.distcoeffs)
        self.timestamp=timestamp
        pass

    @property  # Evaluate reprojection_error on demand by accessing the .reprojection_error property
    def reprojection_error(self):
        if self._reprojection_error is None:
            self._reprojection_error=evaluate_reprojection_error(self.calibration_images, self.cameramatrix, self.distcoeffs)
            pass
        
        return self._reprojection_error
        
    
    def xmlrepr(self,xmldocu,element,defunits=None):

        # Write sequence of snic:calibration_image tags
        for image in self.calibration_images:
            imageel=xmldocu.addelement(element,"snic:calibration_image")
            image.xmlrepr(xmldocu,imageel)
            pass
        
        if self.sizex is not None:
            xmldocu.addsimpleelement(element,"snic:sizex",(self.sizex,))
            pass

        if self.sizey is not None:
            xmldocu.addsimpleelement(element,"snic:sizey",(self.sizey,))
            pass

        # ... followed by camera matrix
        if self.cameramatrix is not None:
            cammtxel=xmldocu.addelement(element,"snic:cameramatrix")
            dc_value.arrayvalue(self.cameramatrix).xmlrepr(xmldocu,cammtxel)
            pass

        # ... followed by distortion coefficients
        if self.distcoeffs is not None:
            distcoeffsel=xmldocu.addelement(element,"snic:distcoeffs")
            dc_value.arrayvalue(self.distcoeffs).xmlrepr(xmldocu,distcoeffsel)
            pass
        
        if self.reprojection_error is not None:
            xmldocu.addsimpleelement(element,"snic:reprojection_error",(self.reprojection_error,))
            pass
        
        if self.timestamp is not None:
            xmldocu.addsimpleelement(element,"snic:calculationtimestamp",(self.timestamp,))
            pass


        pass

    @classmethod
    def fromxml(cls,xmldocu,element,defunits=None,xml_attribute=None):
        calibration_images=[]
        sizex=None
        sizey=None
        cameramatrix=None
        distcoeffs=None
        reprojection_error=None
        timestamp=None

        for child in xmldocu.children(element):
            if xmldocu.tag_is(child,"snic:calibration_image"):
                image=calibration_image.fromxml(xmldocu,child)
                calibration_images.append(image)
                pass
            elif xmldocu.tag_is(child,"snic:sizex"):
                sizex=int(xmldocu.gettext(child))
                pass
            elif xmldocu.tag_is(child,"snic:sizey"):
                sizey=int(xmldocu.gettext(child))
                pass
            elif xmldocu.tag_is(child,"snic:cameramatrix"):
                matrixvalue=dc_value.arrayvalue.fromxml(xmldocu,child)
                cameramatrix=matrixvalue.value()
                pass
            elif xmldocu.tag_is(child,"snic:distcoeffs"):
                distvalue=dc_value.arrayvalue.fromxml(xmldocu,child)
                distcoeffs=distvalue.value()
                pass
            elif xmldocu.tag_is(child,"snic:reprojection_error"):
                reprojection_error=float(xmldocu.gettext(child))
                pass
            elif xmldocu.tag_is(child,"snic:calculationtimestamp"):
                timestamp=xmldocu.gettext(child)
                pass
            else: 
                raise ValueError("Unknown tag in intrinsic calibration: %s" % (xmldocu.get_tag(child)))
            
            pass

            
        return cls(calibration_images,sizex,sizey,cameramatrix,distcoeffs,reprojection_error,timestamp)

    

    def value(self):
        return (copy.copy(self.cameramatrix),copy.copy(self.distcoeffs))

 
    @classmethod
    def fromfile(cls,filename):
        calibdoc=xmldoc.xmldoc.loadfile(filename,nsmap=cls.nsmap)
        
        calibration=cls.fromxml(calibdoc,calibdoc.getroot())
        
        return calibration

    def savetofile(self,filename):
        calibdoc=xmldoc.xmldoc.newdoc(maintagname="snic:intrinsic_calibration",nsmap=self.nsmap)
        self.xmlrepr(calibdoc,calibdoc.getroot())
        calibdoc.set_href(href=filename)
        calibdoc.close()
        pass
    pass



def evaluate_intrinsic_calibration(calibration_image_list,sizex,sizey):
    object_points_list=[]
    image_points_list=[]

    for calibration_image in calibration_image_list:
        object_points_list.append(calibration_image.object_points)
        image_points_list.append(calibration_image.image_points)
        pass
    
    (retval,cameramatrix,distcoeffs, rvecs, tvecs) = cv2.calibrateCamera(object_points_list, image_points_list, (sizex,sizey))

    return (cameramatrix,distcoeffs)


def evaluate_reprojection_error(calibration_image_list, cameramatrix, distcoeffs):
    # Re-projection error based on the camera calibration

    reprojection_error = 0
    for calibration_image in calibration_image_list:
        rvec, tvec, inliers = cv2.solvePnPRansac(calibration_image.object_points, calibration_image.image_points, cameramatrix, distcoeffs) 

        projected_image_points, jacobian = cv2.projectPoints(calibration_image.object_points, rvec, tvec, cameramatrix, distcoeffs)

        arr1 = calibration_image.image_points.reshape(-1,2)
        arr2 = projected_image_points.reshape(-1,2)
        err = np.sqrt(np.sum((arr1[:,0]-arr2[:,0])**2+(arr1[:,1]-arr2[:,1])**2)/arr1.shape[0]) # Root mean square error
        reprojection_error+=err

    # Average error over the calibration images (some may have different grid sizes with points in different locations)

    reprojection_error/=len(calibration_image_list)

    return reprojection_error


#def redistort_image(img, sizex, sizey, cameramatrix, distcoeffs):
#    u, v = np.meshgrid(np.arange(sizex), np.arange(sizey))
#    u = u.flatten()
#    v = v.flatten()
#    img_points = np.vstack([u, v])
#
#    k1, k2, p1, p2, k3 = distcoeffs[0]
#
#    pass




class undistortion_processor(object):
    calibration=None  # intrinsic calibration
    insizex=None
    insizey=None
    outsizex=None
    outsizey=None
    interpolation=None # cv2.INTER_NEAREST, cv2.INTER_LINEAR, cv2.INTER_CUBIC, or cv2.INTER_LANCZOS4

    mapx=None
    mapy=None
    
    newcameramatrix=None
    
    def __init__(self,calibration,insizex,insizey,outsizex,outsizey,interpolation,alpha):
        # Alpha parameter is supposed to control cropping of undistorted image
        # but seems to be broken for alpha!=1.0 (?)
        self.calibration=calibration
        self.insizex=insizex
        self.insizey=insizey
        self.outsizex=outsizex
        self.outsizey=outsizey
        self.interpolation=interpolation
        
        (self.newcameramatrix,retval)=cv2.getOptimalNewCameraMatrix(self.calibration.cameramatrix,self.calibration.distcoeffs,(insizex,insizey),alpha,(outsizex,outsizey))

        (self.mapx,self.mapy)=cv2.initUndistortRectifyMap(self.calibration.cameramatrix,self.calibration.distcoeffs,np.identity(3),self.newcameramatrix,(outsizex,outsizey),cv2.CV_16SC2)#cv2.CV_32FC1)
        
        pass
    def undistort_numpy_image(self,numpy_image):
        
        immax = np.nanmax(numpy_image)
        immin = np.nanmin(numpy_image)

        integer_image=float_to_uint(numpy_image,immin,immax,bits=16)
        integer_image_undistorted=cv2.remap(integer_image,self.mapx,self.mapy,self.interpolation)
        image_undistorted=uint_to_float(integer_image_undistorted, immin, immax, bits=16,dtype=numpy_image.dtype)
        image_undistorted[integer_image_undistorted==0] = np.nan

        return image_undistorted
    pass


