import sys 

import numpy as np
import cv2


from . import intrinsic_calibration

#import dg_eval # now imported when used
#import dg_metadata as dgm # now imported when used

class ImageProjectionModel(object):
    calibration=None   # spatialnde.intrinsic_calibration.intrinsic_calibration
    #new_camera_mtx_cache=None # Dictionary by alpha, image_size_x, image_size_y of new camera matrix
    new_camera_mtx=None
    undistortion_already_applied=None

    #def get_new_camera_mtx(self,alpha,image_size_x,image_size_y):
    #    if (alpha,image_size_x,image_size_y) in self.new_camera_mtx_cache:
    #        return self.new_camera_mtx_cache[(alpha,image_size_x,image_size_y)]
    #
    #    (new_camera_mtx,validpixroi) = cv2.getOptimalNewCameraMatrix(cameraMatrix=self.calibration.cameramatrix, distCoeffs=self.calibration.distcoeffs,imageSize=(self.calibration.sizex, self.calibration.sizey), alpha=alpha, newImgSize=(image_size_x,image_size_y))
    #    self.new_camera_mtx_cache[(alpha,image_size_x,image_size_y)]=new_camera_mtx
    #
    #    return new_camera_mtx
    
    def __init__(self,**kwargs):
        
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
        
            setattr(self,kwarg,kwargs[kwarg])
            pass
        pass

    def evaluaterelativepose(self,part,dataset):
        # Determine orientation of camera relative to part based on dataset
        # Returns affine transformation matrix that takes points in
        # part space and returns points in camera coordinates 
        # NOTE: self.new_camera_mtx is in OpenCV coordinates (x right, y down, z forward)
        # NOTE: dataset.landmarkpixelcoords is in pixels, numbered from upper left hand corner
        
        # Returns dictionary of projected landmark pixel coordinates (numbered from upper left hand corner)
        # and matrix that will convert object coordinates (multiplied on the right) 
        # to camera coordinates in OpenGL coordinates (x right, y up, z backward) 

        assert(dataset.p==2) # image data must have 2 coordinates per point

        imagepoints=[]
        objectpoints=[]
        for landmarkname in dataset.landmarkpixelcoords:
            # NOTE: dataset.landmarkpixelcoords is in pixels, numbered from upper left hand corner
            if part.landmarks.has_landmark(landmarkname):
                #sys.stderr.write("Landmark %s: 2D coords: %s; 3D coords: %s\n" % (landmarkname,str(dataset.landmarkpixelcoords[landmarkname]),str(part.landmarks.lookup3d_objcoords(part,landmarkname))))
                imagepoints.append(dataset.landmarkpixelcoords[landmarkname])
                objectpoints.append(part.landmarks.lookup3d_objcoords(part,landmarkname))
                pass
            else: 
                sys.stderr.write("spatialnde.imageprojectionmodel: Unknown landmark: %s\n" % (landmarkname))
                pass
            pass

        if self.undistortion_already_applied:
            distcoeffs=None
            pass
        else:
            distcoeffs=self.calibration.distcoeffs
            pass

        if len(objectpoints) < 4: 
            raise ValueError("Not enough object/image point pairs to solve for pose (%d): minimum 4 points needed" % (len(objectpoints))) 

        #  image points arrays must be Nx1x2
        # and contiguous, per
        # https://stackoverflow.com/questions/44042323/opencv-error-assertion-failed-in-undistort-cpp-at-line-293
        # and https://github.com/opencv/opencv/issues/4943
        objectpoints_array=np.ascontiguousarray(np.array(objectpoints,dtype='f').reshape(len(objectpoints),1,3))
        imagepoints_array=np.ascontiguousarray(np.array(imagepoints,dtype='f').reshape(len(imagepoints),1,2))
        
        ransac_result= cv2.solvePnPRansac(objectpoints_array,
                                          imagepoints_array,
                                          self.new_camera_mtx,
                                          distCoeffs=distcoeffs)
        if len(ransac_result)==3:
            (rvec,tvec,inliers)=ransac_result
            pass
        else:
            (success,rvec,tvec,inliers)=ransac_result
            pass
        #(retval,rvec,tvec) = cv2.solvePnP(np.array(objectpoints,dtype='f'),
        #                                  np.array(imagepoints,dtype='f'),
        #                                  self.new_camera_mtx,
        #                                 distCoeffs=distcoeffs)
        
        #import pdb
        #pdb.set_trace()
        
        # OK. so rvec is a rotation quaternion, tvec is a translation
        # vector.
        # Assemble both into an affine 4x4 transformation matrix
        # See for reference: http://stackoverflow.com/questions/18637494/camera-position-in-world-coordinate-from-cvsolvepnp
        
        # obtain rotation matrix
        (rmat,jac)=cv2.Rodrigues(rvec)
        
        # assemble left side of 4x4 matrix
        Mleft = np.r_[rmat,np.zeros((1,3),dtype='d')]
        # assemble rightmost column of 4x4 matrix
        Mright = np.r_[ tvec.reshape(3), 1.0 ]
        
        # Concatenate left and right
        M=np.c_[Mleft,Mright]
        
        # This matrix M, given points in part space,
        # Should give us points in camera space.

        # BUT!!! The camera space used by OpenCV is
        # oriented with x right, y down, versus the
        # camera model used in 3D modeling (OpenGL, etc.)
        # is oriented with x right, y up.
        # So we need to multiply on the left by a matrix
        # that will do fix this.
        # just [ 1  0  0 0 ]
        #      [ 0 -1  0 0 ]
        #      [ 0  0  1 0 ]
        #      [ 0 0   0 1 ] will do that but is is a mirror, not a rotation!
        # We rotate y and z to maintain a right-handed frame: 
        #    [ 1  0  0  0 ]
        #    [ 0 -1  0  0 ]
        #    [ 0  0 -1  0 ]
        #    [ 0  0  0  1 ]  i.e. negate the 2nd and 3rd rows
        M[1,:]=-M[1,:]
        M[2,:]=-M[2,:]

        # Therefore output image points from M will be with x right, y up. 

        # Evaluate image points from object points using pose array
        projectedlandmarks={}
        for landmarkname in dataset.landmarkpixelcoords:
            if part.landmarks.has_landmark(landmarkname):
                landmark4dcoords=np.inner(M,np.concatenate((part.landmarks.lookup3d_objcoords(part,landmarkname),(1,))))
                # Convert back to opencv coords
                landmark4dcoords[1]=-landmark4dcoords[1]
                landmark4dcoords[2]=-landmark4dcoords[2]
                
                # Normalize
                landmark4dcoords = landmark4dcoords/landmark4dcoords[3]
                # Multiply by camera matrix
                landmark3dcoords = np.inner(self.new_camera_mtx,landmark4dcoords[:3])
                # Normalize

                landmark3dcoords = landmark3dcoords/landmark3dcoords[2]
                

                projectedlandmarks[landmarkname]=landmark3dcoords[:2]
                pass
            pass


        return (projectedlandmarks,M) # Should we also return corresponding input landmark coords?

    @classmethod
    def fromcalibfile(cls,calibfile,input_image_size_x,input_image_size_y,undistortion_already_applied=False,alpha=1.0):
        calibration = intrinsic_calibration.intrinsic_calibration.fromfile(calibfile)
        if not(undistortion_already_applied):
            # Work directly on the original camera matrix
            new_camera_mtx = calibration.cameramatrix

            assert(calibration.sizex==input_image_size_x)
            assert(calibration.sizey==input_image_size_y)
            
            pass
        else:
            # Work on the undistorted camera matrix
            
            (new_camera_mtx,validpixroi) = cv2.getOptimalNewCameraMatrix(cameraMatrix=calibration.cameramatrix, distCoeffs=calibration.distcoeffs,imageSize=(calibration.sizex, calibration.sizey), alpha=alpha, newImgSize=(input_image_size_x,input_image_size_y))
            pass
        
        return cls(calibration=calibration,
                   new_camera_mtx=new_camera_mtx,
                   undistortion_already_applied=undistortion_already_applied)
    
     
    @classmethod
    def fromcalibfiledataguzzler(cls,wfminfo,calibfile,raw=False,alphaoverride=None):
        # Load in data from a camera calibration file
        # to determine focal properties, etc.
        # also given dataguzzler reference waveform
        # indicating geometry
        
        import dg_eval
        import dg_metadata as dgm
        (ndim,DimLen,IniVal,Step,bases)=dg_eval.geom(wfminfo,raw=raw)
        


        if "spatialnde_cameracalib_alpha" in wfminfo.MetaData or alphaoverride is not None:
            alpha=dgm.GetMetaDatumWIDbl(wfminfo,"spatialnde_cameracalib_alpha",1.0)
            if alphaoverride is not None:
                alpha=alphaoverride
                pass
            # Indicates undistortion already applied
            undistortion_already_applied=True
            pass
        else:
            alpha=1.0
            undistortion_already_applied=False
            pass

        return cls.fromcalibfile(calibfile,
                                 DimLen[0],DimLen[1],
                                 undistortion_already_applied=undistortion_already_applied,
                                 alpha=alpha)
    
    
    
    
    pass

