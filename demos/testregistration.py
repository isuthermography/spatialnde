import sys
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

from spatialnde.ndeobj import ndepart
from spatialnde.coordframes import coordframe,concrete_affine
from spatialnde.imageprojectionmodel import ImageProjectionModel
from spatialnde.dataset import SurfaceStructuredGridDataSet
from spatialnde.imageprojection import imageprojection_prepare_float,imageprojection_float,imgbufzero,validitybufzero

from spatialnde.dataguzzler.dg_3d import ndepartparams_from_landmarked3d,blank_uv_from_landmarked3d

#ndepartparams = '[ [{"XUWKIV": [0, 0.470679, 0.778821], "JUSQOF": [0, 0.470432, 0.470844], "LUDRAZ": [0, 0.000597, 0.778506], "LITJID": [0, 0.470547, 0.530562], "KEXQET": [0, 0.159187, 0.618888], "FIPROJ": [0, 0.000135, 0.530768], "XEZHOD": [0, 0.469243, 0.717868], "SOFBEW": [0, 0.00087, 0.716984], "KILVOK": [0, 0.00053, 0.470313]}, {}], { "UV_ScalingParamBySurfaceNum": [[-0.00038232421875000001, -0.00038085937500000001], [0.78300000000000003, 0.78000000000000003]] } ]'




#x3dname="CR03-SPAR-01H_uv.x3d"
landmarked_3d_dgs = "CR03-SPAR-01H_landmarks.dgs"
calibfile="sc6000_1009_calib.cic"

# Dimensions of texture to create
#texwidth=800
#texheight=800

#texwidth=1200
#texheight=1200

# File with an image of the specimen, and marked landmarks, to map onto the 3D object
imgdgs_withlandmarks = "CR03-SPAR-01H_specimen_640x480_landmarks.dgs"  # Generated with dataguzzler_define_object_landmarks.confm4 and dataguzzler_define_landmarks_load_image.py 
imgdgs_chan="Param"
landmarked_3d_texchan = "CR03_SPAR_01H_tex"
landmarked_3d_x3dchan = "CR03_SPAR_01H"

(imgdgs_metadata,imgdgs_wfmdict)=dgf.loadsnapshot(imgdgs_withlandmarks)

objframe=coordframe()


cameraframe=coordframe()

(landmarked3d_metadata,landmarked_3d_wfmdict) = dgf.loadsnapshot(landmarked_3d_dgs)
# Create ndepartparams from landmarked_3d
ndepartparams=ndepartparams_from_landmarked3d(landmarked_3d_wfmdict,[ landmarked_3d_texchan ] )

# If we load the x3d directly, it doesn't have a texture URL set, so
# we can't match the texture URL to the scaling metadata from the snapshot
# (UV_ScalingParamsByTexURL). The copy within the .dgs DOES have a texture
# url set (from dataguzzler_define_landmarks_load_x3d) so
# we load it instead
#testobj = ndepart.fromx3d(objframe,ndepartparams,x3dname)
x3dbuf = BytesIO(dgm.GetMetaDatumWIStr(landmarked_3d_wfmdict[landmarked_3d_x3dchan],"X3DGeom","").encode('utf-8'))
testobj = ndepart.fromx3d(objframe,ndepartparams,x3dbuf)

projmodel = ImageProjectionModel.fromcalibfiledataguzzler(imgdgs_wfmdict[imgdgs_chan],calibfile,alphaoverride=1.0) # specify alphaoverride to indicate that undistortion has already been applied

# Fudge focal length in projmodel
focal_length=750.0 # pixels
projmodel.new_camera_mtx=np.array(((focal_length,0.0,320.0),
                                   (0.0,focal_length,240.0),
                                   (0.0,0.0,1.0)),dtype='d')

dataset = SurfaceStructuredGridDataSet.fromdataguzzlerpixelimage(imgdgs_wfmdict[imgdgs_chan],imgdgs_wfmdict)

# now we should be able to call projmodel.evaluaterelativepose
(projectedlandmarks,Pose) = projmodel.evaluaterelativepose(testobj,dataset)

# Check if pose works
PoseCopy=np.empty(Pose.shape,dtype='d')
PoseCopy[:,:]=Pose
# Fix 2nd and 3rd rows
PoseCopy[1,:]=-Pose[1,:]
PoseCopy[2,:]=-Pose[2,:]

testlandmark="KEXQET"
testcoordsxyz=testobj.implpart.surfaces[0].intrinsicparameterization.eval_xyz_uv(testobj.implpart.surfaces[0],*testobj.landmarks.landmarkdict2d[testlandmark][1:])
testcoords=np.inner(PoseCopy,np.r_[testcoordsxyz,1.0])
testcoords/=testcoords[2] # normalize by z coordinate

testpixelcoords=np.inner(projmodel.new_camera_mtx,testcoords[:3])
print("testpixelcoords=%s" % (str(testpixelcoords)))
print("datasetcoords=%s" % (str(dataset.landmarkpixelcoords[testlandmark])))

#raise ValueError("exit")

# Pose, given object coordinates returns camera coordinates.

# Define relationship between coord frames
PoseXForm = concrete_affine.fromaugmented(objframe,cameraframe,Pose)

# Assume single texture for all surfaces, overwriting existing texture wfm in loaded wfmdict
paramwfm = blank_uv_from_landmarked3d(landmarked_3d_wfmdict,landmarked_3d_texchan)


#imagebuf=np.zeros((texwidth,texheight),dtype='f')
# paramwfm.data and validitybuf are where the image is being mapped to. 
validitybuf=np.zeros(paramwfm.data.shape,dtype='f',order="F")
#angleofincidencebuf=np.zeros(paramwfm.data.shape,dtype='f',order="F")
parameterizationdict=dict([ (id(surface),(None,paramwfm.data.T,validitybuf.T)) for surface in testobj.implpart.surfaces ])


# Imagedat is the incoming image to be mapped
imagedat=dataset.data[:,:]
#imagedat=dataset.data[::10,::10]
imagedata=np.empty((1,imagedat.shape[0],imagedat.shape[1]),dtype='f')

imagedata[0,:,:]=imagedat

# Perform z-buffering
projparams=imageprojection_prepare_float(projmodel,cameraframe,[ testobj ],parameterizationdict,imagedat,2.0,uv_weightingblur_distance=.01) # pixel blur assumed to be 2 pixels; uv weightingblur of 1cm


surf_imgbuf=projparams.surfacelist[0][1][1]
surf_validitybuf=projparams.surfacelist[0][1][2]
surf_angleofincidencefactorbufuv=projparams.surfacelist[0][3]
surf_weightingbuf=projparams.surfacelist[0][2]

validitybufzero(projparams)
imgbufzero(projparams)

imageprojection_float(projparams,imagedata)

# paramwfm.data /= validitybuf  # division now done inside imageprojection_float()

#pl.imshow(imagebuf)

# dataset.landmarkpixelcoords["KEXQET"]= (369,286.5)  in pixels from upper left on dataset
#   testobj.landmarks.landmarkdict2d["KEXQET"] = [0, 0.159187, 0.618888]
# testobj.implpart.surfaces[0].intrinsicparameterization.eval_xyz_uv(testobj.implpart.surfaces[0],.159187,.618888)
#  -> array([ 0.00570404,  0.002     ,  0.31110466]) 
#  Convert to camera coordinates: np.inner(Pose,np.array([ 0.00570404,  0.002     ,  0.31110466,  1.0]))
#  ->  array([ 0.0151601 , -0.05353841, -6.14003631,  1.        ])
#  Invert Y and Z back to OpenCV camera coordinates
#  -> array([ 0.0151601 , 0.05353841, 6.14003631,  1.        ]) 
#  Normalize to x', y' -> array([ 0.00246906,  0.00871956,  1.        ])
#  Multiply by camera matrix
#  np.inner(projmodel.new_camera_mtx,np.array([ 0.00246906,  0.00871956,  1.        ]))
#  -> KILVOK:(416,351)
#  -> KEXQET:(343,329)

dgf.savesnapshot("/tmp/surfacemapped.dgs",landmarked_3d_wfmdict)
