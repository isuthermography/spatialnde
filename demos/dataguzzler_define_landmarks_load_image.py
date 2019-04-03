# Load in an image file into dataguzzler for assignment of landmarks/fiducials
# For use with dataguzzler_define_object_landmarks.confm4
import sys
import os
import numpy as np
import scipy.ndimage

import dataguzzler as dg
#import dg_file as dgf
import dg_comm as dgc
import dg_metadata as dgm

try:
    from cStringIO import StringIO  # python 2.x
    pass
except ImportError:
    from io import StringIO # python 3.x
    pass


from spatialnde.ndeobj import ndepart



cbrows=20
cbcols=20

imagename=sys.argv[1]
(imagepath,imagefile)=os.path.split(imagename)
imagebasename=os.path.splitext(imagefile)[0]

image_array=scipy.ndimage.imread(imagename,mode='F') # read image as 32-bit float



wfmdict={}

paramname="Param"
paramrawname="ParamRaw"
wfmdict[paramrawname]=dg.wfminfo()
wfmdict[paramrawname].Name=paramrawname
wfmdict[paramrawname].ndim=2
wfmdict[paramrawname].dimlen=np.array(image_array.shape[::-1],dtype='i8')
wfmdict[paramrawname].wfmrevision=0
wfmdict[paramrawname].n=np.prod(wfmdict[paramrawname].dimlen)
wfmdict[paramrawname].data=np.zeros(wfmdict[paramrawname].dimlen,dtype='f')
dgm.AddMetaDatumWI(wfmdict[paramrawname],dgm.CreateMetaDatumStr("Coord1","X Position"))
dgm.AddMetaDatumWI(wfmdict[paramrawname],dgm.CreateMetaDatumDbl("IniVal1",0.0))
dgm.AddMetaDatumWI(wfmdict[paramrawname],dgm.CreateMetaDatumStr("Coord2","Y Position"))
dgm.AddMetaDatumWI(wfmdict[paramrawname],dgm.CreateMetaDatumDbl("IniVal2",0.0))

dgm.AddMetaDatumWI(wfmdict[paramrawname],dgm.CreateMetaDatumStr("Units1","pixels"))
dgm.AddMetaDatumWI(wfmdict[paramrawname],dgm.CreateMetaDatumDbl("Step1",1.0))
dgm.AddMetaDatumWI(wfmdict[paramrawname],dgm.CreateMetaDatumStr("Units2","pixels"))
dgm.AddMetaDatumWI(wfmdict[paramrawname],dgm.CreateMetaDatumDbl("Step2",1.0))


#xpos=np.arange(imageshape[0],dtype='d')
#ypos=np.arange(imageshape[1],dtype='d')

#xchecker = (xpos//(imageshape[0]*1.0/(cbcols)) % 2).astype(np.bool)
#ychecker = (ypos//(imageshape[1]*1.0/(cbrows)) % 2).astype(np.bool)



#wfmdict[paramrawname].data[:,:]=xchecker.reshape(imageshape[0],1) ^ ychecker.reshape(1,imageshape[1])  # XOR operator
wfmdict[paramrawname].data[:,:]=image_array[::-1,:].T

#dgf.savesnapshot(os.path.join(x3dpath,x3dbasename+".dgs"),wfmdict)

dgch=dgc.client();

for wfmname in wfmdict:
    dgc.uploadwfm(dgch,wfmdict[wfmname])
    pass

dgc.command(dgch,"WFM:COPY ParamRaw ParamSrc")
dgc.close(dgch)




