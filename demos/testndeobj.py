import sys
import numpy as np

try:
    from cStringIO import StringIO  # python 2.x
    pass
except ImportError:
    from io import StringIO # python 3.x
    pass

from spatialnde.coordframes import concrete_affine
from spatialnde.coordframes import coordframe
from spatialnde.ndeobj import ndepart
from spatialnde.cadpart.rectangular_plate import rectangular_plate
from spatialnde.cadpart.appearance import simple_material
from spatialnde.cadpart.appearance import texture_url
from spatialnde.exporters.x3d import X3DSerialization
from spatialnde.exporters.vrml import VRMLSerialization

if __name__=="__main__":

    # define lab frame
    labframe = coordframe()

    # define camera frame
    camframe = coordframe()

    # define object frame
    objframe = coordframe()


    # Define 90 deg relation and .2 meter offset between lab frame and camera frame:
    # Define camera view -- in reality this would come from
    # identifying particular reference positions within the lab
    concrete_affine.fromtransform(labframe,camframe,
                                  np.array( (( 0.0,1.0,0.0),
                                             (-1.0,0.0,0.0),
                                             (0.0,0.0,1.0)),dtype='d'),
                                  b=np.array( (0.0,0.2,0.0), dtype='d'),
                                  invertible=True)


    # Define object coordinates relative to lab 
    concrete_affine.fromtransform(labframe,objframe,
                                  np.array( (( 1.0/np.sqrt(2.0),-1.0/np.sqrt(2.0),0.0),
                                             (1.0/np.sqrt(2.0),1.0/np.sqrt(2.0),0.0),
                                             (0.0,0.0,1.0)),dtype='d'),
                                  b=np.array( (0.0,0.0,0.4), dtype='d'),
                                  invertible=True)



    testplate=ndepart.from_implpart_generator(objframe,None,lambda cadpartparams: rectangular_plate(length=.1000,  # m
                                                                                                   width=.0300,
                                                                                                   thickness=.01,
                                                                                                   tol=1e-3))

    x3dnamespace=None


    Appear = simple_material.from_color((.6,.01,.01),DefName="coloring")
    
    
    shape_els=[]
    
    #shape_els.extend(testplate.X3DShapes(objframe,x3dnamespace=x3dnamespace))

    motor=ndepart.fromstl(objframe,None,"/usr/local/src/opencascade/data/stl/shape.stl",tol=1e-6,recalcnormals=False,metersperunit=1.0e-3,defaultappearance=Appear)

    #motor_els = motor.X3DShapes(objframe,x3dnamespace=x3dnamespace,appearance=X3DAppear)
    #shape_els.extend(motor_els)

    # DEFINE ALTERNATE APPEARANCE!!!
    moondial_tex=texture_url.from_url("/usr/local/src/h3d/svn/H3DAPI/examples/x3dmodels/moondial/texture.jpg",DefName="Texture") # reference texture from a URL relative to the x3d/vrml we will write
    

    moondial=ndepart.fromx3d(objframe,None,"/usr/local/src/h3d/svn/H3DAPI/examples/x3dmodels/moondial/moondial_orig.x3d",tol=1e-6)
    
    moondial.assign_appearances(moondial_tex)

    x3dwriter=X3DSerialization.tofileorbuffer("objout.x3d",x3dnamespace=x3dnamespace)

    moondial.X3DWrite(x3dwriter,objframe) #,appearance=X3DAppear)
    
    x3dwriter.finish()

    vrmlwriter = VRMLSerialization.tofileorbuffer("objout.wrl")
    moondial.VRMLWrite(vrmlwriter,objframe,UVparameterization=None) #,appearance=VRMLAppear)
    vrmlwriter.finish()
    
    
    import dataguzzler as dg
    import dg_file as dgf
    import dg_metadata as dgm
    import scipy.ndimage


    # Internal reference to texture for dataguzzler:
    dg_moondial_tex=texture_url.from_url("#moondial_tex",DefName="Texture") # reference to dataguzzler channel... really for .dgs file only
    moondial.assign_appearances(dg_moondial_tex)

    #  Generate dataguzzler VRML string
    
    DGVRMLBuf=StringIO()
    dgvrmlwriter = VRMLSerialization.tofileorbuffer(DGVRMLBuf)
    moondial.VRMLWrite(dgvrmlwriter,objframe,UVparameterization=None) # ,appearance=VRMLAppear)
    dgvrmlwriter.finish()

    DGX3DBuf=StringIO()
    dgx3dwriter = X3DSerialization.tofileorbuffer(DGX3DBuf)
    moondial.X3DWrite(dgx3dwriter,objframe,UVparameterization=None) # ,appearance=VRMLAppear)
    dgx3dwriter.finish()

    
    
    wfmdict={}
    wfmdict["moondial"]=dg.wfminfo()
    wfmdict["moondial"].Name="moondial"
    wfmdict["moondial"].dimlen=np.array((),dtype='i8')
    dgm.AddMetaDatumWI(wfmdict["moondial"],dgm.MetaDatum("VRML97Geom",DGVRMLBuf.getvalue()))
    dgm.AddMetaDatumWI(wfmdict["moondial"],dgm.MetaDatum("X3DGeom",DGX3DBuf.getvalue()))

    teximage=scipy.ndimage.imread("/usr/local/src/h3d/svn/H3DAPI/examples/x3dmodels/moondial/texture.jpg",flatten=True).astype(np.float32).T
    wfmdict["moondial_tex"]=dg.wfminfo()
    wfmdict["moondial_tex"].Name="moondial_tex"
    # ***!!!! NOTE: Must adjust contrast on coloring channel
    # in dg_scope to get a nice colormap
    wfmdict["moondial_tex"].ndim=2
    wfmdict["moondial_tex"].dimlen=np.array(teximage.shape,dtype='i8')
    wfmdict["moondial_tex"].n=np.prod(teximage.shape)
    wfmdict["moondial_tex"].data=teximage[:,::-1] # Flip y axis so it appears correct in scope

    dgm.AddMetaDatumWI(wfmdict["moondial"],dgm.MetaDatum("TextureChan_coloring","coloring:0"))
    dgf.savesnapshot("objout.dgs",wfmdict)
   
    pass
    
