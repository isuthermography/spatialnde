# Example of creating an object from vertices and
# writing it to an X3D file
import numpy as np

from spatialnde.coordframes import concrete_affine
from spatialnde.coordframes import coordframe
from spatialnde.ndeobj import ndepart
from spatialnde.cadpart.cadpart import cadpart
from spatialnde.cadpart.polygonalsurface_texcoordparameterization import polygonalsurface_texcoordparameterization
from spatialnde.cadpart.appearance import simple_material
from spatialnde.cadpart.appearance import texture_url
from spatialnde.exporters.x3d import X3DSerialization
from spatialnde.exporters.vrml import VRMLSerialization

if __name__=="__main__":

    # define object frame
    objframe = coordframe()

    x=np.arange(-5.,5.,0.1)
    y=np.arange(-6.,6.,0.15)

    i=np.arange(x.shape[0])
    j=np.arange(y.shape[0])

    # define (y,x) and (j,i) basis
    (yy,xx)=np.meshgrid(y,x,indexing="ij")
    (jj,ii)=np.meshgrid(j,i,indexing="ij")
    
    # evaluate z as a sinc function
    zz=np.sinc(np.sqrt(xx**2+yy**2))

    # vertices indexed by j,i,0..2
    vertices_ji=np.array((xx,yy,zz)).transpose(1,2,0)   

    # vertices must just be n x 3 
    vertices = np.ascontiguousarray(vertices_ji.reshape(y.shape[0]*x.shape[0],3))

    # Need to be able to map from (j,i) to index into vertices
    verticesindex=np.arange(y.shape[0]*x.shape[0]).reshape(j.shape[0],i.shape[0])

    # use i/j coordinates for texture coordinates, but remap them to 0...1
    texcoord_ji = np.array((ii/(i.shape[0]-1.0),jj/(j.shape[0]-1.0))).transpose(1,2,0)   

    # texcoord must just be n x 2
    texcoord = np.ascontiguousarray(texcoord_ji.reshape(j.shape[0]*i.shape[0],2))

    # now define a bunch of (roughly) quadrilateral polygons
    vertexidx_list = []
    texcoordidx_list = []
    for jcnt in range(j.shape[0]-1):
        for icnt in range(i.shape[0]-1):
            # write vertices ccw viewed from +z direction
            vertexidx_list.extend([ verticesindex[jcnt,icnt],
                                    verticesindex[jcnt,icnt+1],
                                    verticesindex[jcnt+1,icnt+1],
                                    verticesindex[jcnt+1,icnt],
                                    -1])

            # Same for texcoord
            texcoordidx_list.extend([ verticesindex[jcnt,icnt],
                                      verticesindex[jcnt,icnt+1],
                                      verticesindex[jcnt+1,icnt+1],
                                      verticesindex[jcnt+1,icnt],
                                      -1])
            pass
        pass

    # convert vertexidx and texcoordidx to numpy arrays
    vertexidx=np.array(vertexidx_list,dtype='i4')
    texcoordidx=np.array(texcoordidx_list,dtype='i4')

    appearance=texture_url.from_url("texture.jpg",DefName="Texture") # reference texture from a URL relative to the x3d/vrml we will write

    # Define function for creating the texture coordinate parameterization
    buildintrinsicparameterization=lambda surf,cadpartparams: polygonalsurface_texcoordparameterization.new(surf,texcoord,texcoordidx,appearance,cadpartparams=cadpartparams)
    
    # Define the function for constructing the cadpart object from an array of vertices
    buildcadpart = lambda cadpartparams: cadpart.fromvertices(vertices,
                                                              vertexidx,
                                                              buildintrinsicparameterization=buildintrinsicparameterization,
                                                              cadpartparams=cadpartparams,
                                                              appearance=appearance)

    # This creates the object trees by calling the buildcardpart() function
    testplate=ndepart.from_implpart_generator(objframe,None,buildcadpart)



    # Create a writer object for generating x3d
    x3dwriter=X3DSerialization.tofileorbuffer("testsinc.x3d",x3dnamespace=None)

    # Write the object tree to the X3D scene
    testplate.X3DWrite(x3dwriter,objframe) #,appearance=X3DAppear)

    # Close the x3d writer
    x3dwriter.finish()

    pass
    
