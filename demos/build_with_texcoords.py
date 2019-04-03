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

    # vertex coordinates in 3D
    vertices=np.array( ((0,0,0),
                        (1,0,0),
                        (1,1,0),
                        (0,1,0)), dtype='d')

    # vertexidx has the indices of the vertices in the first polygon,
    # followed by -1, followed by the indices of the vertices in the
    # 2nd polygon, followed by -1, etc .
    vertexidx = np.array( (0,1,2,3,-1),dtype='i8')

    # texture (surface parameterization) coordinates in 3D
    texcoord=np.array( ((0,0),
                        (1,0),
                        (1,1),
                        (0,1)), dtype='d')
    # texcoordidx has the indices of the texture coordinates
    # for the vertices in the first polygon,
    # followed by -1, followed by the indices of the texture
    # coordinates for the vertices in the
    # 2nd polygon, followed by -1, etc .
    texcoordidx = np.array( (0,1,2,3,-1),dtype='i8')
    

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
    x3dwriter=X3DSerialization.tofileorbuffer("testplate.x3d",x3dnamespace=None)

    # Write the object tree to the X3D scene
    testplate.X3DWrite(x3dwriter,objframe) #,appearance=X3DAppear)

    # Close the x3d writer
    x3dwriter.finish()

    pass
    
