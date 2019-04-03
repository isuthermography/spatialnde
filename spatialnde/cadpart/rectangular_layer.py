import numpy as np

from .cadpart import cadpart
from .polygonalsurface import polygonalsurface
from .polygonalsurface_intrinsicparameterization import polygonalsurface_texcoordparameterization

class rectangular_layer(cadpart):
    # Thin layer 
    length = None # length in x
    width = None # width in y
    x0=None # x coordinate of starting (min-x) corner
    y0=None # y coordinate of starting (min-y) corner
    z=None # z coordinate of layer

    def __init__(self,**kwargs):
        appearance=None
        
        if "appearance" in kwargs:
            appearance=kwargs["appearance"]
            del kwargs["appearance"]
            pass
        
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
        
            setattr(self,kwarg,kwargs[kwarg])
            pass

        super(rectangular_layer,self).__init__()

        vertices=np.array(((self.x0,self.y0,self.z),
                           (self.x0+self.length,self.y0,self.z),
                           (self.x0+self.length,self.y0+self.width,self.z),
                           (self.x0,self.y0+self.width,self.z)),
                          dtype='d')
        normals=np.array(((0.0,0.0,1.0),), # top surface
                         dtype='d')
        
        surfacepoint=np.array((0.0,0.0,0.0),dtype='d')
        frontnormal=np.array((0.0,0.0,1.0),dtype='d')
        backnormal=np.array((0.0,0.0,-1.0),dtype='d')

        #self.refpoints=np.array(((self.length/2.0,self.width/2.0,0.0),) # top surface
        #dtype='d')

        
        vertexidx=np.array((0,1,2,3,-1),  # front surface
                           dtype='i4')
        vertexidx_indices=np.array((0,),  # front surface
                                   dtype='i4')
        
        #self.surfaceids=np.array((0,  # front surface
        #                          1, # back surface
        #                          2, # right surface
        #                          3, # left surface
        #                          4, # up surface
        #                          5), # down surface
        #                        dtype='d')

        texcoord = np.array( (((0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)),),dtype='d')
        texcoordidx = np.array((0,1,2,3,-1),dtype='i4')
        
        self.surfaces=[]

        for cnt in range(1):
            self.surfaces.append(polygonalsurface.fromvertices(surfaceid=cnt,
                                                               vertices=vertices,
                                                               vertexidx_indices=vertexidx_indices,
                                                               vertexidx=vertexidx[(cnt*5):((cnt+1)*5)],
                                                               normals=normals[cnt:(cnt+1),:],
                                                               vertexnormals=None,
                                                               cadpartparams=None,
                                                               appearance=appearance,
                                                               solid=False,
                                                               buildintrinsicparameterization=lambda surf,params: polygonalsurface_texcoordparameterization.new(surf,texcoord,texcoordidx,params), # NOTE: texcoord cannot be modified between loop iterations without adding another layer of lambdas here. 
                                                               tol=self.tol))
            
            pass
        pass

    pass
