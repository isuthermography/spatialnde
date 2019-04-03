import numpy as np

from .cadpart import cadpart
from .polygonalsurface import polygonalsurface


class rectangular_plate(cadpart):
    length = None # length in x
    width = None # width in y
    thickness = None # Thickness in z

    

    def __init__(self,**kwargs):
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
        
            setattr(self,kwarg,kwargs[kwarg])
            pass

        super(rectangular_plate,self).__init__()

        vertices=np.array(((0.0,0.0,0.0),
                           (self.length,0.0,0.0),
                           (self.length,self.width,0.0),
                           (0.0,self.width,0.0),
                           (0.0,0.0,self.thickness),
                           (self.length,0.0,self.thickness),
                           (self.length,self.width,self.thickness),
                           (0.0,self.width,self.thickness)),
                          dtype='d')
        normals=np.array(((0.0,0.0,1.0), # top surface
                          (0.0,0.0,-1.0), # back surface
                          (1.0,0.0,0.0), # right surface
                          (-1.0,0.0,0.0), # left surface
                          (0.0,1.0,0.0), # up surface
                          (0.0,-1.0,0.0)), # down surface
                         dtype='d')
        
        backsurfacepoint=np.array((0.0,0.0,0.0),dtype='d')
        frontsurfacepoint=np.array((0.0,0.0,self.thickness),dtype='d')
        frontnormal=np.array((0.0,0.0,1.0),dtype='d')
        backnormal=np.array((0.0,0.0,-1.0),dtype='d')

        refpoints=np.array(((self.length/2.0,self.width/2.0,self.thickness), # top surface
                            (self.length/2.0,self.width/2.0,0.0), # back surface
                            (self.length,self.width/2.0,self.thickness/2.0), # right surface
                            (0.0,self.width/2.0,self.thickness/2.0), # left surface
                            (self.length/2.0,self.width,self.thickness/2.0), # up surface
                            (self.length/2.0,0.0,self.thickness/2.0)), # down surface
                        dtype='d')
        
        vertexidx=np.array((4,5,6,7,-1,  # front surface
                            3,2,1,0,-1, # back surface
                            1,2,6,5,-1, # right surface
                            4,7,3,0,-1, # left surface
                            3,7,6,2,-1, # up surface
                            0,1,5,4,-1), # down surface
                           dtype='i4')
        vertexidx_indices=np.array((0,5,10,15,20,15),dtype='i4')
        numvertices=np.array((4,4,4,4,4,4),dtype='i4')
        #self.surfaceids=np.array((0,  # front surface
        #                          1, # back surface
        #                          2, # right surface
        #                          3, # left surface
        #                          4, # up surface
        #                          5), # down surface
        #                        dtype='d')
        self.surfaces=[]

        for cnt in range(6):
            self.surfaces.append(polygonalsurface.fromvertices(surfaceid=cnt,
                                                               vertices=vertices,
                                                               vertexidx_indices=np.array((0,),dtype=np.int32),
                                                               vertexidx=vertexidx[(cnt*5):((cnt+1)*5)],
                                                               numvertices=numvertices,
                                                               normals=normals[cnt:(cnt+1),:],
                                                               normalidx=None,
                                                               cadpartparams=None,
                                                               tol=self.tol))
            
            pass
        pass

    pass
