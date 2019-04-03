import sys
import traceback

import numpy as np

# Concept: Define coordinate frames,
# Then define relations between these frames.
# For now relations are concrete, but in future
# they might include uncertainty.

#
# Note: So long as only a single path from any frame a to any other
#       frame b exists, then that path completely defines the relationship
#       Any relation in that path can be updated, which would then
#       implicitly update the relationship.
#          Once more than one path is possible, there is a question of
#       consistency. If one relation is updated then any parallel
#       relations might need to be updated too.
# Rule: Should only update a relation if there is no other path
#       from one side to the other side. That way, the update
#       cannot affect consistency. 
#


class coordrelation(object):
    frame_u=None  # coordframe u object of this coordinate relation
    frame_u_id=None
    frame_v=None  # coordframe v object of this coordinate relation
    frame_v_id=None

    def frame_from_id(self,frame_id):
        if frame_id==self.frame_u_id:
            return self.frame_u
        elif frame_id == self.frame_v_id:
            return self.frame_v
        else:
            raise ValueError("coordrelation does not involve frame with id %d" % (frame_id))
        pass
    pass

class concrete_affine(coordrelation):
    # Represent frame_u->frame_v as
    #  x_v = A x_u + b
    # or in augmented matrix form 
    # [ x_v ; 1 ] = [ A   b  ; 0   1 ][ x_u ; 1 ] 

    M = None  # augmented matrix, eats u on right, gives v
    Minv = None # augmented inverse, or None if not invertible; eats v on right, gives u

    def __init__(self,**kwargs):
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
        
            setattr(self,kwarg,kwargs[kwarg])
            pass

        assert(self.frame_u is not None)
        self.frame_u_id=id(self.frame_u)

        assert(self.frame_v is not None)
        self.frame_v_id=id(self.frame_v)

        self.frame_u.relations[self.frame_v_id]=self
        self.frame_v.relations[self.frame_u_id]=self
        
        
        pass

    @classmethod
    def fromaugmented(cls,frame_u,frame_v,M,invertible=True):
        # Define transform that eats vectors in frame u, gives vectors in frame v
        if invertible:
            Minv=np.linalg.inv(M)
            pass
        else:
            Minv=None
            pass
        
        return cls(frame_u=frame_u,frame_v=frame_v,M=M,Minv=Minv)

    @classmethod
    def fromtransform(cls,frame_u,frame_v,A,b=None,invertible=True):
        if b is None:
            b=np.zeros((3,),dtype='d')
            pass

        # assemble left side
        Mleft = np.r_[A,np.zeros((1,3),dtype='d')]
        # Assemble right side
        Mright = np.r_[ b, 1.0 ]
        # Concatenate left and right
        M=np.c_[Mleft,Mright]

        return cls.fromaugmented(frame_u,frame_v,M,invertible=invertible)

    def gettransform(self,source_id):
        if source_id==self.frame_u_id:
            return self.M
            pass
        elif source_id==self.frame_v_id:
            return self.Minv
            pass
        else:
            raise ValueError("gettransform(): source_id is not u or v of transform")
        pass

    
    def destframe(self,source_id):
        if source_id==self.frame_u_id:
            return self.frame_v
        elif source_id==self.frame_v_id:
            return self.frame_u
        else:
            raise ValueError("destframe(): source_id is not u or v of transform")
        pass
    

    def applytransform(self,source_id,coords,cartesianaxis=-1,inhomogeneous_vector=False):
        # source_id is id of the source coordinate frame
        # coords is an array with an axis representing input coordinates.
        # The input coordinate axis much either be length 3,
        # or length 4 (representing homogenous coordinates)
        # it is selected by the "cartesianaxis" parameter, which
        # defaults to -1 (the last axis)
        # if the input coordinate axis is length 3 and
        # inhomogenous_vector==True, then the transformation
        # will be done "Vector style" i.e. with the 4th coordinate
        # interpreted as 0 rather than 1
        
        
        if cartesianaxis < 0:
            cartesianaxis += len(coords.shape)
            pass
                
        if coords.shape[cartesianaxis]==3:
            #coords=np.c_[coords,np.ones((coords.shape[0],1),dtype='d')]
            washomogeneous=False
            pass
        else:
            assert(coords.shape[cartesianaxis]==4)
            washomogeneous=True
            pass

        if washomogeneous: 
            if source_id==self.frame_u_id:
                # Axis rearrangement
                # Tensordot leaves operating axis out front
                # so if axes were ABCDEF, where D was cartesian (cartesianaxis=3)
                # then after tensordot we have
                # GABCEF where G is transformed cartesian
                # So we transpose with 1,2,3,0,4,5
                
                res=np.tensordot(self.M,coords,axes=((1,),(cartesianaxis,))).transpose(np.concatenate((np.arange(1,cartesianaxis+1),0,np.arange(cartesianaxis+1,len(coords.shape)))))
                pass
            elif source_id==self.frame_v_id:
                res=np.tensordot(self.Minv,coords,axes=((1,),(cartesianaxis,))).transpose(np.concatenate((np.arange(1,cartesianaxis+1),0,np.arange(cartesianaxis+1,len(coords.shape)))))
                pass
            else:
                raise ValueError("applytransform(): source_id is not u or v of transform")
            pass
        else:
            # not washomogeneous
            # Apply M or Minv in pieces:
            # M = [ M11 M12 ] [ coords                  ]
            #     [ M21 M22 ] [ ~inhomogeneous_vector    ] 
            if source_id==self.frame_u_id:
                res=np.tensordot(self.M[:3,:3],coords,axes=((1,),(cartesianaxis,)))
                # Evaluate normalization factor
                normfactor=np.tensordot(self.M[3,:3],coords,axes=((0,),(cartesianaxis,)))
                # Axis 0 is now the cartesian axis
                if not inhomogeneous_vector:
                    # Add in [ M12 ] to axis 0 
                    res += self.M[:3,3].reshape(*([4]+[1]*(len(coords.shape)-1)))
                    normfactor += self.M[3,3]
                    pass

                # normalize
                res /= normfactor.reshape([1]+coords.shape)

                # transpose cartesian axis back where it belongs
                res=res.transpose(np.concatenate((np.arange(1,cartesianaxis+1),0,np.arange(cartesianaxis+1,len(coords.shape)))))                
                pass
            elif source_id==self.frame_v_id:
                res=np.tensordot(self.Minv[:3,:3],coords,axes=((1,),(cartesianaxis,)))
                # Evaluate normalization factor
                normfactor=np.tensordot(self.Minv[3,:3],coords,axes=((0,),(cartesianaxis,)))
                # Axis 0 is now the cartesian axis
                if not inhomogeneous_vector:
                    # Add in [ M12 ] to axis 0 
                    res += self.Minv[:3,3].reshape(*([4]+[1]*(len(coords.shape)-1)))
                    normfactor += self.Minv[3,3]
                    pass

                # normalize
                res /= normfactor.reshape([1]+coords.shape)
                
                # transpose cartesian axis back where it belongs
                res=res.transpose(np.concatenate((np.arange(1,cartesianaxis+1),0,np.arange(cartesianaxis+1,len(coords.shape)))))                

                pass
            else:
                raise ValueError("applytransform(): source_id is not u or v of transform")
            
            ## normalize homogeneous results
            #res=res[:,0:3]/res[:,3:4] # Return 3D normalized coordinates
            pass

        return res
    
    
    pass


    
class coordframe(object):
    # Coordinate frame defined through its relations to
    # other coordinate frames
    
    id=None # Unique ID for this coordinate frame
    relations=None # Dictionary, by id, of relations to other coordinate frames

    
    def __init__(self):
        self.id=id(self)
        self.relations={}
        pass


    def _find_paths(self, other, path=[]):
        # Loosely based on https://www.python.org/doc/essays/graphs/
        if self is other:
            return [ path ]
        
        paths = []
        for frame_id in self.relations:
            if frame_id not in path:

                dest=self.relations[frame_id].frame_from_id(frame_id)
                newpaths = dest._find_paths(other, path + [ dest.id ] )
                for newpath in newpaths:
                    paths.append(newpath)
                    pass
                pass
        return paths
    
    
    def transformto(self,other,coords,cartesianaxis=-1,inhomogeneous_vector=False):
        # Other is another coordinate frame
        # coords is an array with an axis representing input coordinates.
        # The input coordinate axis much either be length 3,
        # or length 4 (representing homogenous coordinates)
        # it is selected by the "cartesianaxis" parameter, which
        # defaults to -1 (the last axis)
        # if the input coordinate axis is length 3 and
        # inhomogenous_vector==True, then the transformation
        # will be done "Vector style" i.e. with the 4th coordinate
        # interpreted as 0 rather than 1
        
        # Walk relations graph to find path
        paths=self._find_paths(other)

        if len(paths) < 1:
            raise ValueError("Unable to find relation path to output frame")
        elif len(paths) > 1:
            sys.stderr.write("WARNING: Multiple transformation paths to output frame. Using first path only; traceback:\n")
            traceback.print_stack(sys.stderr)
            pass

        path=paths[0]  # path is a list of frame ids

        curframe=self
        curcoords=coords
        for pathentry in path:
            curcoords=curframe.relations[pathentry].applytransform(curframe.id,curcoords,cartesianaxis=cartesianaxis,inhomgeneous_vector=inhomogeneous_vector)
            curframe=curframe.relations[pathentry].destframe(curframe.id)
            pass
        assert(curframe==other)
        return curcoords


    def transformationmatrix(self,other):
        # Create transformation matrix that will convert coordinates in our frame into
        # coordinates in other frame
        # Other is another coordinate frame
        # coords is an nx3 or nx4 array of cartesian or homogenious coordinates
        
        # Walk relations graph to find path
        paths=self._find_paths(other)

        if len(paths) < 1:
            raise ValueError("Unable to find relation path to output frame")
        elif len(paths) > 1:
            sys.stderr.write("WARNING: Multiple transformation paths to output frame. Using first path only; traceback:\n")
            traceback.print_stack(sys.stderr)
            pass

        path=paths[0]  # path is a list of frame ids

        curframe=self
        curx=np.eye(4,dtype='d') # Identity matrix
        for pathentry in path:
            
            curx=np.dot(curframe.relations[pathentry].gettransform(curframe.id),curx)
            curframe=curframe.relations[pathentry].destframe(curframe.id)
            
            pass
        assert(curframe==other)
        return curx
    
    pass

