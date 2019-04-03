import sys
import json

class point_landmarks(object):
    landmarkdict2d=None  # Dictionary by name of (parameterization_region_key,u,v) tuples, with reference
                         # to canonicalparameterization
    landmarkdict3d=None  # Dictionary by name of (x,y,z) tuples, 
    
    def __init__(self,**kwargs):

        self.landmarkdict={}
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
        
            setattr(self,kwarg,kwargs[kwarg])
            pass
        pass

    @classmethod
    def fromparams(cls,landmarksparams):
        if landmarksparams is not None:
            (landmarkdict2d,landmarkdict3d)=landmarksparams
            pass
        else:
            landmarkdict2d=None
            landmarkdict3d=None
            pass
        
        return cls(landmarkdict2d=landmarkdict2d,landmarkdict3d=landmarkdict3d)

    def toparams(self):
        return json.dumps((self.landmarkdict2d,self.landmarkdict3d))

    def addlandmark3d(self,name,xyz):
        if name in self.landmarkdict2d:
            del self.landmarkdict2d[name]
            pass
        
        self.landmarkdict3d[name]=xyz
        pass

    def has_landmark(self,name):
        return name in self.landmarkdict2d or name in self.landmarkdict3d
    

    def lookup3d_objcoords(self,ndepart,name):
        if name in self.landmarkdict3d:
            return self.landmarkdict3d[name]
        if name in self.landmarkdict2d:
            #if ndepart.uvparam is not None:
            #    raise NotImplementedError("Need to implement landmark lookup in non-canonical parameterization")
            #else:
                # parameterization_region_key is in this case is assumed
                # to be the surface index within the cadpart
                sys.stderr.write("spatialnde.landparks.py: WARNING: Performing lookup assuming surface index is zero\n")
                # Surface index used to be the surfacenum, but
                # we will reorg how we define surfaces anyway.
                # So right now this is hardwired to 0.
                #
                # landmarkdict2d[name][0] is now the texture url,
                # at least for landmarks_params generated from dg_3d.ndepartparams_from_landmarked3d() 
                #implpartsurface=ndepart.implpart.surfaces[self.landmarkdict2d[name][0]]
                implpartsurface=ndepart.implpart.surfaces[0]
                try :
                    retval=implpartsurface.intrinsicparameterization.eval_xyz_uv(implpartsurface,self.landmarkdict2d[name][1],self.landmarkdict2d[name][2])
                    pass
                except ValueError:
                    raise ValueError("Could not find (X,Y,Z) coordinates for landmark %s at (u,v)=(%g,%g)" % (name,self.landmarkdict2d[name][1],self.landmarkdict2d[name][2]))
                return retval
            #pass

        raise KeyError(name)

    def lookup3d(self,ndepart,coordframe,name):
        objcoords=self.lookup3d_objcoords(ndepart,name)
        if coordframe is not None:
            
            coords=ndepart.frame.transformto(coordframe,objcoords)
            return coords
        return objcoords
    pass
