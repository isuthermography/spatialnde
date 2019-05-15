import sys
import json
import collections

import numpy as np

from .loaders.x3d import x3d_indexedfaceset_loader
from .loaders.x3d import x3d_indexedfaceset
from .loaders.stl import loadstl

from .polygonalsurface_planarparameterization import identify_surfaces
from .polygonalsurface_planarparameterization import uniquevertices
from .polygonalsurface_planarparameterization import polygonalsurface_planarparameterization
from .polygonalsurface import polygonalsurface
from .polygonalsurface import find_vertexidx_indices
from .polygonalsurface import calcnormals

def vecnorm(a,axis=-1):
    
    # compatibility for old versions of Numpy
    return np.sqrt(np.sum(a**2.0,axis))


class cadpart(object):
    #  class for a part
    # that we are doing NDE testing on
    surfaces=None # List or array of surfaces defined per frame
    implpartparams=None # string parameters of extrinsic info about how the
                       # CAD object has been interpreted... e.g for the STL laoder,
                       # JSON string identifying surfaces that come from particular polygons 


    # !!!*** Need to define and be able to evaluate surface curvature at intersection points
    # Need to identify u and v coordinates of intersection points
    # From u and v coordinates, need to be able to evaluate surface curvature.
    
                    
    tol=None # tolerance


    
    def __init__(self,**kwargs):
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
        
            setattr(self,kwarg,kwargs[kwarg])
            pass        
        pass

    def adjacent(self,ray1_intersect_info,ray2_intersect_info):
        # determine if the intersections of ray1 and ray2 with the
        # object are adjacent, or if there is some substantive step between them
        # *** Not yet implemented
        raise NotImplementedError("cadpart/adjacent")
        pass


    def ray_intersections(self,localpoints,localdirecs):
        # rays is a set of n rays, each of which is specified by 6 coordinates, source point then unit direction vector, frame is the coordframe of those rays

        n=localpoints.shape[0]
        assert(n==localdirecs.shape[0])

        firstintersectpoints=np.empty((n,3),dtype='d')
        firstintersectpoints.fill(np.inf)

        firstintersectdists=np.empty((n,),dtype='d')
        firstintersectdists.fill(np.inf)

        firstintersectnormals=np.zeros((n,3),dtype='d')

        for surface in self.surfaces:
            surface.ray_intersections(self,localpoints,localdirecs,firstintersectpoints,firstintersectnormals,firstintersectdists)


            pass

        return (firstintersectpoints,firstintersectnormals,firstintersectdists)

        

    @classmethod
    def fromx3d(cls,filename,implpartparams=None,tol=1e-4,recalcnormals=False,defaultappearance=None):
        loader=x3d_indexedfaceset_loader.parsefile(filename)
        
        # For now, treat everything as a part
        surfaces=[]
        cnt=0
        for x3dshape in loader.surfaces:
            x3dsurface=x3dshape.geometry
            

            if x3dsurface is None:
                continue
            
            #textureurl=None
            #if x3dshape.appearance is not None and x3dshape.appearance.texture is not None and hasattr(x3dshape.appearance.texture,"url") and x3dshape.appearance.texture.url is not None:
            #    textureurl=x3dshape.appearance.texture.url
            #    pass

            # Meshed surface
            if isinstance(x3dsurface,x3d_indexedfaceset):
                # use implpartparams as intrinsicparameterizationparams
                surface=polygonalsurface.fromx3d(cnt,x3dsurface,x3dshape.appearance,implpartparams,tol=tol,defaultappearance=defaultappearance)
                pass
            else:
                # NURBS surface
                raise ValueError("FromX3D: Unknown geometry class %s" % (x3dsurface.__class__.__name__))
            

            surfaces.append(surface)
            cnt+=1
        return cls(tol=tol,implpartparams=implpartparams,surfaces=surfaces)

    @classmethod
    def fromvertices(cls,vertices,vertexidx,buildintrinsicparameterization=None,cadpartparams=None,appearance=None,solid=True,tol=1e-6,principal_curvatures=None,curvature_tangent_axes=None):
        surfaces=[]

        (vertexidx_indices,numvertices,missing_final_terminator) = find_vertexidx_indices(vertexidx)
        normals=calcnormals(vertices,vertexidx_indices,vertexidx,numvertices)
        
        surface=polygonalsurface.fromvertices(0,vertices,vertexidx_indices,vertexidx,numvertices,normals,None,buildintrinsicparameterization=buildintrinsicparameterization,cadpartparams=cadpartparams,appearance=appearance,solid=solid,tol=tol,principal_curvatures=principal_curvatures,curvature_tangent_axes=curvature_tangent_axes)

        surfaces.append(surface)
        return cls(tol=tol,implpartparams=cadpartparams,surfaces=surfaces)
        


    @classmethod
    def fromstl(cls,filename,implpartparams=None,tol=1e-4,recalcnormals=False,metersperunit=1.0,defaultappearance=None):
        facets=loadstl(filename)

        nfacets=facets.shape[0]

        facets[:,1:,:]*=metersperunit  # convert distances to meters
        # !!!*** Should use calcnormals() here 
        if recalcnormals or (facets[0,0,:]==np.zeros(3,dtype='d')).all():
            V = facets[:,2,:]-facets[:,1,:]
            W = facets[:,3,:]-facets[:,1,:]
            N = np.cross(V,W)
            N = N/vecnorm(N,axis=1).reshape(N.shape[0],1) # normalize

            facets[:,0,:]=N  # Store in facets
            pass

        # identify surfaces, either from parameter (JSON string)
        # or by algorithm 
        if implpartparams is None:
            (list_of_surface_facets,facetsbyvertexnum,uniquevertexcoords)=identify_surfaces(facets,tol)
            list_of_facets_surface_params = [ (facetnumarray,None) for facetnumarray in list_of_surface_facets ]
            pass
        else:
            non_numpy_list_of_facets_surface_params = json.loads(implpartparams)
            list_of_facets_surface_params =  [ (np.array(facetnums,dtype='d'),this_surface_params) for (facetnums,this_surface_params) in non_numpy_list_of_facets_surface_params ]

            # identify unique vertices
            (vertexnumbers,uniquevertexcoords) = uniquevertices(facets[:,1:,:].reshape(nfacets*3,3),tol)
            facetsbyvertexnum = vertexnumbers.reshape(nfacets,3) # This is a nfaces x 3 array of integers that gives unique vertex number of each of the three vertices of each triangle. 
            
            pass
        
        
        # m=3  # STL files are limited to triangles
        
        surfaces=[]
        for cnt in range(len(list_of_facets_surface_params)):
            (thesefacets,thesesurfparams)=list_of_facets_surface_params[cnt]
            vertexids = facetsbyvertexnum[thesefacets,:]
            vertexidx=np.concatenate((vertexids,np.ones((vertexids.shape[0],1),dtype=np.int32)*-1),axis=1).reshape(vertexids.shape[0]*4)
            vertexidx_indices=np.arange(0,vertexids.shape[0]*4,4,dtype=np.uint32)
            numvertices=3*np.ones(vertexidx_indices.shape[0],dtype=np.uint32)

            buildintrinsicparameterization=lambda surf,params: polygonalsurface_planarparameterization.new(surf,params)
            
            surface = polygonalsurface.fromvertices(cnt,
                                                    uniquevertexcoords,
                                                    vertexidx_indices,
                                                    vertexidx,
                                                    numvertices,
                                                    facets[thesefacets,0,:],
                                                    None,
                                                    buildintrinsicparameterization=buildintrinsicparameterization,
                                                    cadpartparams={"intrinsicparameterization":  thesesurfparams},
                                                    tol=tol,
                                                    appearance=defaultappearance)
            
            
            surfaces.append(surface)
            pass

        if implpartparams is None:
            # Create JSON definition of our parameterization
            non_numpy_list_of_facets_surface_params = []

            for cnt in range(len(list_of_facets_surface_params)):
                this_surface_params=surfaces[cnt].intrinsicparameterization.intrinsicparameterizationparams
                non_numpy_list_of_facets_surface_params.append((tuple(list_of_facets_surface_params[cnt][0].astype(object)),this_surface_params))
                pass
            
            implpartparams=json.dumps(non_numpy_list_of_facets_surface_params)
            pass
        
        
        
        return cls(tol=tol,implpartparams=implpartparams,surfaces=surfaces)

    @classmethod
    def fromDMobject(cls, srcobject, implpartparams=None, tol=1e-4, recalcnormals=False, metersperunit=1.0, defaultappearance=None):

        facets = np.empty((len(srcobject.faces[0].triangles), 4, 3))
        for i in range(0, len(srcobject.faces)):
            for j in range(0, len(srcobject.faces[i].triangles)):
                tri = srcobject.faces[i].triangles[j]
                facets[j, 0, :] = tri.facenormal
                facets[j, 1, :] = tri[0].point
                facets[j, 2, :] = tri[1].point
                facets[j, 3, :] = tri[2].point

        nfacets = facets.shape[0]

        facets[:, 1:, :] *= metersperunit  # convert distances to meters
        # !!!*** Should use calcnormals() here
        if recalcnormals or (facets[0, 0, :] == np.zeros(3, dtype='d')).all():
            V = facets[:, 2, :] - facets[:, 1, :]
            W = facets[:, 3, :] - facets[:, 1, :]
            N = np.cross(V, W)
            N = N / vecnorm(N, axis=1).reshape(N.shape[0], 1)  # normalize

            facets[:, 0, :] = N  # Store in facets
            pass

        # identify surfaces, either from parameter (JSON string)
        # or by algorithm
        if implpartparams is None:
            (list_of_surface_facets, facetsbyvertexnum, uniquevertexcoords) = identify_surfaces(facets, tol)
            list_of_facets_surface_params = [(facetnumarray, None) for facetnumarray in list_of_surface_facets]
            pass
        else:
            non_numpy_list_of_facets_surface_params = json.loads(implpartparams)
            list_of_facets_surface_params = [(np.array(facetnums, dtype='d'), this_surface_params) for
                                             (facetnums, this_surface_params) in
                                             non_numpy_list_of_facets_surface_params]

            # identify unique vertices
            (vertexnumbers, uniquevertexcoords) = uniquevertices(facets[:, 1:, :].reshape(nfacets * 3, 3), tol)
            facetsbyvertexnum = vertexnumbers.reshape(nfacets,
                                                      3)  # This is a nfaces x 3 array of integers that gives unique vertex number of each of the three vertices of each triangle.

            pass

        # m=3  # DMObjects are limited to triangles

        surfaces = []
        for cnt in range(len(list_of_facets_surface_params)):
            (thesefacets, thesesurfparams) = list_of_facets_surface_params[cnt]
            vertexids = facetsbyvertexnum[thesefacets, :]
            vertexidx = np.concatenate((vertexids, np.ones((vertexids.shape[0], 1), dtype=np.int32) * -1),
                                       axis=1).reshape(vertexids.shape[0] * 4)
            vertexidx_indices = np.arange(0, vertexids.shape[0] * 4, 4, dtype=np.uint32)
            numvertices = 3 * np.ones(vertexidx_indices.shape[0], dtype=np.uint32)

            buildintrinsicparameterization = lambda surf, params: polygonalsurface_planarparameterization.new(surf,
                                                                                                              params)

            surface = polygonalsurface.fromvertices(cnt,
                                                    uniquevertexcoords,
                                                    vertexidx_indices,
                                                    vertexidx,
                                                    numvertices,
                                                    facets[thesefacets, 0, :],
                                                    None,
                                                    buildintrinsicparameterization=buildintrinsicparameterization,
                                                    cadpartparams={"intrinsicparameterization": thesesurfparams},
                                                    tol=tol,
                                                    appearance=defaultappearance)

            surfaces.append(surface)
            pass

        if implpartparams is None:
            # Create JSON definition of our parameterization
            non_numpy_list_of_facets_surface_params = []

            for cnt in range(len(list_of_facets_surface_params)):
                this_surface_params = surfaces[cnt].intrinsicparameterization.intrinsicparameterizationparams
                non_numpy_list_of_facets_surface_params.append((tuple(list_of_facets_surface_params[cnt][0].astype(object)), this_surface_params))
                pass

            implpartparams = json.dumps(non_numpy_list_of_facets_surface_params)
            pass

        return cls(tol=tol, implpartparams=implpartparams, surfaces=surfaces)

    def X3DWrite(self,serializer,ourframe,destframe,UVparameterization=None):
        # return a list of lxml X3D <Shape> elements with texture coordinates
        # given according to UVparameterization.
        #
        # If UVparameterization is a single parameterization object, then
        # the shapes will have a single texcoord mapping.
        #
        # If UVparameterization is a list then the shapes
        # will have multitexture mappings
        #
        # x3dnamespace gives the X3D namespace to use for XML elements.
        # Likely values are:
        #   None  (For namespaceless X3D)
        #   x3dnamespace="http://www.web3d.org/specifications/x3d-namespace"
        #   x3dnamespace="http://www.web3d.org/specifications/x3d"
        #   x3dnamespace="http://www.web3d.org/specifications/x3d-3.2.xsd"
        # See also
        # http://www.web3d.org/specifications/x3d/PotentialX3dNamespace.html
        #

        
        #if not isinstance(appearances,collections.Sequence):
        #    appearance=appearances
        #    appearances=[appearance] * len(self.surfaces)  # repeat a single appearance
        #    pass

                
        shapes=[]
        cnt=0
        for surface in self.surfaces:

            
            surface.X3DWrite(serializer,ourframe,destframe,UVparameterization) #x3dnamespace=x3dnamespace,appearance=appearances[cnt])
            

            # shapes.append(shape)
            cnt+=1
            pass
        
        #return shapes
        pass

    def VRMLWrite(self,serializer,ourframe,destframe,UVparameterization=None): #,appearances=None):
        # return a list of VRML Shape { }  elements with texture coordinates
        # given according to UVparameterization.
        #
        # If UVparameterization is a single parameterization object, then
        # the shapes will have a single texcoord mapping.
        #
        # If UVparameterization is a list it is an error because
        # VRML does not support  multitexture mappings
        #
        # Appearance is instance of class VRMLAppearance or sequence of class VRMLAppearance (one per surface)
        
        #shapes=[]
        cnt=0

        #if not isinstance(appearances,collections.Sequence):
        #    appearance=appearances
        #    appearances=[appearance] * len(self.surfaces)  # repeat a single appearance
        #    pass
        
        
        for surface in self.surfaces:

            surface.VRMLWrite(serializer,ourframe,destframe,UVparameterization) #,appearance=appearances[cnt])
            
            
            #shapes.append(shape)
            cnt+=1
            pass
        
        #return shapes
        pass

    def assign_appearances(self,appearances):
        
        if not isinstance(appearances,collections.Sequence):
            # Just a single appearance... broadcast it over all surfaces
            appearances = [appearances]*len(self.surfaces)
            pass

        for cnt in range(len(self.surfaces)):
            self.surfaces[cnt].assign_appearance(appearances[cnt])
            pass
        pass
            
            

    
    pass

