import sys
import copy
import itertools

from lxml import etree

import numpy as np
try: 
    from numpy.linalg import svd
except ImportError:
    from scipy.linalg import svd
    pass
#import scipy
#import scipy.linalg

#try:
#    from cStringIO import StringIO  # python 2.x
#    pass
#except ImportError:
#    from io import StringIO # python 3.x
#    pass

from io import BytesIO

try: 
    from collections.abc import Sequence # python3
    pass
except ImportError:
    from collections import Sequence # python 2.x
    pass

from .surface import surface
from .polygonalsurface_planarparameterization import polygonalsurface_planarparameterization
from .polygonalsurface_texcoordparameterization import polygonalsurface_texcoordparameterization
from .appearance import simple_material
from .appearance import texture_url

#from ..ray import ray_to_plane_distance
#from .loaders.x3d import ReadX3DIndexedVertexField
from ..geometry import vertex_all_in_set


use_geom_accel=True

if use_geom_accel:
    from ..geometry_accel import polygon_intersects_box_3d
    pass
else:
    from ..geometry import polygon_intersects_box_3d
    pass

use_ps_accel=True
if use_ps_accel:
    from .polygonalsurface_accel import enclosed_or_intersecting_polygons_3d_accel
    
    pass

# plotted=False

def vecnorm(a,axis=-1):
    
    # compatibility for old versions of Numpy
    return np.sqrt(np.sum(a**2.0,axis))

def rms(a, axis=-1):

    if axis < 0:
        axis+=len(a.shape)
        pass
    
    return np.sqrt(np.sum(a**2.0,axis)/a.shape[axis])


def enclosed_or_intersecting_polygons_3d(polypool,vertexpool,vertices,vertexidx_indices,vertexidx,numvertices,inplanemats,facetnormals,box_v0,box_v1):
    #res=np.zeros(polypool.shape[0],dtype=np.bool)
    res=[]
    num_fully_enclosed=0
    
    vertexset=frozenset(vertexpool)
    for poolidx in np.arange(len(polypool)):
        idx=polypool[poolidx]

        if idx < 0:
            # masked out polygon
            continue
        
        firstidx=vertexidx_indices[idx]
        
        polygon_fully_enclosed = vertex_all_in_set(vertexidx[firstidx:(firstidx+numvertices[idx])],vertexset)  
        if polygon_fully_enclosed:
            res.append(idx)
            polypool[poolidx] = -1 # mask out polygon
            num_fully_enclosed+=1
            pass
        
        # if it's fully enclosed, nothing else need look at at, so we filter it here from the broader sibling pool
        if not polygon_fully_enclosed:
            polygon_vertices = vertices[vertexidx[firstidx:(firstidx+numvertices[idx])],:]
            # does it intersect?
            if polygon_intersects_box_3d(box_v0,box_v1,polygon_vertices,inplanemats[idx,:,:],facetnormals[idx,:]):
                res.append(idx)
                # Don't filter it out in this case because it must
                # intersect with a sibiling too 
                pass
            pass
        
        pass
    # return updated pool as array
    return (num_fully_enclosed,np.array(res,dtype=np.int32))


def calcnormals(vertices,vertexidx_indices,vertexidx,numvertices):
    # need to recalc normals
    
    #for polycnt in range(self.vertexidx_indices.shape[0]):
    #    firstidx=self.vertexidx_indices[cnt]
    #    termidx=firstidx+self.numvertices[cnt]
    
    
    
    V=vertices[vertexidx[vertexidx_indices+1],:]-vertices[vertexidx[vertexidx_indices+0],:]
    V=V/vecnorm(V,axis=1).reshape(V.shape[0],1)
            
    W=vertices[vertexidx[vertexidx_indices+2],:]-vertices[vertexidx[vertexidx_indices+0],:]
    W=W/vecnorm(W,axis=1).reshape(W.shape[0],1)
    
    N=np.cross(V,W)

    # If vector from 0th to 1st vertex and vector from 0th to 2nnd
    # vertex are too close to parallel, find another vertex
    
    min_cross_product = 1e-3
    tooparallel=np.where(vecnorm(N,axis=1) < min_cross_product)[0]

    for polycnt in tooparallel:
        third_vertex=3
        
        firstidx=vertexidx_indices[polycnt]
        
        while vecnorm(N[polycnt,:]) < min_cross_product:
            if third_vertex < numvertices[polycnt]:
                #  Try another vertex
                Wnew=vertices[vertexidx[firstidx+third_vertex],:]-vertices[vertexidx[firstidx],:]
            elif third_vertex == numvertices[polycnt]:
                # last possible try
                # In this case use vector from element 2 to element 1
                Wnew=vertices[vertexidx[firstidx+2],:]-vertices[vertexidx[firstidx+1],:]
                pass
            else:                        
                raise ValueError("Facet %d has extreme aspect ratio: Error calculating normal" % (polycnt))
            
            Wnew=Wnew/vecnorm(Wnew)
            N[polycnt,:]=np.cross(V[polycnt,:],Wnew)
            third_vertex+=1
            pass
        pass
    
    normals = N/vecnorm(N,axis=1).reshape(N.shape[0],1) # normalize
    
    
    return normals 



def FindBasis(coords,centroid,roughnormal):
    # coords should be n x 3 array of 3D coordinates
    # centroid should be mean 
    
    coordvals=(coords-centroid[np.newaxis,:]).T

    # calculate SVD
    (U,s,Vt)=svd(coordvals,full_matrices=True,compute_uv=True)
    # extract columns for 2d coordinate basis vectors
    # want columns x and y that correspond to the largest two
    # singular values and z that corresponds to the smallest singular value
    
    # We also want the x column cross the y column to give
    # the outward normal
    xcolindex=0
    ycolindex=1
    zcolindex=2
    
    
    # First, select colums to ensure zcolindex
    # corresponds to minimum singular value (normal direction)
    if abs(s[0]) < abs(s[1]) and abs(s[0]) < abs(s[2]):
        # element 0 is smallest s.v.
        xcolindex=2
        zcolindex=0
        pass
    if abs(s[1]) < abs(s[2]) and abs(s[1]) < abs(s[0]):
        # element 1 is smallest s.v.
        ycolindex=2
        zcolindex=1
        pass
    
    # Second, check to see if xcol cross ycol is in the
    # normal direction
    if np.dot(np.cross(U[:,xcolindex],U[:,ycolindex]),roughnormal) < 0:
        # x cross y is in wrong direction
        temp=xcolindex
        xcolindex=ycolindex
        ycolindex=temp
        pass

    To2D=U[:,np.array((xcolindex,ycolindex))].T # 2x3... Rows of To2D are x and y basis vectors, respectively

    # Get facet normal from the 3rd element of U,
    # flipping it to make sure it is pointing outward
    if np.inner(U[:,zcolindex],roughnormal) >= 0.0:
        normal = U[:,zcolindex]
        pass
    else:
        normal = -U[:,zcolindex]
        pass
    return (To2D,normal,np.array((s[xcolindex],s[ycolindex],s[zcolindex]),dtype='d'))


def find_vertexidx_indices(VertexIndex):

    # VertexIndex is an array of vertex indices, with -1 terminators
    
    # need to determine numvertices, terminatorindices, missing_final_terminator
    polynum=0
    vertexnum=0

    vertexidx_indices=[]
    
    numvertices=[]

    for cnt in range(VertexIndex.shape[0]):
        if vertexnum==0:
            # First vertex in a polygon
            vertexidx_indices.append(cnt)
            pass
        if VertexIndex[cnt]==-1:
            numvertices.append(vertexnum)
            polynum+=1
            vertexnum=0
            #terminatorindices.append(cnt)
            continue
        vertexnum+=1
        pass
    
    if vertexnum > 0:
        # if coordIndex did not end with -1, need to
        # do cleanup here
        numvertices.append(vertexnum)
        polynum+=1
        missing_final_terminator=True
        pass
    else:
        missing_final_terminator=False
        pass
    numpolys = polynum

    assert(numpolys==len(vertexidx_indices))
    
       
    return (np.array(vertexidx_indices,dtype=np.uint32),np.array(numvertices,dtype=np.uint32),missing_final_terminator)

class polygonalsurface(surface):
    # The rule in this implementation is that
    # surfaces are collections of polygon, should be
    # contiguous and the normal of any member polygon
    # should not exceed 80 degress away from the direction
    # of the mean of the normals of all polygons

    # ***!!! NOTE: This will be very memory-ineffcient if
    # there are small percentage of polygons with a LOT of vertices,
    # e.g. something with a cross section slice, where the slice
    # is a single polygon.
    # When converting to C++ probably
    # makes sense to change structure of vertexids field
    # to use separate field to identify polygon boundaries
    # rather than having a rectangular matrix
    
    # n polygons each with up to m vertices from p vertices total
    # q unique normals
    # Polygons presumed to be convex
    surfaceid = None # index number of this surface within this object

    
    vertices=None  # px3 array of vertex coordinates... must call invalidate() when modified
    vertexidx_indices=None # n integer array of indices into vertexidx. Identifies first entry corresponding to a particular polygon... must call invalidate() when modified
    vertexidx=None # integer array of vertex ids.... Each polygon terminated with -1 ... must call invalidate() when modified
    numvertices=None # integer array indicating number of vertices in each polygon
    #normals=None  # nx3 element array indicating nominal normal direction per facet... Should be auto-set by loader

    # vertexnormals... obsolete 
    normals = None # optional qx3 array indicating normal direction 
    normalidx = None # integer array of normal ids, for the normals to each vertex in each polygon, each polygon terminated by -1... or None, in which case we just have a single normal per polygon and q must match # of polygons
    
    boxes=None # Array of bounding boxes. Each box is identified by integers:
               #  8 box children, or -1, identified as indices into
               # this array, then an index (or -1) into boxpolys array
               # So this is a 9 x num_children array
               #
               # Note that current boxing is inefficient because polygons on the boundary of big boxes tend to get tested a lot
               # auto-set by buildboxes
    boxpolys=None # array of polygon indices for each leaf box,
                  # each terminated by -1
    boxcoords=None # Array of box coordinates, same structure as boxes, 6 values per box: 6 box bounds, Box bounds are (minx,miny,minz,maxx,maxy,maxz)

    refpoints=None # nx3 array of reference points (center points) -- auto-set by buildprojinfo()
    maxradius=None # n  maximum radius of vector from ref point -- autoset by buildprojinfo()
    # projection info -- used by imageprojection.pyx
    facetnormals=None   # This is probably redundant with normals, above, but the requirements may be a bit different, so 
                        # we keep it for the time being. This is calculated by SVD as perpendicular to the elements
                        # of inplanemats... autoset by buildprojinfo()
    refsizes=None       # autoset by buildprojinfo()
    inplanemats=None     # autoset by buildprojinfo() ... numpolys by 2 by 3 array of two orthogonal in-plane vectors in each polygon. The first vector cross the second vector should be in the facetnormal direction. 

    
    intrinsicparameterization=None
    parameterizations=None # list of additional parameterizations (not including intrinsicparameterization)
    cadpartparams=None
    solid=None # If True, this surface is part of a solid, so the back side does not need to be rendered 



    # Curvature information
    
    principal_curvatures = None # Optional per-vertex principal curvatures
    curvature_tangent_axes = None # Optional tangent directions corresponding to principal curvatures, per-vertex x 2 x 3.... the 2 is which vector, corresponding to the principal curvature, 3 corresponds to the component of the tangent vector.


    # Connectivity information
    polynum_by_vertex=None # object (list) array indexed by vertexnum of list of polynums
    adjacent_vertices=None # object (list) array indexed by vertexnum of list of adjacent vertexnums
    # NOTE WARNING: Avoid using sets as they are unordered and might not interate over in a consistent order. 


    
    
    # outofplanevec=None # outofplanevec is really facet normal
    #inplane2texcoords=None # use intrinsicparameterization.inplane2texcoords  or parameterization.inplane2texcoords

    # appearance=None  # cadpart.appearance object. Note: This is actually declared in surface.py
    
    tol=None

    def __init__(self,**kwargs):
        self.solid=True
        self.parameterizations=[]
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
        
            setattr(self,kwarg,kwargs[kwarg])
            pass

        pass



    def _eval_xyz_polygonuv(self,parameterization,polynum,u,v):
        # Evaluate (x,y,z) coordinates given a particular polygon
        # and the (u,v) coordinates
        #
        # Don't call this directly... use parameterization.eval_xy_uv()
        # instead, which figures out which polygon for you. 
        #
        # See also scope_coin3d.cpp:Get3DCoordsGivenTexPolygonCoords

        if parameterization is None:
            parameterization=self.intrinsicparameterization
            pass
        
        
        # (centroid,s,xcolindex,ycolindex,To2D, AijMat,AijMatInv)=parameterization._determine_tex_xform(polysurf,polysurf_polynum)
        #refsize=np.sqrt(s[xcolindex]**2.0 + s[ycolindex]**2.0)
        centroid=self.refpoints[polynum]
        # refsize=self.refsizes[polynum]
        To2D = self.inplanemats[polynum,:,:]
        AijMat = parameterization.inplane2texcoords[polynum,:,:]
        AijMatInv = parameterization.texcoords2inplane[polynum,:,:]

        #import pdb
        #pdb.set_trace()
        

        TexXYExt=np.array(((u-parameterization.lowerleft_meaningfulunits[0])/parameterization.meaningfulunits_per_texcoord[0],(v-parameterization.lowerleft_meaningfulunits[1])/parameterization.meaningfulunits_per_texcoord[1],1.0),dtype='d')

        # Note Capital UV represent a different (not u,v!) parameterization
        # of the in-plane 3D space of this facet. 
        TexUVExt = np.inner(AijMatInv,TexXYExt)
        TexUVExt /= TexUVExt[2] # normalize inhomogeneous coordinates

        # These coordinates of this facet are relative to its centroid
        TexUV = TexUVExt[:2]

        # Get 3D coordinates relative to centroid
        Tex3D = np.inner(To2D.T,TexUV)

        # Add to centroid
        Tex3DObject = Tex3D + centroid   # This gives 3D coordinates of the point in object coordinates 
        
        return Tex3DObject 

    
    def _buildbox(self,boxlist,boxcoordlist,polys,boxpolys,vertexpool,cnt,depth,minx,miny,minz,maxx,maxy,maxz):

        #import pdb
        #pdb.set_trace()
        # polys is an index array
        # boxpolys is accumulated array of where
        # polys is true. 
        thisbox=np.ones(9,dtype='i4')*-1

        #if 266241 in polys:
        #    firstidx=self.vertexidx_indices[266241]
        #    sys.stderr.write("got 266241... vertices=%s\n" % (str(self.vertices[self.vertexidx[firstidx:(firstidx+self.numvertices[266241])],:])))
        #    pass


        box_v0=np.array((minx,miny,minz),dtype='d')
        box_v1=np.array((maxx,maxy,maxz),dtype='d')
        
        # FIXME: Current algorithm is slow for lookups because we only box up
        # fully-enclosed polygons, so polygons which cross high-level boundaries
        # get tested A LOT!!!
        #
        # The fix for this is to identify intersecting 
        # polygons as well as fully enclosed polygons.
        # But these cannot be pulled from the list immediately
        # because they might intersect a sibling.
        #
        # Instead, if they intersect or are contained in this level
        # They should be flushed down to the levels underneath.
        # The key missing link is a box-polygon intersection algorithm
        # such as this one: 
        # https://github.com/erich666/GraphicsGems/blob/master/gemsv/ch7-2/pcube.c
        # (implementation of the above currently completed!)

        
        # filter down polys according to what is in this box
        if depth != 0: # all pass for depth = 0

            ## for debugging
            #debugpolys=copy.copy(polys)
            
            if use_ps_accel:
                (num_fully_enclosed,ourpolys)=enclosed_or_intersecting_polygons_3d_accel(polys,self.vertices,self.vertexidx_indices,self.vertexidx,self.numvertices,self.inplanemats,self.facetnormals,box_v0,box_v1)
                ourvertexpool=None


                pass
            else:
                vertices=self.vertices[vertexpool,:]
                
                
                # Filter down vertexpool from parent
                ourvertexpool=vertexpool[np.where(((vertices[:,0] >= minx) &
                                                   (vertices[:,0] <= maxx) &
                                                   (vertices[:,1] >= miny) &
                                                   (vertices[:,1] <= maxy) &
                                                   (vertices[:,2] >= minz) &
                                                   (vertices[:,2] <= maxz)))[0]]

                (num_fully_enclosed,ourpolys)=enclosed_or_intersecting_polygons_3d(polys,ourvertexpool,self.vertices,self.vertexidx_indices,self.vertexidx,self.numvertices,self.inplanemats,self.facetnormals,box_v0,box_v1)
                pass

            ## Filter down vertexpool from parent (debugging)
            if use_ps_accel:
                allvertexpool=np.where(((self.vertices[:,0] >= minx) &
                                        (self.vertices[:,0] <= maxx) &
                                        (self.vertices[:,1] >= miny) &
                                        (self.vertices[:,1] <= maxy) &
                                        (self.vertices[:,2] >= minz) &
                                        (self.vertices[:,2] <= maxz)))[0]
                pass
            else:
                allvertexpool=ourvertexpool
                pass
            ##debugpolys=np.arange(self.vertexidx_indices.shape[0])
            
            #(junk,allpolys)=enclosed_or_intersecting_polygons_3d(debugpolys,allvertexpool,self.vertices,self.vertexidx_indices,self.vertexidx,self.numvertices,self.inplanemats,self.facetnormals,box_v0,box_v1)
            #assert((ourpolys==allpolys).all())


            pass
        else:
            ourpolys=polys
            num_fully_enclosed=ourpolys.shape[0]
            ourvertexpool=vertexpool
            pass
        
        
        
        boxlist.append(thisbox)
        boxcoordlist.append((minx,miny,minz,maxx,maxy,maxz))

        newcnt=cnt+1

        if num_fully_enclosed > 10 and depth <= 18:
            # split up box
            distx=maxx-minx
            disty=maxy-miny
            distz=maxz-minz
            eps=1e-4*np.sqrt(distx**2 + disty**2 + distz**2)

            thisbox[0]=newcnt
            newcnt=self._buildbox(boxlist,boxcoordlist,ourpolys,boxpolys,ourvertexpool,newcnt,depth+1,minx,miny,minz,minx+distx/2.0+eps,miny+disty/2.0+eps,minz+distz/2.0+eps)
            thisbox[1]=newcnt
            newcnt=self._buildbox(boxlist,boxcoordlist,ourpolys,boxpolys,ourvertexpool,newcnt,depth+1,minx+distx/2.0-eps,miny,minz,maxx,miny+disty/2.0+eps,minz+distz/2.0+eps)
            thisbox[2]=newcnt
            newcnt=self._buildbox(boxlist,boxcoordlist,ourpolys,boxpolys,ourvertexpool,newcnt,depth+1,minx,miny+disty/2.0-eps,minz,minx+distx/2.0+eps,maxy,minz+distz/2.0+eps)
            thisbox[3]=newcnt
            newcnt=self._buildbox(boxlist,boxcoordlist,ourpolys,boxpolys,ourvertexpool,newcnt,depth+1,minx+distx/2.0-eps,miny+disty/2.0-eps,minz,maxx,maxy,minz+distz/2.0+eps)
            thisbox[4]=newcnt
            newcnt=self._buildbox(boxlist,boxcoordlist,ourpolys,boxpolys,ourvertexpool,newcnt,depth+1,minx,miny,minz+distz/2.0-eps,minx+distx/2.0+eps,miny+disty/2.0+eps,maxz)
            thisbox[5]=newcnt
            newcnt=self._buildbox(boxlist,boxcoordlist,ourpolys,boxpolys,ourvertexpool,newcnt,depth+1,minx+distx/2.0-eps,miny,minz+distz/2.0-eps,maxx,miny+disty/2.0+eps,maxz)
            thisbox[6]=newcnt
            newcnt=self._buildbox(boxlist,boxcoordlist,ourpolys,boxpolys,ourvertexpool,newcnt,depth+1,minx,miny+disty/2.0-eps,minz+distz/2.0-eps,minx+distx/2.0+eps,maxy,maxz)
            thisbox[7]=newcnt
            newcnt=self._buildbox(boxlist,boxcoordlist,ourpolys,boxpolys,ourvertexpool,newcnt,depth+1,minx+distx/2.0-eps,miny+disty/2.0-eps,minz+distz/2.0-eps,maxx,maxy,maxz)
            pass
        #
        else:
            # This is a leaf node
            # Record our polygons... These are those which are
            # fully enclosed or intersecting.
            # The index where they start is boxpolys[8]


            thisbox[8]=len(boxpolys)
            boxpolys.extend(ourpolys[ourpolys >= 0])
            boxpolys.append(-1)
            pass
        
        return newcnt

    def invalidateboxes(self):
        self.boxes=None
        self.boxcoords=None
        self.boxpolys=None
        pass
    

    def buildboxes(self):
        if self.boxes is not None:
            return
        
        self.buildprojinfo() # need inplanemats

        boxlist=[]
        boxcoordlist=[]
        # polys: indices of polygons which are inside this box         
        polys=np.arange(self.vertexidx_indices.shape[0],dtype=np.int32)  # -1 will be used as a flag to indicate this poly is no longer present

        if use_ps_accel:
            vertexpool=None # not needed
            pass
        else:
            vertexpool=np.arange(self.vertices.shape[0],dtype=np.uint32)
            pass

                      
        boxpolys=[]
        boxcnt=0

        (minx,miny,minz)=np.min(self.vertices,axis=0)
        (maxx,maxy,maxz)=np.max(self.vertices,axis=0)

        distx=maxx-minx
        disty=maxy-miny
        distz=maxz-minz
        eps=1e-4*np.sqrt(distx**2 + disty**2 + distz**2)

        self._buildbox(boxlist,boxcoordlist,polys,boxpolys,vertexpool,0,0,minx-eps,miny-eps,minz-eps,maxx+eps,maxy+eps,maxz+eps)
        
        self.boxes=np.array(boxlist,dtype='i4')
        self.boxcoords=np.array(boxcoordlist,dtype='d')
        self.boxpolys=np.array(boxpolys,dtype='i4')
        
        pass

    def eval_all_curvatures_from_mesh(self):
        # Evaluate and return curvatures at all vertices.
        # Does not assign to curvature variables... you must
        # call assigncurvature() to do that!

        curvatures = np.empty((self.vertices.shape[0],2),dtype='d')
        curvature_tangent_axes = np.empty((self.vertices.shape[0],2,3),dtype='d')

        for vertnum in range(self.vertices.shape[0]):
            (curvatures[vertnum,:],curvature_tangent_axes[vertnum,:,:])=self._evalcurvature(vertnum)
            pass

        return (curvatures,curvature_tangent_axes)
    
    def add_curvatures_to_surface_from_mesh(self):
        (curvatures,curvature_tangent_axes)=self.eval_all_curvatures_from_mesh()
        self.assign_curvatures(curvatures,curvature_tangent_axes)
        pass

    
    def _evalcurvature(self,vertexnum):
        # Evaluate local curvature at a vertex, based on
        # the local mesh
        

        if self.adjacent_vertices is None:
            self.buildconnectivity()
            pass

        #if vertexnum==519 or vertexnum==520 or vertexnum==9943:
        #    import pdb
        #    pdb.set_trace()


        #if vertexnum==8304:
        #    import pdb
        #    pdb.set_trace()

        
        vertexidxs=list(set([vertexnum]) | set(self.adjacent_vertices[vertexnum]))
        # Find best-fit plane for local vertices
        centroid=np.mean(self.vertices[vertexidxs,:],0)

        coordvals = (self.vertices[vertexidxs,:]-centroid[np.newaxis,:]).T # coordvals is the coordinates relative to centroid, 3 x numpoints

        # get roughnormal from 1st facet
        firstpolynum=self.polynum_by_vertex[vertexnum][0]
        if self.normalidx is None:
            roughnormal=self.normals[firstpolynum]
            pass
        else:
            firstidx=self.vertexidx_indices[firstpolynum]
            numvertices=self.numvertices[firstpolynum]

            polyvertices=self.vertexidx[firstidx:(firstidx+numvertices)]         
            polyvertexidx=np.where(polyvertices==vertexnum)[0][0]
            
            roughnormal=self.normals[self.normalidx[firstidx+polyvertexidx],:]
            pass
        
        (To2D,basisnormal,s)=FindBasis(self.vertices[vertexidxs,:],centroid,roughnormal)
        # We now have a basis. Fit local coords to a polynomial,
        # then evaluate the polynomial curvature using the usual methods
        #
        # To2D is 2x3 and rows are x and y basis vectors

        # Vertex2D coords are 2D coords relative to vertex we are interested in
        # ... 2 rows, representing u,v directions, by n columns, representing vertices

        # broadervertexidxs includes vertices two hops away
        broadervertexidxs=set(vertexidxs)
        for nearvertexnum in vertexidxs:
            broadervertexidxs.update(self.adjacent_vertices[nearvertexnum])
            pass
        broadervertexidxs=list(broadervertexidxs) # convert set to list for consistent ordering

        Vertex2DCoords=np.inner(To2D,self.vertices[broadervertexidxs,:]-self.vertices[vertexnum,:])
        VertexZCoords=np.inner(basisnormal,self.vertices[broadervertexidxs,:]-self.vertices[vertexnum,:])

        # obtain all combinations of powers xpow==0..2, ypow==0..2
        xypows = list(itertools.product(range(3),range(3)))
        # each element of xypows is an (xpow,ypow)

        polyfitmat = np.zeros((VertexZCoords.shape[0],len(xypows)),dtype='d')
        for col in range(len(xypows)):
            (xpow,ypow)=xypows[col]
            polyfitmat[:,col]=Vertex2DCoords[0,:]**xpow*Vertex2DCoords[1,:]**ypow
            
            pass

        # Perform surface fitting
        (xypowcoeffs,residuals,rank,s)=np.linalg.lstsq(polyfitmat,VertexZCoords)
        
        # Now surface approximated by sum(xypowcoeffs*(x**xpow)*(y**ypow))
        # e.g. xpows=np.array([xpow for (xpow,ypow) in xypows ])
        #      ypows=np.array([ypow for (xpow,ypow) in xypows ])
        # Recon = np.sum(xypowcoeffs[:,np.newaxis]*Vertex2DCoords[0:1,:]**xpows[:,np.newaxis]*Vertex2DCoords[1:2,:]**ypows[:,np.newaxis],axis=0)
        #global plotted
        
        if False: # not plotted:
            xpows=np.array([xpow for (xpow,ypow) in xypows ])
            ypows=np.array([ypow for (xpow,ypow) in xypows ])
            Recon = np.sum(xypowcoeffs[:,np.newaxis]*Vertex2DCoords[0:1,:]**xpows[:,np.newaxis]*Vertex2DCoords[1:2,:]**ypows[:,np.newaxis],axis=0)

            import matplotlib.pyplot as pl
            from mpl_toolkits.mplot3d import Axes3D
            
            fig=pl.figure()
            ax=fig.add_subplot(111,projection='3d')
            ax.scatter(Vertex2DCoords[0,:],Vertex2DCoords[1,:],VertexZCoords,c=np.array(((0.0,0.0,1.0),)))
            ax.scatter(Vertex2DCoords[0,:],Vertex2DCoords[1,:],Recon,c=np.array(((0.0,1.0,0.0),)))
            
            pl.show()
            pl.drawnow()
            #ax.plot_surface(Vertex2DCoords[0,np.newaxis,:],Vertex2DCoords[1,:
            plotted=True
            pass
        
        # If we wanted more accurate, we could probably iterate the
        # fitting process with an adjusted normal, etc.
        # Sigma u is the derivative

        # This quadratic fit is a parametric surface
        # X(u,v) where X has x, y, and z coefficients
        # w.o.l.o.g, let x = u and y=v

        # We evaluate at the origin (relative to our chosen vertex) in the To2D frame
        
        # Sigma_u is dX/du
        Sigma_u = np.array((1.0,0.0,xypowcoeffs[xypows.index((1,0))]),dtype='d')
        # Sigma_v is dX/dv
        Sigma_v = np.array((0.0,1.0,xypowcoeffs[xypows.index((0,1))]),dtype='d')
        # Sigma_uu is d^2X/du^2
        Sigma_uu = np.array((0.0,0.0,2.0*xypowcoeffs[xypows.index((2,0))]),dtype='d')
        # Sigma_vv is d^2X/dv^2
        Sigma_vv = np.array((0.0,0.0,2.0*xypowcoeffs[xypows.index((0,2))]),dtype='d')
        # Sigma_uv is d^2X/dudv
        Sigma_uv = np.array((0.0,0.0,xypowcoeffs[xypows.index((1,1))]),dtype='d')

        # Now follow curvature evaluation process
        # see addcurvatures.py / evalcurvature.py for more details
        
        # First fundamental form... (metric tensor coeffs)
        E = np.inner(Sigma_u,Sigma_u)
        F = np.inner(Sigma_u,Sigma_v)
        G = np.inner(Sigma_v,Sigma_v)

        # Second fundamental form:
        # ... The normal in the To2D frame is (0,0,1)
        NormalInTo2DFrame = np.array((0,0,1.0),dtype='d')
        L = np.inner(Sigma_uu,NormalInTo2DFrame)
        M = np.inner(Sigma_uv,NormalInTo2DFrame)
        N = np.inner(Sigma_vv,NormalInTo2DFrame)

        # Evaluate shape operator
        S = (1.0/(E*G-F**2.0)) * np.array(((L*G-M*F,M*G-N*F),
                                           (M*E-L*F,N*E-M*F)),dtype='d') #  (see wikipedia + https://math.stackexchange.com/questions/536563/trouble-computing-the-shape-operator) http://mathworld.wolfram.com/WeingartenEquations.html

        # S should be symmetric since basis is orthogonal
        # ... evaluate asymmetry
        asymmetry = S[1,0]-S[0,1]

        if asymmetry > .35*rms(np.ravel(S)):
            sys.stderr.write("polygonalsurface.py:_evalcurvature: Warning: Shape matrix is substantally asymmetric. vertexnum = %d; S = %s\n" % (vertexnum,str(S)))
            pass
        
        # correct asymmetry
        S[1,0] -= asymmetry/2.0
        S[0,1] += asymmetry/2.0

        (curvatures,evects) = np.linalg.eig(S)
        curvaturetangents=np.dot(To2D.T,evects)

        # We don't want the eigenframe to be mirrored relative to the (U,V)
        # frame, for consistency in interpreting positive vs. negative curvature.
        # ...  so if the dot/inner product of (UxV) with (TANGENT0xTANGENT1)
        # is negative, that indicates mirroring 
        # Negating one of the eigenvectors will un-mirror it. 

        if np.inner(np.cross(To2D[0,:],To2D[1,:]),np.cross(curvaturetangents[:,0],curvaturetangents[:,1])) < 0.0:
            curvaturetangents[:,0]=-curvaturetangents[:,0]
            pass

        return (curvatures,curvaturetangents.transpose(1,0))
        

    def invalidateconnectivity(self):
        self.polynum_by_vertex=None
        self.adjacent_vertices=None
        pass

    def buildconnectivity(self):
        # build polynum_by_vertex

        if self.polynum_by_vertex is not None:
            return

        
        self.polynum_by_vertex=np.empty(self.vertices.shape[0],dtype='O')
        for vertcnt in range(self.vertices.shape[0]):
            self.polynum_by_vertex[vertcnt]=[]
            pass

        for polycnt in range(self.vertexidx_indices.shape[0]):
            firstidx=self.vertexidx_indices[polycnt]
            for vertexnum in self.vertexidx[firstidx:(firstidx+self.numvertices[polycnt])]:
                self.polynum_by_vertex[vertexnum].append(polycnt)
                pass
            pass

        # build adjacent_vertices array
        self.adjacent_vertices=np.empty(self.vertices.shape[0],dtype='O')
        for vertcnt in range(self.vertices.shape[0]):
            vertexidxset=set([])
            for adjacent_poly in self.polynum_by_vertex[vertcnt]:
                firstidx=self.vertexidx_indices[adjacent_poly]
                numvertices=self.numvertices[adjacent_poly]
                polyvertices=self.vertexidx[firstidx:(firstidx+numvertices)]         
                polyvertexidx=np.where(polyvertices==vertcnt)[0][0]
                nextpolyvertex = polyvertices[(polyvertexidx + 1) % numvertices]
                prevpolyvertex = polyvertices[(polyvertexidx+numvertices-1) % numvertices]
                vertexidxset.add(nextpolyvertex)
                vertexidxset.add(prevpolyvertex)
                pass
            self.adjacent_vertices[vertcnt]=list(vertexidxset)

            pass
        pass

    def invalidateprojinfo(self):
        self.refpoints=None
        self.facetnormals = None
        self.inplanemats=None
        self.outofplanevec=None

        self.refsizes=None
        self.maxradius=None

        self.intrinsicparameterization.invalidateprojinfo()
        
        for parameterization in self.parameterizations:
            parameterization.invalidateprojinfo()
            pass

        pass
        
    
    def buildprojinfo(self):
        # see also polygonalsurface_texcoordparameterization.py/_determine_tex_xform()
        # and following steps in polygonalsurface_texcoordparameterization.py:buildprojinfo()

        # See also scope_coin3d.cpp:DetermineTexXform
        # See also polygonalsurface.py: buildprojinfo()

        # NOTE: AijMat is 2x3, but AijMatInv is 3x3 for historical reasons 

        # OK. So we extract the numpoints points. These are in
        # 3-space and nominally on a plane. How to flatten them
        # into 2D?
        # 1. Subtract out the centroid (best-fit plane
        #    will pass through centroid)
        #    (Also subtract same location from ObjectSpace intersection
        #    point coordinates, from above)
        # 2. Assemble point vectors into a matrix of vectors from centroid
        # 3. Apply SVD. Resulting two basis vectors corresponding to
        #    largest-magnitude two singular values will span the plane
        # 4. Re-evaluate point vectors and intersection location in
        #    terms of these basis vectors. Now our point coordinates
        #    and intersection coordinates are in 2-space. 

        if self.refpoints is not None:
            return # already built
        
        numpolys=self.vertexidx_indices.shape[0]
        
        self.refpoints=np.zeros((numpolys,3),dtype='d')  # refpoints are centroids of the polygons
        self.facetnormals = np.zeros((numpolys,3),dtype='d')
        self.inplanemats=np.zeros((numpolys,2,3),dtype='d')
        self.outofplanevec=np.zeros((numpolys,2,3),dtype='d')

        self.refsizes=np.zeros(numpolys,dtype='d')
        self.maxradius=np.zeros(numpolys,dtype='d')
        
        for polynum in range(numpolys):
            numpoints=self.numvertices[polynum] #np.count_nonzero(self.vertexids[polynum,:] >= 0)
            firstidx=self.vertexidx_indices[polynum]
            centroid = np.mean(self.vertices[self.vertexidx[firstidx:(firstidx+numpoints)],:],axis=0)
            self.refpoints[polynum,:]=centroid

            if self.normalidx is None:
                # single normal per facet
                normalvec = self.normals[polynum]
                pass
            else:
                # normal per vertex
                vertexnormals = self.normals[self.normalidx[firstidx:(firstidx+numpoints)],:]
                # Average the normals from the vertices
                normalvec = np.mean(self.normals,axis=0)
                normalvec /= np.linalg.norm(normalvec)  # renormalize
                pass

            (To2D,basisnormal,s)=FindBasis(self.vertices[self.vertexidx[firstidx:(firstidx+numpoints)],:],centroid,normalvec)
            
            self.facetnormals[polynum,:]=basisnormal
            self.inplanemats[polynum,:,:]=To2D

            
            self.refsizes[polynum]=np.sqrt(s[0]**2.0 + s[1]**2.0)
            self.maxradius[polynum] = np.max(vecnorm(self.vertices[self.vertexidx[firstidx:(firstidx+numpoints)],:]-centroid.reshape(1,3),axis=1),axis=0)
            
            
            pass

        if self.intrinsicparameterization is not None:
            self.intrinsicparameterization.buildprojinfo(self)
            pass
        
        for parameterization in self.parameterizations:
            parameterization.buildprojinfo(self)
            pass
        
        
        pass

        
    
    @classmethod
    def fromvertices(cls,surfaceid,vertices,vertexidx_indices,vertexidx,numvertices,normals,normalidx,buildintrinsicparameterization=None,cadpartparams=None,appearance=None,solid=True,tol=1e-6,principal_curvatures=None,curvature_tangent_axes=None):
        
    
        n=vertexidx_indices.shape[0]
        # vertexids is (number of facets, max_vertices_per_facet)
        
        # Evaluate refpoints from average of polygon boundary locations
        # refpoints now built by buildprojinfo()
        refpoints=np.zeros((n,3),dtype='d')
        navg = np.zeros(n,dtype=np.int32)
        #m=vertexids.shape[1]
        # Evaluate maxradius
        maxradius = np.zeros(n,dtype='d')

        for cnt in range(n):
            firstidx=vertexidx_indices[cnt]
            followingidx=firstidx+numvertices[cnt]
            #haveanothervertex=vertexids[:,cnt] >= 0
            #refpoints[haveanothervertex] += vertices[vertexids[haveanothervertex,cnt],:]
            #navg[haveanothervertex] += 1
            refpoints[cnt,:] = np.mean(vertices[vertexidx[firstidx:followingidx],:],axis=0)
            maxradius[cnt] = np.max(vecnorm(refpoints[cnt:(cnt+1),:]-vertices[vertexidx[firstidx:followingidx],:],axis=1))
            pass
        
        # refpoints /= navg.reshape(n,1)

        
        #for cnt in range(m):
        #    haveanothervertex=vertexids[:,cnt] >= 0
        #    radius=vecnorm(refpoints[haveanothervertex]-vertices[vertexids[haveanothervertex,cnt],:],axis=1)
        #    biggerradius=np.where(radius > maxradius[haveanothervertex])[0]
        #    maxradius[haveanothervertex][biggerradius]=radius[biggerradius]
        #    pass


        # refpoints NOT provided.... build by buildprojinfo()
        retval=  cls(surfaceid=surfaceid,
                     vertices=vertices,
                     maxradius=maxradius,
                     vertexidx_indices=vertexidx_indices,
                     vertexidx=vertexidx,
                     numvertices=numvertices,
                     normals=normals,
                     normalidx=normalidx,
                     cadpartparams=cadpartparams,
                     appearance=appearance,
                     solid=solid,
                     tol=tol,
                     principal_curvatures=principal_curvatures,
                     curvature_tangent_axes=curvature_tangent_axes)

        #if buildintrinsicparameterization is None:
        #    buildintrinsicparameterization=polygonalsurface_planarparameterization.new
        #    pass

        if buildintrinsicparameterization is not None: 
            retval.intrinsicparameterization = buildintrinsicparameterization(retval,cadpartparams)
            # retval.intrinsciparameterizationparams=retval.intrinsicparameterization.cadpartparams
            pass
        
        retval.buildboxes()
        retval.buildprojinfo()
        
        return retval
        
    
    def ray_intersections_OBSOLETE(self,localpoints,localdirecs,firstintersectpoints=None,firstintersectnormals=None,firstintersectdists=None):
        # NOTE: Untested and not a very good implementation...
        # See imageprojection.pyx for something better
        
        # rays is a set of n rays, each of which is specified by 6 coordinates, source point then unit direction vector, 

        # n here is number of rays, m is number of polygons
        n=localpoints.shape[0]
        m=self.refpoints.shape[0]


        
        dists=ray_to_plane_distance(localpoints,localdirecs,self.refpoints,self.normals,self.normalidx)
        # dists is (number of rays,number of polygons) matrix of distances
        # points will be (number of rays, number of polygons, 3) array of
        # intersection point coordinates


        # intersectpoints is (numrays, numpolygons, 3) 
        intersectpoints = localpoints.reshape(n,1,3) + dists*localdirecs.reshape(n,1,3)

        # Mark intersections behind us as at infinity -- i.e. no intersection
        dists[dists < 0] = np.Inf

        # radii is (numrays, numpolygons)
        radii = vecnorm(intersectpoints - self.refpoints.reshape(1,m,3),axis=2)

        # Mark intersection outside maxradius at infinity -- i.e. no interaction
        dists[radii > self.maxradius.reshape(1,m)+self.tol]=np.inf
        

        # for each ray, determine if intersectpoint is inside the polygon
        # Check by making sure the point is laying on the same
        # side of each polygon segment when traversed
        # ... do this by taking cross product of polygon border segment
        # with vector from point to one of the border segment vertices.
        #
        #    ----------
        #    |  x     |
        #    |        /
        #    \       /
        #     \     /
        #      \---/
        #  Cross product of clockwise segment with vector from point to
        # segment endpoint is always out-of-screen
        # or consistent sign in general case, except when point
        # is on segment or on extension of segment.

        # Center point of polygon is refpoint.
        # So we push each vertex tolerance away from refpoint so that
        # we make polygon slightly larger and therefore err on the side of
        # affirming an intersection

        userays = np.isfinite(dists).sum(1) > 0  # count number of finite distances, identify where that count is > 0
        #usepolys = np.isfinite(dists).sum(0) > 0 #

        if firstintersectpoints is None:
            firstintersectpoints=np.empty((n,3),dtype='d')
            firstintersectpoints.fill(np.inf)
            pass

        if firstintersectnormals is None:
            firstintersectnormals=np.zeros((n,3),dtype='d')
            pass

        if firstintersectdists is None:
            
            firstintersectdists=np.empty((n,),dtype='d')
            firstintersectdists.fill(np.inf)
            pass

        
        
        for raycnt in range(np.count_nonzero(userays)):
            raydists=dists[userays,:][raycnt,:]
            polyorder = np.argsort(raydists)

            # trim infinite distances from polyorder
            polyorder=polyorder[:np.count_nonzero(np.isfinite(raydists[polyorder]))]
            for polycnt in range(len(polyorder)):
                # determine if intersectpoints[userays,polyorder[polycnt],:][raycnt,:]
                # is inside polygon # polyorder[polycnt]

                pointcoords=intersectpoints[userays,polyorder[polycnt],:][raycnt,:]
                # sign=None

                firstidx=self.vertexidx_indices[polyorder[polycnt]]
                num_vertices=self.numvertices[polyorder[polycnt]]
                vertices=self.vertexidx[firstidx:(firstidx+num_vertices)]

                verticescoords=self.vertices[vertices[:num_vertices],:] # get vertex coordinates -- (num_vertices x 3)

                # Add tolerance to point_to_edge to avoid missing anything in the seams
                
                
                

                nextverticescoords=np.roll(verticescoords,1,axis=0)
                
                edge_vectors=verticescoords-nextverticescoords
                point_to_edge = verticescoords-pointcoords.reshape(1,3)

                # NOTE: ***!!!FIXME *** Should use normalidx here when appropriate
                
                signs = np.sign(np.inner(np.cross(edge_vectors,point_to_edge),self.normals[polyorder[polycnt]]))   # length num_vertices
                if (0.0 in signs # on boundary.... count as outside because of tolerance trick above
                    or not((signs==signs[0]).all())):
                    # mark distance as infinite i.e. no ray intersection
                    dists[userays,polyorder[polycnt]][raycnt] = np.inf
                    pass
                else:
                    # have a ray intersection... farther distances don't
                    # matter because we already intersected the ray
                    dists[userays,polyorder[(polycnt+1):]][raycnt,:] = np.inf
                    firstintersectpoints[userays,:][raycnt,:]=pointcoords
                    firstintersectnormals[userays,:][raycnt,:]=self.normals[polyorder[polycnt],:]
                    firstintersectdists[userays][raycnt]=dists[userays,polyorder[polycnt]][raycnt]
                    

                    break
                    
                    
                pass
            pass
        
        return (firstintersectpoints,firstintersectnormals,firstintersectdists)


    def _generate_texcoord(self,NSPRE,NSMAP,UVparameterization):
        texcoordbuf=BytesIO()
        texcoord=etree.Element(NSPRE+"TextureCoordinate",nsmap=NSMAP)
            
        for cnt in range(self.vertexidx_indices.shape[0]):
            #numvertices=np.count_nonzero(self.vertexids[cnt,:] >= 0) # How many vertices does this polygon have?

            
            
            #iuvcoords=self.intrinsicparameterization.eval_uv(self,self.vertexids[cnt,:numvertices])
            if UVparameterization is None:
                # UVparameterization of None means use intrinsic values
                UVparameterization = self.intrinsicparameterization
                pass
            # texcoords=UVparameterization.eval_texcoords_from_intrinsic_texcoords(iuvcoords)
            
            texcoords=UVparameterization.eval_texcoord_polygonvertex(self,cnt,range(self.numvertices[cnt]))
            
            
            
            np.savetxt(texcoordbuf,texcoords,delimiter=" ",newline=" ")
            
            pass
        texcoord.attrib["point"]=texcoordbuf.getvalue().decode('utf-8')
        
        return texcoord
    
    def X3DWrite(self,serializer,ourframe=None,destframe=None,UVparameterization=None): # ,x3dnamespace=None,appearance=None):
        shape=etree.Element(serializer.NSPRE+"Shape",nsmap=serializer.NSMAP)
        if self.appearance is not None:
            shape.append(self.appearance.generate_x3d(serializer))
            pass
        
        indexedfaceset=etree.Element(serializer.NSPRE+"IndexedFaceSet",nsmap=serializer.NSMAP)
        
        # Write coordIndex as attribute of IndexedFaceSet
        coordindexbuf=BytesIO()
        for cnt in range(self.vertexidx_indices.shape[0]):
            firstvertex=self.vertexidx_indices[cnt]
            #numvertices=np.count_nonzero(self.vertexids[cnt,:] >= 0)
            np.savetxt(coordindexbuf,self.vertexidx[firstvertex:(firstvertex+self.numvertices[cnt])],fmt="%d",delimiter=" ",newline=" ")
            coordindexbuf.write(b" -1 ") # barrier between polygons
            pass
        indexedfaceset.attrib["coordIndex"]=coordindexbuf.getvalue().decode('utf-8')

        # indexedfaceset.attrib["creaseAngle"]="0.01" # If you don't want creases, set per-vertex normals
        indexedfaceset.attrib["solid"]=str(self.solid).lower()

        # Add <Coordinate> tag
        coordinate=etree.Element(serializer.NSPRE+"Coordinate",nsmap=serializer.NSMAP)

        # Add point= attribute of <Coordinate> element
        coordbuf=BytesIO()
        np.savetxt(coordbuf,ourframe.transformto(destframe,self.vertices),delimiter=" ",newline=" ")
        coordinate.attrib["point"]=coordbuf.getvalue().decode('utf-8')
        indexedfaceset.append(coordinate)

        # Need to generate normals
        normal=etree.Element(serializer.NSPRE+"Normal",nsmap=serializer.NSMAP)
        normalbuf=BytesIO()
        normalindexbuf=BytesIO()
        
        if self.normalidx is None:
            # Single normal per polygon
            assert(self.normals.shape[0]==self.vertexidx_indices.shape[0]) # consistency...
            
            # Evaluate all normals
            np.savetxt(normalbuf,ourframe.transformto(destframe,self.normals,inhomogeneous_vector=True),delimiter=" ",newline=" ")
            
            for cnt in range(self.normals.shape[0]):
                normalindexbuf.write(((str(cnt)+" ") * self.numvertices[cnt]).encode('utf-8')) # Write out normal index, numvertices times
                normalindexbuf.write(b" -1 ") # barrier between polygons
                pass
            
            pass
        else:
            # Seperate normal per vertex
            #assert(self.vertexnormals.shape[0]==self.vertexidx_indices.shape[0]) # consistency...
            normalpos=0
        
            for cnt in range(self.vertexidx_indices.shape[0]):
                #numvertices=np.count_nonzero(self.vertexids[cnt,:] >= 0) # How many vertices does this polygon have?
                firstvertex=self.vertexidx_indices[cnt]
                
                np.savetxt(normalbuf,ourframe.transformto(destframe,self.normals[self.normalidx[firstvertex:(firstvertex+self.numvertices[cnt])],:],inhomogeneous_vector=True),delimiter=" ",newline=" ")

                np.savetxt(normalindexbuf,np.arange(normalpos,normalpos+self.numvertices[cnt],dtype='i32'),fmt="%d",delimiter=" ",newline=" ")
                
                normalpos+=self.numvertices[cnt]
                pass
            
            pass
        indexedfaceset.attrib["normalIndex"] = normalindexbuf.getvalue().decode('utf-8')
        normal.attrib["vector"]=normalbuf.getvalue().decode('utf-8')
        indexedfaceset.append(normal)

        
        
        # Need to generate texCoordIndex attribute
        # and TextureCoordinate element
        # Need to generate normals
        
        # TODO: ***!!! Should be able to pass on tex coord
        # indexing from polygonalsurface_texcoordparameterization,
        # But this will only work with multitextures if
        # all of the parameterizations share the same texcoordindex 

        texcoord=None
        if isinstance(UVparameterization,Sequence):  # User provided a list
            texcoord=etree.Element(serializer.NSPRE+"MultiTextureCoordinate",nsmap=serializer.NSMAP)
            for UVparam in UVparameterization:
                subtexcoord=self._generate_texcoord(serializer.NSPRE,serializer.NSMAP,UVparam)
                texcoord.append(subtexcoord)
                pass
            
            pass
        elif UVparameterization is not None or self.intrinsicparameterization is not None:  # Single parameterization
            texcoord=self._generate_texcoord(serializer.NSPRE,serializer.NSMAP,UVparameterization)
            pass
        if texcoord is not None:
            indexedfaceset.append(texcoord)
            texcoordindexbuf=BytesIO()
            
            texcoordpos=0
            for cnt in range(self.vertexidx_indices.shape[0]):
                #numvertices=np.count_nonzero(self.vertexids[cnt,:] >= 0) # How many vertices does this polygon have?
                np.savetxt(texcoordindexbuf,np.arange(texcoordpos,texcoordpos+self.numvertices[cnt],dtype='i4'),fmt="%d",delimiter=" ",newline=" ")
                texcoordindexbuf.write(b" -1 ")
                texcoordpos+=self.numvertices[cnt]
                pass
            
            indexedfaceset.attrib["texCoordIndex"] = texcoordindexbuf.getvalue().decode('utf-8')
            pass
        
        # add curvature data, if present
        if self.principal_curvatures is not None and self.curvature_tangent_axes is not None:
            principal_curvatures=etree.Element(serializer.spatialnde_NSPRE+"PrincipalCurvatures",nsmap=serializer.NSMAP)
            principcurvbuf = BytesIO()
            np.savetxt(principcurvbuf,self.principal_curvatures,delimiter=" ",newline=" ")
            principal_curvatures.attrib["curvatures"]=principcurvbuf.getvalue().decode('utf-8')
            indexedfaceset.append(principal_curvatures)

            curvature_tangent_axes=etree.Element(serializer.spatialnde_NSPRE+"CurvatureTangentAxes",nsmap=serializer.NSMAP)
            tangaxbuf = BytesIO()
            np.savetxt(tangaxbuf,ourframe.transformto(destframe,self.curvature_tangent_axes,cartesianaxis=1,inhomogeneous_vector=True).ravel(),delimiter=" ",newline=" ")
            curvature_tangent_axes.attrib["axes"]=tangaxbuf.getvalue().decode('utf-8')
            indexedfaceset.append(curvature_tangent_axes)
            
            pass
        
        
        shape.append(indexedfaceset)
        serializer.scene.append(shape)

        return shape

    
    def VRMLWrite(self,serializer,ourframe=None,destframe=None,UVparameterization=None): #,appearance=None):

        serializer.buf.write(b"Shape {\n")
        if self.appearance is not None:
            serializer.buf.write(self.appearance.generate_vrml(serializer))
            pass
        
        serializer.buf.write(b"  geometry IndexedFaceSet {\n")

        # Write coordIndex as attribute of IndexedFaceSet
        serializer.buf.write(b"    coordIndex [\n")
        for cnt in range(self.vertexidx_indices.shape[0]):
            firstidx=self.vertexidx_indices[cnt]
            #numvertices=np.count_nonzero(self.vertexids[cnt,:] >= 0)
            np.savetxt(serializer.buf,self.vertexidx[firstidx:(firstidx+self.numvertices[cnt])],fmt="%d",delimiter=" ",newline=" ")
            serializer.buf.write(b" -1 ") # barrier between polygons
            pass
        serializer.buf.write(b"    ]\n")  # terminate coordIndex

        # indexedfaceset.attrib["solid"]="true"

        # Add <Coordinate> tag
        serializer.buf.write(b"    coord Coordinate {\n")
        
        # Add point= attribute of <Coordinate> element
        serializer.buf.write(b"      point [\n")
        np.savetxt(serializer.buf,ourframe.transformto(destframe,self.vertices),delimiter=" ",newline=" ")
        serializer.buf.write(b"      ]\n") # end point= attribute
        serializer.buf.write(b"    }\n") # end Coordinate

        # Need to generate normals

        serializer.buf.write(b"    normal Normal {\n")
        serializer.buf.write(b"      vector [\n")
        normalindexbuf=BytesIO()
        
        if self.normalidx is None:
            # Single normal per polygon
            assert(self.normals.shape[0]==self.vertexidx_indices.shape[0]) # consistency...
            
            # Evaluate all normals
            np.savetxt(serializer.buf,ourframe.transformto(destframe,self.normals,inhomogeneous_vector=True),delimiter=" ",newline=" ")
            
            for cnt in range(self.normals.shape[0]):
                #numvertices=np.count_nonzero(self.vertexids[cnt,:] >= 0)
                normalindexbuf.write(((str(cnt)+" ") * self.numvertices[cnt]).encode('utf-8')) # Write out normal index, numvertices times
                normalindexbuf.write(b" -1 ") # barrier between polygons
                pass
            
            pass
        else:
            # Seperate normal per vertex
            #assert(self.vertexnormals.shape[0]==self.vertexidx_indices.shape[0]) # consistency...
            normalpos=0
        
            for cnt in range(self.vertexidx_indices.shape[0]):
                #numvertices=np.count_nonzero(self.vertexids[cnt,:] >= 0) # How many vertices does this polygon have?
                firstvertex=self.vertexidx_indices[cnt]
                
                np.savetxt(serializer.buf,ourframe.transformto(destframe,self.normals[self.normalidx[firstvertex:(firstvertex+self.numvertices[cnt])],:],inhomogeneous_vector=True),delimiter=" ",newline=" ")
                
                np.savetxt(normalindexbuf,np.arange(normalpos,normalpos+self.numvertices[cnt],dtype='i32'),fmt="%d",delimiter=" ",newline=" ")
                
                normalpos+=self.numvertices[cnt]
                pass
            
            pass
        serializer.buf.write(b"      ]\n") # close Vector
        serializer.buf.write(b"    }\n") # close Normal
 
        serializer.buf.write(b"    normalIndex [\n")
        serializer.buf.write(normalindexbuf.getvalue())
        serializer.buf.write(b"    ]\n") # close normalIndex        
        
        # Need to generate texCoordIndex attribute
        # and TextureCoordinate element
        # Need to generate normals
        

        if UVparameterization is None:
            # UVparameterization of None means use intrinsic values
            UVparameterization = self.intrinsicparameterization
            pass

        if UVparameterization is not None:
            serializer.buf.write(b"    texCoord TextureCoordinate {\n")
            serializer.buf.write(b"      point [\n")
        
            for cnt in range(self.vertexidx_indices.shape[0]):
                # numvertices=np.count_nonzero(self.vertexids[cnt,:] >= 0) # How many vertices does this polygon have?
                numvertices=self.numvertices[cnt]
            
            

                #iuvcoords=self.intrinsicparameterization.eval_uv(self,self.vertexids[cnt,:numvertices])
                uvcoords=UVparameterization.eval_texcoord_polygonvertex(self,cnt,range(numvertices))
            
            
                np.savetxt(serializer.buf,uvcoords,delimiter=" ",newline=" ")
            
                pass
            serializer.buf.write(b"      ]\n") # close point
            serializer.buf.write(b"    }\n") # close TexCoord

            # Finally, texcoordindex
            serializer.buf.write(b"    texCoordIndex [\n")
            texcoordpos=0
            for cnt in range(self.vertexidx_indices.shape[0]):
                numvertices=self.numvertices[cnt]
                #numvertices=np.count_nonzero(self.vertexids[cnt,:] >= 0) # How many vertices does this polygon have?
                np.savetxt(serializer.buf,np.arange(texcoordpos,texcoordpos+numvertices,dtype='i4'),fmt="%d",delimiter=" ",newline=" ")
                serializer.buf.write(b" -1 ")
                texcoordpos+=numvertices
                pass
        
            serializer.buf.write(b"    ]\n") # close texCoordIndex        
            pass
        serializer.buf.write(b"    solid %s\n" % (str(self.solid).upper()))
        
        serializer.buf.write(b"  }\n") # Close IndexedFaceSet
        serializer.buf.write(b"}\n") # Close Shape

        pass
    

    @classmethod
    def fromx3d(cls,surfaceid,x3dsurface,x3dappearance,cadpartparams=None,tol=1e-6,defaultappearance=None):
        # WARNING: Do not modify x3dsurface subelements
        # after making this call, as the resulting object may
        # (or may not) reference the internal variables!!!
        
        # Meshed surface
        # Available elements
        # x3dsurface.coord
        # x3dsurface.normal
        # x3dsurface.texCoord
        # x3dsurface.texCoordIndex
        # x3dsurface.ccw
        # x3dsurface.coordIndex
        # x3dsurface.normalIndex
        # x3dsurface.normalPerVertex
        # x3dsurface.transform

        
        # go through coordIndex and count number of vertices in each
        # element

            
        (vertexidx_indices,numvertices,missing_final_terminator)=find_vertexidx_indices(x3dsurface.coordIndex)
        if (missing_final_terminator):
            vertexidx=np.concatenate(x3dsurface.coordIndex,np.array((-1,),dtype=np.int32))
            pass
        else:
            vertexidx=x3dsurface.coordIndex
            pass
        

        

        texcoord=x3dsurface.texCoord # may be None
        texcoordidx=x3dsurface.texCoordIndex

        # texcoord should define
        # the intrinsic parameterization. 

        # transformedcoord=x3dsurface.coord
        
        transformedcoord=np.inner(x3dsurface.transform[:3,:3],x3dsurface.coord).T + x3dsurface.transform[:3,3].reshape(1,3)

        # We will need vertices to be C-contiguous later, so let us
        # make it that way now
        vertices=np.empty(transformedcoord.shape,dtype='d')
        vertices[:,:]=transformedcoord

        normals=calcnormals(vertices,vertexidx_indices,vertexidx,numvertices)
        normalidx=None
        
        
        # x3dappearance.texture.url
        if x3dappearance is not None and x3dappearance.texture is not None and x3dappearance.texture.url is not None:
            appearance = texture_url.from_url(x3dappearance.texture.url)
            pass
        elif x3dappearance is not None and x3dappearance.material is not None:
            # !!!*** fixme: Should be able to map in textures, more than just diffuse color, etc. 
            appearance = simple_material.from_color(x3dappearance.material.diffuseColor)
            pass
        else:
            appearance = defaultappearance
            pass

        principal_curvatures=x3dsurface.principal_curvatures
        curvature_tangent_axes=x3dsurface.curvature_tangent_axes

        if texcoord is not None:
            buildintrinsicparameterization=lambda surf,params: polygonalsurface_texcoordparameterization.new(surf,texcoord,texcoordidx,appearance,params)
            pass
        else:
            buildintrinsicparameterization=None  # may default to planarparameterization
            pass
        
        surface = cls.fromvertices(surfaceid,
                                   vertices,
                                   vertexidx_indices,
                                   vertexidx,
                                   numvertices,
                                   normals,
                                   normalidx,
                                   buildintrinsicparameterization=buildintrinsicparameterization,
                                   cadpartparams=cadpartparams,
                                   appearance=appearance,
                                   solid=x3dsurface.solid,
                                   tol=tol,
                                   principal_curvatures=principal_curvatures,
                                   curvature_tangent_axes=curvature_tangent_axes)
        return surface
    
    def assign_appearance(self,appearance):
        self.appearance=appearance
        pass
    
    def assign_curvatures(self,curvatures,curvature_tangent_axes):
        # Note that assignment uses a reference to the provided arrays
        self.principal_curvatures=curvatures
        self.curvature_tangent_axes=curvature_tangent_axes
        pass

    def clear_curvatures(self):
        self.principal_curvatures=None
        self.curvature_tangent_axes=None
        pass

    pass


