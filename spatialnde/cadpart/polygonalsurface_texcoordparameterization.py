import sys
import copy
import json

import numpy as np
import numpy.linalg    
#import scipy
#import scipy.linalg


from ..geometry import vertex_all_in_set

use_geom_accel=True

if use_geom_accel:
    from ..geometry_accel import point_in_polygon_2d
    from ..geometry_accel import polygon_intersects_box_2d
    from ..geometry_accel import box_inside_polygon_2d
else:
    from ..geometry import point_in_polygon_2d
    from ..geometry import polygon_intersects_box_2d
    from ..geometry import box_inside_polygon_2d
    pass

use_pstp_accel=True

if use_pstp_accel:
    from .polygonalsurface_texcoordparameterization_accel import enclosed_or_intersecting_polygons_2d_accel
    from .polygonalsurface_texcoordparameterization_accel import _evaluate_curvature
    from .polygonalsurface_texcoordparameterization_accel import _identify_polynum_uv
    #from .polygonalsurface_texcoordparameterization_accel import _test_polynum_uv
    from .polygonalsurface_texcoordparameterization_accel import _linelength_avgcurvature
    from .polygonalsurface_texcoordparameterization_accel import _linelength_avgcurvature_meshbased
    from .polygonalsurface_texcoordparameterization_accel import _linelength_avgcurvature_mirroredbox_meshbased
    pass


def vecnorm(a,axis=-1):
    
    # compatibility for old versions of Numpy
    return np.sqrt(np.sum(a**2.0,axis))

def vecnormkeepshape(a,axis=-1):
    if axis < 0:
        axis+=len(a.shape)
        pass

    retarray = np.sqrt(np.sum(a**2.0,axis))

    # Insert an axis of length one where we summed
    # So if orig shape is (4, 8, 9 ,3, 2)
    # and we summed over axis 3,
    # result is (4, 8, 9, 2)
    #
    # We reshape result to (4, 8, 9, 1, 2)
    return retarray.reshape(*( a.shape[:axis] + (1,) + a.shape[(axis+1):]))

    


def enclosed_or_intersecting_polygons_2d(polypool,vertexpool,texcoord,vertexidx_indices,texcoordidx,numvertices,texcoordredundant_numcopies,texcoordredundant_polystartindexes,texcoordredundant_polystartpolynum,texcoordredundant_texcoordidx,minu,minv,maxu,maxv):
    #res=copy.copy(polypool)
    #res=np.zeros(polypool.shape[0],dtype=np.bool)
    res=[]
    
    vertexset=frozenset(vertexpool)
    num_fully_enclosed=0

    for poolidx in np.arange(len(polypool)):
        idx=polypool[poolidx]
        if idx < 0:
            # masked out polygon
            continue

        if idx < vertexidx_indices.shape[0]:
            firstidx=vertexidx_indices[idx]
            vertexidxs=texcoordidx[firstidx:(firstidx+numvertices[idx])]
            polygon_fully_enclosed = vertex_all_in_set(vertexidxs,vertexset)
            pass
        else:
            # redundant texcoords
            firstidx=texcoordredundant_polystartindexes[idx-vertexidx_indices.shape[0]]
            polynum=texcoordredundant_polystartpolynum[idx-vertexidx_indices.shape[0]]
            vertexidxs = texcoordredundant_texcoordidx[firstidx:(firstidx+numvertices[polynum])]
            polygon_fully_enclosed = vertex_all_in_set(vertexidxs,vertexset)
            pass

        if polygon_fully_enclosed:
            res.append(idx)
            # if it's fully enclosed, nothing else need look at at, so we filter it here from the broader sibling pool
            polypool[poolidx] = -1 # mask out polygon
            num_fully_enclosed +=1
            pass
        

        # ***Could optimize further here by checking if all of the vertices of a
        # polygon are beyond each bound of the box. 

        if not polygon_fully_enclosed:
            box_v0=np.array((minu,minv),dtype='d')
            box_v1=np.array((maxu,maxv),dtype='d')
            
            # does it intersect?
            if polygon_intersects_box_2d(box_v0,
                                         box_v1,
                                         texcoord[vertexidxs,:]):
                res.append(idx)
                # Don't filter it out in this case because it must
                # intersect with a sibiling too 
                pass

            # What if the box is entirely inside the polygon?
            # Then we should return it also

            elif box_inside_polygon_2d(box_v0,box_v1,texcoord[vertexidxs,:]):
                res.append(idx)
                # Don't filter it out in this case because it must
                # intersect with sibilings too 
                pass
            pass
        
        pass
    # return updated pool as array
    return (num_fully_enclosed,np.array(res,dtype=np.int32))


class polygonalsurface_texcoordparameterization(object):
    
    cadpartparams=None
    lowerleft_meaningfulunits=None # Coordinates of the lower left corner of the texture, in meaningful units.... these are the coordinates of the lower left corner of the lower left pixel
    meaningfulunits_per_texcoord=None

    # These next members are based on the polygons from the surface
    # They should probably be modified to index on which surface, so that
    # multiple surfaces parameterized into a single space
    # can share this parameterization object

    # p is number of independent texture coordinate vertices
    
    texcoord=None # p x 3 texture coordinate array.... Only valid for valid vertices.
    texcoordidx=None # indexes into texcoord for each polygon, terminated by -1, must be structured identically to polygonalsurface.vertexidx, and indices can be identified with polygonalsurface.vertexidx_indices and polygonalsurface.numvertices

    # Support for redundant texture mappings 
    texcoordredundant_firstpolynum=None  # indexed by which polygon, gives index into texcoordredundant_polystartindexes
    texcoordredundant_numcopies=None # indexed by which polygon, gives number of redundant copies of the polygon, i.e. number of elements in texcoordredundant_polystartindexes that correspond to this polygon
                                          
    texcoordredundant_polystartindexes=None # single index indicates which redundant copy of which polygon. Value gives starting index into texcoordredundant_texcoordidx (# of indices into texcoordredundant_texcoordidx is given by polysurf.numvertices[polynum].  The length of this vector is the total of texcoordredundant_numcopies.
    texcoordredundant_polystartpolynum=None # indexed like texcoordredundant_polystartindexes, but gives the original polynum ( < vertexidx_indices.shape[0] )
                                          
    texcoordredundant_texcoordidx=None # this is a vector of indices into the texcoord array. Entries for the same physical polygon (facet) are separated by -2. Entries for different polygons are separated by -1. If there is a polygon that has no redundant texcoords then there will still an extra -1 in the vector. 
     

    ### ***!!! Will need redundant copies of these too... or
    ### More precisely, they need to be as big as the total of
    ### original + redundant polygons. 
    inplane2texcoords = None # 2x3 matrix that will convert in-plane (based on polysurf.inplanemats)
                             # coordinates to texture coordinates
    texcoords2inplane = None # 3x3 matrix that will convert texture coordinates
                             # to in-plane (polysurf.inplanemats) coords. 







    boxes = None # Array of bounding boxes. Each box is identified by integers:
                 #  4 box children, or -1, identified as indices into
                 # this array, then an index into boxpolys array
                 # then the number of entries in boxpolys array
                 # So this is a num_boxes x 6 array
                 #
                 # Note that current boxing is inefficient because polygons on the boundary of big boxes tend to get tested a lot
                 # auto-set by buildboxes
    boxpolys=None # array of polygon indices for each leaf box,
                  # each terminated by -1
    boxcoords=None # Array of box coordinates, same structure as boxes, 4 values per box: 4 box bounds, Box bounds are (minu,minv,maxu,maxv)



    def _buildbox(self,polysurf,boxlist,boxcoordlist,polys,boxpolys,vertexpool,cnt,depth,minu,minv,maxu,maxv):

        #import pdb
        # pdb.set_trace()
        # polys is a truth value array
        # boxpolys is accumulated array of where
        # polys is true. 
        thisbox=np.ones(6,dtype='i4')*-1

        
        #fullvertexpool=np.where((self.texcoord[:,0] >= minu) & (self.texcoord[:,0] <= maxu) & (self.texcoord[:,1] >= minv) & (self.texcoord[:,1] <= maxv))[0]
        #assert((fullvertexpool==ourvertexpool).all())
        
        # filter down polys according to what is in this box
        if depth != 0: # all pass for depth = 0
            #debugpolys = copy.copy(polys)
            
            if use_pstp_accel:
                (num_fully_enclosed,ourpolys)=enclosed_or_intersecting_polygons_2d_accel(polys,self.texcoord,polysurf.vertexidx_indices,self.texcoordidx,polysurf.numvertices,self.texcoordredundant_numcopies,self.texcoordredundant_polystartindexes,self.texcoordredundant_polystartpolynum,self.texcoordredundant_texcoordidx,minu,minv,maxu,maxv)

                ourvertexpool=None
                pass
            else:
                vertices=self.texcoord[vertexpool,:]
                
                # Filter down vertexpool from parent
                ourvertexpool=vertexpool[np.where((vertices[:,0] >= minu) &
                                                  (vertices[:,0] <= maxu) &
                                                  (vertices[:,1] >= minv) &
                                                  (vertices[:,1] <= maxv))[0]]
                
                (num_fully_enclosed,ourpolys)=enclosed_or_intersecting_polygons_2d(polys,ourvertexpool,self.texcoord,polysurf.vertexidx_indices,self.texcoordidx,polysurf.numvertices,self.texcoordredundant_numcopies,self.texcoordredundant_polystartindexes,self.texcoordredundant_polystartpolynum,self.texcoordredundant_texcoordidx,minu,minv,maxu,maxv)
                pass
            ##debugpolys=np.arange(polysurf.vertexidx_indices.shape[0])
            #allvertexpool=np.where((self.texcoord[:,0] >= minu) &
            #                       (self.texcoord[:,0] <= maxu) &
            #                       (self.texcoord[:,1] >= minv) &
            #                       (self.texcoord[:,1] <= maxv))[0]
            #(junk,allpolys)=enclosed_or_intersecting_polygons_2d(debugpolys,allvertexpool,self.texcoord,polysurf.vertexidx_indices,self.texcoordidx,polysurf.numvertices,self.texcoordredundant_numcopies,self.texcoordredundant_polystartindexes,self.texcoordredundant_polystartpolynum,self.texcoordredundant_texcoordidx,minu,minv,maxu,maxv)
            #assert((ourpolys==allpolys).all())
            
            pass
        else:
            ourpolys=polys
            num_fully_enclosed=ourpolys.shape[0]
            ourvertexpool=vertexpool
            pass
        
        
        
        boxlist.append(thisbox)
        boxcoordlist.append((minu,minv,maxu,maxv))

        newcnt=cnt+1

        if num_fully_enclosed > 6 and depth <= 18:
            # split up box
            distu=maxu-minu
            distv=maxv-minv
            eps=1e-4*np.sqrt(distu**2 + distv**2)

            thisbox[0]=newcnt
            newcnt=self._buildbox(polysurf,boxlist,boxcoordlist,ourpolys,boxpolys,ourvertexpool,newcnt,depth+1,minu,minv,minu+distu/2.0+eps,minv+distv/2.0+eps)
            thisbox[1]=newcnt
            newcnt=self._buildbox(polysurf,boxlist,boxcoordlist,ourpolys,boxpolys,ourvertexpool,newcnt,depth+1,minu+distu/2.0-eps,minv,maxu,minv+distv/2.0+eps)
            thisbox[2]=newcnt
            newcnt=self._buildbox(polysurf,boxlist,boxcoordlist,ourpolys,boxpolys,ourvertexpool,newcnt,depth+1,minu,minv+distv/2.0-eps,minu+distu/2.0+eps,maxv)
            thisbox[3]=newcnt
            newcnt=self._buildbox(polysurf,boxlist,boxcoordlist,ourpolys,boxpolys,ourvertexpool,newcnt,depth+1,minu+distu/2.0-eps,minv+distv/2.0-eps,maxu,maxv)
            pass
        #
        else:
            # This is a leaf node
            # Record our polygons... these are those which are fully enclosed
            # or intersecting us, the index where they start is boxpolys[4]

            goodpolys=ourpolys >= 0
            thisbox[4]=len(boxpolys)
            thisbox[5]=np.count_nonzero(goodpolys)
            boxpolys.extend(ourpolys[goodpolys])
            boxpolys.append(-1)
            pass
        
        return newcnt

    

    def invalidateboxes(self):
        self.boxes=None
        self.boxcoords=None
        self.boxpolys=None
        pass

    def buildboxes(self,polysurf):
        if self.boxes is not None:
            return # already built
        
        boxlist=[]
        boxcoordlist=[]

        numpolys=polysurf.vertexidx_indices.shape[0]
        if self.texcoordredundant_numcopies is not None:
            numpolys += self.texcoordredundant_numcopies.sum()
            pass
        
        # polys: indices of polygons which are inside this box 

        polys=np.arange(numpolys,dtype=np.int32)  # -1 will be used as a flag to indicate this poly is no longer present
        if use_pstp_accel:
            vertexpool=None # not used in accelerated mode
            pass
        else:
            vertexpool=np.arange(self.texcoord.shape[0],dtype=np.uint32) # indexes of texture vertices within a given box
            pass
        
        boxpolys=[]
        boxcnt=0

        minu=0.0
        maxu=1.0

        minv=0.0
        maxv=1.0

        #(minx,miny,minz)=np.min(self.vertices,axis=0)
        #(maxx,maxy,maxz)=np.max(self.vertices,axis=0)

        distu=maxu-minu
        distv=maxv-minv

        eps=1e-4*np.sqrt(distu**2 + distv**2)

        self._buildbox(polysurf,boxlist,boxcoordlist,polys,boxpolys,vertexpool,0,0,minu-eps,minv-eps,maxu+eps,maxv+eps)
        
        self.boxes=np.array(boxlist,dtype='i4')
        self.boxcoords=np.array(boxcoordlist,dtype='d')
        self.boxpolys=np.array(boxpolys,dtype='i4')
        
        pass

    
    #def eval_uv(self,polysurf,polysurf_vertexids):
    #    coords=polysurf.vertices[polysurf_vertexids,:] # n x 3 matrix
    #    relcoords=coords-self.centercoords.reshape(1,3)
    #
    #    return np.array((np.inner(relcoords,self.u),np.inner(relcoords,self.v)),dtype='d').T


    # def _determine_tex_xform(self,polysurf,polysurf_polynum):
    #     # See also scope_coin3d.cpp:DetermineTexXform
    #    # See also polygonalsurface.py: buildprojinfo()
    # 
    #     # NOTE: AijMat is 2x3, but AijMatInv is 3x3 for historical reasons 
    #
    #     # OK. So we extract the numpoints points. These are in
    #     # 3-space and nominally on a plane. How to flatten them
    #     # into 2D?
    #     # 1. Subtract out the centroid (best-fit plane
    #     #    will pass through centroid)
    #     #    (Also subtract same location from ObjectSpace intersection
    #     #    point coordinates, from above)
    #     # 2. Assemble point vectors into a matrix of vectors from centroid
    #     # 3. Apply SVD. Resulting two basis vectors corresponding to
    #     #    largest-magnitude two singular values will span the plane
    #     # 4. Re-evaluate point vectors and intersection location in
    #     #    terms of these basis vectors. Now our point coordinates
    #     #    and intersection coordinates are in 2-space. 
    #     # 5. Evaluate a transform
    #     #    [ A11 A12 A13 ][ x ] = [ tex ]
    #     #    [ A21 A22 A23 ][ y ] = [ tey ]
    #     #    [ 0    0  1   ][ 1 ] = [  1  ]
    #     #    or
    #     #    [ x  y  1  0  0  0 ] [ A11 ] = [ tex ]
    #     #    [ 0  0  0  x  y  1 ] [ A12 ] = [ tey ]
    #     #                         [ A13 ]
    #     #                         [ A21 ]
    #     #                         [ A22 ]
    #     #                         [ A23 ]
    #     # With rows repeated for each point.
    #     # Solve for Axx values from the known coordinates
    #     # Then substitute the 2D intersection coordinates as (x,y)
    #     # and multiply to get (tex,tey), the desired texture coordinates.

    #     numpoints=np.count_nonzero(polysurf.vertexids[polysurf_polynum,:] >= 0)
        
    #     centroid = np.mean(polysurf.vertices[polysurf.vertexids[polysurf_polynum,:numpoints],:],axis=0)
    #     coordvals = (polysurf.vertices[polysurf.vertexids[polysurf_polynum,:numpoints],:]-centroid.reshape(1,3)).T  # coordvals is the coordinates relative to centroid, 3 x numpoints

    #     texcoordvals = self.texcoord[polysurf_polynum,:numpoints].T # texcoordvals is the texture coordinates, 2 rows by numpoints cols... # Note that textures are in range 0...1 by convention
    #
    #     # calculate SVD
    #     (U,s,Vt)=scipy.linalg.svd(coordvals,full_matrices=True,compute_uv=True)
    #
    #     # extract columns for 2d coordinate basis vectors
    #     # want columns that correspond to the largest two
    #     # singular values
    #     xcolindex=0
    #     ycolindex=1
    #
    #     if abs(s[0]) < abs(s[1]) and abs(s[0]) < abs(s[2]):
    #         # element 0 is smallest s.v.
    #         xcolindex=2
    #         pass
    #     if abs(s[1]) < abs(s[2]) and abs(s[1]) < abs(s[0]):
    #         # element 1 is smallest s.v.
    #         ycolindex=2
    #         pass
    #
    #     To2D=U[:,np.array((xcolindex,ycolindex))].T # 2x3... Rows of To2D are x and y basis vectors, respectively
    #
    #     coordvals2d = np.dot(To2D,coordvals) # 2 rows by numpoints cols... in 2D basis relative to centroid
    #
    #     TexXformMtx=np.zeros((2*numpoints,6),dtype='d')
    #     TexXformMtx[:(2*numpoints):2,0]=coordvals2d[0,:] # assign 'x' elements
    #     TexXformMtx[:(2*numpoints):2,1]=coordvals2d[1,:] # assign 'y' elements
    #     TexXformMtx[:(2*numpoints):2,2]=1   # assign '1' entries
    #     TexXformMtx[1:(2*numpoints):2,3]=coordvals2d[0,:] # assign 'x' elements
    #     TexXformMtx[1:(2*numpoints):2,4]=coordvals2d[1,:] # assign 'y' elements
    #     TexXformMtx[1:(2*numpoints):2,5]=1   # assign '1' entries
    #     
    #     TexCoordVec=np.zeros((2*numpoints),dtype='d')
    #     TexCoordVec[:(2*numpoints):2] = texcoordvals[0,:]  # assign tex
    #     TexCoordVec[1:(2*numpoints):2] = texcoordvals[1,:] # assign tey
    #
    #     (AijVals,residuals,rank,lstsq_s) = np.linalg.lstsq(TexXformMtx,TexCoordVec)
    #     AijMat=AijVals.reshape(2,3) # reshape to 2x3
    #     AijMatExt = np.concatenate((AijMat,np.array((0.0,0.0,1.0),dtype='d').reshape(1,3)),axis=0) # Add 0.0, 0.0, 1.0 row to bottom of matrix
    #    
    #     AijMatInv=np.linalg.inv(AijMatExt)
    #
    #     return (centroid,s,xcolindex,ycolindex,To2D, AijMat,AijMatInv)
    
    def eval_texcoord_polygonvertex(self,polysurf,polysurf_polynum,polysurf_vertexnum):
        # Can supply vectors as polysurf_polynum and/or polysurf_vertexnum

        #texcoords = self.texcoord[polysurf_polynum,polysurf_vertexnum,:]
        firstidx=polysurf.vertexidx_indices[polysurf_polynum]
        
        texcoords = self.texcoord[self.texcoordidx[firstidx+polysurf_vertexnum],:]
        

        return texcoords


    def invalidateprojinfo(self):
        self.inplane2texcoords = None
        self.texcoords2inplane = None
        pass
        
    
    def buildprojinfo(self,polysurf):
        # See also scope_coin3d.cpp:DetermineTexXform
        # see also polygonalsurface_intrinsicparameterization.py/_determine_tex_xform()
        
        
        # and preceding steps in polygonalsurface.py:buildprojinfo()
        # 5. Evaluate a transform
        #    [ A11 A12 A13 ][ x ] = [ tex ]
        #    [ A21 A22 A23 ][ y ] = [ tey ]
        #    [ 0    0  1   ][ 1 ] = [  1  ]
        #    or
        #    [ x  y  1  0  0  0 ] [ A11 ] = [ tex ]
        #    [ 0  0  0  x  y  1 ] [ A12 ] = [ tey ]
        #                         [ A13 ]
        #                         [ A21 ]
        #                         [ A22 ]
        #                         [ A23 ]
        # With rows repeated for each point.
        # Solve for Axx values from the known coordinates
        # Then substitute the 2D intersection coordinates as (x,y)
        # and multiply to get (tex,tey), the desired texture coordinates.

        if self.inplane2texcoords is not None:
            return # already built
        
        numpolys=polysurf.vertexidx_indices.shape[0]

        self.inplane2texcoords = np.zeros((numpolys,2,3),dtype='d')
        self.texcoords2inplane = np.zeros((numpolys,3,3),dtype='d')

        for polynum in range(numpolys):
            firstidx=polysurf.vertexidx_indices[polynum]

            numpoints=polysurf.numvertices[polynum]

            centroid = polysurf.refpoints[polynum,:]
            coordvals = (polysurf.vertices[polysurf.vertexidx[firstidx:(firstidx+numpoints)],:]-centroid.reshape(1,3)).T # coordvals is the coordinates relative to centroid, 3 x numpoints
            
            To2D = polysurf.inplanemats[polynum,:,:]
            
            coordvals2d = np.dot(To2D,coordvals) # 2 rows by numpoints cols... in 2D basis relative to centroid

            
            texcoordvals = self.texcoord[self.texcoordidx[firstidx:(firstidx+numpoints)],:].T # texcoordvals is the texture coordinates, 2 rows by numpoints cols... # Note that textures are in range 0...1 by convention
            
            
            
            TexXformMtx=np.zeros((2*numpoints,6),dtype='d')
            TexXformMtx[:(2*numpoints):2,0]=coordvals2d[0,:] # assign 'x' elements
            TexXformMtx[:(2*numpoints):2,1]=coordvals2d[1,:] # assign 'y' elements
            TexXformMtx[:(2*numpoints):2,2]=1   # assign '1' entries
            TexXformMtx[1:(2*numpoints):2,3]=coordvals2d[0,:] # assign 'x' elements
            TexXformMtx[1:(2*numpoints):2,4]=coordvals2d[1,:] # assign 'y' elements
            TexXformMtx[1:(2*numpoints):2,5]=1   # assign '1' entries
            
            TexCoordVec=np.zeros((2*numpoints),dtype='d')
            TexCoordVec[:(2*numpoints):2] = texcoordvals[0,:]  # assign tex
            TexCoordVec[1:(2*numpoints):2] = texcoordvals[1,:] # assign tey


            (AijVals,residuals,rank,lstsq_s) = np.linalg.lstsq(TexXformMtx,TexCoordVec,rcond=-1)
            AijMat=AijVals.reshape(2,3) # reshape to 2x3
            AijMatExt = np.concatenate((AijMat,np.array((0.0,0.0,1.0),dtype='d').reshape(1,3)),axis=0) # Add 0.0, 0.0, 1.0 row to bottom of matrix

            # NOTE: Possible bug: This matrix inversion (next line) will
            # fail if the polygon has zero area in texture space due to
            # (for example) limited precision in writing down the
            # texture coordinates in the data file.
            #
            # Not sure what to do in this case...
            
            AijMatInv=np.linalg.inv(AijMatExt)
            # Assign AijMat 
            self.inplane2texcoords[polynum,:,:]=AijMat
            self.texcoords2inplane[polynum,:,:]=AijMatInv
            pass

        pass


    def _evaluate_curvature(self,polysurf,polynum,u,v):
        # Evaluate the curvature, within polygon # polynum
        # at (u,v) coordinates...  (u,v) in texture coordinate
        # range [0...1]
        # ... C accelerated version available
    
                                
        if polynum >= polysurf.vertexidx_indices.shape[0]:
            # This polynum corresponds to a redundant texture
            polysurf_polynum=self.texcoordredundant_polystartpolynum[polynum]
            pass
        else:
            polysurf_polynum=polynum
            pass
        
        To2D=polysurf.inplanemats[polysurf_polynum,:,:] # To2D is 2x3
        #AijMat=self.inplane2texcoords[polynum,:,:]
        AijMatInv=self.texcoords2inplane[polynum,:,:]
        
        # Note Capital UV represent the texture  parameterization
        # of the in-plane 3D space of this facet. 
        TexUVExt = np.inner(AijMatInv,np.array((u,v,1.0)))
        TexUVExt /= TexUVExt[2] # normalize inhomogeneous coordinates
        
        # These coordinates of this (u,v) of this facet are relative to its centroid,
        # and are in terms of the basis vectors in To2D
        TexUV = TexUVExt[:2]
        
        # Get 3D coordinates relative to centroid
        Tex3D = np.inner(To2D.T,TexUV)
        
        # Need to evaluate 3D vertex coords, relative to centroid,
        # Use them to weight the vertex curvatures
        # according to distance from our point.
        
        centroid = polysurf.refpoints[polysurf_polynum,:] # Centroied in 3d coords
        firstidx=polysurf.vertexidx_indices[polysurf_polynum]
        numpoints=polysurf.numvertices[polysurf_polynum]
        
        # Check to see if we have curvatures at all vertices:
        if np.isnan(polysurf.principal_curvatures[polysurf.vertexidx[firstidx:(firstidx+numpoints)],0]).any():
            # abort if we are missing a curvature
            #curvmat[vcnt,ucnt,:,:]=np.NaN
            return np.array(((np.NaN,np.NaN),(np.NaN,np.NaN)),dtype='d')
        
        # For this facet, the 3D coords of the vertices are
        coordvals = (polysurf.vertices[polysurf.vertexidx[firstidx:(firstidx+numpoints)],:]-centroid.reshape(1,3)).T # coordvals is the coordinates relative to centroid, 3 x numpoints

        # Now coordvals is 3 x numvertices, coordinates of the vertices
        # relative to centroid
        # Tex3D is 3 vector, coordinates of our (u,v) location
        # relative to centroid.
        #
        # Perform weighted average
        dists = vecnorm(Tex3D.reshape(3,1) - coordvals,axis=0)
        eps = np.max(dists)/10000.0 # small number, so we don't divide by 0 
                
        rawweights=1.0/(dists+eps)
        totalweights=np.sum(rawweights)
        
        weights=rawweights/totalweights
        
        ## The 2D coords of the vertices are 
        #coordvals2d = np.dot(To2D,coordvals) # 2 rows by numpoints cols... in 2D basis relative to centroid
        
        # Likewise 2D coords of the curvature_tangent_axes
        CTA_2D = np.inner(To2D,polysurf.curvature_tangent_axes[polysurf.vertexidx[firstidx:(firstidx+numpoints)],:,:]).transpose(1,0,2) # Transpose to keep broadcast axis to the left. Pre-transpose axes lengths are: 2 (2D axes) by # of vertices by 2 (principal curvature)
        # CTA_2D axes: # of vertices by 2 (2D axes) by 2 (principal curvature)

        # Normalize curvature_tangent_axes (should be unit length)
        CTA_2D /= vecnormkeepshape(CTA_2D,1) # Axis is axis 0 because it came from To2D
        # Construct curvature matrices ...
        # Need to construct V*K*V', broadcasting over which vertex
        curvmatrices=np.einsum('...ij,...j,...jk->...ik', CTA_2D,polysurf.principal_curvatures[polysurf.vertexidx[firstidx:(firstidx+numpoints)],:],CTA_2D.transpose(0,2,1)) # result is # of vertices by 2x2 curvature matrix

        # Weighting of vertices relative to our point (u,v)
        weightedcurvmatrices = weights.reshape(numpoints,1,1)*curvmatrices

        # meancurvmatrix (weighted average)
        meancurvmatrix = weightedcurvmatrices.sum(axis=0)
        
        # meancurvmatrix is a 2x2 which should be close to symmetric
        asymmetry = meancurvmatrix[1,0]-meancurvmatrix[0,1]
        
        if abs(asymmetry) > 0.1*np.linalg.norm(meancurvmatrix):
            sys.stderr.write("_evaluate_curvature: WARNING Large asymmetry in mean curvature matrix at (u,v) = (%g,%g). Matrix = %s\n" % (u,v,str(meancurvmatrix)))
            pass
        
        # correct asymmetry
        meancurvmatrix[1,0] -= asymmetry/2.0
        meancurvmatrix[0,1] += asymmetry/2.0
        
        ## Determine principal curvatures
        #(princcurvs,evects) = np.linalg.eig(meancurvmatrix)
        
        # curvtangentaxes3d = np.dot(To2D.T,evects)
        # 
        # # We don't want the eigenframe to be mirrored relative to the (U,V)
        # # frame, for consistency in interpreting positive vs. negative curvature.
        # # ...  so if the dot/inner product of (UxV) with (TANGENT0xTANGENT1)
        # is negative, that indicates mirroring 
        # Negating one of the eigenvectors will un-mirror it. 
        #if np.inner(np.cross(To2D[0,:],To2D[1,:]),np.cross(curvtangentaxes[:,0],curvtangentaxes[:,1])) < 0.0:
        #    curvtangentaxes3d[:,0]=-curvtangentaxes3d[:,0]
        #    evects[:,0]=-evects[:,0]
        #    pass
        
        
        ## Then assign curvature
        #principal_curvatures[vcnt,ucnt,:]=princcurvs
        #curvature_tangent_axes[vcnt,ucnt,:,:]=curvtangentaxes.T
        
        # meancurvmatrix is in the polysurf.inplanemats (i.e. To2D)
        # orthonormal basis, with units of meters
        
        # We want to return 2D vectors in the AijMat i.e. self.inplane2texcoords) frame
        # with units of texcoords.
        # i.e. AijMat * meancurv * inv(AijMat)
        #
        # NOTE: AijMat is in inhomogeneous coordinates
        # ... meancurvmatrix represents vectors, not coords
        # so use only first two rows
        return np.dot(self.inplane2texcoords[polynum,:,:2],np.dot(meancurvmatrix,self.texcoords2inplane[polynum,:2,:2]))
        
    # For testing, run x3d_add_curvature.py <_UV.x3d> <.brep> <_UV_curvature.x3d>
    # then curvmat=obj.implpart.surfaces[0].intrinsicparameterization.interpolate_curvature(obj.implpart.surfaces[0],100,100)
    # pl.imshow(curvmat[::-1,:,0,0])  # Gives curvature along X
    # pl.imshow(curvmat[::-1,:,1,1])  # Gives curvature along Y
    # pl.imshow(curvmat[::-1,:,0,1])  # Gives curvature along XY
    def interpolate_curvature(self,polysurf,ny,nx):
        # Interpolate the surface curvature on a pixel grid
        # Grid size is nx by ny going from texcoords 0..1

        if self.boxes is None:
            self.buildboxes(polysurf)
            pass
        
        
        #principal_curvatures=np.empty((ny,nx,2),dtype='f')
        #principal_curvatures[::]=np.NaN
        #curvature_tangent_axes=np.empty((ny,nx,2,3),dtype='f')
        curvmat=np.empty((ny,nx,2,2),dtype='f')
        #curvature_tangent_axes[::]=np.NaN

        # For each pixel, need to figure out which polygon
        # contains that pixel

        polynum=None
        
        for vcnt in range(ny):

            for ucnt in range(nx):
                # ucoord_pixels=u*(nx)-0.5
                # vcoord_pixels=v*(ny)-0.5
                # so u = (ucoord_pixels + 0.5)/nx
                # so v = (vcoord_pixels + 0.5)/ny

                # these u,v in [0,1] range
                u = (ucnt+0.5)/nx
                v = (vcnt+0.5)/ny

                # Try to find which polygon... use result from last attempt as a candidate
                if use_pstp_accel:
                    polynum=_identify_polynum_uv(polysurf.vertexidx_indices,
                                                 polysurf.numvertices,
                                                 self.texcoordidx,
                                                 self.texcoordredundant_polystartindexes,
                                                 self.texcoordredundant_polystartpolynum,
                                                 self.texcoordredundant_texcoordidx,
                                                 self.texcoord,
                                                 self.boxes,
                                                 self.boxpolys,
                                                 self.boxcoords,
                                                 u,v,candidate_polynum=polynum)
                        
                    pass
                else:
                    polynum=self._identify_polynum_uv(polysurf,u,v,candidate_polynum=polynum)
                    pass
                if polynum is None:
                    # Abort if this pixel is outside all polygons
                    curvmat[vcnt,ucnt,:,:]=np.NaN
                    continue

                #if vcnt==50 and ucnt==50:
                #    import pdb
                #    pdb.set_trace()
                #    pass

                if use_pstp_accel: # use Cython accelerated version
                    curvmat[vcnt,ucnt,:,:]=_evaluate_curvature(polysurf.vertexidx_indices,
                                                               polysurf.numvertices,
                                                               polysurf.vertexidx,
                                                               polysurf.vertices,
                                                               polysurf.refpoints,
                                                               self.texcoordredundant_polystartpolynum,
                                                               polysurf.inplanemats,
                                                               self.inplane2texcoords,
                                                               self.texcoords2inplane,
                                                               polysurf.principal_curvatures,
                                                               polysurf.curvature_tangent_axes,
                                                               polynum,u,v)
                    pass
                else:
                    curvmat[vcnt,ucnt,:,:]=self._evaluate_curvature(polysurf,polynum,u,v)
                    pass
                
                
                
                pass

            pass
        
        return curvmat

    


    def _evaluate_stepsize(self,polysurf,polynum,nv,nu):
        # Evaluate the step sizes within polygon # polynum
    
        
        #if polynum >= polysurf.vertexidx_indices.shape[0]:
        #    # This polynum corresponds to a redundant texture
        #    polysurf_polynum=self.texcoordredundant_polystartpolynum[polynum]
        #    pass
        #else:
        #    polysurf_polynum=polynum
        #    pass

        # surf.intrinsicparameterization.inplane2texcoords[something,:,:]
        # is a 2x3 matrix that when multiplied on the right by in-polygon-plane
        # inhomogeneous 2D coordinates relative to polygon centroid, gives (u,v)
        # i.e. left 2x2 gives [ du/dX du/dY ; dv/dX dv/dY ] where X,Y
        # is the polygons arbitrary 2D frame
        #
        # We need to operate on the inverse to get [ dX/du dX/dv ; dY/du dY/dv ]
        # because X and Y are measured in physical distance units so we can
        # get the physical distance corresponding to a pixel in the
        # parameterization. 
        #
        # The original Matrix (extended) is
        #  [  du/dX du/dY  uoffs ]
        #  [  dv/dX dv/dY  voffs ]
        #  [    0    0       1   ]
        #
        #  Represent this as [ Ared offs ]
        #                    [ 0     1   ]
        #  Inverse is then   [ Aredinv   offsinv   ]
        #                    [ 0         1         ]
        #
        # Aredinv = [ dX/du dX/dv ]
        #           [ dY/du dY/dv ]
        # Check: Original*Inverse =  [ Ared*Aredinv    0     ] 
        #                            [ 0               1     ]
        #   Where Ared*Aredinv is eye(2) by definition
        # So we can get the inverse we are looking for from
        # the upper 2x2 of surf.intrinsicparameterization.texcoord2inplane
        #
        # ... So the net result is we multiply texcoord2inplane by [ 1.0/num_u_steps;0;0 ]
        # and then take the vector norm to get the physical step size corresponding
        # to one pixel in u for that polygon.
        # Likewise we multiply texcoord2inplane by [ 0;1.0/num_v_steps;0 ]
        # and then take the vector norm to get the physical step size corresponding
        # to one pixel in u for that polygon.
        #
        # ... Or more efficiently, we take the first two rows of the first column,
        # divide by num_u_steps, and take the norm to get the u step size (m)
        # and take the first two rows of the 2nd column, divide by num_v_steps
        # and take the norm to get the v step size (m)
        
        dRdu=np.linalg.norm(self.texcoords2inplane[polynum,:2,0]/nu)
        dRdv=np.linalg.norm(self.texcoords2inplane[polynum,:2,1]/nv)
        return (dRdu,dRdv)

    
    def interpolate_stepsizes(self,polysurf,ny,nx):
        # Interpolate the surface dv and du physical stepsizes on a pixel grid
        # Grid size is nx by ny going from texcoords 0..1

        if self.boxes is None:
            self.buildboxes(polysurf)
            pass
        
        
        stepsizemat=np.empty((ny,nx,2),dtype='f')

        polynum=None
        
        for vcnt in range(ny):

            for ucnt in range(nx):
                # ucoord_pixels=u*(nx)-0.5
                # vcoord_pixels=v*(ny)-0.5
                # so u = (ucoord_pixels + 0.5)/nx
                # so v = (vcoord_pixels + 0.5)/ny

                # these u,v in [0,1] range
                u = (ucnt+0.5)/nx
                v = (vcnt+0.5)/ny

                # Try to find which polygon... use result from last attempt as a candidate
                if use_pstp_accel:
                    polynum=_identify_polynum_uv(polysurf.vertexidx_indices,
                                                 polysurf.numvertices,
                                                 self.texcoordidx,
                                                 self.texcoordredundant_polystartindexes,
                                                 self.texcoordredundant_polystartpolynum,
                                                 self.texcoordredundant_texcoordidx,
                                                 self.texcoord,
                                                 self.boxes,
                                                 self.boxpolys,
                                                 self.boxcoords,
                                                 u,v,candidate_polynum=polynum)
                    pass
                else:
                    polynum=self._identify_polynum_uv(polysurf,u,v,candidate_polynum=polynum)
                    pass
                
                if polynum is None:
                    # Abort if this pixel is outside all polygons
                    stepsizemat[vcnt,ucnt,:]=np.NaN
                    continue

                stepsizemat[vcnt,ucnt,:]=self._evaluate_stepsize(polysurf,polynum,ny,nx)
                pass

            pass
        
        return stepsizemat

    
    def _test_polynum_uv(self,polysurf,u,v,polynum):
        """This takes [0,1] u,v coordinates """
        # Acclerated Cython version available
        
        # Find whether our (u,v) coordinate is inside this polygon
        #// ... is our (u,v) point inside this polygon?
        if polynum < polysurf.vertexidx_indices.shape[0]:
            firstidx=polysurf.vertexidx_indices[polynum]
            numpoints=polysurf.numvertices[polynum]
            vertexidxs=self.texcoordidx[firstidx:(firstidx+numpoints)]
            
            pass
        else:
            firstidx=self.texcoordredundant_polystartindexes[polynum-polysurf.vertexidx_indices.shape[0]]
            polysurf_polynum=self.texcoordredundant_polystartpolynum[polynum-polysurf.vertexidx_indices.shape[0]]
            vertexidxs = self.texcoordredundant_texcoordidx[firstidx:(firstidx+polysurf.numvertices[polysurf_polynum])]
            
            
            pass

        vertices = self.texcoord[vertexidxs,:]
        
        return point_in_polygon_2d(vertices-np.array(((u,v),),dtype='d'))
    
    def _identify_polynum_uv(self,polysurf,u,v,candidate_polynum=None):
        """This takes [0,1] u,v coordinates """
        # NOTE: Accelerated Cython version available
        
        if self.boxes is None:
            self.buildboxes(polysurf)
            pass


        if candidate_polynum is not None:
            # Check provided polygon to see if it fits
            point_in_poly=self._test_polynum_uv(polysurf,u,v,candidate_polynum) 
            if point_in_poly:
                return candidate_polynum
            pass
        
        
        # Search boxes to find candidate polygons
        candidatepolys=[]
        boxes_to_test=[]
        boxes_to_test.append(0)
        while len(boxes_to_test) > 0:
            curbox=boxes_to_test.pop()
            

            if (u >= self.boxcoords[curbox,0] and
                v >= self.boxcoords[curbox,1] and
                u <= self.boxcoords[curbox,2] and
                v <= self.boxcoords[curbox,3]):
                # we are inside box
                children=self.boxes[curbox,:4]
                if children[0] >= 0:
                    boxes_to_test.extend(children)
                    pass
                if self.boxes[curbox,4] >= 0:
                    polys_in_curbox_idx=self.boxes[curbox,4]
                    polys_in_curbox_num=self.boxes[curbox,5]
                    candidatepolys.extend(self.boxpolys[polys_in_curbox_idx:(polys_in_curbox_idx+polys_in_curbox_num)])
                    pass
                    
                pass
            
            pass
        
        
        for polynum in candidatepolys:

            point_in_poly=self._test_polynum_uv(polysurf,u,v,polynum) 
            if point_in_poly:
                return polynum
            pass
        

        # If we got this far, the search failed! */
        return None
    
    def linelength_avgcurvature_meshbased(self,
                                          polysurf,
                                          curvmats,
                                          stepsizearray,
                                          param_lowerleft_meaningfulunits_u,
                                          param_lowerleft_meaningfulunits_v,
                                          param_stepsize_u,
                                          param_stepsize_v,
                                          du,dv,u1,v1,u2,v2):
        # param_lowerleft_meaningfulunits_u is the u coordinate of the 
        # lower-left corner of the lower left element of curvmats and 
        # stepsize array. Likewise param_stepsize_u and _v are the 
        # step sizes of the stepsize arrya. 
        # NOTE: Unlike the non-meshbased routines, the meshbased 
        # linelength routines take scaled (i.e. meaningfulunits) 
        # u,v coordinates
        
        #u1_unscaled = (u1-self.lowerleft_meaningfulunits[0])/self.meaningfulunits_per_texcoord[0]
        #v1_unscaled = (v1-self.lowerleft_meaningfulunits[1])/self.meaningfulunits_per_texcoord[1]
        
        #u2_unscaled = (u2-self.lowerleft_meaningfulunits[0])/self.meaningfulunits_per_texcoord[0]
        #v2_unscaled = (v2-self.lowerleft_meaningfulunits[1])/self.meaningfulunits_per_texcoord[1]
        
        return _linelength_avgcurvature_meshbased(curvmats,
                                                  stepsizearray,
                                                  param_lowerleft_meaningfulunits_u,
                                                  param_lowerleft_meaningfulunits_v,
                                                  param_stepsize_u,
                                                  param_stepsize_v,
                                                  du,
                                                  dv,
                                                  u1,
                                                  v1,
                                                  u2,
                                                  v2)



    def linelength_avgcurvature_mirroredbox_meshbased(self,
                                                      polysurf,
                                                      curvmats,
                                                      stepsizearray,
                                                      param_lowerleft_meaningfulunits_u,
                                                      param_lowerleft_meaningfulunits_v,
                                                      param_stepsize_u,
                                                      param_stepsize_v,
                                                      boxu1,boxv1,boxu2,boxv2,du,dv,u1,v1,u2,v2):
        # param_lowerleft_meaningfulunits_u is the u coordinate of the 
        # lower-left corner of the lower left element of curvmats and 
        # stepsize array. Likewise param_stepsize_u and _v are the 
        # step sizes of the stepsize arrya. 
        # NOTE: Unlike the non-meshbased routines, the meshbased 
        # linelength routines take scaled (i.e. meaningfulunits) 
        # u,v coordinates

        #boxu1_unscaled = (boxu1-self.lowerleft_meaningfulunits[0])/self.meaningfulunits_per_texcoord[0]
        #boxv1_unscaled = (boxv1-self.lowerleft_meaningfulunits[1])/self.meaningfulunits_per_texcoord[1]

        #boxu2_unscaled = (boxu2-self.lowerleft_meaningfulunits[0])/self.meaningfulunits_per_texcoord[0]
        #boxv2_unscaled = (boxv2-self.lowerleft_meaningfulunits[1])/self.meaningfulunits_per_texcoord[1]
        
        #u1_unscaled = (u1-self.lowerleft_meaningfulunits[0])/self.meaningfulunits_per_texcoord[0]
        #v1_unscaled = (v1-self.lowerleft_meaningfulunits[1])/self.meaningfulunits_per_texcoord[1]
        
        #u2_unscaled = (u2-self.lowerleft_meaningfulunits[0])/self.meaningfulunits_per_texcoord[0]
        #v2_unscaled = (v2-self.lowerleft_meaningfulunits[1])/self.meaningfulunits_per_texcoord[1]
        
        return _linelength_avgcurvature_mirroredbox_meshbased(curvmats,
                                                              stepsizearray,
                                                              param_lowerleft_meaningfulunits_u,
                                                              param_lowerleft_meaningfulunits_v,
                                                              param_stepsize_u,
                                                              param_stepsize_v,
                                                              boxu1,
                                                              boxv1,
                                                              boxu2,
                                                              boxv2,
                                                              du,
                                                              dv,
                                                              u1,
                                                              v1,
                                                              u2,
                                                              v2)
        


    def linelength_avgcurvature(self,polysurf,du,dv,u1,v1,u2,v2):
        """NOTE: This takes meaningfully scaled u,v coordinates"""

        if self.boxes is None:
            self.buildboxes(polysurf)
            pass
        
        u1_unscaled = (u1-self.lowerleft_meaningfulunits[0])/self.meaningfulunits_per_texcoord[0]
        v1_unscaled = (v1-self.lowerleft_meaningfulunits[1])/self.meaningfulunits_per_texcoord[1]

        u2_unscaled = (u2-self.lowerleft_meaningfulunits[0])/self.meaningfulunits_per_texcoord[0]
        v2_unscaled = (v2-self.lowerleft_meaningfulunits[1])/self.meaningfulunits_per_texcoord[1]
        
        
        
        return _linelength_avgcurvature(
            polysurf.vertexidx_indices,
            polysurf.numvertices,
            polysurf.vertexidx,
            polysurf.vertices,
            polysurf.refpoints,
            self.texcoordidx,
            self.texcoordredundant_polystartindexes,
            self.texcoordredundant_polystartpolynum,
            self.texcoordredundant_texcoordidx,
            self.texcoord,
            self.boxes,
            self.boxpolys,
            self.boxcoords,
            polysurf.inplanemats,        
            self.inplane2texcoords,
            self.texcoords2inplane,
            polysurf.principal_curvatures,
            polysurf.curvature_tangent_axes,
            du/self.meaningfulunits_per_texcoord[0],
            dv/self.meaningfulunits_per_texcoord[1],
            u1_unscaled,v1_unscaled,u2_unscaled,v2_unscaled)

    
    def eval_xyz_uv(self,polysurf,u,v):
        """NOTE: This takes meaningfully scaled u,v coordinates"""
        # See also scope_coin3d.cpp:Get3DCoordsGivenTexCoords()

        if self.boxes is None:
            self.buildboxes(polysurf)
            pass
        
        
        # Evaluate x,y,z coordinates given (u,v) coordinates


        # Convert meaningfully scaled (u,v) to 0-1 range
        TexXY=np.array(((u-self.lowerleft_meaningfulunits[0])/self.meaningfulunits_per_texcoord[0],(v-self.lowerleft_meaningfulunits[1])/self.meaningfulunits_per_texcoord[1]),dtype='d')

        if use_pstp_accel:
            polynum=_identify_polynum_uv(polysurf.vertexidx_indices,
                                         polysurf.numvertices,
                                         self.texcoordidx,
                                         self.texcoordredundant_polystartindexes,
                                         self.texcoordredundant_polystartpolynum,
                                         self.texcoordredundant_texcoordidx,
                                         self.texcoord,
                                         self.boxes,
                                         self.boxpolys,
                                         self.boxcoords,
                                         TexXY[0],TexXY[1])
            pass
        else:
            polynum = self._identify_polynum_uv(polysurf,TexXY[0],TexXY[1])
            pass
        

        if polynum is None: 
            raise ValueError("polygonalsurface_intrinsicparameterization: polygonalsurface_texcoordparameterization.eval_uv_xyz(): failed to find polygon for (U,V) point\n")

        return polysurf._eval_xyz_polygonuv(self,polynum,u,v)


    #def eval_uv_xyz(self,polysurf,xyz):
    #    
    #    pass
    
    
    
    def __init__(self,**kwargs):
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
            
            setattr(self,kwarg,kwargs[kwarg])
            pass
        pass

    # intrinsicparameterizationparams
    # are ( (lowerleft_meaningfulunits_u, lowerleft_meaningfulunits_v), meaningfulunits_per_texcoord_u, meaningfulunits_per_texcoord_v )
    @classmethod
    def new(cls,polysurface,texcoord,texcoordidx,appearance,cadpartparams=None):
        # appearance is a spatialnde.cadpart.appearance.vrml_x3d_appearance subclass
        if cadpartparams is not None and "UV_ScalingParamsByTexURL" in cadpartparams and hasattr(appearance,"texture_url"):
            scalingparams = cadpartparams["UV_ScalingParamsByTexURL"]
            (lowerleft_meaningfulunits,meaningfulunits_per_texcoord) = scalingparams[appearance.texture_url]
            pass
        elif cadpartparams is not None and "UV_ScalingParamsBySurfaceNum" in cadpartparams:
            # Obsolete structure
            # I think it is safe to remove this
            assert(0)  # (if not would be constantly crashing here) 
            scalingparams = cadpartparams["UV_ScalingParamsBySurfaceNum"]
            (lowerleft_meaningfulunits,meaningfulunits_per_texcoord) = scalingparams[str(polysurface.surfaceid)]
            pass
        else:
            lowerleft_meaningfulunits=(0.0,0.0)
            meaningfulunits_per_texcoord=(1.0,1.0)
            pass
        
        return cls(texcoord=texcoord,
                   texcoordidx=texcoordidx,
                   lowerleft_meaningfulunits=lowerleft_meaningfulunits,
                   meaningfulunits_per_texcoord=meaningfulunits_per_texcoord,
                   cadpartparams=cadpartparams)
# intrinsicparameterizationparams=intrinsicparameterizationparams)


    
    pass


