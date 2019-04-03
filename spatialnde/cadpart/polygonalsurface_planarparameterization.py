import sys
import copy
import json

import numpy as np
import numpy.linalg    
#import scipy
#import scipy.linalg


def vecnorm(a,axis=-1):
    
    # compatibility for old versions of Numpy
    return np.sqrt(np.sum(a**2.0,axis))


def uniquevertices(vertices,tolerance):
    nv = vertices.shape[0]
    
    vertexdists = vecnorm(vertices.reshape(nv,1,3)-vertices.reshape(1,nv,3),axis=2)
    
    vertexnumbers=np.empty(nv,dtype='i4')
    
    
    
    # assign an index to each vertex
    uniquevertexcoords=[]
    # uniquevertexsets=[]
    uniquecnt=0
    
    vertexcnt=0
    
    vertexmask = np.zeros(nv,dtype=np.bool)  # We set an element here to True once we assign it to a unique vertex
    
    while vertexcnt < nv:
        if vertexmask[vertexcnt]:
            # already handled
            vertexcnt+=1
            continue
        
        # Find other vertices within tolerance
        
        withintolerance = vertexdists[vertexcnt,:] <= tolerance
        setmask = withintolerance & (~vertexmask)
        setmask[vertexcnt] = True
        #vertexset=frozenset(np.where(setmask))
        
        if np.count_nonzero(withintolerance & vertexmask):
            sys.stderr.write("WARNING: Tolerance spheres overlap (is tolerance too big or too small?) \n")
            pass
        
        
        vertexmask |= withintolerance
        vertexcoords = vertices[vertexcnt,:]
        vertexnumbers[setmask] = uniquecnt
        
        
        
        uniquevertexcoords.append(vertexcoords)
        #uniquevertexsets.append(vertexset)
        uniquecnt+=1
        
        vertexcnt+=1
        pass

    
    
    return (vertexnumbers,np.array(uniquevertexcoords,dtype='d'))

def findadjacentfacets(facetsbyvertexnum):
    n=facetsbyvertexnum.shape[0]


    ## Find edges
    #edges1 = facetsbyvertexnum[:,0:2]
    #edges2 = facetsbyvertexnum[:,1:3]
    #edges3 = 
    
    # Now try to find which facets are adjacent by shared vertices
    # create n x n matrix of truth values for vertex i matches vertex i
    facetvertexmatch_ii=(facetsbyvertexnum[:,0].reshape(n,1)-facetsbyvertexnum[:,0].reshape(1,n)).astype(np.bool) 

    facetvertexmatch_ij=(facetsbyvertexnum[:,0].reshape(n,1)-facetsbyvertexnum[:,1].reshape(1,n)).astype(np.bool) 

    facetvertexmatch_ik=(facetsbyvertexnum[:,0].reshape(n,1)-facetsbyvertexnum[:,2].reshape(1,n)).astype(np.bool) 

    facetvertexmatch_ji=(facetsbyvertexnum[:,1].reshape(n,1)-facetsbyvertexnum[:,0].reshape(1,n)).astype(np.bool) 

    facetvertexmatch_jj=(facetsbyvertexnum[:,1].reshape(n,1)-facetsbyvertexnum[:,1].reshape(1,n)).astype(np.bool) 

    facetvertexmatch_jk=(facetsbyvertexnum[:,1].reshape(n,1)-facetsbyvertexnum[:,2].reshape(1,n)).astype(np.bool) 

    facetvertexmatch_ki=(facetsbyvertexnum[:,2].reshape(n,1)-facetsbyvertexnum[:,0].reshape(1,n)).astype(np.bool) 

    facetvertexmatch_kj=(facetsbyvertexnum[:,2].reshape(n,1)-facetsbyvertexnum[:,1].reshape(1,n)).astype(np.bool) 

    facetvertexmatch_kk=(facetsbyvertexnum[:,2].reshape(n,1)-facetsbyvertexnum[:,2].reshape(1,n)).astype(np.bool) 


    have_common_edge=((facetvertexmatch_ii & (facetvertexmatch_jk | facetvertexmatch_kj))
                      | (facetvertexmatch_jj & (facetvertexmatch_ik | facetvertexmatch_ki))
                      | (facetvertexmatch_kk & (facetvertexmatch_ij | facetvertexmatch_ji))
                      | (facetvertexmatch_ij & (facetvertexmatch_jk | facetvertexmatch_kj)) 
                      | (facetvertexmatch_ij & (facetvertexmatch_ik | facetvertexmatch_ki))
                      | (facetvertexmatch_ji & (facetvertexmatch_jk | facetvertexmatch_kj))
                      | (facetvertexmatch_ji & (facetvertexmatch_ik | facetvertexmatch_ki))
                      | (facetvertexmatch_ik & (facetvertexmatch_kj | facetvertexmatch_jk))
                      | (facetvertexmatch_ik & (facetvertexmatch_ij | facetvertexmatch_ji))
                      | (facetvertexmatch_ki & (facetvertexmatch_kj | facetvertexmatch_jk))
                      | (facetvertexmatch_ki & (facetvertexmatch_ij | facetvertexmatch_ji))
                      | (facetvertexmatch_jk & (facetvertexmatch_ki | facetvertexmatch_ik))
                      | (facetvertexmatch_jk & (facetvertexmatch_ji | facetvertexmatch_ij))
                      | (facetvertexmatch_kj & (facetvertexmatch_ki | facetvertexmatch_ik))
                      | (facetvertexmatch_kj & (facetvertexmatch_ji | facetvertexmatch_ij))
                      | (facetvertexmatch_ij & facetvertexmatch_ji)
                      | (facetvertexmatch_jk & facetvertexmatch_kj)
                      | (facetvertexmatch_ik & facetvertexmatch_ki))
      

    
    #facetvertexmatches = np.zeros((n,n),dtype=np.uint32)
    
    #facetvertexmatches += facetvertexmatch_ii.astype(np.uint32)
    #facetvertexmatches += facetvertexmatch_ij.astype(np.uint32)
    #facetvertexmatches += facetvertexmatch_ik.astype(np.uint32)
    #facetvertexmatches += facetvertexmatch_ji.astype(np.uint32)
    #facetvertexmatches += facetvertexmatch_jj.astype(np.uint32)
    #facetvertexmatches += facetvertexmatch_jk.astype(np.uint32)
    #facetvertexmatches += facetvertexmatch_ki.astype(np.uint32)
    #facetvertexmatches += facetvertexmatch_kj.astype(np.uint32)
    #facetvertexmatches += facetvertexmatch_kk.astype(np.uint32)

    #adjacent = (facetvertexmatches == 2)
    #identity = (facetvertexmatches == 3)

    #errors = (facetvertexmatches > 3)
    #if errors.any():
    #    raise ValueError("Triangles smaller than tolerance!")
    #    pass
    #
    ## identity should only apply on the diagonal 
    #(identity_row,identity_col)=np.where(identity)
    #if (identity_row != identity_col).any():
    #    raise ValueError("Identical triangles!")

    return have_common_edge
    



def identify_surfaces(facets,tolerance):
    # Break facet array (n_triangles x 4 x 3) array
    # index #1 represents 0 for normal, then 1,2,3 for vertices
    # ... down into specific surfaces

    # Fairly dumb algorithm, not super scalable

    # identify unique vertices (within tolerance)

    n=facets.shape[0]
    vertices = facets[:,1:,:].reshape(n*3,3)


    # identify unique vertices
    (vertexnumbers,uniquevertexcoords) = uniquevertices(vertices,tolerance)

    facetsbyvertexnum = vertexnumbers.reshape(n,3) # This is a nfaces x 3 array of integers that gives unique vertex number of each of the three vertices of each triangle.

    adjacent=findadjacentfacets(facetsbyvertexnum)
    # adjacent is an nxn bool array indicating adjacency
    # FYI this algorithm assumes that triangle edges are not subdivided

    used_mask=np.zeros(n,dtype=np.bool)

    facetnum=0
    surfaces=[]

    while facetnum < n:
        if used_mask[facetnum]:
            facetnum+=1
            continue

        surface_mask=np.zeros(n,dtype=np.bool)

        surface_mask[facetnum]=True
        surface_normal = facets[facetnum,0,:]


        # try to add adjacent facets
        addedfacets=True

        while addedfacets:
            addedfacets=False

            surface_adjacent_facet_mask = adjacent[surface_mask,:][0,:] & (~used_mask)
            surface_adjacent_facets = np.where(surface_adjacent_facet_mask)[0]

            # Compute inner product of adjacent surface normals relative to our surface normal
            normal_innerprod = np.inner(facets[surface_adjacent_facets,0,:],surface_normal)

            # normals should not be > 80 deg away from surface_normal
            # i.e. innerprod > cos(80 deg)
            facets_use = normal_innerprod > np.cos(80.0*np.pi/180.0)

            if np.count_nonzero(facets_use) > 0:
                facets_use_nums = surface_adjacent_facets[facets_use]

                # Add these facets to the surface
                surface_mask[facets_use_nums]=True
                addedfacets=True

                # mark these facets as used
                used_mask[facets_use_nums] = True

                pass
            pass


        surfaces.append(np.where(surface_mask)[0])


        facetnum+=1
        pass

    # surfaces is now a list of integer arrays indicating
    # which facets are in each surface
    
    return (surfaces,facetsbyvertexnum,uniquevertexcoords)


class polygonalsurface_planarparameterization(object):
    
    intrinsicparameterizationparams=None
    u = None # 3 element array with unit vector in the u direction for this surface
    v = None # 3 element array with unit vector in the v direction for this surface
    meannormal=None # 3 element array representing nominal normal vector for the surface
    centercoords=None # 3 element array representing a nominal center coordinate for the surface
    
    #def eval_uv(self,polysurf,polysurf_vertexids):
    #    coords=polysurf.vertices[polysurf_vertexids,:] # n x 3 matrix
    #    relcoords=coords-self.centercoords.reshape(1,3)
    #    
    #    return np.array((np.inner(relcoords,self.u),np.inner(relcoords,self.v)),dtype='d').T

    def eval_xyz_uv(self):
        # Need to be able to evaluate (x,y,z) coords given
        # (u,v) coordinates within this parameterization
        raise NotImplementedError()

    def buildprojinfo(self,polysurf):
        sys.stderr.write("polygonalsurface_intrinsicparameterization: buildprojinfo() not implemented. Will not be able to perform projection calculations\n")
        pass
    
    def eval_texcoord_polygonvertex(self,polysurf,polysurf_polynum,polysurf_vertexnum):
        # Can supply vectors as polysurf_polynum and/or polysurf_vertexnum

        firstidx=polysurf.vertexidx_indices[polysurf_polynum]
        
        
        vertexids = polysurf.vertexidx[firstidx+polysurf_vertexnum]
        coords=polysurf.vertices[vertexids.ravel(),:] # n x 3 matrix
        relcoords=coords-self.centercoords.reshape(1,3)

        res_shape = [2]
        if not np.isscalar(polysurf_vertexnum):
            res_shape.insert(0,len(polysurf_vertexnum))
            pass
        if not np.isscalar(polysurf_polynum):
            res_shape.insert(0,len(polysurf_polynum))
            pass
        return np.array((np.inner(relcoords,self.u),np.inner(relcoords,self.v)),dtype='d').T.reshape(*res_shape)

        
    
    def __init__(self,**kwargs):
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
            
            setattr(self,kwarg,kwargs[kwarg])
            pass
        pass

    @classmethod
    def new(cls,polysurface,cadpartparams=None):

        intrinsicparameterizationparams=None
        if cadpartparams is not None and "intrinsicparameterization" in cadpartparams:
            intrinsicparameterizationparams=cadpartparams["intrinsicparameterization"]
            pass
        
        
        if intrinsicparameterizationparams is None:
            # evaluate surfacecentercoords
            # Get list of unique vertex ides used by this surface
            centercoords=np.mean(polysurface.vertices,axis=0) # Get mean location

            
            meannormal=np.mean(polysurface.normals,axis=0)
            # normalize
            meannormal /= np.linalg.norm(meannormal)
            
            # evaluate u -- use in-surfacemeannormal-plane component
            # of vector between first two vertices of first polygon
            
            # Get u
            firsttwovertices = polysurface.vertexidx[:2]
            verticesdirec = polysurface.vertices[firsttwovertices[1],:]-polysurface.vertices[firsttwovertices[0],:]

            # normalize vector
            verticesdirec /= np.linalg.norm(verticesdirec)

            # subtract component parallel to surfacemeannormal
            verticesdirec = verticesdirec - np.inner(verticesdirec,meannormal)*meannormal
            
            # renormalize
            verticesdirec /= np.linalg.norm(verticesdirec)
            
            u=verticesdirec

            # Get v so that u cross v = meannormal
            # meannormal cross u = v
            v=np.cross(meannormal,u)
            v /= np.linalg.norm(verticesdirec) # ensure normalized

            # Create JSON definition of our parameterization
            paramsdict = { "centercoords": tuple(centercoords.astype(object)),
                           "meannormal": tuple(meannormal.astype(object)),
                           "u": tuple(u.astype(object)),
                           "v": tuple(v.astype(object)),
            }
            intrinsicparameterizationparams=copy.deepcopy(paramsdict)
            pass
        else:
            # intrinsicparameterizationparams given
            paramsdict=copy.deepcopy(intrinsicparameterizationparams)
            
            centercoords=np.array(paramsdict["centercoords"],dtype='d')
            meannormal=np.array(paramsdict["meannormal"],dtype='d')
            u=np.array(paramsdict["u"],dtype='d')
            v=np.array(paramsdict["v"],dtype='d')
            
            pass
        
        return cls(centercoords=centercoords,
                   meannormal=meannormal,
                   u=u,
                   v=v,
                   intrinsicparameterizationparams=intrinsicparameterizationparams)
    
        
    
    
    pass


