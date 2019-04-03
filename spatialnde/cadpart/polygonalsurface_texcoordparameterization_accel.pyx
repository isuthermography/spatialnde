import numpy as np
cimport numpy as np
cimport cython

cimport libc.math

from libc.stdint cimport int32_t,uint32_t
#from libc.stdlib cimport calloc,malloc,free
from libc.stdio cimport fprintf,stderr


cdef extern from "../spatialnde_ctypes.h":
    ctypedef float float32_t
    ctypedef double float64_t
    pass

cdef extern from "../vecops.h":
    pass

cdef extern from "../geometry_accel_ops.h":
    pass

cdef extern from "polygonalsurface_texcoordparameterization_accel_ops.h":
    size_t enclosed_or_intersecting_polygons_2d_c(int32_t *polypool,
					          size_t polypoollen,
					          float64_t *texcoord,
					          uint32_t *vertexidx_indices,
					          uint32_t num_3d_polygons,
					          int32_t *texcoordidx,
					          uint32_t *numvertices,
					          uint32_t *texcoordredundant_numcopies,
					          uint32_t *texcoordredundant_polystartindexes,
                                                  uint32_t *texcoordredundant_polystartpolynum,
					          int32_t *texcoordredundant_texcoordidx,
					          float64_t minu,
					          float64_t minv,
					          float64_t maxu,
					          float64_t maxv,
					          int32_t *retpolys,
					          uint32_t *num_fully_enclosed) nogil

    int test_polynum_uv_c(uint32_t *vertexidx_indices,
			  size_t num_3d_polys,
			  uint32_t *numvertices,
			  int32_t *texcoordidx,
			  uint32_t *texcoordredundant_polystartindexes,
			  uint32_t *texcoordredundant_polystartpolynum,
			  int32_t *texcoordredundant_texcoordidx,
			  float64_t *texcoord,
			  float64_t u,
			  float64_t v,
			  size_t polynum) nogil

    int32_t identify_polynum_uv_c(uint32_t *vertexidx_indices,
				  size_t num_3d_polys,
				  uint32_t *numvertices,
				  int32_t *texcoordidx,
				  uint32_t *texcoordredundant_polystartindexes,
				  uint32_t *texcoordredundant_polystartpolynum,
				  int32_t *texcoordredundant_texcoordidx,
				  float64_t *texcoord,
				  int32_t *boxes,
				  int32_t *boxpolys,
				  float64_t *boxcoords,
				  float64_t u,
				  float64_t v,
				  int32_t candidate_polynum,
				  int32_t *boxpool_mem,
				  size_t boxpool_mem_nelem) nogil
    
    void evaluate_curvature_c(uint32_t *vertexidx_indices,
			 size_t num_3d_polys,
			  uint32_t *numvertices,
			  int32_t *vertexidx,
			  float64_t *vertices,
			  float64_t *refpoints,
			  uint32_t *texcoordredundant_polystartpolynum,
			  float64_t *inplanemats,
			  float64_t *inplane2texcoords,
			  float64_t *texcoord2inplane,
			  float64_t *principal_curvatures,
			  float64_t *curvature_tangent_axes,
			  size_t polynum,float64_t u, float64_t v,
			  float64_t *curvmatout) nogil

    void linelength_avgcurvature_c_array(
        uint32_t *vertexidx_indices,
        size_t num_3d_polys,
        uint32_t *numvertices,
        int32_t *vertexidx,
        float64_t *vertices,
        float64_t *refpoints,
        int32_t *texcoordidx,
        uint32_t *texcoordredundant_polystartindexes,
        uint32_t *texcoordredundant_polystartpolynum,
        int32_t *texcoordredundant_texcoordidx,
        float64_t *texcoord,
        int32_t *boxes,
	size_t nboxes,
        int32_t *boxpolys,
        float64_t *boxcoords,
        float64_t *inplanemats,        
        float64_t *inplane2texcoords,
        float64_t *texcoords2inplane,
        float64_t *principal_curvatures,
        float64_t *curvature_tangent_axes,
        float64_t du,float64_t dv,
	float64_t *u1,
	float64_t *v1,
	float64_t *u2,
	float64_t *v2,
	size_t numlines,
	float64_t *linelengthout,
	float64_t *avgcurvatureout) nogil

    void linelength_avgcurvature_meshbased_c_array(
        float32_t *curvmats, 
        float32_t *stepsizearray,
        size_t nv, 
        size_t nu, # number of u steps in curvmats or stepsizearray
        float64_t param_lowerleft_meaningfulunits_u,
        float64_t param_lowerleft_meaningfulunits_v,
        float64_t param_stepsize_u, # /* meaningful unit stepsize for curvmats/stepsizearray (nominal value), used to measure u1,v1,etc. */
        float64_t param_stepsize_v,
	float64_t du, # nominal resolution in u, meaning desired step size
	float64_t dv, # nominal resolution in v, desired step size, in texcoord [0-1] units
	float64_t *u1,
	float64_t *v1,
	float64_t *u2,
	float64_t *v2,
	size_t numlines,
	float64_t *linelengthout,
	float64_t *avgcurvatureout,
        float64_t *avgcrosscurvatureout) nogil
    
    void linelength_avgcurvature_mirroredbox_meshbased_c_array(
        float32_t *curvmats, # indexed c-style (v,u,2,2)  ... must cover entire range of [0,1] texture coordinates
        float32_t *stepsizearray,
        size_t nv, # number of v steps in curvmats or stepsizearray
        size_t nu, # number of u steps in curvmats or stepsizearray
        float64_t param_lowerleft_meaningfulunits_u,
        float64_t param_lowerleft_meaningfulunits_v,
        float64_t param_stepsize_u, #/* meaningful unit stepsize for curvmats/stepsizearray (nominal value), used to measure u1,v1,etc. */
        float64_t param_stepsize_v,
        float64_t boxu1,
        float64_t boxv1,
        float64_t boxu2,
        float64_t boxv2,
	float64_t du, # nominal resolution in u, meaning desired step size
	float64_t dv, # nominal resolution in v, desired step size, in texcoord [0-1] units
	float64_t *u1,
	float64_t *v1,
	float64_t *u2,
	float64_t *v2,
	size_t numlines,
	float64_t *linelengthout,
	float64_t *avgcurvatureout,
        float64_t *avgcrosscurvatureout,
        float64_t *thetaout) nogil
    
    pass


@cython.boundscheck(False)
def enclosed_or_intersecting_polygons_2d_accel(
        np.ndarray[np.int32_t,ndim=1,mode="c"] polypool,
        np.ndarray[np.float64_t,ndim=2,mode="c"] texcoord,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] vertexidx_indices,
        np.ndarray[np.int32_t,ndim=1,mode="c"] texcoordidx,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] numvertices,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] texcoordredundant_numcopies,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] texcoordredundant_polystartindexes,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] texcoordredundant_polystartpolynum,

        np.ndarray[np.int32_t,ndim=1,mode="c"] texcoordredundant_texcoordidx,
        float64_t minu,
        float64_t minv,
        float64_t maxu,
        float64_t maxv):

    cdef size_t num_returned_polys

    cdef np.ndarray[np.int32_t,ndim=1,mode="c"] retpolys
    cdef np.ndarray[np.uint32_t,ndim=1] num_fully_enclosed

    retpolys=np.empty(polypool.shape[0],dtype=np.int32)
    num_fully_enclosed=np.zeros(1,dtype=np.uint32)
    
    num_returned_polys = enclosed_or_intersecting_polygons_2d_c(<int32_t *>polypool.data,polypool.shape[0],<float64_t *>texcoord.data,<uint32_t *>vertexidx_indices.data,vertexidx_indices.shape[0],<int32_t *>texcoordidx.data,<uint32_t *>numvertices.data,<uint32_t *>texcoordredundant_numcopies.data,<uint32_t *>texcoordredundant_polystartindexes.data,<uint32_t *>texcoordredundant_polystartpolynum.data,<int32_t *>texcoordredundant_texcoordidx.data,minu,minv,maxu,maxv,<int32_t *>retpolys.data,<uint32_t *>num_fully_enclosed.data)
    
    return (num_fully_enclosed[0],retpolys[:num_returned_polys])

def _test_polynum_uv(
        np.ndarray[np.uint32_t,ndim=1,mode="c"] vertexidx_indices,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] numvertices,
        np.ndarray[np.int32_t,ndim=1,mode="c"] texcoordidx,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] texcoordredundant_polystartindexes,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] texcoordredundant_polystartpolynum,
        np.ndarray[np.int32_t,ndim=1,mode="c"] texcoordredundant_texcoordidx,
        np.ndarray[np.float64_t,ndim=2,mode="c"] texcoord,
        u,v,polynum):
    """This takes [0,1] u,v coordinates """

    return (test_polynum_uv_c(<uint32_t *>vertexidx_indices.data,
                              vertexidx_indices.shape[0],
                              <uint32_t *>numvertices.data,
                              <int32_t *>texcoordidx.data,
                              <uint32_t *>texcoordredundant_polystartindexes.data,
                              <uint32_t *>texcoordredundant_polystartpolynum.data,
                              <int32_t *>texcoordredundant_texcoordidx.data,
                              <float64_t *>texcoord.data,
                              u,v,polynum) != 0)


def _identify_polynum_uv(
        np.ndarray[np.uint32_t,ndim=1,mode="c"] vertexidx_indices,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] numvertices,
        np.ndarray[np.int32_t,ndim=1,mode="c"] texcoordidx,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] texcoordredundant_polystartindexes,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] texcoordredundant_polystartpolynum,
        np.ndarray[np.int32_t,ndim=1,mode="c"] texcoordredundant_texcoordidx,
        np.ndarray[np.float64_t,ndim=2,mode="c"] texcoord,
        np.ndarray[np.int32_t,ndim=2,mode="c"] boxes,
        np.ndarray[np.int32_t,ndim=1,mode="c"] boxpolys,
        np.ndarray[np.float64_t,ndim=2,mode="c"] boxcoords,
        u,v,candidate_polynum=None):

    cdef np.ndarray[np.int32_t,ndim=1,mode="c"] boxpool_mem


    if candidate_polynum is None:
        candidate_polynum=-1
        pass

    boxpool_mem=np.empty(boxes.shape[0],dtype=np.int32) # large enough by definition! (includes all possible boxes)

    polynum = identify_polynum_uv_c(<uint32_t *>vertexidx_indices.data,
                                    vertexidx_indices.shape[0],
                                    <uint32_t *>numvertices.data,
                                    <int32_t *>texcoordidx.data,
                                    <uint32_t *>texcoordredundant_polystartindexes.data,
                                    <uint32_t *>texcoordredundant_polystartpolynum.data,
                                    <int32_t *>texcoordredundant_texcoordidx.data,
                                    <float64_t *>texcoord.data,
                                    <int32_t *>boxes.data,
                                    <int32_t *>boxpolys.data,
                                    <float64_t *>boxcoords.data,
                                    u,v,candidate_polynum,
                                    <int32_t *>boxpool_mem.data,
                                    boxpool_mem.shape[0])
    if polynum < 0:
        return None

    return polynum


def _evaluate_curvature(
        np.ndarray[np.uint32_t,ndim=1,mode="c"] vertexidx_indices,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] numvertices,
        np.ndarray[np.int32_t,ndim=1,mode="c"] vertexidx,
        np.ndarray[np.float64_t,ndim=2,mode="c"] vertices,
        np.ndarray[np.float64_t,ndim=2,mode="c"] refpoints,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] texcoordredundant_polystartpolynum,
        np.ndarray[np.float64_t,ndim=3,mode="c"] inplanemats,
        np.ndarray[np.float64_t,ndim=3,mode="c"] inplane2texcoords,
        np.ndarray[np.float64_t,ndim=3,mode="c"] texcoord2inplane,
        np.ndarray[np.float64_t,ndim=2,mode="c"] principal_curvatures,
        np.ndarray[np.float64_t,ndim=3,mode="c"] curvature_tangent_axes,
        polynum, u, v):

    cdef np.ndarray[np.float64_t,ndim=2,mode="c"] curvmatout

    curvmatout=np.empty((2,2),dtype=np.float64)

    evaluate_curvature_c(<uint32_t *>vertexidx_indices.data,
                         vertexidx_indices.shape[0],
                         <uint32_t *>numvertices.data,
                         <int32_t *>vertexidx.data,
                         <float64_t *>vertices.data,
                         <float64_t *>refpoints.data,
                         <np.uint32_t *>texcoordredundant_polystartpolynum.data,
                         <float64_t *>inplanemats.data,
                         <float64_t *>inplane2texcoords.data,
                         <float64_t *>texcoord2inplane.data,
                         <float64_t *>principal_curvatures.data,
                         <float64_t *>curvature_tangent_axes.data,
                         polynum,u,v,
                         <float64_t *>curvmatout.data)
    return curvmatout


    
def _linelength_avgcurvature(
        np.ndarray[np.uint32_t,ndim=1,mode="c"] vertexidx_indices,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] numvertices,
        np.ndarray[np.int32_t,ndim=1,mode="c"] vertexidx,
        np.ndarray[np.float64_t,ndim=2,mode="c"] vertices,
        np.ndarray[np.float64_t,ndim=2,mode="c"] refpoints,
        np.ndarray[np.int32_t,ndim=1,mode="c"] texcoordidx,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] texcoordredundant_polystartindexes,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] texcoordredundant_polystartpolynum,
        np.ndarray[np.int32_t,ndim=1,mode="c"] texcoordredundant_texcoordidx,
        np.ndarray[np.float64_t,ndim=2,mode="c"] texcoord,
        np.ndarray[np.int32_t,ndim=2,mode="c"] boxes,
        np.ndarray[np.int32_t,ndim=1,mode="c"] boxpolys,
        np.ndarray[np.float64_t,ndim=2,mode="c"] boxcoords,
        np.ndarray[np.float64_t,ndim=3,mode="c"] inplanemats,        
        np.ndarray[np.float64_t,ndim=3,mode="c"] inplane2texcoords,
        np.ndarray[np.float64_t,ndim=3,mode="c"] texcoords2inplane,
        np.ndarray[np.float64_t,ndim=2,mode="c"] principal_curvatures,
        np.ndarray[np.float64_t,ndim=3,mode="c"] curvature_tangent_axes,
        du,dv,
        u1,v1,u2,v2):

    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] u1_unwrapped
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] v1_unwrapped
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] u2_unwrapped
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] v2_unwrapped
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] linelengthout
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] avgcurvatureout
    
    if not isinstance(u1,np.ndarray):
        u1_unwrapped=np.empty(1,dtype=np.float64)
        v1_unwrapped=np.empty(1,dtype=np.float64)
        u2_unwrapped=np.empty(1,dtype=np.float64)
        v2_unwrapped=np.empty(1,dtype=np.float64)
        linelengthout=np.empty(1,dtype=np.float64)
        avgcurvatureout=np.empty(1,dtype=np.float64)
        n=1;
        
        u1_unwrapped[0]=u1
        v1_unwrapped[0]=v1
        u2_unwrapped[0]=u2
        v2_unwrapped[0]=v2

        resultshape=None
        pass
    elif len(u1.shape) > 1:
        if (u1.shape != v1.shape or
            u1.shape != u2.shape or
            u1.shape != v2.shape or
            v1.shape != u2.shape or
            v1.shape != v2.shape or
            u2.shape != v2.shape):
            
            shapes=np.array((u1.shape,v1.shape,u2.shape,v2.shape),dtype=np.uint32)
            maxsize = np.max(shapes,axis=0)
            resultshape=maxsize
            zerosbroadcast = np.zeros(maxsize,dtype=np.float64)
            u1=u1+zerosbroadcast
            v1=v1+zerosbroadcast
            u2=u2+zerosbroadcast
            v2=v2+zerosbroadcast
            pass
        else:
            resultshape=u1.shape
            pass
        u1_unwrapped=u1.ravel()
        v1_unwrapped=v1.ravel()
        u2_unwrapped=u2.ravel()
        v2_unwrapped=v2.ravel()

        n=u1_unwrapped.shape[0]

        linelengthout=np.empty(n,dtype=np.float64)
        avgcurvatureout=np.empty(n,dtype=np.float64)
        pass
    else:
        # single axis
        n=u1.shape[0]
        resultshape=(n,)
        
        u1_unwrapped=u1
        v1_unwrapped=v1
        u2_unwrapped=u2
        v2_unwrapped=v2

        linelengthout=np.empty(n,dtype=np.float64)
        avgcurvatureout=np.empty(n,dtype=np.float64)

        pass
    
    
    linelength_avgcurvature_c_array(
        <uint32_t *> vertexidx_indices.data,
        vertexidx_indices.shape[0],
        <uint32_t *> numvertices.data,
        <int32_t *> vertexidx.data,
        <float64_t *>vertices.data,
        <float64_t *>refpoints.data,
        <int32_t *> texcoordidx.data,
        <uint32_t *> texcoordredundant_polystartindexes.data,
        <uint32_t *> texcoordredundant_polystartpolynum.data,
        <int32_t *> texcoordredundant_texcoordidx.data,
        <float64_t *> texcoord.data,
        <int32_t *> boxes.data,
        boxes.shape[0],
        <int32_t *> boxpolys.data,
        <float64_t *> boxcoords.data,
        <float64_t *> inplanemats.data,        
        <float64_t *> inplane2texcoords.data,
        <float64_t *> texcoords2inplane.data,
        <float64_t *> principal_curvatures.data,
        <float64_t *> curvature_tangent_axes.data,
        du,dv,
        <float64_t *>u1_unwrapped.data,
        <float64_t *>v1_unwrapped.data,
        <float64_t *>u2_unwrapped.data,
        <float64_t *>v2_unwrapped.data,
        n,
        <float64_t *>linelengthout.data,
        <float64_t *>avgcurvatureout.data)

    if resultshape is None:
        return (linelengthout[0],avgcurvatureout[0])
    else:
        return (linelengthout.reshape(*resultshape),
                avgcurvatureout.reshape(*resultshape))
    pass





def _linelength_avgcurvature_meshbased(
        np.ndarray[np.float32_t,ndim=4,mode="c"] curvmats,
        np.ndarray[np.float32_t,ndim=3,mode="c"] stepsizearray,
        param_lowerleft_meaningfulunits_u,
        param_lowerleft_meaningfulunits_v,
        param_stepsize_u, # /* meaningful unit stepsize for curvmats/stepsizearray (nominal value), used to measure u1,v1,etc. */
        param_stepsize_v,
        du,dv,
        u1,v1,u2,v2):

    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] u1_unwrapped
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] v1_unwrapped
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] u2_unwrapped
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] v2_unwrapped
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] linelengthout
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] avgcurvatureout
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] avgcrosscurvatureout
    
    if not isinstance(u1,np.ndarray):
        u1_unwrapped=np.empty(1,dtype=np.float64)
        v1_unwrapped=np.empty(1,dtype=np.float64)
        u2_unwrapped=np.empty(1,dtype=np.float64)
        v2_unwrapped=np.empty(1,dtype=np.float64)
        linelengthout=np.empty(1,dtype=np.float64)
        avgcurvatureout=np.empty(1,dtype=np.float64)
        avgcrosscurvatureout=np.empty(1,dtype=np.float64)
        n=1;
        
        u1_unwrapped[0]=u1
        v1_unwrapped[0]=v1
        u2_unwrapped[0]=u2
        v2_unwrapped[0]=v2

        resultshape=None
        pass
    elif len(u1.shape) > 1:
        if (u1.shape != v1.shape or
            u1.shape != u2.shape or
            u1.shape != v2.shape or
            v1.shape != u2.shape or
            v1.shape != v2.shape or
            u2.shape != v2.shape):
            
            shapes=np.array((u1.shape,v1.shape,u2.shape,v2.shape),dtype=np.uint32)
            maxsize = np.max(shapes,axis=0)
            resultshape=maxsize
            zerosbroadcast = np.zeros(maxsize,dtype=np.float64)
            u1=u1+zerosbroadcast
            v1=v1+zerosbroadcast
            u2=u2+zerosbroadcast
            v2=v2+zerosbroadcast
            pass
        else:
            resultshape=u1.shape
            pass
        u1_unwrapped=u1.ravel()
        v1_unwrapped=v1.ravel()
        u2_unwrapped=u2.ravel()
        v2_unwrapped=v2.ravel()

        n=u1_unwrapped.shape[0]

        linelengthout=np.empty(n,dtype=np.float64)
        avgcurvatureout=np.empty(n,dtype=np.float64)
        avgcrosscurvatureout=np.empty(n,dtype=np.float64)
        pass
    else:
        # single axis
        n=u1.shape[0]
        resultshape=(n,)
        
        u1_unwrapped=u1
        v1_unwrapped=v1
        u2_unwrapped=u2
        v2_unwrapped=v2

        linelengthout=np.empty(n,dtype=np.float64)
        avgcurvatureout=np.empty(n,dtype=np.float64)
        avgcrosscurvatureout=np.empty(n,dtype=np.float64)

        pass
    
    assert(curvmats.shape[0]==stepsizearray.shape[0])
    assert(curvmats.shape[1]==stepsizearray.shape[1])
    assert(curvmats.shape[2]==2)
    assert(stepsizearray.shape[2]==2)
    assert(curvmats.shape[3]==2)

    linelength_avgcurvature_meshbased_c_array(
        <float32_t *>curvmats.data,
        <float32_t *>stepsizearray.data,
        curvmats.shape[0], curvmats.shape[1],
        param_lowerleft_meaningfulunits_u,
        param_lowerleft_meaningfulunits_v,
        param_stepsize_u, #/* meaningful unit stepsize for curvmats/stepsizearray (nominal value), used to measure u1,v1,etc. */
        param_stepsize_v,
        du, dv,
        <float64_t *>u1_unwrapped.data,
        <float64_t *>v1_unwrapped.data,
        <float64_t *>u2_unwrapped.data,
        <float64_t *>v2_unwrapped.data,
        n,
        <float64_t *>linelengthout.data,
        <float64_t *>avgcurvatureout.data,
        <float64_t *>avgcrosscurvatureout.data)
    
    if resultshape is None:
        return (linelengthout[0],avgcurvatureout[0],avgcrosscurvatureout[0])
    else:
        return (linelengthout.reshape(*resultshape),
                avgcurvatureout.reshape(*resultshape),
                avgcrosscurvatureout.reshape(*resultshape))
    pass





def _linelength_avgcurvature_mirroredbox_meshbased(
        np.ndarray[np.float32_t,ndim=4,mode="c"] curvmats,
        np.ndarray[np.float32_t,ndim=3,mode="c"] stepsizearray,
        param_lowerleft_meaningfulunits_u,
	param_lowerleft_meaningfulunits_v,
        param_stepsize_u, # /* meaningful unit stepsize for curvmats/stepsizearray (nominal value), used to measure u1,v1,etc. */
        param_stepsize_v,
        boxu1,boxv1,boxu2,boxv2,
        du,dv,
        u1,v1,u2,v2):

    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] u1_unwrapped
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] v1_unwrapped
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] u2_unwrapped
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] v2_unwrapped
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] linelengthout
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] avgcurvatureout
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] avgcrosscurvatureout
    cdef np.ndarray[np.float64_t,ndim=1,mode="c"] thetaout
    
    if not isinstance(u1,np.ndarray):
        u1_unwrapped=np.empty(1,dtype=np.float64)
        v1_unwrapped=np.empty(1,dtype=np.float64)
        u2_unwrapped=np.empty(1,dtype=np.float64)
        v2_unwrapped=np.empty(1,dtype=np.float64)
        linelengthout=np.empty(1,dtype=np.float64)
        avgcurvatureout=np.empty(1,dtype=np.float64)
        avgcrosscurvatureout=np.empty(1,dtype=np.float64)
        thetaout=np.empty(1,dtype=np.float64)
        n=1;
        
        u1_unwrapped[0]=u1
        v1_unwrapped[0]=v1
        u2_unwrapped[0]=u2
        v2_unwrapped[0]=v2

        resultshape=None
        pass
    elif len(u1.shape) > 1:
        if (u1.shape != v1.shape or
            u1.shape != u2.shape or
            u1.shape != v2.shape or
            v1.shape != u2.shape or
            v1.shape != v2.shape or
            u2.shape != v2.shape):
            
            shapes=np.array((u1.shape,v1.shape,u2.shape,v2.shape),dtype=np.uint32)
            maxsize = np.max(shapes,axis=0)
            resultshape=maxsize
            zerosbroadcast = np.zeros(maxsize,dtype=np.float64)
            u1=u1+zerosbroadcast
            v1=v1+zerosbroadcast
            u2=u2+zerosbroadcast
            v2=v2+zerosbroadcast
            pass
        else:
            resultshape=u1.shape
            pass
        u1_unwrapped=u1.ravel()
        v1_unwrapped=v1.ravel()
        u2_unwrapped=u2.ravel()
        v2_unwrapped=v2.ravel()

        n=u1_unwrapped.shape[0]

        linelengthout=np.empty(n,dtype=np.float64)
        avgcurvatureout=np.empty(n,dtype=np.float64)
        avgcrosscurvatureout=np.empty(n,dtype=np.float64)
        thetaout=np.empty(n,dtype=np.float64)
        pass
    else:
        # single axis
        n=u1.shape[0]
        resultshape=(n,)
        
        u1_unwrapped=u1
        v1_unwrapped=v1
        u2_unwrapped=u2
        v2_unwrapped=v2

        linelengthout=np.empty(n,dtype=np.float64)
        avgcurvatureout=np.empty(n,dtype=np.float64)
        avgcrosscurvatureout=np.empty(n,dtype=np.float64)
        thetaout=np.empty(n,dtype=np.float64)

        pass
    
    assert(curvmats.shape[0]==stepsizearray.shape[0])
    assert(curvmats.shape[1]==stepsizearray.shape[1])
    assert(curvmats.shape[2]==2)
    assert(stepsizearray.shape[2]==2)
    assert(curvmats.shape[3]==2)

    linelength_avgcurvature_mirroredbox_meshbased_c_array(
        <float32_t *>curvmats.data,
        <float32_t *>stepsizearray.data,
        curvmats.shape[0], curvmats.shape[1],
        param_lowerleft_meaningfulunits_u,
        param_lowerleft_meaningfulunits_v,
        param_stepsize_u, #/* meaningful unit stepsize for curvmats/stepsizearray (nominal value), used to measure u1,v1,etc. */
        param_stepsize_v,
        boxu1,boxv1,boxu2,boxv2,
        du, dv,
        <float64_t *>u1_unwrapped.data,
        <float64_t *>v1_unwrapped.data,
        <float64_t *>u2_unwrapped.data,
        <float64_t *>v2_unwrapped.data,
        n,
        <float64_t *>linelengthout.data,
        <float64_t *>avgcurvatureout.data,
        <float64_t *>avgcrosscurvatureout.data,
        <float64_t *>thetaout.data)
    
    if resultshape is None:
        return (linelengthout[0],avgcurvatureout[0],avgcrosscurvatureout[0],thetaout[0])
    else:
        return (linelengthout.reshape(*resultshape),
                avgcurvatureout.reshape(*resultshape),
                avgcrosscurvatureout.reshape(*resultshape),
                thetaout.reshape(*resultshape))
    pass

