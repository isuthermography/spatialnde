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

cdef extern from "polygonalsurface_accel_ops.h":
    int enclosed_or_intersecting_polygons_3d_c(int32_t *polypool,
					       size_t polypoollen,
					       float64_t *vertices,
					       uint32_t *vertexidx_indices,
					       int32_t num_polygons,
					       int32_t *vertexidx,
					       uint32_t *numvertices,
					       float64_t *inplanemats,
					       float64_t *facetnormals,
					       float64_t *box_v0,
					       float64_t *box_v1,
					       int32_t *retpolys,
					       uint32_t *num_fully_enclosed) nogil
    pass



@cython.boundscheck(False)
def enclosed_or_intersecting_polygons_3d_accel(
        np.ndarray[np.int32_t,ndim=1,mode="c"] polypool,
        np.ndarray[np.float64_t,ndim=2,mode="c"] vertices,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] vertexidx_indices,
        np.ndarray[np.int32_t,ndim=1,mode="c"] vertexidx,
        np.ndarray[np.uint32_t,ndim=1,mode="c"] numvertices,
        np.ndarray[np.float64_t,ndim=3,mode="c"] inplanemats,
        np.ndarray[np.float64_t,ndim=2,mode="c"] facetnormals,
        np.ndarray[np.float64_t,ndim=1,mode="c"] box_v0,
        np.ndarray[np.float64_t,ndim=1,mode="c"] box_v1):

    cdef size_t num_returned_polys

    cdef np.ndarray[np.int32_t,ndim=1,mode="c"] retpolys
    cdef np.ndarray[np.uint32_t,ndim=1] num_fully_enclosed

    retpolys=np.empty(polypool.shape[0],dtype=np.int32)
    num_fully_enclosed=np.zeros(1,dtype=np.uint32)
    
    num_returned_polys = enclosed_or_intersecting_polygons_3d_c(<int32_t *>polypool.data,polypool.shape[0],<float64_t *>vertices.data,<uint32_t *>vertexidx_indices.data,vertexidx_indices.shape[0],<int32_t *>vertexidx.data,<uint32_t *>numvertices.data,<float64_t *>inplanemats.data,<float64_t *>facetnormals.data,<float64_t *>box_v0.data,<float64_t *>box_v1.data,<int32_t *>retpolys.data,<uint32_t *>num_fully_enclosed.data)
    
    return (num_fully_enclosed[0],retpolys[:num_returned_polys])
