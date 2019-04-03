import numpy as np
cimport numpy as np
cimport cython

cimport libc.math

from libc.stdint cimport int32_t,uint32_t
#from libc.stdlib cimport calloc,malloc,free
from libc.stdio cimport fprintf,stderr

cdef extern from "spatialnde_ctypes.h":
    ctypedef float float32_t
    ctypedef double float64_t
    pass

cdef extern from "vecops.h":
    pass

cdef extern from "geometry_accel_ops.h":
    int point_in_polygon_2d_c(float64_t *vertices_rel_point,size_t numvertices) nogil 
    int polygon_intersects_box_3d_c(float64_t *box_v0, float64_t *box_v1, float64_t *vertices, int32_t *optional_vertexidx,size_t nvertices, float64_t *inplanemat, float64_t *facetnormal) nogil
    int polygon_intersects_box_2d_c(float64_t *box_v0,float64_t *box_v1,float64_t *vertices,int32_t *optional_vertexidx,size_t numvertices)
    int box_inside_polygon_2d_c(float64_t *box_v0,float64_t *box_v1,float64_t *vertices,int32_t *optional_vertexidx,size_t num_vertices)
    
    pass
    




@cython.boundscheck(False)
def polygon_intersects_box_2d(
        np.ndarray[np.float64_t,ndim=1,mode="c"] box_v0,
        np.ndarray[np.float64_t,ndim=1,mode="c"] box_v1,
        np.ndarray[np.float64_t,ndim=2,mode="c"] vertices):
    cdef size_t num_vertices

    num_vertices=vertices.shape[0]
    
    return polygon_intersects_box_2d_c(<float64_t *>box_v0.data,<float64_t *>box_v1.data, <float64_t *>vertices.data, NULL,num_vertices) != 0

@cython.boundscheck(False)
def polygon_intersects_box_3d(
        np.ndarray[np.float64_t,ndim=1,mode="c"] box_v0,
        np.ndarray[np.float64_t,ndim=1,mode="c"] box_v1,
        np.ndarray[np.float64_t,ndim=2,mode="c"] vertices,
        np.ndarray[np.float64_t,ndim=2,mode="c"] inplanemat,
        np.ndarray[np.float64_t,ndim=1,mode="c"] facetnormal):
    cdef size_t num_vertices

    num_vertices=vertices.shape[0]
    
    return polygon_intersects_box_3d_c(<float64_t *>box_v0.data,<float64_t *>box_v1.data, <float64_t *>vertices.data, NULL, num_vertices, <float64_t *>inplanemat.data, <float64_t *>facetnormal.data) != 0




@cython.boundscheck(False)
def point_in_polygon_2d(
        np.ndarray[np.float64_t,ndim=2,mode="c"] vertices_rel_point):
    cdef size_t num_vertices

    num_vertices=vertices_rel_point.shape[0]
    
    return point_in_polygon_2d_c(<float64_t *>vertices_rel_point.data,num_vertices) != 0
    
def box_inside_polygon_2d(
        np.ndarray[np.float64_t,ndim=1,mode="c"] box_v0,
        np.ndarray[np.float64_t,ndim=1,mode="c"] box_v1,
        np.ndarray[np.float64_t,ndim=2,mode="c"] vertices):

    return (box_inside_polygon_2d_c(<float64_t *>box_v0.data,<float64_t *>box_v1.data,<float64_t *>vertices.data,NULL,vertices.shape[0]) != 0)

