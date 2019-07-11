import sys
import time
try: 
    from numpy.linalg import svd
except ImportError:
    from scipy.linalg import svd
    pass

import numpy as np
cimport numpy as np
cimport cython

cimport libc.math

from libc.stdint cimport int32_t,uint32_t,int64_t,uint64_t
from libc.stdlib cimport calloc,malloc,free,exit
from libc.stdio cimport fprintf,stderr

from scipy.ndimage import binary_opening

# from cython.parallel import prange,parallel


cdef extern from "spatialnde_ctypes.h":
    ctypedef float float32_t
    ctypedef double float64_t
    pass


cdef extern from "vecops.h":
    void multmat23vec(float64_t *mat,float64_t *vec,float64_t *out) nogil
    void multmatvec4(float64_t *mat,float64_t *vec,float64_t *out) nogil
    void multmatvec2(float64_t *mat,float64_t *vec,float64_t *out) nogil
    float64_t dotvecvec3(float64_t *vec1,float64_t *vec2) nogil
    float64_t scalevec3(float64_t coeff,float64_t *vec1,float64_t *out) nogil
    float64_t scalevec2(float64_t coeff,float64_t *vec1,float64_t *out) nogil
    void subvecvec3(float64_t *vec1,float64_t *vec2,float64_t *out) nogil
    void addvecscaledvec3(float64_t *vec1,float64_t coeff, float64_t *vec2,float64_t *out) nogil
    void normalize_wcoord4(float64_t *vec) nogil
    double to_unit_vector4(float64_t *vec) nogil
    double to_unit_vector3(float64_t *vec) nogil
    pass

cdef extern from "imageprojection_ops.h":
    ctypedef uint32_t atomicpixel_t # actually _Atomic

    float32_t *atomicbuf_sanitize(atomicpixel_t *buf,size_t n) nogil
    atomicpixel_t *atomicbuf_prepare(float32_t *buf,size_t n) nogil
    float32_t atomicpixel_load(atomicpixel_t *var) nogil
    float32_t atomicpixel_nonatomicload(atomicpixel_t *var) nogil
    atomicpixel_store(atomicpixel_t *var,float32_t val) nogil
    atomicpixel_nonatomicstore(atomicpixel_t *var,float32_t val) nogil
    void atomicpixel_accumulate(atomicpixel_t *var,float toadd) nogil
    void atomicbuf_zero(atomicpixel_t *buf,size_t n) nogil
    void normalize_atomic_buffer(atomicpixel_t *tonormalize,atomicpixel_t *normbuf,size_t n) nogil
    void atomicbuffer_dropthreshold(atomicpixel_t *buf,atomicpixel_t *thresholdtestbuf,float32_t threshold,size_t n) nogil

    
    struct projectionsurface:
        atomicpixel_t *imgbuf  # volatile imgbuf must be big enough for all frames of uv parameterization... must be normalized (divided) by validitybuf
        atomicpixel_t *weightingbuf # volatile u-v space weighting coefficients
        atomicpixel_t *validitybuf # volatile validitybuf should be big enough for one frame of uv parameterization
        atomicpixel_t *angle_of_incidence_factor_buf_uv # volatile... filled out during prepare. Needs to be normalized by validitybuf just like imgbuf. 
        size_t nx
        size_t ny

        float64_t u_meaningfulunits_per_texcoord
        float64_t v_meaningfulunits_per_texcoord
        
        float64_t *tosurfaceframe  # 4x4
        float64_t *vecfromsurfaceframe #3x3
        float64_t *boxcoords
        size_t nboxes
        int32_t *boxes
        int32_t *boxpolys
        uint32_t *vertexidx_indices
        int32_t *vertexidx
        uint32_t *numvertices
 	#size_t max_vertices
        size_t npolys
        float64_t *vertices

        # !!!*** Need to generate facetnormals, refpoints, inplanemats, and inplane2texcoords!!!***
        float64_t *facetnormals
        float64_t *refpoints
        float64_t *inplanemats
        float64_t *inplane2texcoords  # 2x3 projective matrix
        float64_t *maxradius # n_polys

        pass
    void iterate_diffuse_mismatch(float32_t *upperbound,
    	                          float32_t *value,
                                  float32_t *extra_buf,
				  double delta_tbar,
				  uint64_t ny,
				  uint64_t nx,
				  double dybarsq,
				  double dxbarsq,
				  uint64_t ntbar) nogil

    float64_t ray_to_plane_distance(float64_t *starting_point,float64_t *ray_direc, float64_t *planepoint, float64_t *planenormal) nogil
    void ray_to_plane_raydirec_shift(float64_t *starting_point,float64_t *ray_direc, float64_t *planepoint, float64_t *planenormal,float64_t *ray_direc_deriv,float64_t *deriv) nogil
    int ray_intersects_polygon(float64_t *vertexarray,uint32_t *vertexidx_indices,int32_t *vertexidx,uint32_t *numvertices,float64_t *normal_vector,  float64_t *refpoints, float64_t *inplanemats, float64_t *starting_point, float64_t *ray_direc, float64_t ray_to_plane_dist,int include_edges) nogil

    void find_ray_intersections(projectionsurface *surf,size_t surfacecnt,size_t boxnum,float64_t *starting_point,float64_t *ray_direc,float32_t *zbufpt,uint32_t *surfaceidpt,uint32_t *facetidpt) nogil
    
    int ray_box_intersection(float64_t min_dist, float64_t *boxcoords, float64_t *starting_point, float64_t *ray_direc) nogil
    
    void projecttoimgbuf(float32_t pixelval,float32_t pixelweighting,float64_t *uvcoords,float64_t *uvcoords_deriv_horiz, float64_t *uvcoords_deriv_vert,atomicpixel_t *imgbuf, atomicpixel_t *weightingbuf,atomicpixel_t *validitybuf,size_t framenum,size_t imgbuf_nx,size_t imgbuf_ny,float64_t min_radius_uv_pixels,float64_t min_radius_source_pixels,float64_t bandwidth_fraction) nogil

    void evaluate_zbuffer(projectionsurface *surfacearray,
		          size_t src_nx,
		          size_t src_ny,
		          size_t numsurfaces,
		          float64_t *cam_mtx_Ut,
		          float64_t *cam_mtx_Vtsubinv,
		          float64_t s1,float64_t s2,
		          float64_t V31,float64_t V32,
		          float32_t *imagedata_zbuffer,
		          uint32_t *imagedata_surfaceid,
		          uint32_t *imagedata_facetid,
                          float32_t *imagedata_angleofincidence,
                          float32_t *imagedata_angleofincidence_weighting,
                          #float32_t *imagedata_projectzup,
                          #float32_t *imagedata_projectzdown,
                          #float32_t *imagedata_projectzleft,
                          #float32_t *imagedata_projectzright
                          float32_t *imagedata_horiz_ray_intersect_shift_deriv_cam_z,
                          float32_t *imagedata_vert_ray_intersect_shift_deriv_cam_z,
                          float32_t *imagedata_uvcoords) nogil
    void project_image(projectionsurface *surfacearray,
                       float32_t *framedata,
                       size_t framecnt,
		       size_t src_nx,
		       size_t src_ny,
		       float64_t *cam_mtx_Ut,
		       float64_t *cam_mtx_Vtsubinv,
		       float64_t s1,float64_t s2,
		       float64_t V31,float64_t V32,
		       float32_t *imagedata_zbuffer,
		       uint32_t *imagedata_surfaceid,
		       uint32_t *imagedata_facetid,
                       float32_t *imagedata_weighting,
                       int debug) nogil
    
    pass

pass


cdef class imageprojection_params:
    cdef public int nsurfaces
    cdef projectionsurface *surfacearray

    cdef public object surfacelist  # Python list containing (surface,parameterizationdictentry,weightingbuf,angleofincidencefactorbuf_in_uv_frame), for purposes of reference counting
    
    
    # parameters of the projection
    cdef public size_t src_nx
    cdef public size_t src_ny
    

    cdef public np.ndarray cam_mtx_Ut
    cdef public np.ndarray cam_mtx_Vtsub
    cdef public np.ndarray cam_mtx_Vtsubinv

    cdef public object tosurfaceframes  # Store list of tosurfaceframe matrices here too for purposes of reference counting
    cdef public object vecfromsurfaceframes  # Store list of vecfromsurfaceframe matrices here too for purposes of reference counting


    cdef public float64_t s1
    cdef public float64_t s2
    cdef public float64_t V31
    cdef public float64_t V32


    cdef public np.ndarray imagedata_zbuffer
    cdef public np.ndarray imagedata_surfaceid
    cdef public np.ndarray imagedata_facetid 
    cdef public np.ndarray imagedata_angleofincidence 
    cdef public np.ndarray imagedata_angleofincidence_weighting
    cdef public np.ndarray imagedata_weighting # Weighting coefficients in camera frame
    cdef public np.ndarray imagedata_horiz_ray_intersect_shift_deriv_cam_z 
    cdef public np.ndarray imagedata_vert_ray_intersect_shift_deriv_cam_z 
    cdef public np.ndarray imagedata_uvcoords
    cdef public object uv_weightingbufs


    def __init__(self, size_t nsurfaces):
        self.surfacelist=[]
        pass
    
    def __cinit__(self, size_t nsurfaces):
        self.surfacearray=<projectionsurface *>calloc(sizeof(projectionsurface)*nsurfaces,1)
        self.nsurfaces=nsurfaces
        pass

    def __dealloc__(self):
        free(self.surfacearray)
        pass
    
    

    pass
    
# Required library initializations
np.import_array()



#def diffuse_mismatch(mismatcharray,mismatchvalue,nomismatchvalue,diffusiondistance_horizpixels,diffusiondistance_vertpixels):
    # Given a 2D boolean mismatch array and a desired diffusion
    # distance (measured in horizontal and vertical pixels),
    # Simulate a diffusion process that is initialized to
    # nomismatchvalue, but with mismatchvalue overlaid as a Dirichlet
    # (specified value) boundary condition. The diffusion
    # equation is then iterated to diffuse the specified diffusion distance.

def diffuse_mismatch(np.ndarray[np.float32_t,mode="c",ndim=2] upper_bound,initialvalue,diffusiondistance_horizpixels,diffusiondistance_vertpixels):
    cdef np.ndarray[np.float32_t,mode="c",ndim=2] value
    cdef np.ndarray[np.float32_t,mode="c",ndim=2] extra_buf
    
    # Given a spatially distributed lower bound, and an initial value
    # or initial value field, and a desired diffusion 
    # distance (measured in horizontal and vertical pixels),
    # Simulate a diffusion process that keeps the field
    # less than or equal to the upper_bound everywhere at each timestep, 
    # (specified value) boundary condition. The diffusion
    # equation is then iterated to diffuse the specified diffusion distance.

    # ***!!!! Ideally diffuse_mismatch would use the local stepsizes at each point in the parameterization, but it doesn't now!!!

    #
    # This assumes that the edges of the array should be listed as mismatches.
    # Edge behavior may be incorrect if they are not. 
    #
    # The diffusion distance is defined such that one-dimensional
    # diffusion of a pulse gives that half-distance for a 1/e
    # reduction from peak, 
    # i.e. exp(-x^2/4alphat) = exp(-1)
    # i.e. -x^2/4alphat = -1
    # i.e. x^2 = 4alphat
    # i.e. t = x^2/4alpha
    # Where x (x0 from here on) is the specified diffusion distance.

    # ... Nondimensionalization diffusion equation:
    # alpha del^2 T = dT/dt, alpha in m^2/s
    # let xbar = x/x0  -> dxbar = dx/x0
    # let tbar = t*alpha/x0^2  -> dtbar = dt*alpha/x0^2
    # del = d/dx = d/(x0dxbar)
    # let delbar = d/dxbar = x0*del
    # Diffusion eq. is now
    # alpha/x0^2 delbar^2 T = dT/(dtbar*x0^2/alpha) = (alpha/x0^2) * dT/dtbar
    # delbar^2 T = dT/dtbar
    # Where xbar=x/x0 and tbar = t*alpha/x0^2.
    # Per above, t should go from 0 up to x0^2/4alpha,
    # so tbar should go from 0 to 1/4

    # per https://me.ucsb.edu/~moehlis/APC591/tutorials/tutorial5/node3.html
    # simple finite-difference forward iteration only convergent
    # when delta_tbar/delta_xbar^2 <= 0.5
    #
    # ... Given x and x0 in pixels, delta_xbar = 1/x0
    # Therefore delta_tbar <= 0.5 * delta_xbar^2
    # or delta_tbar <= 0.5/x0^2

    biggest_diffusiondistance=max(diffusiondistance_horizpixels,diffusiondistance_vertpixels)
    
    delta_tbar = 0.25/biggest_diffusiondistance**2  # half required value
    ntbar_desired = 0.25/delta_tbar

    # require at least 8 iterations so a single iteration gained or lost
    # is not a big deal.... If we need more efficiency here
    # we could probably size the iterations to be an even multiple
    # naturally...
    if ntbar_desired < 8:
        delta_tbar = delta_tbar * ntbar_desired/8.0
        ntbar_desired = 8.0
        pass
    ntbar = int(round(ntbar_desired))

    tbar = np.arange(ntbar+1,dtype='d')*delta_tbar
    value=np.zeros((upper_bound.shape[0],upper_bound.shape[1]),dtype='f')+initialvalue
    extra_buf=np.zeros((upper_bound.shape[0],upper_bound.shape[1]),dtype='f')

    dxbar = 1.0/diffusiondistance_horizpixels
    dybar = 1.0/diffusiondistance_vertpixels
    # perform differences  (dx = 1.0 pixels)

    #    assert(value.shape[0]==upper_bound.shape[0])
    
    iterate_diffuse_mismatch(<np.float32_t *>upper_bound.data,
                             <np.float32_t *>value.data,
                             <np.float32_t *>extra_buf.data,
			     delta_tbar,
			     value.shape[0],
			     value.shape[1],
			     dybar**2,
			     dxbar**2,
			     ntbar)
    #for tcnt in range(ntbar):
    #    exceedance_zone = value > upper_bound
    #    value[ exceedance_zone ]=upper_bound[ exceedance_zone ]  # Assign specified value condition where the value exceeds the upper bound
    #    
    #
    #    # Iterate time 
    #    # T_t+ = T_t- + dt*( (T_x+ -2T + T_x-)/dx^2 + (T_y+ -2T + T_y-)/dy^2 )
    #    value[1:-1,1:-1]=value[1:-1,1:-1] + delta_tbar * ( ((value[1:-1,2:] - 2*value[1:-1,1:-1] + value[1:-1,:-2])/dxbar**2) + ((value[2:,1:-1] - 2*value[1:-1,1:-1] + value[:-2,1:-1])/dybar**2) )
    #    pass
    return value


@cython.boundscheck(False)
def imageprojection_prepare_float(projectionmodel,cameraframe,partlist,parameterizationdict,NaNmap,edgeblurriness_pixels,uv_weightingblur_distance):
    # NaNmap was imageshape...
    
    # np.ndarray[np.float64_t,ndim=2,mode="c"]
    # cameraframe is a coordframe object
    # part list is a list of part objects, e.g. from an assembly
    # parameterizationdict is a dictionary, indexed by surface object id,
    # of tuples (parameterization,imagebuf,validitybuf), where parameterization
    # is the parameterization object or None to represent the
    # intrinsic parameterization of that surface
    # edgeblurriness_pixels is rough slight overestimate of the blurring of edges in
    # the image to be projected, in pixels (used to reduce weighting
    # of possibly erroneous data near an edge)
    # uv_weightingblur_distance is distance for the weighting blur of edges in
    # (u,v) space (used to reduce weighting
    # of possibly erroneous or incompatible data near an edge or
    # shadow. For thermography it should be a lateral thermal diffusion length
    # For simple optical projection 0 is fine. 
    
    # imagedata is top-row-first

    # Also parameterizations of multiple surfaces that end up
    # in the same space should be unified.
    #
    #
    # So the parameterization structure could have a place
    # to drop an image to be filled out, or something.
    # and multiple parameterizations can refer to the same
    # thing
    # In c++ we can use std::shared_ptr for reference counting.

    # We can't cvWarpPerspective from a distorted image, and
    # if we undistorted first, the transform matrix might
    # need to be modified (?)
    #
    # So at this point, we require the undistortion to
    # have already been applied
    assert(projectionmodel.undistortion_already_applied)
    
    
    cdef np.ndarray[np.float64_t,ndim=2,mode="c"] new_camera_mtx = projectionmodel.new_camera_mtx
    # cdef np.ndarray[np.float64_t,ndim=1,mode="c"] distcoeffs = projectionmodel.calibration.distcoeffs

    # Could probably speed things up for multi-frame by
    # doing all of the perspective calculations only once...


    # Count surfaces
    cdef size_t nsurfaces=0
    cdef size_t imgbuf_nx,imgbuf_ny
    # loop through assembly
    for ndepart in partlist:
        
        # loop through surfaces of part
        
        for surface in ndepart.implpart.surfaces:
            if id(surface) in parameterizationdict:
                nsurfaces+=1
                pass
            pass
        pass
    
    params=imageprojection_params(nsurfaces)
    params.tosurfaceframes=[]
    params.vecfromsurfaceframes=[]
    
    cdef np.ndarray[np.float32_t,mode="c",ndim=2] imgbuf
    cdef atomicpixel_t *imgbuf_ptr
    cdef np.ndarray[np.float32_t,mode="c",ndim=2] weightingbuf
    cdef np.ndarray[np.float32_t,mode="c",ndim=2] validitybuf
    cdef atomicpixel_t *validitybuf_ptr
    cdef np.ndarray[np.float32_t,mode="c",ndim=2] angle_of_incidence_factor_buf_uv
    cdef atomicpixel_t *angle_of_incidence_factor_buf_uv_ptr
    cdef np.ndarray[np.float64_t,mode="c",ndim=2] tosurfaceframe
    cdef np.ndarray[np.float64_t,mode="c",ndim=2] vecfromsurfaceframe

    cdef np.ndarray[np.float64_t,mode="c",ndim=2] boxcoords
    cdef np.ndarray[np.int32_t,mode="c",ndim=2] boxes
    cdef np.ndarray[np.int32_t,mode="c",ndim=1] boxpolys
    cdef np.ndarray[np.uint32_t,mode="c",ndim=1] vertexidx_indices
    cdef np.ndarray[np.int32_t,mode="c",ndim=1] vertexidx
    cdef np.ndarray[np.uint32_t,mode="c",ndim=1] numvertices
    cdef np.ndarray[np.float64_t,mode="c",ndim=2] vertices

    cdef np.ndarray[np.float64_t,mode="c",ndim=2] facetnormals
    cdef np.ndarray[np.float64_t,mode="c",ndim=2] refpoints
    cdef np.ndarray[np.float64_t,mode="c",ndim=3] inplanemats
    #cdef np.ndarray[np.float64_t,mode="c",ndim=3] outofplanevec # outofplanevec is same as facetnormal
    
    cdef np.ndarray[np.float64_t,mode="c",ndim=3] inplane2texcoords
    cdef np.ndarray[np.float64_t,mode="c",ndim=1] maxradius
    
    # Fill out surfacearray
    cdef size_t surfacecnt=0
    # loop through assembly
    for ndepart in partlist:
        
        # loop through surfaces of part
        
        for surface in ndepart.implpart.surfaces:
            if id(surface) in parameterizationdict:

                (parameterization,imgbuf,validitybuf)=parameterizationdict[id(surface)] # used to include angleofincidencebuf

                #weightingbuf=np.ones(imgbuf.shape,dtype='f')  # weightingbuf generated below after z-buffering
                angle_of_incidence_factor_buf_uv=np.empty((validitybuf.shape[0],validitybuf.shape[1]),dtype='f')                  
                params.surfacelist.append([surface,parameterizationdict[id(surface)],None,angle_of_incidence_factor_buf_uv ])  # Here we keep reference to surface and parameterizationdict and (once assigned) weightingbuf entries so its stuff doesn't get free'd
                
                tosurfaceframe = (cameraframe.transformationmatrix(ndepart.frame)).astype(np.float64) # 4x4 matrix that will transform camera coordinates into part coordinates 
                params.tosurfaceframes.append(tosurfaceframe) # ensure it doesn't get freed

                assert(tosurfaceframe.shape[0]==4 and tosurfaceframe.shape[1]==4)

                # tosurfaceframe matrix should be strictly a rotation and translation 
                assert(tosurfaceframe[3,0]==0.0)
                assert(tosurfaceframe[3,1]==0.0)
                assert(tosurfaceframe[3,2]==0.0)
                assert(tosurfaceframe[3,3]==1.0)

                vecfromsurfaceframe = np.ascontiguousarray(tosurfaceframe[:3,:3].transpose()) # transpose of upper left 3x3 of tosurfaceframe
                params.vecfromsurfaceframes.append(vecfromsurfaceframe) # Keep variable around by object reference

                # check that transpose is inverse
                assert(np.linalg.norm((np.dot(tosurfaceframe[:3,:3],vecfromsurfaceframe)-np.eye(3)).ravel()) < 1e-5)
                
                boxcoords=surface.boxcoords
                boxes=surface.boxes
                boxpolys=surface.boxpolys
                vertexidx_indices=surface.vertexidx_indices
                vertexidx=surface.vertexidx
                numvertices=surface.numvertices
                vertices=surface.vertices

                facetnormals=surface.facetnormals
                refpoints=surface.refpoints
                inplanemats=surface.inplanemats
                #outofplanevec=surface.outofplanevec
                
                

                if parameterization is None:
                    parameterization=surface.intrinsicparameterization
                    pass
                
                inplane2texcoords=parameterization.inplane2texcoords

                maxradius=surface.maxradius

                nx=imgbuf.shape[1]
                ny=imgbuf.shape[0]
                assert(nx==validitybuf.shape[1])
                assert(ny==validitybuf.shape[0])
                #assert(nx==angleofincidencebuf.shape[1])
                #assert(ny==angleofincidencebuf.shape[0])
                params.surfacearray[surfacecnt].nx=nx
                params.surfacearray[surfacecnt].ny=ny

                params.surfacearray[surfacecnt].u_meaningfulunits_per_texcoord=parameterization.meaningfulunits_per_texcoord[0]
                params.surfacearray[surfacecnt].v_meaningfulunits_per_texcoord=parameterization.meaningfulunits_per_texcoord[1]
                
                # Here we are hardwired into the meshed representation
                #
                # We want to get away from this for C++:
                #
                # Specifically, rather than copying key parameters here,
                # we should get a pointer to the implpart internal
                # method that can do this calculation for us.
                # Then we can call that method on those surfaces
                # rather than using the hardcoded engine here.
                #    That is more practical in C++ because in the
                # Cython environment, calling Python code is very
                # expensive by comparison

                # clear out imgbuf and validitybuf
                #imgbuf[::]=0.0;
                imgbuf_ptr=<atomicpixel_t *>imgbuf.data
                atomicbuf_zero(imgbuf_ptr,imgbuf.size)
		
                #validitybuf[::]=0.0;
                validitybuf_ptr=<atomicpixel_t *>validitybuf.data
                atomicbuf_zero(validitybuf_ptr,validitybuf.size)

                angle_of_incidence_factor_buf_uv_ptr=<atomicpixel_t *>angle_of_incidence_factor_buf_uv.data
                atomicbuf_zero(angle_of_incidence_factor_buf_uv_ptr,angle_of_incidence_factor_buf_uv.size)
                
                params.surfacearray[surfacecnt].imgbuf=imgbuf_ptr
                
                #params.surfacearray[surfacecnt].weightingbuf=<float32_t *>weightingbuf.data
                params.surfacearray[surfacecnt].weightingbuf=NULL
                params.surfacearray[surfacecnt].validitybuf=validitybuf_ptr
                params.surfacearray[surfacecnt].angle_of_incidence_factor_buf_uv=angle_of_incidence_factor_buf_uv_ptr
                params.surfacearray[surfacecnt].tosurfaceframe=<float64_t *>tosurfaceframe.data
                params.surfacearray[surfacecnt].vecfromsurfaceframe=<float64_t *>vecfromsurfaceframe.data

                params.surfacearray[surfacecnt].boxcoords=<float64_t *>boxcoords.data
                params.surfacearray[surfacecnt].nboxes=boxcoords.shape[0]
                params.surfacearray[surfacecnt].boxes=<int32_t *>boxes.data
                params.surfacearray[surfacecnt].boxpolys=<int32_t *>boxpolys.data
                params.surfacearray[surfacecnt].vertexidx_indices=<uint32_t *>vertexidx_indices.data
                params.surfacearray[surfacecnt].vertexidx=<int32_t *>vertexidx.data
                params.surfacearray[surfacecnt].numvertices=<uint32_t *>numvertices.data
 		#params.surfacearray[surfacecnt].max_vertices=vertexids.shape[1]
                params.surfacearray[surfacecnt].npolys=vertexidx_indices.shape[0]
                params.surfacearray[surfacecnt].vertices=<float64_t *>vertices.data

                params.surfacearray[surfacecnt].facetnormals=<float64_t *>facetnormals.data
                params.surfacearray[surfacecnt].refpoints=<float64_t *>refpoints.data
                params.surfacearray[surfacecnt].inplanemats=<float64_t *>inplanemats.data
                #surfacearray[surfacecnt].outofplanevec=<float64_t *>outofplanevec.data
                params.surfacearray[surfacecnt].inplane2texcoords=<float64_t *>inplane2texcoords.data
                params.surfacearray[surfacecnt].maxradius=<float64_t *>maxradius.data
                

                surfacecnt+=1
                pass
            pass
        pass

    assert(surfacecnt==nsurfaces)
    
    # cdef size_t numsurfaces=surfacecnt

    imageshape=NaNmap.shape
    
    assert(len(imageshape)==2)

    # Image is presumed to be stored raster style C indexed, so coordinates are (y,x)
    params.src_ny=imageshape[0]
    params.src_nx=imageshape[1]

    
        
    # define z-buffer, currently at infinity
    params.imagedata_zbuffer = (np.ones((imageshape[0],imageshape[1]),dtype='f')*np.inf).astype(np.float32)
    # define array of which surface id
    params.imagedata_surfaceid = np.zeros((imageshape[0],imageshape[1]),dtype='u4')
    # define array of which facet id
    params.imagedata_facetid = np.zeros((imageshape[0],imageshape[1]),dtype='u4')
    # define array of angles of incidence
    params.imagedata_angleofincidence = np.zeros((imageshape[0],imageshape[1]),dtype='f')
    # define array of weightings of camera image
    params.imagedata_angleofincidence_weighting = np.ones((imageshape[0],imageshape[1]),dtype='f')
    

    # Define array of z derivatives with respect to adjacent rays 
    params.imagedata_horiz_ray_intersect_shift_deriv_cam_z = np.zeros((imageshape[0],imageshape[1]),dtype='f')
    params.imagedata_vert_ray_intersect_shift_deriv_cam_z = np.zeros((imageshape[0],imageshape[1]),dtype='f')

    # Place to store the (u,v) coordinates of the pixels
    params.imagedata_uvcoords=np.zeros((imageshape[0],imageshape[1],2),dtype='f')
    
    
    ## loop through assembly
    #for ndepart in partlist:
    #    
    #    # loop through surfaces of part
    #
    #    for surface in ndepart.implpart.surfaces:
    #        if id(surface) in parameterizationdict:
    #            (parameterization,imgbuf,validitybuf,angleofincidencebuf)=parameterizationdict[id(surface)]

                    # Find intersection point from (u,v):
                    # [u]                                      [X]
                    # [v]   =  new_camera_matrix*tocameraframe*[Y]
                    # [1]                                      [Z]
                    #                                          [1]
                    #
                    #surfaces.append(projectionsurface(parameterization=parameterization,
                    #                                  imgbuf=imgbuf,
                    #                                  validitybuf=validitybuf,
                    #                                  tocameraframe=tocameraframe,
                    #                                  tosurfaceframe=tosurfaceframe,
                    #                                  boxcoords=surface.boxcoords,
                    #boxpolys=surface.boxpolys,
                    #                                  vertexids=surface.vertexids,
                    #vertices=surface.vertices
                    #))
                    #pass
                #pass
            #pass
            
    # See doc/Projection.pdf
    (U,S,Vt) = svd(new_camera_mtx[:2,:3])
    params.cam_mtx_Ut=np.ascontiguousarray(U.T)
    params.cam_mtx_Vtsub=np.zeros((2,2),dtype='d')
    params.cam_mtx_Vtsub[:,:]=Vt[0:2,0:2]
    params.cam_mtx_Vtsubinv=np.zeros((2,2),dtype='d')
    params.cam_mtx_Vtsubinv[:,:]=np.linalg.inv(params.cam_mtx_Vtsub)
    params.s1=S[0]
    params.s2=S[1]
    params.V31=Vt[0,2]
    params.V32=Vt[1,2]
    zbufstart=time.time()

    cdef projectionsurface *surfacearray=params.surfacearray

    cdef size_t src_nx=params.src_nx
    cdef size_t src_ny=params.src_ny
    cdef np.ndarray[np.float64_t,ndim=2,mode="c"] cam_mtx_Ut = params.cam_mtx_Ut
    cdef np.ndarray[np.float64_t,ndim=2,mode="c"] cam_mtx_Vtsubinv = params.cam_mtx_Vtsubinv

    cdef float64_t s1=params.s1
    cdef float64_t s2=params.s2
    cdef float64_t V31=params.V31
    cdef float64_t V32=params.V32
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_zbuffer=params.imagedata_zbuffer
    cdef np.ndarray[np.uint32_t, ndim=2, mode="c"] imagedata_surfaceid=params.imagedata_surfaceid
    cdef np.ndarray[np.uint32_t, ndim=2, mode="c"] imagedata_facetid=params.imagedata_facetid
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_angleofincidence=params.imagedata_angleofincidence
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_angleofincidence_weighting=params.imagedata_angleofincidence_weighting
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_horiz_ray_intersect_shift_deriv_cam_z=params.imagedata_horiz_ray_intersect_shift_deriv_cam_z
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_vert_ray_intersect_shift_deriv_cam_z=params.imagedata_vert_ray_intersect_shift_deriv_cam_z
    cdef np.ndarray[np.float32_t, ndim=3, mode="c"] imagedata_uvcoords=params.imagedata_uvcoords

    
    with nogil:
        evaluate_zbuffer(surfacearray,
                         src_nx,
                         src_ny,
                         nsurfaces,
                         <float64_t *>cam_mtx_Ut.data,
                         <float64_t *>cam_mtx_Vtsubinv.data,
                         params.s1,params.s2,
                         params.V31,params.V32,
                         <float32_t *>imagedata_zbuffer.data, 
                         <uint32_t *>imagedata_surfaceid.data,
                         <uint32_t *>imagedata_facetid.data,
                         <float32_t *>imagedata_angleofincidence.data,
                         <float32_t *>imagedata_angleofincidence_weighting.data, # initial weighting is just 1.0... initial validitybuf will be from projection operation alone. This will get filled in with the desired angle-of-incidence weighting, but that waiting will NOT affect the validitybuf output
                         <float32_t *>imagedata_horiz_ray_intersect_shift_deriv_cam_z.data,
                         <float32_t *>imagedata_vert_ray_intersect_shift_deriv_cam_z.data,
                         <float32_t *>imagedata_uvcoords.data)
        pass

    # NOTE: evaluate_zbuffer drops baseline projection validity into validitybuf
    
    print("zbuffer calc time: %f s" % (time.time()-zbufstart))
    print("num facets: %d" % (np.count_nonzero(imagedata_facetid)))

    ## Debugging: disable z prediction
    #for surfacecnt in range(nsurfaces):
    #    surf_imgbuf=params.surfacelist[surfacecnt][1][1]
    #    weightingbuf = (np.array(surf_imgbuf.copy()*0.0,dtype='f')+1.0).astype(np.float32)
    #
    #
    #    # store weightingbuf in params.surfacelist, which will keep the objects reference-counted
    #    params.surfacelist[surfacecnt][2]=weightingbuf
    #
    #
    #    # store C pointer as well
    #    params.surfacearray[surfacecnt].weightingbuf=<float32_t *>weightingbuf.data
    #params.imagedata_weighting=np.ones((src_ny,src_nx),dtype='f')
    #cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_weighting0=params.imagedata_weighting
    #
    #
    #return params

    # Find z prediction mismatches
    zmismatch=np.zeros((imageshape[0],imageshape[1]),dtype=np.bool)


    # Evaluate z predictions to right
    # Mark mismatch when angular prediction error exceeds 30deg.
    # Could also use distance prediction error, but this will
    # only occur without angular prediction error at steep angles
    # of incidence, which are windowed out separately
 
 
    fx=new_camera_mtx[0,0] # x focal length of camera
    fy=new_camera_mtx[1,1] # x focal length of camera
    
    
    #zpredict_right=zbuffer + imagedata_horiz_ray_intersect_shift_deriv_cam_z
    predict_leftright_angle = np.arctan2(imagedata_horiz_ray_intersect_shift_deriv_cam_z[1:-1,1:-1],imagedata_zbuffer[1:-1,1:-1]/fx)  # atan(dz/(z/fx)) where the z/fx is the dx corresponding to one pixel at the detected plane
    #predict_leftright_dist = np.sqrt(1.0+imagedata_horiz_ray_intersect_shift_deriv_cam_z**2.0)
    actual_angle_right=np.arctan2(imagedata_zbuffer[1:-1,1:-1]-imagedata_zbuffer[1:-1,2:],imagedata_zbuffer[1:-1,1:-1]/fx)
    actual_angle_left=np.arctan2(imagedata_zbuffer[1:-1,:-2]-imagedata_zbuffer[1:-1,1:-1],imagedata_zbuffer[1:-1,1:-1]/fx)  # atan(dz/(z/fx)) where the z/fx is the dx corresponding to one pixel at the detected plane

    #actual_dist_right=np.sqrt(1.0+(zbuffer[:,2:]-zbuffer[:,1:-1])**2.0)
    #actual_dist_left=np.sqrt(1.0+(zbuffer[:,1:-1]-zbuffer[:,:2])**2.0)

    zmismatch[1:-1,1:-1] = zmismatch[1:-1,1:-1] | (~(np.abs(predict_leftright_angle-actual_angle_right) < np.pi/5.0))  # NOTE: ~(a < 30deg) is used rather than a >= 30 deg so that NaNs will always show as mismatches 

    zmismatch[1:-1,1:-1] = zmismatch[1:-1,1:-1] | (~(np.abs(predict_leftright_angle-actual_angle_left) < np.pi/5.0))


    # Same for up/down (image y axis)
    # Note that the sign of this presumes that going down is an index increase
    # (raster scan order)
    predict_updown_angle = np.arctan2(imagedata_vert_ray_intersect_shift_deriv_cam_z[1:-1,1:-1],imagedata_zbuffer[1:-1,1:-1]/fx)
    actual_angle_down=np.arctan2(imagedata_zbuffer[1:-1,1:-1]-imagedata_zbuffer[2:,1:-1],imagedata_zbuffer[1:-1,1:-1]/fy)
    actual_angle_up=np.arctan2(imagedata_zbuffer[:-2,1:-1]-imagedata_zbuffer[1:-1,1:-1],imagedata_zbuffer[1:-1,1:-1]/fy)

    zmismatch[1:-1,1:-1] = zmismatch[1:-1,1:-1] | (~(np.abs(predict_updown_angle-actual_angle_down) < np.pi/5.0))

    zmismatch[1:-1,1:-1] = zmismatch[1:-1,1:-1] | (~(np.abs(predict_updown_angle-actual_angle_up) < np.pi/5.0))
    

    # edges always mismatch
    zmismatch[0,:]=1.0
    zmismatch[:,0]=1.0
    zmismatch[-1,:]=1.0
    zmismatch[:,-1]=1.0

    zmismatch[np.isnan(NaNmap)]=1.0
 
    # perform binary opening to eliminate isolated
    # mismatch pixels
    zmismatch=binary_opening(zmismatch);


    # edges still always mismatch
    zmismatch[0,:]=1.0
    zmismatch[:,0]=1.0
    zmismatch[-1,:]=1.0
    zmismatch[:,-1]=1.0

    # Ok... We now have boolean mismatch picture in camera coordinates.

    # Need to diffuse the mismatch per specified distance x0... and store
    # in camera image weighting factor
    #    params.imagedata_weighting=params.imagedata_weighting*diffuse_mismatch(zmismatch,0.01,1.0,edgeblurriness_pixels,edgeblurriness_pixels)
    # angle-of-incidence weighting is now done in (u,v) space
    # so here we just do the mismatch dffusion

    # We store the diffused mismatch in imagedata_weighting
    upper_bound=np.ones(zmismatch.shape,dtype='f')
    upper_bound[zmismatch]=0.01
    #params.imagedata_weighting=diffuse_mismatch(zmismatch,0.01,1.0,edgeblurriness_pixels,edgeblurriness_pixels)
    params.imagedata_weighting=diffuse_mismatch(upper_bound,1.0,edgeblurriness_pixels,edgeblurriness_pixels)
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_weighting=params.imagedata_weighting
    
    # Also translate the mismatch onto the u,v parameterizations
    # for diffusion-based weighting
    #
    # so here we use no weighting, but the input array is imagedata_weighting (the result of diffusing zmismatch)

    # Perform a projection of the diffused zmismatch in order to do this.
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] unit_imagedata_weighting=np.ones((src_ny,src_nx),dtype='f')
    

    
    for surfacecnt in range(nsurfaces):
        #params.surfacearray[surfacecnt].angle_of_incidence_factor_buf_uv

        # normalize angle_of_incidence_factor_buf_uv by validitybuf
        #surf_validitybuf=params.surfacelist[surfacecnt][1][2]
        #surf_angleofincidencefactorbufuv=params.surfacelist[surfacecnt][3]
        #surf_angleofincidencefactorbufuv /=  surf_validitybuf
        imgbuf_nx=params.surfacearray[surfacecnt].nx
        imgbuf_ny=params.surfacearray[surfacecnt].ny

        
        

        normalize_atomic_buffer(params.surfacearray[surfacecnt].angle_of_incidence_factor_buf_uv,params.surfacearray[surfacecnt].validitybuf,imgbuf_nx*imgbuf_ny)
        
        
        # Anything with a small validity is by definition invalid... this also gets rid of NaN's and Inf's
        #surf_angleofincidencefactorbufuv[surf_validitybuf < 0.01] = 0.0; 
        atomicbuffer_dropthreshold(params.surfacearray[surfacecnt].angle_of_incidence_factor_buf_uv,params.surfacearray[surfacecnt].validitybuf,0.01,imgbuf_nx*imgbuf_ny)
        

        # clear out validitybuf before we reproject on top
        #surf_validitybuf[::]=0.0
        atomicbuf_zero(params.surfacearray[surfacecnt].validitybuf,imgbuf_ny*imgbuf_nx)
        pass

    
    with nogil: 
        project_image(params.surfacearray,
                      <float32_t *>imagedata_weighting.data, # Project our camera image weighting field onto the parameterization
                      0, # NOTE: If we provide a nonzero framecnt then both imgbuf and validitybuf need to be multi-frame
                      src_nx,
                      src_ny,
                      <float64_t *>cam_mtx_Ut.data,
                      <float64_t *>cam_mtx_Vtsubinv.data,
                      s1,s2,
                      V31,V32,
                      <float32_t *>imagedata_zbuffer.data,
                      <uint32_t *>imagedata_surfaceid.data,
                      <uint32_t *>imagedata_facetid.data,
                      <float32_t *>unit_imagedata_weighting.data,0)  # unit_imagedata_weighting is just 1.0 everywhere

        for surfacecnt in range(nsurfaces):
            # Sanitize the data that project_image
            # just wrote to atomic spaces, so it is OK to access them as floats
            imgbuf_nx=params.surfacearray[surfacecnt].nx
            imgbuf_ny=params.surfacearray[surfacecnt].ny
            atomicbuf_sanitize(params.surfacearray[surfacecnt].imgbuf,imgbuf_ny*imgbuf_nx)
            atomicbuf_sanitize(params.surfacearray[surfacecnt].validitybuf,imgbuf_ny*imgbuf_nx)
            atomicbuf_sanitize(params.surfacearray[surfacecnt].angle_of_incidence_factor_buf_uv,imgbuf_ny*imgbuf_nx)
            # weightingbuf hasn't been set yet (should be NULL)
        
        pass

    # now the imgbufs contain projections of the diffused zmismatch validity regions
    
    #return params,zmismatch,params.surfacelist[0][1][1],params.surfacelist[0][1][2],params.surfacelist[0][3]
    # loop through surfaces 
    for surfacecnt in range(nsurfaces):
        # ... Here we want to treat everything that wasn't projected
        # as a mismatch
        
        #(imgbuf_ny,imgbuf_nx)=params.surfacearray[surfacecnt].imgbuf.shape
        #imgbuf_ny=params.surfacearray[surfacecnt].ny
        #imgbuf_nx=params.surfacearray[surfacecnt].nx

        surf_imgbuf=params.surfacelist[surfacecnt][1][1]
        surf_validitybuf=params.surfacelist[surfacecnt][1][2]
        surf_angleofincidencefactorbufuv=params.surfacelist[surfacecnt][3]
        surf_parameterization=params.surfacelist[surfacecnt][1][0]

        #imgbuf_ptr=params.surfacearray[surfacecnt].imgbuf
        #validitybuf_ptr=params.surfacearray[surfacecnt].validitybuf
        #angle_of_incidence_factor_buf_uv_ptr = params.surfacearray[surfacecnt].angle_of_incidence_factor_buf_uv

        if surf_parameterization is None:
            surf_surface = params.surfacelist[surfacecnt][0]
            surf_parameterization=surf_surface.intrinsicparameterization
            pass

        
        
        (imgbuf_ny,imgbuf_nx)=surf_imgbuf.shape

        # apply validity normalization  # NO LONGER NECESSARY
        #surf_imgbuf = surf_imgbuf/surf_validitybuf

        # imgbuf is now the projection of the diffused mismatch validity
        # i.e. near 0 = mismatch (invalid), 1=no mismatch (valid)


        # surf_imgbuf will be our diffusion upper_bound
        
        # Anything with a small validity is by definition invalid... this also gets rid of NaN's and Inf's
        surf_imgbuf[surf_validitybuf < 0.01] = 0.0; 
        #surf_imgbuf[np.isnan[surf_imgbuf]]=0.0
        #surf_imgbuf[np.isinf[surf_imgbuf]]=0.0

        # blur it in (u,v) space
        if uv_weightingblur_distance > 0.0:
            # ***!!!! Ideally diffuse_mismatch would use the local stepsizes at each point in the parameterization!!!
            weightingbuf = diffuse_mismatch(surf_imgbuf,1.0,uv_weightingblur_distance*imgbuf_nx/surf_parameterization.meaningfulunits_per_texcoord[0],uv_weightingblur_distance*imgbuf_ny/surf_parameterization.meaningfulunits_per_texcoord[1]).astype(np.float32)
            pass
        else:
            weightingbuf=surf_imgbuf.copy()
            pass
        
        # use angleofincidencefactorbufuv weighting as the nomismatchvalue for diffusion
        # so that low weights due to steep angles/changing curvature will be blurred out. 

        normalized_angleofincidencefactor_uv = surf_angleofincidencefactorbufuv #/surf_validitybuf
        # Anything with a small validity is by definition invalid... this also gets rid of NaN's and Inf's
        #normalized_angleofincidencefactor_uv[surf_validitybuf < 0.01] = 0.0; 
        #normalized_angleofincidencefactor_uv[np.isnan(normalized_angleofincidencefactor_uv)]=0.0
        #normalized_angleofincidencefactor_uv[np.isinf(normalized_angleofincidencefactor_uv)]=0.0

        
        weightingbuf *=  normalized_angleofincidencefactor_uv 

        # store weightingbuf in params.surfacelist, which will keep the objects reference-counted
        params.surfacelist[surfacecnt][2]=weightingbuf


        # Make it legit to read weightingbuf as atomic pixels
        # store C pointer as well
        params.surfacearray[surfacecnt].weightingbuf=<atomicpixel_t *>atomicbuf_prepare(<float32_t *>weightingbuf.data,imgbuf_ny*imgbuf_nx)

        # Now that we've extracted baseline validity, clear the validitybuf for projection use. PROBABLY UNNECESSARY BECAUSE USER DOES THIS AND WE atomicbuf_prepare in imageprojection_float()
        
        #surf_validitybuf[::]=0.0
        validitybuf_ptr = params.surfacearray[surfacecnt].validitybuf
        atomicbuf_zero(validitybuf_ptr,imgbuf_ny*imgbuf_nx)
        
        # clear imgbuf for projection use. 
        #surf_imgbuf=params.surfacelist[surfacecnt][1][1]
        #surf_imgbuf[::]=0.0        
        imgbuf_ptr = params.surfacearray[surfacecnt].imgbuf
        atomicbuf_zero(imgbuf_ptr,imgbuf_ny*imgbuf_nx)
        
        
        #zmismatch_uv = np.zeros((imgbuf_ny,imgbuf_nx),dtype='bool')
        #ucoords_pixels = np.round((imagedata_uvcoords[:,:,0][zmismatch])*imgbuf_nx - 0.5) # -0.5 corresponds to lower-left in unscaled (u,v) being at -0.5 pixel
        
        #vcoords_pixels = np.round((imagedata_uvcoords[:,:,1][zmismatch])*imgbuf_ny - 0.5)
        ## Eliminate pixels that don't project anywhere
        #ucoords_pixels=ucoords_pixels[~np.isnan(ucoords_pixels)]
        #vcoords_pixels=vcoords_pixels[~np.isnan(vcoords_pixels)]
        
        ## bound coordinate values and round to integer type
        #ucoords_pixels[ucoords_pixels < 0]=0.0
        #ucoords_pixels[ucoords_pixels >= imgbuf_nx]=imgbuf_nx-1
        #ucoords_pixels_int=ucoords_pixels.astype(np.uint32)
        #
        #vcoords_pixels[vcoords_pixels < 0]=0.0
        #vcoords_pixels[vcoords_pixels >= imgbuf_ny]=imgbuf_ny-1
        #vcoords_pixels_int=vcoords_pixels.astype(np.uint32)
                
        ## Now use these u and v coordinates to assign mismatch into zmismatch_uv
        #zmismatch_uv[vcoords_pixels_int,ucoords_pixels_int]=True

        ## also define z mismatch wherever there isn't a valid projection, from the validitybuf
        #thisvaliditybuf=params.surfacelist[surfacecnt][1][2]      
        #
        #thisangleofincidencefactorbufuv=params.surfacelist[surfacecnt][3]      
        
        ##zmismatch_uv[thisvaliditybuf < 0.1] = True
        ## use angleofincidencefactorbufuv weighting as the nomismatchvalue for diffusion
        ## so that low weights due to steep angles/changing curvature will be blurred out. 

        #normalized_angleofincidencefactor_uv = thisangleofincidencefactorbufuv/thisvaliditybuf
        #normalized_angleofincidencefactor_uv[np.isnan(normalized_angleofincidencefactor_uv)]=0.0
        #normalized_angleofincidencefactor_uv[np.isinf(normalized_angleofincidencefactor_uv)]=0.0
        
        #mismatchvalue=normalized_angleofincidencefactor_uv[zmismatch_uv]*0.01
        #
        #nomismatchvalue=normalized_angleofincidencefactor_uv
        
        ## Now diffuse the mismatch in (u,v) space
        ## minimal weighting here is .05 because typically weighting at edges will also be attenuated by camera domain weighting
        

        

        pass

    #from matplotlib import pyplot as pl
    #pl.imshow(imagedata_zbuffer)
    #pl.show()

    return params #,zmismatch


def imgbufzero(imageprojection_params params):
    cdef size_t imgbuf_nx
    cdef size_t imgbuf_ny
    cdef atomicpixel_t *imgbuf_ptr
    
    for surfacecnt in range(params.nsurfaces):
        imgbuf_nx = params.surfacearray[surfacecnt].nx
        imgbuf_ny = params.surfacearray[surfacecnt].ny
        imgbuf_ptr = params.surfacearray[surfacecnt].imgbuf
        atomicbuf_zero(imgbuf_ptr,imgbuf_nx*imgbuf_ny)
        pass
    
    pass

def validitybufzero(imageprojection_params params):
    cdef size_t imgbuf_nx
    cdef size_t imgbuf_ny
    cdef atomicpixel_t *validitybuf_ptr
    
    for surfacecnt in range(params.nsurfaces):
        imgbuf_nx = params.surfacearray[surfacecnt].nx
        imgbuf_ny = params.surfacearray[surfacecnt].ny
        validitybuf_ptr = params.surfacearray[surfacecnt].validitybuf
        atomicbuf_zero(validitybuf_ptr,imgbuf_nx*imgbuf_ny)
        pass
    
    pass


def imageprojection_float(imageprojection_params params not None,imagedata):
    # !!!*** NOTE: Currently assumes single surface (!)

    cdef size_t framecnt
    cdef size_t imgbuf_nx
    cdef size_t imgbuf_ny
    cdef np.ndarray[np.float32_t,ndim=2,mode="c"] framedata
    cdef int surfacecnt
    
    assert(imagedata.ndim >= 2)

    # imagedata is presumed to be an array of images with the
    # last indices y, and x respectively, stored as a raster scan
    # unwrap all but last two indices
    if len(imagedata.shape) > 2:
        imagedata_unwrap=imagedata.reshape(np.prod(imagedata.shape[:-2]),imagedata.shape[-2],imagedata.shape[-1])
        pass
    else:
        imagedata_unwrap=imagedata.reshape(1,imagedata.shape[0],imagedata.shape[1])
        pass

    
    cdef size_t src_nx=params.src_nx
    cdef size_t src_ny=params.src_ny
    cdef np.ndarray[np.float64_t,ndim=2,mode="c"] cam_mtx_Ut = params.cam_mtx_Ut
    cdef np.ndarray[np.float64_t,ndim=2,mode="c"] cam_mtx_Vtsubinv = params.cam_mtx_Vtsubinv
    cdef float64_t s1=params.s1
    cdef float64_t s2=params.s2
    cdef float64_t V31=params.V31
    cdef float64_t V32=params.V32
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_zbuffer=params.imagedata_zbuffer
    cdef np.ndarray[np.uint32_t, ndim=2, mode="c"] imagedata_surfaceid=params.imagedata_surfaceid
    cdef np.ndarray[np.uint32_t, ndim=2, mode="c"] imagedata_facetid=params.imagedata_facetid
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_weighting=params.imagedata_weighting

    for framecnt in range(imagedata_unwrap.shape[0]):
        # !!!*** NOTE: 
        # We should really expand current image to power of 2 size
        # for mip-mapping to accommodate (perhaps rare) scenarios 
        # where a lot of image
        # pixels are being mapped down to a few parameterization pixels

        # This would come up if we zoomed the camera way in, but forgot
        # to upscale the parameterization image. 
        
        #import resource
        #sys.stderr.write("pt 6a memusage: %s\n" % (str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))

        framedata=imagedata[framecnt,:,:]



        projectstart=time.time()
        #print("projection start time=%s\n" % (str(projectstart)))

    
        with nogil: 
                
            project_image(params.surfacearray,
                          <float32_t *>framedata.data,
                          0, # NOTE: If we provide a nonzero framecnt then both imgbuf and validitybuf need to be multi-frame
                          src_nx,
                          src_ny,
                          <float64_t *>cam_mtx_Ut.data,
                          <float64_t *>cam_mtx_Vtsubinv.data,
                          s1,s2,
                          V31,V32,
                          <float32_t *>imagedata_zbuffer.data,
                          <uint32_t *>imagedata_surfaceid.data,
                          <uint32_t *>imagedata_facetid.data,
                          <float32_t *>imagedata_weighting.data,1)
            for surfacecnt in range(params.nsurfaces):
                imgbuf_nx = params.surfacearray[surfacecnt].nx
                imgbuf_ny = params.surfacearray[surfacecnt].ny
                imgbuf_ptr = params.surfacearray[surfacecnt].imgbuf
                validitybuf_ptr = params.surfacearray[surfacecnt].validitybuf
                normalize_atomic_buffer(imgbuf_ptr,validitybuf_ptr,imgbuf_nx*imgbuf_ny)

                # Sanitize the data that project_image
                # just wrote to atomic spaces, so it is OK to access them as floats
                atomicbuf_sanitize(imgbuf_ptr,imgbuf_nx*imgbuf_ny)
                atomicbuf_sanitize(validitybuf_ptr,imgbuf_nx*imgbuf_ny)
                
                pass
            pass

        
        #new_imgbuf=params.surfacelist[0][1][1]
        #new_validitybuf=params.surfacelist[0][1][2]

        ## factor out the validity
        #new_imgbuf[:,:]=new_imgbuf/new_validitybuf
        
        #print("project calc time: %f s" % (time.time()-projectstart))

        #print(np.count_nonzero(np.isfinite(params.imagedata_zbuffer)))
        
        #if not projectionmodel.undistortion_already_applied:
        #    # Need to perform undistortion
        #    undistortion_alpha=1.0
        #
        #    processor=intrinsic_calibration.undistortion_processor(projectionmodel.calibration,projectionmodel.calibration.sizex,projectionmodel.calibration.sizey,int(projectionmodel.calibration.sizex*undistortion_grow_factor),int(projectionmodel.calibration.sizey*undistortion_grow_factor),cv2.INTER_LINEAR,1.0)
        #    
        #    undistorted_image=processor.undistort_numpy_image(imagedata)
        #    pass

        # cvWarpPerspective operates only on rectangles for (somewhat) obvious reasons.
        # So what we have to do is find each polygon, define a rectangular bounding box on the source image,
        # warp the bounding box onto the destination,
        # then project the individual vertices, and copy what is enclosed
        # by the lines between the vertices. 
        # ...
        # Or we could manually transform coordinates of each pixel on the
        # destination to extract the corresponding pixel value from the
        # source. 
        
        pass
        
    #cdef cvMat cv_camera_mtx 
    #cdef cvMat cv_dist_coeffs # = cv::Mat(1,distcoeffs.shape[0],cv::CV_64F,distcoeffs.data)
    #cv_camera_mtx= cvMat(4,4,CV_64F, new_camera_mtx.data)
	
    pass
