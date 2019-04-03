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

    void evaluate_zbuffer_orthographic(projectionsurface *surfacearray,
				       size_t src_nx,
				       size_t src_ny,
				       size_t numsurfaces,
				       float64_t *orthographic_proj_matrix,
				       float32_t *imagedata_zbuffer,
				       uint32_t *imagedata_surfaceid,
				       uint32_t *imagedata_facetid,
				       float32_t *imagedata_angleofincidence,
				       float32_t *imagedata_angleofincidence_weighting,
				       float32_t *imagedata_uvcoords) nogil

    void project_orthographic(projectionsurface *surfacearray,
			      float32_t *framedata,
			      size_t framecnt,
			      size_t src_nx,
			      size_t src_ny,
			      float64_t *orthographic_proj_matrix,
			      float32_t *imagedata_zbuffer,
			      uint32_t *imagedata_surfaceid,
			      uint32_t *imagedata_facetid,
			      float32_t *imagedata_weighting) nogil
    
    pass

pass


cdef class imageorthographic_params:
    cdef public int nsurfaces
    cdef projectionsurface *surfacearray

    cdef public object surfacelist  # Python list containing (surface,parameterizationdictentry), for purposes of reference counting
    
    
    # parameters of the projection
    cdef public size_t src_nx
    cdef public size_t src_ny
    


    cdef public object tosurfaceframes  # Store list of tosurfaceframe matrices here too for purposes of reference counting
    cdef public object vecfromsurfaceframes  # Store list of vecfromsurfaceframe matrices here too for purposes of reference counting


    #cdef public float64_t s1
    #cdef public float64_t s2
    #cdef public float64_t V31
    #cdef public float64_t V32 
    cdef public np.ndarray orthographic_proj_matrix

    cdef public np.ndarray imagedata_zbuffer
    cdef public np.ndarray imagedata_surfaceid
    cdef public np.ndarray imagedata_facetid 
    #cdef public np.ndarray imagedata_angleofincidence 
    #cdef public np.ndarray imagedata_angleofincidence_weighting
    #cdef public np.ndarray imagedata_weighting # Weighting coefficients in camera frame
    #cdef public np.ndarray imagedata_horiz_ray_intersect_shift_deriv_cam_z 
    #cdef public np.ndarray imagedata_vert_ray_intersect_shift_deriv_cam_z 
    #cdef public np.ndarray imagedata_uvcoords
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



@cython.boundscheck(False)
def imageorthographic_prepare_float(np.ndarray[np.float64_t,ndim=2,mode="c"] orthographic_proj_matrix,cameraframe,partlist,parameterizationdict,imageshape):
    
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
    
    params=imageorthographic_params(nsurfaces)
    params.tosurfaceframes=[]
    params.vecfromsurfaceframes=[]
    
    cdef np.ndarray[np.float32_t,mode="c",ndim=2] imgbuf
    cdef atomicpixel_t *imgbuf_ptr
    #cdef np.ndarray[np.float32_t,mode="c",ndim=2] weightingbuf
    cdef np.ndarray[np.float32_t,mode="c",ndim=2] validitybuf
    cdef atomicpixel_t *validitybuf_ptr
    #cdef np.ndarray[np.float32_t,mode="c",ndim=2] angle_of_incidence_factor_buf_uv
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

		#weightingbuf=np.ones(imgbuf.shape,dtype='f')  # weightingbuf just all-ones
		
                #angle_of_incidence_factor_buf_uv=np.zeros((validitybuf.shape[0],validitybuf.shape[1]),dtype='f')                  
                params.surfacelist.append([surface,parameterizationdict[id(surface)] ])  # Here we keep reference to surface and parameterizationdict and weightingbuf
                
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

                
                params.surfacearray[surfacecnt].imgbuf=imgbuf_ptr
                
                params.surfacearray[surfacecnt].weightingbuf=<atomicpixel_t *>NULL   # weightingbuf.data
                params.surfacearray[surfacecnt].validitybuf=validitybuf_ptr
                params.surfacearray[surfacecnt].angle_of_incidence_factor_buf_uv=NULL # <float32_t *>angle_of_incidence_factor_buf_uv.data
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

    
    assert(len(imageshape)==2)

    # Image is presumed to be stored raster style C indexed, so coordinates are (y,x)
    params.src_ny=imageshape[0]
    params.src_nx=imageshape[1]
    params.orthographic_proj_matrix=orthographic_proj_matrix
    
        
    # define z-buffer, currently at infinity
    params.imagedata_zbuffer = (np.ones((imageshape[0],imageshape[1]),dtype='f')*np.inf).astype(np.float32)
    # define array of which surface id
    params.imagedata_surfaceid = np.zeros((imageshape[0],imageshape[1]),dtype='u4')
    # define array of which facet id
    params.imagedata_facetid = np.zeros((imageshape[0],imageshape[1]),dtype='u4')
    # define array of angles of incidence
    #params.imagedata_angleofincidence = np.zeros((imageshape[0],imageshape[1]),dtype='f')
    # define array of weightings of camera image
    #params.imagedata_angleofincidence_weighting = np.ones((imageshape[0],imageshape[1]),dtype='f')
    

    # Define array of z derivatives with respect to adjacent rays 
    #params.imagedata_horiz_ray_intersect_shift_deriv_cam_z = np.zeros((imageshape[0],imageshape[1]),dtype='f')
    #params.imagedata_vert_ray_intersect_shift_deriv_cam_z = np.zeros((imageshape[0],imageshape[1]),dtype='f')

    ## Place to store the (u,v) coordinates of the pixels
    #params.imagedata_uvcoords=np.zeros((imageshape[0],imageshape[1],2),dtype='f')
    
    
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
            
    ## See doc/Projection.pdf
    #(U,S,Vt) = svd(new_camera_mtx[:2,:3])
    zbufstart=time.time()

    cdef projectionsurface *surfacearray=params.surfacearray

    cdef size_t src_nx=params.src_nx
    cdef size_t src_ny=params.src_ny

    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_zbuffer=params.imagedata_zbuffer
    cdef np.ndarray[np.uint32_t, ndim=2, mode="c"] imagedata_surfaceid=params.imagedata_surfaceid
    cdef np.ndarray[np.uint32_t, ndim=2, mode="c"] imagedata_facetid=params.imagedata_facetid
    #cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_angleofincidence=params.imagedata_angleofincidence
    #cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_angleofincidence_weighting=params.imagedata_angleofincidence_weighting
    #cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_horiz_ray_intersect_shift_deriv_cam_z=params.imagedata_horiz_ray_intersect_shift_deriv_cam_z
    #cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_vert_ray_intersect_shift_deriv_cam_z=params.imagedata_vert_ray_intersect_shift_deriv_cam_z
    #cdef np.ndarray[np.float32_t, ndim=3, mode="c"] imagedata_uvcoords=params.imagedata_uvcoords

    
    with nogil:
        evaluate_zbuffer_orthographic(surfacearray,
                         src_nx,
                         src_ny,
                         nsurfaces,
                         <float64_t *>orthographic_proj_matrix.data,
                         <float32_t *>imagedata_zbuffer.data, 
                         <uint32_t *>imagedata_surfaceid.data,
                         <uint32_t *>imagedata_facetid.data,
                         NULL, #<float32_t *>imagedata_angleofincidence.data,
                         NULL, #<float32_t *>imagedata_angleofincidence_weighting.data, # initial weighting is just 1.0... initial validitybuf will be from projection operation alone. This will get filled in with the desired angle-of-incidence weighting, but that waiting will NOT affect the validitybuf output
                         #<float32_t *>imagedata_horiz_ray_intersect_shift_deriv_cam_z.data,
                         #<float32_t *>imagedata_vert_ray_intersect_shift_deriv_cam_z.data,
                         #<float32_t *>imagedata_uvcoords.data)
			 NULL)
        pass

    # NOTE: evaluate_zbuffer drops baseline projection validity into validitybuf
    
    print("zbuffer calc time: %f s" % (time.time()-zbufstart))



    return params #,zmismatch


def imageorthographic_float(imageorthographic_params params not None,imagedata,np.ndarray[np.float32_t,ndim=2,mode="c"] validitydata):

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
    #cdef float64_t s1=params.s1
    #cdef float64_t s2=params.s2
    #cdef float64_t V31=params.V31
    #cdef float64_t V32=params.V32
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] imagedata_zbuffer=params.imagedata_zbuffer
    cdef np.ndarray[np.uint32_t, ndim=2, mode="c"] imagedata_surfaceid=params.imagedata_surfaceid
    cdef np.ndarray[np.uint32_t, ndim=2, mode="c"] imagedata_facetid=params.imagedata_facetid
    cdef np.ndarray[np.float64_t, ndim=2, mode="c"] orthographic_proj_matrix=params.orthographic_proj_matrix

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

        framedata=np.ascontiguousarray(imagedata_unwrap[framecnt,:,:])



        projectstart=time.time()
        #print("projection start time=%s\n" % (str(projectstart)))

    
        with nogil: 
            project_orthographic(params.surfacearray,
                                 <float32_t *>framedata.data,
                                 0, # NOTE: If we provide a nonzero framecnt then both imgbuf and validitybuf need to be multi-frame
                                 src_nx,
                                 src_ny,
                                 <float64_t *>orthographic_proj_matrix.data,
                                 <float32_t *>imagedata_zbuffer.data,
                                 <uint32_t *>imagedata_surfaceid.data,
                                 <uint32_t *>imagedata_facetid.data,
                                 <float32_t *>validitydata.data)
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
