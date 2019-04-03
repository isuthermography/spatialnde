#include <omp.h>

#ifdef _MSC_VER
#define IPOPS_INLINE  __inline
#include <float.h>
#define IPOPS_ISINF(x) (!_finite(x))
#define IPOPS_ISNAN(x) (_isnan(x))

// MSVC STILL doesn't support C11
typedef uint32_t atomicpixel_t;
#else
//#include <pthread.h>
#include <fenv.h>
#include <stdatomic.h>
#define IPOPS_INLINE  inline
#define IPOPS_ISINF isinf
#define IPOPS_ISNAN isnan

typedef _Atomic uint32_t atomicpixel_t; /* underlying data is actually a float,
					   but we use this type because it permits
					   us to use atomic operations. YOU WILL
					   GET BAD DATA IF YOU ACCESS THIS DIRECTLY!!!
					*/

#endif


static IPOPS_INLINE float32_t *atomicbuf_sanitize(atomicpixel_t *buf,size_t n)
{
  /* sanitize atomic data in buf so it is OK to read it as floats */
  int64_t cnt;
  union cvt {
    float32_t *fptr;
    atomicpixel_t *iptr;
  } cvtptr;

  assert(sizeof(atomicpixel_t)==4);
  assert(sizeof(float32_t)==4);
  assert(sizeof(uint32_t)==4);

  cvtptr.iptr = buf;
  
#pragma omp parallel default(shared) private(cnt) //,n)
#pragma omp for
  for (cnt=0;cnt < n;cnt++) {
    union {
      char cbuf[4];
      float32_t fbuf;
      uint32_t ibuf;
    } cvtbuf;
    memcpy(&cvtbuf.ibuf,&cvtptr.iptr[cnt],4);
    memcpy(&cvtptr.fptr[cnt],&cvtbuf.fbuf,4);    
  }
  return cvtptr.fptr;
}


static IPOPS_INLINE atomicpixel_t *atomicbuf_prepare(float32_t *buf,size_t n)
{
  /* prepare float data in buf so it is OK to read it as atomics */
  int64_t cnt;
  assert(sizeof(atomicpixel_t)==4);
  assert(sizeof(float32_t)==4);
  assert(sizeof(uint32_t)==4);
  
  union cvt {
    float32_t *fptr;
    atomicpixel_t *iptr;
  } cvtptr;
  
  cvtptr.fptr = buf;

#pragma omp parallel default(shared) private(cnt)//,n)
#pragma omp for
  for (cnt=0;cnt < n;cnt++) {
    union {
      char cbuf[4];
      float32_t fbuf;
      atomicpixel_t ibuf;
    } cvtbuf;
    memcpy(&cvtbuf.fbuf,&cvtptr.fptr[cnt],4);
    memcpy(&cvtptr.iptr[cnt],&cvtbuf.ibuf,4);    
  }
  return cvtptr.iptr;
}


static IPOPS_INLINE float32_t atomicpixel_load(volatile atomicpixel_t *var)
{
  // Use of a union like this is legal in C11, even under the strictest
  // aliasing rules
  union {
    uint32_t intval;
    float32_t floatval;
  } pixval;
  pixval.intval=atomic_load_explicit(var,memory_order_acquire);//,memory_order_consume);

  return pixval.floatval;
}

static IPOPS_INLINE float32_t atomicpixel_nonatomicload(atomicpixel_t *var)
{
  // Use of a union like this is legal in C11, even under the strictest
  // aliasing rules
  union {
    uint32_t intval;
    float32_t floatval;
  } pixval;
  pixval.intval=*var;

  return pixval.floatval;
}

static IPOPS_INLINE void atomicpixel_store(volatile atomicpixel_t *var,float32_t val)
{
  // Use of a union like this is legal in C11, even under the strictest
  // aliasing rules
  union {
    uint32_t intval;
    float32_t floatval;
  } pixval;
  pixval.floatval=val;
  atomic_store_explicit(var,pixval.intval,memory_order_release);
  
}
static IPOPS_INLINE void atomicpixel_nonatomicstore(atomicpixel_t *var,float32_t val)
{
  // Use of a union like this is legal in C11, even under the strictest
  // aliasing rules
  union {
    uint32_t intval;
    float32_t floatval;
  } pixval;
  pixval.floatval=val;
  *var=pixval.intval;
  
}

//pthread_mutex_t accumulatemutex=PTHREAD_MUTEX_INITIALIZER;

static IPOPS_INLINE void atomicpixel_accumulate(volatile atomicpixel_t *var,float toadd)
{
  // Use of a union like this is legal in C11, even under the strictest
  // aliasing rules
  union {
    uint32_t intval;
    float32_t floatval;
    char workbuf[4];
  } oldvalue,newvalue; // ,workvalue;

  //  pthread_mutex_lock(&accumulatemutex);

  
  //oldvalue.floatval=atomicpixel_load(var);
  oldvalue.intval=atomic_load_explicit(var,memory_order_acquire);//,memory_order_consume);
  
  do {
    //memcpy(workvalue.workbuf,&oldvalue.intval,4);
    newvalue.floatval=oldvalue.floatval+toadd;
    //workvalue.floatval+=toadd;
    //memcpy(&newvalue.intval,&workvalue.workbuf,4);
  } while (!atomic_compare_exchange_strong_explicit(var,&oldvalue.intval,newvalue.intval,memory_order_seq_cst,memory_order_acquire)); //,memory_order_consume));


  //  pthread_mutex_unlock(&accumulatemutex);

}

static IPOPS_INLINE void atomicbuf_zero(atomicpixel_t *buf,size_t n)
{
  /* sanitize atomic data in buf so it is OK to read it as floats */
  int64_t cnt;
  assert(sizeof(atomicpixel_t)==4);
  assert(sizeof(float32_t)==4);
  
#pragma omp parallel default(shared) private(cnt) //,n)
#pragma omp for
  for (cnt=0;cnt < n;cnt++) {
    atomicpixel_nonatomicstore(&buf[cnt],0.0);
    
  }
}

static IPOPS_INLINE void normalize_atomic_buffer(atomicpixel_t *tonormalize,atomicpixel_t *normbuf,size_t n)
{
  /* Divide each element in tonormalize by the corresponding element in normbuf */
  int64_t cnt;
  assert(sizeof(atomicpixel_t)==4);
  assert(sizeof(float32_t)==4);
  
#pragma omp parallel default(shared) private(cnt)//,n)
#pragma omp for
  for (cnt=0;cnt < n;cnt++) {
    float32_t value,normalization;
    value=atomicpixel_nonatomicload(&tonormalize[cnt]);
    normalization=atomicpixel_nonatomicload(&normbuf[cnt]);

    value/=normalization;
    atomicpixel_nonatomicstore(&tonormalize[cnt],value);
    
  }
}


static IPOPS_INLINE void atomicbuffer_dropthreshold(atomicpixel_t *buf,atomicpixel_t *thresholdtestbuf,float32_t threshold,size_t n)
{
  /* Set each element in buf corresponding to an element of thresholdtestbuf that is less than threshold to 0 */
  int64_t cnt;
  assert(sizeof(atomicpixel_t)==4);
  assert(sizeof(float32_t)==4);
  
#pragma omp parallel default(shared) private(cnt)//,n)
#pragma omp for
  for (cnt=0;cnt < n;cnt++) {
    float32_t thresholdtest;
    thresholdtest=atomicpixel_nonatomicload(&thresholdtestbuf[cnt]);
    if (thresholdtest < threshold) {
      //bufval=atomicpixel_nonatomicload(&buf[cnt]);
      atomicpixel_nonatomicstore(&buf[cnt],0.0);
    }
    
  }
}



// C++: This should be implemented as an abstract interface
// that can either use a meshed representation (like this one)
// or a NURBS representation, depending on the surface object.
// ... Might also call OpenCL to do the heavy lifting

//#include <math.h>


///* Right now, OpenCV is used solely for matrix inversion... */
//#include <opencv2/core/version.hpp>
//#if (CV_MAJOR_VERSION >= 3)
//#include "opencv2/core/fast_math.hpp" // Temporary workaround for https://github.com/opencv/opencv/issues/6585
//#endif


//#include <opencv/cv.h>

/* splat characteristics for projecttoimgbuf */
#define MIN_RADIUS_UV_PIXELS 1.5
#define MIN_RADIUS_SRC_PIXELS 1.5
//#define BANDWIDTH_FRACTION 0.7
//#define BANDWIDTH_FRACTION 0.25
#define BANDWIDTH_FRACTION 0.4

struct projectionsurface {
  atomicpixel_t *imgbuf;
  atomicpixel_t *weightingbuf;
  atomicpixel_t *validitybuf;
  atomicpixel_t *angle_of_incidence_factor_buf_uv;
  size_t nx;
  size_t ny;
  
  float64_t u_meaningfulunits_per_texcoord;
  float64_t v_meaningfulunits_per_texcoord;
  
  float64_t *tocameraframe; // 4x4
  float64_t *tosurfaceframe; // 4x4
  float64_t *vecfromsurfaceframe; // 3x3
  float64_t *boxcoords; // nboxes x 6
  size_t nboxes;
  int32_t *boxes; // nboxes x 9
  int32_t *boxpolys;
  uint32_t *vertexidx_indices; // npolys 
  int32_t *vertexidx; // indexed by vertexidx_indices, -1 terminator after each polygon
  uint32_t *numvertices; // npolys 
  size_t max_vertices;
  //size_t nvertices;
  size_t npolys;
  float64_t *vertices; // nvertices x 3

  float64_t *facetnormals;
  float64_t *refpoints;
  float64_t *inplanemats;
  //float64_t *outofplanevec;
  float64_t *inplane2texcoords; // 2x3 projective matrix
  float64_t *maxradius; // n_polys
};



static IPOPS_INLINE void iterate_diffuse_mismatch(float32_t *upper_bound,
						  float32_t *value,
						  float32_t *extra_buf,
						  double delta_tbar,
						  uint64_t ny,
						  uint64_t nx,
						  double dybarsq,
						  double dxbarsq,
						  uint64_t ntbar)
{
  uint64_t tcnt;
  int64_t ycnt; // must be signed for msvc compatibility
  float32_t *ping_buf;
  float32_t *pong_buf; 
  float32_t *temp_ptr;
  
  if (ntbar % 2) ntbar++; // ntbar must be even so our result ends up
  // in value, not in extra_buf

  // interchange ping_buf and pong_buf every iteration
  ping_buf = value;
  pong_buf = extra_buf; 
    
  for (tcnt=0;tcnt < ntbar;tcnt++,temp_ptr=ping_buf,ping_buf=pong_buf,pong_buf=temp_ptr) {

    //    exceedance_zone = value > upper_bound
    //    value[ exceedance_zone ]=upper_bound[ exceedance_zone ]  # Assign specified value condition where the value exceeds the upper bound

#pragma omp parallel default(shared) private(ycnt)//,ny)
    {
#pragma omp for
      for (ycnt=0;ycnt < ny;ycnt++) {
	uint64_t xcnt;
	
	for(xcnt=0;xcnt < nx;xcnt++) {
	  if (ping_buf[nx*ycnt + xcnt] > upper_bound[nx*ycnt + xcnt]) {
	    ping_buf[nx*ycnt + xcnt] = upper_bound[nx*ycnt + xcnt];
	  }
	}
      }
    }


#pragma omp parallel default(shared) private(ycnt)//,ny)
    {
    
      // Iterate time 
      // T_t+ = T_t- + dt*( (T_x+ -2T + T_x-)/dx^2 + (T_y+ -2T + T_y-)/dy^2 )
#pragma omp for
      for (ycnt=1;ycnt < ny-1;ycnt++) {
	uint64_t xcnt;
	for(xcnt=1;xcnt < nx-1;xcnt++) {
	  // value[1:-1,1:-1]=value[1:-1,1:-1] + delta_tbar * ( ((value[1:-1,2:] - 2*value[1:-1,1:-1] + value[1:-1,:-2])/dxbar**2) + ((value[2:,1:-1] - 2*value[1:-1,1:-1] + value[:-2,1:-1])/dybar**2) )
	  pong_buf[nx*ycnt + xcnt] = ping_buf[nx*ycnt + xcnt] + delta_tbar * ( ((ping_buf[nx*ycnt + xcnt+1] - 2*ping_buf[nx*ycnt + xcnt] + ping_buf[nx*ycnt + xcnt-1])/dxbarsq) + ((ping_buf[nx*(ycnt+1) + xcnt] - 2*ping_buf[nx*ycnt + xcnt] + ping_buf[nx*(ycnt-1) + xcnt])/dybarsq));
	}
      }
    
    }
  }
}


static IPOPS_INLINE int ray_box_intersection(float64_t *boxcoords, float64_t *starting_point, float64_t *ray_direc)
{
  /* Slab method: Look at distance t along the ray where we first 
     intersect each slab. That's tnear. Look at the distance t 
     along the ray where we last intersect each slab. That's tfar. 
     Look at the largest tnear value. If that's greater than the 
     smallest tfar vluae, the ray misses the box. 
     Also special cases if the ray is parallel to an axis */
  /* See http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm */
  //double tnear=infnan(-ERANGE); /* -infinity */
  //double tfar=infnan(ERANGE); /* infinity */
  /*** NEED MSVC FIX ***/

#ifdef _MSC_VER
  uint64_t infinity_val=0x000000007FF00000l;
  uint64_t neginfinity_val=0x00000000FFF00000l;
  double tnear,tfar;

  tnear = *(float64_t *)&neginfinity_val; /* -infinity */
  tfar = *(float64_t *)&infinity_val; /* infinity */
#else 
  double tnear=-1.0/0.0; /* -infinity */
  double tfar=1.0/0.0; /* infinity */
#endif
  {
    double curnear,curfar,temp;
    
    int index; 
    int max_index;

    for (index=0;index < 3;index++) {
      /* index indexes starting_point and ray_direc, 
	 and finds the min value from box_coords */
      max_index=index+3; /* finds the max value from box_coords */
      /* Start with X=xmin and X=xmax planes */
      if (ray_direc[index] == 0.0) {
	/* Ray is normal to given axis, parallel to these planes */
	
	/* if origin not between planes, does not intersect */
	if (starting_point[index] < boxcoords[index] || starting_point[index] > boxcoords[max_index]) {
	  return FALSE;
	}
	
      } else {
	curnear=(boxcoords[index]-starting_point[index])/ray_direc[index]; /* distance to reach x value of boxcoords[0] */
	curfar=(boxcoords[max_index]-starting_point[index])/ray_direc[index];
	
	/* did we get curnear and curfar in the correct order? */
	if (curfar < curnear) {
	  /* swap */
	  temp=curfar;
	  curfar=curnear;
	  curnear=temp;
	}
	if (curnear > tnear) {
	  tnear=curnear; /* largest tnear value */
	}
	if (curfar < tfar) {
	  tfar=curfar; /* smallest tfar value */
	}
	if (tnear > tfar) {
	  /* missed box */
	  return FALSE;
	}
	if (tfar < 0.0) {
	  /* box is behind */
	  return FALSE;
	}
      }
    }
    return TRUE;
  }
}


static IPOPS_INLINE float64_t ray_to_plane_distance(float64_t *starting_point,float64_t *ray_direc, float64_t *planepoint, float64_t *planenormal)
{
  // starting point and planepoint should be 3-vectors (or normalized 4-vectors)
  // (planepoint - raypoint) dot planenormal / (raydirec dot planenormal)
  // returns -1.0 if ray is parallel to plane
  
  float64_t pointdiff[3];
  float64_t denominator;
  

  subvecvec3(planepoint,starting_point,pointdiff);
  denominator=dotvecvec3(ray_direc,planenormal);
  if (denominator==0.0) return -1.0;
  return dotvecvec3(pointdiff,planenormal)/denominator;
  
}


static IPOPS_INLINE void ray_to_plane_raydirec_shift(float64_t *starting_point,float64_t *ray_direc, float64_t *planepoint, float64_t *planenormal,float64_t *ray_direc_deriv,float64_t *deriv)
  /* Calculate change in ray_to_plane intersection point with respect to ray_direc (change in normalized vector) */ 
{
  // starting point and planepoint should be 3-vectors (or normalized 4-vectors)
  // (planepoint - raypoint) dot planenormal / (raydirec dot planenormal)
  // returns -1.0 if ray is parallel to plane
  
  float64_t pointdiff[3],firstterm[3],secterm[3],diff[3];
  float64_t ray_direc_dot_planenormal;

  
  // Regular evaluation of intersection point
  subvecvec3(planepoint,starting_point,pointdiff);
  ray_direc_dot_planenormal=dotvecvec3(ray_direc,planenormal);
  //intersectpoint = startingpoint + dotvecvec3(pointdiff,planenormal)/ray_direc_dot_planenormal * ray_direc;

  // We need to calculate dintersectpoint/dray_direc:
  // pointdiff is constant, planenormal is constant
  //intersectpoint = startingpoint + dotvecvec3(pointdiff,planenormal)* (ray_direc/(ray_direc dot planenormal)

  // deriv intersectpoint = dotvecvec3(pointdiff,planenormal)* ( deriv ray_direc * (ray_direc dot planenormal) - ((deriv ray_direc) dot planenormal) * ray_direc)/(ray_direc dot planenormal)^2 

  // deriv intersectpoint = (dotvecvec3(pointdiff,planenormal)/(ray_direc_dot_planenormal)^2) * ( deriv ray_direc * (ray_direc dot planenormal) - ((deriv ray_direc) dot planenormal) * ray_direc)

  scalevec3(ray_direc_dot_planenormal,ray_direc_deriv,firstterm);
  scalevec3(dotvecvec3(ray_direc_deriv,planenormal),ray_direc,secterm);
  subvecvec3(firstterm,secterm,diff);
  
  scalevec3(dotvecvec3(pointdiff,planenormal)/(ray_direc_dot_planenormal*ray_direc_dot_planenormal),diff,deriv); // store result in deriv
  
  
  
  // you should multiply return value by discrete shift (factor for ray_direc_deriv)
  
  // (but this is 1.0 for a one pixel shift if factors are correct on the way in)

}





static IPOPS_INLINE void ray_to_plane_raypos_shift(float64_t *starting_point,float64_t *ray_direc, float64_t *planepoint, float64_t *planenormal,float64_t *starting_point_direc_deriv,float64_t *deriv)
  /* Calculate change in ray_to_plane intersection point with respect to ray_direc (change in normalized vector) */ 
{
  // starting point and planepoint should be 3-vectors (or normalized 4-vectors)
  // (planepoint - raypoint) dot planenormal / (raydirec dot planenormal)
  // returns -1.0 if ray is parallel to plane
  
  float64_t pointdiff[3];
  float64_t ray_direc_dot_planenormal;

  
  // Regular evaluation of intersection point
  subvecvec3(planepoint,starting_point,pointdiff);
  ray_direc_dot_planenormal=dotvecvec3(ray_direc,planenormal);
  //intersectpoint = startingpoint + dotvecvec3(planepoint-starting_point,planenormal)/ray_direc_dot_planenormal * ray_direc;

  // We need to calculate dintersectpoint/dstartingpoint:
  // ray_direc is constant, planenormal is constant
  
  //intersectpoint = startingpoint + dotvecvec3(planepoint-starting_point,planenormal)*(ray_direc/(ray_direc dot planenormal))
  //intersectpoint = startingpoint + dotvecvec3(planepoint,planenormal)*(ray_direc/(ray_direc dot planenormal)) - dotvecvec3(starting_point,planenormal)*(ray_direc/(ray_direc dot planenormal))

  // deriv intersectpoint = deriv startingpoint  - dotvecvec3(deriv starting_point,planenormal)*(ray_direc/(ray_direc dot planenormal))

  addvecscaledvec3(starting_point_direc_deriv,-dotvecvec3(starting_point_direc_deriv,planenormal)/ray_direc_dot_planenormal,ray_direc,deriv);
  
  // you should multiply return value by discrete shift (factor for ray_direc_deriv)
  
  // (but this is 1.0 for a one pixel shift if factors are correct on the way in)

}




static IPOPS_INLINE int ray_intersects_polygon(float64_t *vertexarray,
					 int32_t *vertexids,
					 uint32_t numvertices,
					 float64_t *normal_vector, // unit vector
					 float64_t *refpoint,
					 float64_t *inplanemat, // 2x3 matrix: two orthogonal unit vectors normal to normal_vector
					 float64_t *starting_point,
					 float64_t *ray_direc,
					 float64_t maxradius,
					 float64_t ray_to_plane_dist,
					 int include_edges,int trace)
{
  float64_t intersectpoint[3];
  float64_t u1,v1,u2,v2;
  float64_t windingnum=0.0;
  size_t vertnum,nextvertnum;
  float64_t vec1[3];
  float64_t vec2[3];
  float64_t intersectdist[3];

  float64_t magn1,magn2,det,cosparam;
  
  addvecscaledvec3(starting_point,ray_to_plane_dist,ray_direc,intersectpoint);
  
  subvecvec3(intersectpoint,refpoint,intersectdist);
  
  /* checking whether intersection point within maximum outer
     radius of polygon */

  if (trace) {
    fprintf(stderr,"intersectpoint=%g,%g,%g; refpoint=%g,%g,%g\n",intersectpoint[0],intersectpoint[1],intersectpoint[2],refpoint[0],refpoint[1],refpoint[2]);
  
    fprintf(stderr,"intersectdist=%g; maxradius=%g\n",normvec3(intersectdist),maxradius);

    {
      float64_t intersect_outofplane=dotvecvec3(intersectdist,normal_vector);
      fprintf(stderr,"out of plane component = %g\n",intersect_outofplane);
    //if (fabs(intersect_outofplane) > (0.3*normvec3(intersectdist))) {
    //  *((char *)0) = 0;  // force segfault
    // }
    }

  }
  
  if (normvec3(intersectdist) <= maxradius) {

    if (trace) fprintf(stderr,"within outer radius\n");

  
    /* Within outer radius */
	  
	  
    // for the ray, determine if intersectpoint is inside the polygon

    // Apply winding number algorithm.
    // This algorithm is selected -- in its most simple form --
    // because it is so  simple and robust in the case of the
    // intersect point being on or near the edge. It may well
    // be much slower than optimal.
    //
    // Should probably implement a faster algorithm then drop
    // down to this for the special cases.

    // See Hormann and Agathos, The point in polygon problem
    // for arbitrary polygons, Computational Geometry 20(3) 131-144 (2001)
    // http://dx.doi.org/10.1016/S0925-7721(01)00012-8
    // https://pdfs.semanticscholar.org/e90b/d8865ddb7c7af2b159d413115050d8e5d297.pdf
    
    // Winding number is sum over segments of
    // acos((point_to_vertex1 dot point_to_vertex2)/(magn(point_to_vertex1)*magn(point_to_vertex_2))) * sign(det([ point_to_vertex1  point_to_vertex2 ]))
    // where sign(det) is really: What is the sign of the z
    // component of (point_to_vertex1 cross point_to_vertex2)
    //
    // Special cases: magn(point_to_vertex1)==0 or
    // magn_point_to_vertex2   -> point is on edge
    // det([ point_to_vertex1  point_to_vertex2 ]) = 0 -> point may be on edge
    //

    for (vertnum=0;vertnum < numvertices && vertexids[vertnum] >= 0;vertnum++) {
      nextvertnum=vertnum+1;
      if (vertexids[nextvertnum] < 0 || nextvertnum >= numvertices) nextvertnum=0;

      // calculate (thisvertex - intersectionpoint) -> vec1
      subvecvec3(&vertexarray[vertexids[vertnum]*3],intersectpoint,vec1);
      
      // calculate (nextvertex - intersectionpoint) -> vec2
      subvecvec3(&vertexarray[vertexids[nextvertnum]*3],intersectpoint,vec2);

      // Project points into 2-space:
      u1=dotvecvec3(vec1,inplanemat+0);
      v1=dotvecvec3(vec1,inplanemat+3);

      u2=dotvecvec3(vec2,inplanemat+0);
      v2=dotvecvec3(vec2,inplanemat+3);

      magn1 = sqrt(u1*u1 + v1*v1);
      magn2 = sqrt(u2*u2 + v2*v2);

      if (magn1==0.0 || magn2==0.0) {
	return include_edges; 
      }

      /* Normalize vectors */
      u1 /= magn1;
      v1 /= magn1;

      u2 /= magn2;
      v2 /= magn2;
      
      det=(u1*v2-u2*v1); // matrix determinant


      cosparam=(u1*u2 + v1*v2); //  /(magn1*magn2);

      if (cosparam < -1.0) {
	cosparam=-1.0; // Shouldn't be possible...just in case of weird roundoff
      }
      if (cosparam > 1.0) {
	cosparam=1.0; // Shouldn't be possible...just in case of weird roundoff
      }

      if (det > 0) {	
	windingnum += acos(cosparam);
      } else if (det < 0) {
	windingnum -= acos(cosparam);
      } else {
	/* det==0  */

	/* Vectors parallel or anti-parallel */
	
	if (cosparam > 0.9) {
	  // Vectors parallel. We are OUTSIDE. Do Nothing */
	}
	else if (cosparam < -0.9) {
	  // Vectors anti-parallel. We are ON EDGE */	  
	  return include_edges;
	} else {
	  assert(0); /* Should only be able to get cosparam = +/- 1.0 if det==0.0 */
	}
	
      }
      
      
      
    }
    
    windingnum=fabs(windingnum)*(1.0/(2.0*M_PI)); // divide out radians to number of winds; don't care about clockwise vs. ccw
    
    if (windingnum > 0.999 && windingnum < 1.001) {
      //fprintf(stderr,"Winding number ~1.0\n");
      return TRUE;  /* almost exactly one loop */
    }
    if (windingnum < 0.001) {
      //fprintf(stderr,"Winding number ~0.0\n");
      return FALSE;
      
    }
    fprintf(stderr,"imageprojection_ops.h/ray_intersects_polygon Got weird winding number of %e; assuming inaccurate calculation on polygon edge\n",windingnum);
    // Could also be self intersecting polygon 
    return include_edges; 
    
  }
  return FALSE;
}

static IPOPS_INLINE int box_contains_polygon(struct projectionsurface *surf,size_t boxnum,size_t polynum,int trace) /* used for debugging only */
{
  int retval=0;
  size_t cnt,subbox;
  
  if (surf->boxes[boxnum*9 + 8]!=-1) {
    /* Have index into boxpolys array: This box contains polygons */
    for (cnt=surf->boxes[boxnum*9 + 8];surf->boxpolys[cnt] >= 0;cnt++) {
      /* Got a polygon */
      if (polynum==surf->boxpolys[cnt]) {
	if (trace) {
	  fprintf(stderr,"Box %d directly contains polygon %d\n",(int)boxnum,(int)polynum);
	}
	return TRUE;
      }
    }
  }
  /* This box may be sub-divided into 8 */
  for (subbox=0;subbox < 8; subbox++) {
    if (surf->boxes[boxnum*9+subbox] >= 0) {
      retval = retval || box_contains_polygon(surf,surf->boxes[boxnum*9 + subbox],polynum,FALSE);
    }
  }
  if (retval && trace) {
    fprintf(stderr,"Box %d indirectly contains polygon %d\n",(int)boxnum,(int)polynum);
  } else if (trace) {
    fprintf(stderr,"Box %d does not contain polygon %d\n",(int)boxnum,(int)polynum);   
  }
  return retval;
  
}


static IPOPS_INLINE int find_ray_intersections(struct projectionsurface *surf,
					  size_t surfacecnt,
					  size_t boxnum,
					  float64_t *starting_point,
					  float64_t *ray_direc,
					  float32_t *zbufpt,
					  uint32_t *surfaceidpt,
					  uint32_t *facetidpt,int trace)
/* returns nonzero of at least one nearer intersection was found */
{
  size_t cnt,subbox;
  float64_t dist;
  size_t firstidx;
  size_t polynum;
  int retval=FALSE;
  
  if (trace) fprintf(stderr,"find_ray_intersections(boxnum=%d)\n",(int)boxnum);
  
  if (!ray_box_intersection(&surf->boxcoords[boxnum*6],starting_point,ray_direc)) {
    /* Ray does not pass through our box. We can not possibly have an intersection */
    if (trace) {
      fprintf(stderr,"find_ray_intersections(): Ray does not intersect box %d\n",(int)boxnum);
      //box_contains_polygon(surf,boxnum,17281,TRUE);
    }
    return FALSE;  
  }
  //if (trace) {
  //  box_contains_polygon(surf,boxnum,17281,TRUE);
  //}

  
  if (surf->boxes[boxnum*9 + 8]!=-1) {
    /* Have index into boxpolys array: This box contains polygons */
    for (cnt=surf->boxes[boxnum*9 + 8];surf->boxpolys[cnt] >= 0;cnt++) {
      /* Got a polygon */
      polynum=surf->boxpolys[cnt];

      firstidx=surf->vertexidx_indices[polynum];
      
      dist = ray_to_plane_distance(starting_point,ray_direc, &surf->vertices[3*surf->vertexidx[firstidx + 0]], &surf->facetnormals[polynum*3]);
      
      if (trace) fprintf(stderr,"polynum=%d; zbufdist=%g\n",(int)polynum,dist);
      
      if (dist < *zbufpt) {
	
	if (ray_intersects_polygon(surf->vertices,&surf->vertexidx[firstidx],surf->numvertices[polynum],surf->facetnormals+3*polynum,surf->refpoints+3*polynum,surf->inplanemats+6*polynum,starting_point,ray_direc,surf->maxradius[polynum],dist,TRUE,trace)) {
	  // Will give us a truth value for whether we actually
	  // intersected this polygon.

	  if (trace) fprintf(stderr,"ray_intersects_polygon; zbufdist=%g\n",dist);
	  
	  // If so, record the distance into the z-buffer,
	  // (if we are closest)
	  // the surface count and facet ids into their
	  // proper places, 

	  *zbufpt=dist;
	  *surfaceidpt=surfacecnt;
	  *facetidpt=polynum;

	  retval = retval || TRUE;
	} else{
	  if (trace) fprintf(stderr,"ray does not intersect polygon\n");
	  
	}
	
      }
    }
  }
  
  /* This box may be sub-divided into 8 */
  for (subbox=0;subbox < 8; subbox++) {
    if (surf->boxes[boxnum*9+subbox] >= 0) {
      retval = find_ray_intersections(surf,surfacecnt,surf->boxes[boxnum*9 + subbox],starting_point,ray_direc,zbufpt,surfaceidpt,facetidpt,trace) || retval;  // Note order of OR is important because if it were first, retval could mask find_ray_intersections from execution
    }
  }
  return retval;
}



#define FRIN_STACKSIZE 100

static IPOPS_INLINE int find_ray_intersections_nonrecursive(struct projectionsurface *surf,
							    size_t surfacecnt,
							    size_t boxnum_first,
							    float64_t *starting_point,
							    float64_t *ray_direc,
							    float32_t *zbufpt,
							    uint32_t *surfaceidpt,
							    uint32_t *facetidpt,int trace)
/* returns nonzero if at least one nearer intersection was found */
{
  
  size_t cnt;
  size_t boxnum_stack[FRIN_STACKSIZE];
  size_t subbox;
  size_t boxnum;
  size_t stackentries=0;
  float64_t dist;
  size_t firstidx;
  size_t polynum;
  int retval=FALSE;

  /* Push given box onto stack */
  
  boxnum_stack[stackentries]=boxnum_first;
  stackentries++;
  
  while (stackentries > 0) {
    boxnum=boxnum_stack[stackentries-1];
    if (trace) fprintf(stderr,"find_ray_intersections_nonrecursive(boxnum=%d)\n",(int)boxnum);
    
    if (!ray_box_intersection(&surf->boxcoords[boxnum*6],starting_point,ray_direc)) {
      /* Ray does not pass through our box. We can not possibly have an intersection */
      if (trace) {
	fprintf(stderr,"find_ray_intersections(): Ray does not intersect box %d\n",(int)boxnum);
	//box_contains_polygon(surf,boxnum,17281,TRUE);
      }

      // Pop this entry off the stack
      stackentries--;

      // loop back
      continue;
    }
    //if (trace) {
    //  box_contains_polygon(surf,boxnum,17281,TRUE);
    //}
    
    
    if (surf->boxes[boxnum*9 + 8]!=-1) {
      /* Have index into boxpolys array: This box contains polygons */
      for (cnt=surf->boxes[boxnum*9 + 8];surf->boxpolys[cnt] >= 0;cnt++) {
	/* Got a polygon */
	polynum=surf->boxpolys[cnt];
	
	firstidx=surf->vertexidx_indices[polynum];
	
	dist = ray_to_plane_distance(starting_point,ray_direc, &surf->vertices[3*surf->vertexidx[firstidx + 0]], &surf->facetnormals[polynum*3]);
	
	if (trace) fprintf(stderr,"polynum=%d; zbufdist=%g\n",(int)polynum,dist);
	
	if (dist < *zbufpt) {
	  
	  if (ray_intersects_polygon(surf->vertices,&surf->vertexidx[firstidx],surf->numvertices[polynum],surf->facetnormals+3*polynum,surf->refpoints+3*polynum,surf->inplanemats+6*polynum,starting_point,ray_direc,surf->maxradius[polynum],dist,TRUE,trace)) {
	    // Will give us a truth value for whether we actually
	    // intersected this polygon.
	    
	    if (trace) fprintf(stderr,"ray_intersects_polygon; zbufdist=%g\n",dist);
	    
	    // If so, record the distance into the z-buffer,
	    // (if we are closest)
	    // the surface count and facet ids into their
	    // proper places, 
	    
	    *zbufpt=dist;
	    *surfaceidpt=surfacecnt;
	    *facetidpt=polynum;
	    
	    retval = retval || TRUE;
	  } else{
	    if (trace) fprintf(stderr,"ray does not intersect polygon\n");
	    
	  }
	  
	}
      }
    }

    /* Pop this box off of the stack */
    stackentries--;

    
    /* This box may be sub-divided into 8 */
    /* Push subboxes onto stack */
    for (subbox=0;subbox < 8; subbox++) {
      if (surf->boxes[boxnum*9+subbox] >= 0) {
	assert(stackentries < FRIN_STACKSIZE); /* check for stack overflow */
	
	boxnum_stack[stackentries]=surf->boxes[boxnum*9+subbox];	
	stackentries++;

      }
    }
  }
  return retval;
}



static IPOPS_INLINE void projecttoimgbuf(float32_t pixelval,float32_t pixelweighting,float64_t *uvcoords,float64_t *uvcoords_deriv_horiz, float64_t *uvcoords_deriv_vert,volatile atomicpixel_t *imgbuf, volatile atomicpixel_t *weightingbuf,volatile atomicpixel_t *validitybuf,size_t framenum,size_t imgbuf_nx,size_t imgbuf_ny,float64_t min_radius_uv_pixels,float64_t min_radius_src_pixels,float64_t bandwidth_fraction)
//size_t max_projection_pixels,uint32 *pixelnumbuf,float64_t *pixelvecbuf,float64_t *pixelmatbuf,
{

  long long arraywidth,arrayheight;
  int arrayx0,arrayy0;
  float64_t projecthalfwidth,projecthalfheight;
  float64_t newwidth,newheight;

  float64_t uvcoords0_pixels=uvcoords[0]*(imgbuf_nx) - 0.5;
  float64_t uvcoords1_pixels=uvcoords[1]*(imgbuf_ny) - 0.5;

  //CvMat *jacobian,*jacinv;
  float64_t jacobian[4],jacinv[4],detinv;

  float64_t weightingfactor;
    
  size_t xcnt,ycnt;

  float64_t r2_uv,r2_src,coeff,cosval,cosparam;  //sincparam
  float64_t pos[2],pos_frac[2],srcpos[2];
  //float64_t angle_of_incidence_factor;

  if (IPOPS_ISNAN(pixelval)) {
    return; /* never project NaN */
  }

  // Ignore anything at extreme angles of incidence
  //if (angle_of_incidence > 3*M_PI/8) return;

  //jacobian=cvCreateMat(2,2,CV_64F);

  //jacobian->data.db[0]=uvcoords_deriv_horiz[0]; 
  //jacobian->data.db[1]=uvcoords_deriv_vert[0];
  //jacobian->data.db[2]=uvcoords_deriv_horiz[1]; 
  //jacobian->data.db[3]=uvcoords_deriv_vert[1]; 
  jacobian[0]=uvcoords_deriv_horiz[0]; 
  jacobian[1]=uvcoords_deriv_vert[0];
  jacobian[2]=uvcoords_deriv_horiz[1]; 
  jacobian[3]=uvcoords_deriv_vert[1]; 
  //jacinv=cvCreateMat(2,2,CV_64F);
  
  //cvInvert(jacobian,jacinv,CV_LU);

  // 2x2 matrix inverse by determinant
  // method:
  // inv([ a b ; c d ]) = (1.0/(ad-bc)) * [d -b ; -c a ]
  detinv=1.0/(jacobian[0]*jacobian[3]-jacobian[1]*jacobian[2]);
  jacinv[0] = detinv*jacobian[3];
  jacinv[1] = -detinv*jacobian[1];
  jacinv[2] = -detinv*jacobian[2];
  jacinv[3] = detinv*jacobian[0];

  //// Define factor by which we de-emphasize data at larger angles of incidence
  // (moved into evaluate_zbuffer)
  //angle_of_incidence_factor = cos(angle_of_incidence * (M_PI/2)/(3*M_PI/8));

  
  projecthalfwidth=min_radius_uv_pixels;  // texture coordinates are relative to image size, still 
  projecthalfheight=min_radius_uv_pixels;

  newwidth=uvcoords_deriv_horiz[0]*min_radius_src_pixels* imgbuf_nx; // mul by imgbuf_nx to convert from tex coordinate to tex pixels
  if (newwidth > projecthalfwidth) projecthalfwidth=newwidth;
  newheight=uvcoords_deriv_horiz[1]*min_radius_src_pixels*imgbuf_ny;
  if (newheight > projecthalfheight) projecthalfheight=newheight;

  newwidth=uvcoords_deriv_vert[0]*min_radius_src_pixels*imgbuf_nx;
  if (newwidth > projecthalfwidth) projecthalfwidth=newwidth;
  newheight=uvcoords_deriv_vert[1]*min_radius_src_pixels*imgbuf_ny;
  if (newheight > projecthalfheight) projecthalfheight=newheight;
  
  arraywidth = (size_t) (projecthalfwidth*2+1);
  arrayheight= (size_t) (projecthalfheight*2+1);
  arrayx0 = (int)(uvcoords0_pixels-projecthalfwidth);
  arrayy0 = (int)(uvcoords1_pixels-projecthalfheight);
  
  if (arrayx0 < 0) arrayx0=0;
  if (arrayy0 < 0) arrayy0=0;
  if (arrayx0 + arraywidth >= imgbuf_nx) arraywidth = imgbuf_nx-arrayx0-1;
  if (arrayy0 + arrayheight >= imgbuf_ny) arrayheight = imgbuf_ny-arrayy0-1;
  if (arraywidth < 0) arraywidth=0;
  if (arrayheight < 0) arrayheight=0;

  //pixelcnt=0;
  for (ycnt=0;ycnt < arrayheight;ycnt++) {
    for (xcnt=0;xcnt < arraywidth;xcnt++) {
      // xcnt+arrayx0, ycnt+arrayy0 are the indices into the pixel image
      // xcnt+arrayx0=0 corresponds to the center of the leftmost pixel
      // The left edge of this pixel should map to u=0.0,u_pixels=0.0
      // this left edge corresponds to xcnt+arrayx0=-0.5
      // That gives pos[0] = 0.0 = -0.5 -uvcoords0_pixels = -0.5 -uvcoords[0]*(imgbuf_nx) + 0.5
      //                 pos[0] = 0.0 = -0.5 - 0.0 + 0.5
      //
      // xcnt+arrayx0=imgbuf_nx-1 corresponds to the rightmost pixel
      // the right edge of thix pixel should map to u=1.0
      // This right edge corresponds to xcnt+arrayx0 = imgbuf_nx - 0.5
      // that gives pos[0] = imgbuf_nx-0.5 - uvcoords0_pixels
      //            pos[0] = imgbuf_nx-0.5 - uvcoords[0]*(imgbuf_nx) + 0.5
      //            pos[0] = imgbuf_nx-0.5 - (imgbuf_nx) + 0.5
      //            pos[0] = -0.5  + 0.5  = 0 ... so we are right on target
      
      pos[0]=xcnt+arrayx0 - uvcoords0_pixels;
      pos[1]=ycnt+arrayy0 - uvcoords1_pixels;
      
      r2_uv = pos[0]*pos[0] + pos[1]*pos[1];

      // pos is in pixels so far, but jacobian/jacimg are in terms
      // of the unity width of the UV parameterization. Scale it down

      pos_frac[0]=pos[0]/imgbuf_nx;
      pos_frac[1]=pos[1]/imgbuf_ny;
      
      //multmatvec2(jacinv->data.db,pos_frac,srcpos);
      multmatvec2(jacinv,pos_frac,srcpos);
      r2_src = srcpos[0]*srcpos[0]+srcpos[1]*srcpos[1];
      //fprintf(stderr,"r_uv=%g; min_radius_uv = %g\n",sqrt(r2_uv),fabs(min_radius_uv_pixels));
      //fprintf(stderr,"r_src=%g; min_radius_src = %g\n",sqrt(r2_src),fabs(min_radius_src_pixels));
      if (r2_uv <= min_radius_uv_pixels*min_radius_uv_pixels ||
	  r2_src <= min_radius_src_pixels*min_radius_src_pixels) {
	/* Include this point */
	// Forget complicated bandlimited interpolation
	// Instead project 2D generalized circular sinc function
	// from source into UV space

	
	//if (pixelcnt >= max_projection_pixels) {
	//  /* too many pixels */
	//  goto fail; 
	//}

	weightingfactor=1.0;
	if (weightingbuf) {
	  //weightingfactor=weightingbuf[(arrayy0+ycnt)*imgbuf_nx + (arrayx0+xcnt)];
	  weightingfactor=atomicpixel_load(&weightingbuf[(arrayy0+ycnt)*imgbuf_nx + (arrayx0+xcnt)]);
	}
	
	// Generalized 2D circular sinc: based on
	// http://www.ebyte.it/library/docs/math07/SincN.html
	// Eq. 18: For 2D, n=2 
	// sinc(n,x) = gamma(1+n/2)*J_(n/2)(x)/(x/2)^(n/2)
	// at n=2
	// sinc(2,x) = gamma(2)*J1(x)/(x/2)
	// where gamma(2) is 1.0 and J1 is the Bessel function of the first kind
	// WARNING: This fails at x=0
	//
	//sincparam=sqrt(r2_src)*bandwidth_fraction;
	//if (sincparam < 1e-6) {
	//  // include weighting effects of both the pixel in the camera
	//  // image (pixelweighting) and the u,v space weighting (weightingbuf)
	//  coeff=pixelweighting*weightingfactor; //angle_of_incidence_factor; // effectively 1.0, not including angle of incidence correction
	//} else{
	//  // include weighting effects of both the pixel in the camera
	//  // image (pixelweighting) and the u,v space weighting (weightingbuf)
	//  coeff=pixelweighting*weightingfactor * j1(sincparam)/(sincparam/2.0);  //angle_of_incidence_factor * j1(sincparam)/(sincparam/2.0);
	//}

	// Replacement -- Just use raised Cosine
	// 1+cos(sqrt(r2_src)*bandwidth_fraction)
	cosparam=sqrt(r2_src)*bandwidth_fraction*M_PI/M_SQRT2; // 
	if (cosparam > M_PI) {
	  cosval=0.0;
	} else {
	  cosval=0.5+0.5*cos(cosparam);
	}
	coeff=pixelweighting*weightingfactor*cosval;
	
	//fprintf(stderr,"imgbuf[%d]+=%g\n",framenum*imgbuf_ny*imgbuf_nx + (arrayy0+ycnt)*imgbuf_nx + (arrayx0+xcnt),coeff*pixelval);

	// ***!!!!***!!!  THIS IS WHERE THE PARALLEL CPU THREAD TRIP OVER EACH OTHER
	// SHOULD WRITE WITH ATOMIC OPERATIONS
	if (imgbuf) {
	  //imgbuf[framenum*imgbuf_ny*imgbuf_nx + (arrayy0+ycnt)*imgbuf_nx + (arrayx0+xcnt)] += coeff * pixelval;
	  atomicpixel_accumulate(&imgbuf[framenum*imgbuf_ny*imgbuf_nx + (arrayy0+ycnt)*imgbuf_nx + (arrayx0+xcnt)], coeff*pixelval);
	  
	}
	//validitybuf[framenum*imgbuf_ny*imgbuf_nx + (arrayy0+ycnt)*imgbuf_nx + (arrayx0+xcnt)] += coeff;
	atomicpixel_accumulate(&validitybuf[framenum*imgbuf_ny*imgbuf_nx + (arrayy0+ycnt)*imgbuf_nx + (arrayx0+xcnt)],coeff);
	
	//if (angleofincidencebuf) angleofincidencebuf[framenum*imgbuf_ny*imgbuf_nx + (arrayy0+ycnt)*imgbuf_nx + (arrayx0+xcnt)] += coeff * angle_of_incidence;
	//pixelnumbuf[pixelcnt]=ycnt*arraywidth + xcnt;
	//pixelcnt++;
	      
      }
    }
  }
  
  //numpixels=pixelcnt;

  /*
  // Form inner product matrix 
  for (pixelcnt=0;pixelcnt < numpixels;pixelcnt++) {
    xcnt=pixelnumbuf[pixelcnt] % arraywidth;
    ycnt=pixelnumbuf[pixelcnt] / arraywidth;

    pos[0]=xcnt+arrayx0 - uvcoords0_pixels;
    pos[1]=ycnt+arrayy0 - uvcoords1_pixels;

    // pos is in pixels so far, but jacobian/jacimg are in terms
    // of the width of the UV image. Scale it down
    pos_frac[0]=pos[0]/imgbuf_nx;
    pos_frac[1]=pos[1]/imgbuf_ny;

    multmatvec2(jacinv->data.dbl,pos_frac,srcpos);
    // srcpos indicates how many pixels this point is off of the
    // source image

    src_r=sqrt(srcpos[0]*srcpos[0]+srcpos[1]*srcpos[1]);
	       
    
    //// inner product of sinc function centered at pos, parameter
    //// bandwidth_fraction*pixel_pos, with sinc function
    //// centered at (uvcoords0_pixels,uvcoords1_pixels)
    //pixelmatbuf[pixelcnt]=bandwidth_fraction*bandwidth_fraction*(improj_sinc(bandwidth_fraction * pos[0]) + improj_sinc 


    
									        
    // srcpos is now (horizontal, vertical) shift in pixels in source image
    */
  //cvReleaseMat(&jacobian);
  //cvReleaseMat(&jacinv);

}

static IPOPS_INLINE void basic_raycalcs(struct projectionsurface *surfacearray,
					size_t surfaceid,
					float64_t *cam_mtx_Ut,
					float64_t *cam_mtx_Vtsubinv,
					float64_t s1,
					float64_t s2,
					float64_t V31,
					float64_t V32,
					size_t src_pixelx,
					size_t src_pixely,
					/* outputs ... */
					double *dc1_du0, /* also a flag for whether we need derivs */
					double *dc2_du0,
					double *dc1_dv0,
					double *dc2_dv0,
					double *rayvecobj, /* 4-vector */
					double *rayvecfactor)
{
  float64_t q,r,c1,c2;  
  double dq_du0,dr_du0,dq_dv0,dr_dv0;
  float64_t rayvec[4];
  
  //Ut is 2x2
  //q=cam_mtx_Ut[0,0]*xcnt + cam_mtx_Ut[0,1]*ycnt;
  //r=cam_mtx_Ut[1,0]*xcnt + cam_mtx_Ut[1,1]*ycnt;
  q=cam_mtx_Ut[0]*src_pixelx + cam_mtx_Ut[1]*src_pixely;
  r=cam_mtx_Ut[2+0]*src_pixelx + cam_mtx_Ut[2+1]*src_pixely;
  
  //c1=cam_mtx_Vtsubinv[0,0]*(q/s1-V31) + cam_mtx_Vtsubinv[0,1]*(r/s2-V32);
  //c2=cam_mtx_Vtsubinv[1,0]*(q/s1-V31) + cam_mtx_Vtsubinv[1,1]*(r/s2-V32);
  // Vtsubinv is 2x2
  c1=cam_mtx_Vtsubinv[0]*(q/s1-V31) + cam_mtx_Vtsubinv[1]*(r/s2-V32);
  c2=cam_mtx_Vtsubinv[2+0]*(q/s1-V31) + cam_mtx_Vtsubinv[2+1]*(r/s2-V32);


  if (dc1_du0) {
    // also store dervatives
    
    //dq_du0=cam_mtx_Ut[0,0];
    //dr_du0=cam_mtx_Ut[1,0];
    //dq_dv0=cam_mtx_Ut[0,1];
    //dr_dv0=cam_mtx_Ut[1,1];
    
    // Ut is 2x2
    dq_du0=cam_mtx_Ut[0];
    dr_du0=cam_mtx_Ut[2+0];
    dq_dv0=cam_mtx_Ut[1];
    dr_dv0=cam_mtx_Ut[2+1];
    
    
    //dc1_du0 = cam_mtx_Vtsubinv[0,0]*(dq_du0/s1) + cam_mtx_Vtsubinv[0,1]*(dr_du0/s2);
    //dc2_du0 = cam_mtx_Vtsubinv[1,0]*(dq_du0/s1) + cam_mtx_Vtsubinv[1,1]*(dr_du0/s2);
    
    //dc1_dv0 = cam_mtx_Vtsubinv[0,0]*(dq_dv0/s1) + cam_mtx_Vtsubinv[0,1]*(dr_dv0/s2);
    //dc2_dv0 = cam_mtx_Vtsubinv[1,0]*(dq_dv0/s1) + cam_mtx_Vtsubinv[1,1]*(dr_dv0/s2);
    // Vtsubinv is 2x2
    *dc1_du0 = cam_mtx_Vtsubinv[0]*(dq_du0/s1) + cam_mtx_Vtsubinv[1]*(dr_du0/s2);
    *dc2_du0 = cam_mtx_Vtsubinv[2+0]*(dq_du0/s1) + cam_mtx_Vtsubinv[2+1]*(dr_du0/s2);
    
    *dc1_dv0 = cam_mtx_Vtsubinv[0]*(dq_dv0/s1) + cam_mtx_Vtsubinv[1]*(dr_dv0/s2);
    *dc2_dv0 = cam_mtx_Vtsubinv[2+0]*(dq_dv0/s1) + cam_mtx_Vtsubinv[2+1]*(dr_dv0/s2);
  }
    
  // OK... We have a ray leaving (0,0,0) in the direction
  // of (unnormalized) (c1,c2,1) in camera coordinates
	
  rayvec[0]=c1;
  rayvec[1]=c2;
  rayvec[2]=1.0;
  rayvec[3]=0.0; // # "vector" has 0 w-component

  /* apply correction for camera matrix being in OpenCV coords, but
     our transforms being in OpenGL coords */
  /* We negate rayvec[1] and rayvec[2] rather than messing with tosurfaceframe. */
  
  rayvec[1]=-rayvec[1];
  rayvec[2]=-rayvec[2];
  
  multmatvec4(surfacearray[surfaceid].tosurfaceframe,rayvec,rayvecobj);
  *rayvecfactor = to_unit_vector4(rayvecobj);

}

static IPOPS_INLINE void deriv_raycalcs(struct projectionsurface *surfacearray,
					float64_t *focalpointobj,
					float64_t *rayvecobj,
					double rayvecfactor,
					size_t surfaceid,
					size_t facetid,
					double dc1_du0,
					double dc2_du0,
					double dc1_dv0,
					double dc2_dv0,
					/* outputs */
					float64_t *horiz_ray_intersect_shift_deriv,
					float64_t *vert_ray_intersect_shift_deriv)
{
  float64_t drayvec_du0[4];
  float64_t drayvec_du0obj[4];
  float64_t drayvec_dv0[4];
  float64_t drayvec_dv0obj[4];
  size_t firstidx;
  float64_t horiz_ray_direc_shift[3];
  float64_t vert_ray_direc_shift[3];
  
  /* Note minus signs inserted in elements 1 and 2 to compensate for OpenGL vs. OpenCV and needing to negate columns of tosurfaceframe  */
      
  drayvec_du0[0]=dc1_du0*rayvecfactor;
  drayvec_du0[1]=-dc2_du0*rayvecfactor;
  drayvec_du0[2]=0.0; // -1.0*rayvecfactor;
  drayvec_du0[3]=0.0;
  multmatvec4(surfacearray[surfaceid].tosurfaceframe,drayvec_du0,drayvec_du0obj);
	
  drayvec_dv0[0]=dc1_dv0*rayvecfactor;
  drayvec_dv0[1]=-dc2_dv0*rayvecfactor;
  drayvec_dv0[2]=0.0; // -1.0*rayvecfactor;
  drayvec_dv0[3]=0.0;
  multmatvec4(surfacearray[surfaceid].tosurfaceframe,drayvec_dv0,drayvec_dv0obj);
  
  firstidx=surfacearray[surfaceid].vertexidx_indices[facetid];
  
  
  // Figure out what happens as we shift to adjacent pixel
  // in image being projected.
  // ray_direc + ray_direc_deriv_horiz -> horizshifted_ray_direc
  
  // Consider derivative of
  // vector normalization: deriv (a/|a|)
  //  = (deriv(a)*|a| - (deriv(|a|)  a)/|a|^2
  //  = deriv(a)/|a| - (a/|a|) deriv(|a|)/|a|
  // since a is already normalized, |a| = 1.0
  // so = deriv(a) - a deriv|a|
  // where |a| = sqrt(x*2 + y*2 + z*2)
  // so deriv |a| = d|a|/dwhatever = 0.5|a|^(-1)*2x dx/dwhatever + 0.5|a|^(-1)*2y dy/dwhatever + 0.5*|a|^(-1)*2z dz/dwhatever
  // deriv |a| = (1/|a|) x dx/dwhatever + y dy/dwhatever + z dz/dwhatever
  // deriv |a| = (1/|a|) (a dot deriv(a))
  // So full derivative = deriv(a) - a ((1/|a|) (a dot deriv(a)))
  //      = deriv(a) - a (a dot deriv(a))
  
  
  // Therefore,
  // horiz_ray_direc_shift = ray_direc_deriv_horiz - ray_direc * (dot(ray_direc,ray_direc_deriv_horiz)
  addvecscaledvec3(drayvec_du0obj,-dotvecvec3(rayvecobj,drayvec_du0obj),rayvecobj,horiz_ray_direc_shift);
	
  // horiz_ray_direct_shift is the derivative of the ray unit vector with respect to a horizontal pixel shift, in object coordinates
	
  // Determine intersect position shift in object coordinates with respect to a unit horizontal pixel shift 
  ray_to_plane_raydirec_shift(focalpointobj,rayvecobj,&surfacearray[surfaceid].vertices[3*surfacearray[surfaceid].vertexidx[firstidx + 0]], &surfacearray[surfaceid].facetnormals[facetid*3],horiz_ray_direc_shift,horiz_ray_intersect_shift_deriv);

  // Do same for vert, store vert_ray_intersect_shift_deriv 
  addvecscaledvec3(drayvec_dv0obj,-dotvecvec3(rayvecobj,drayvec_dv0obj),rayvecobj,vert_ray_direc_shift);
	
  // vert_ray_direct_shift is the derivative of the ray unit vector with respect to a vertical pixel shift, in object coordinates
	
  // Determine intersect position shift in object coordinates with respect to unit a vertical pixel shift 
  ray_to_plane_raydirec_shift(focalpointobj,rayvecobj,&surfacearray[surfaceid].vertices[3*surfacearray[surfaceid].vertexidx[firstidx + 0]], &surfacearray[surfaceid].facetnormals[facetid*3],vert_ray_direc_shift,vert_ray_intersect_shift_deriv);
  
}




static IPOPS_INLINE void orthographic_deriv_raycalcs(struct projectionsurface *surfacearray,
						     float64_t *raysrclocobj,
						     float64_t *rayvecobj,
						     size_t surfaceid,
						     size_t facetid,
						     float64_t *horiz_raysrc_shift,
						     float64_t *vert_raysrc_shift,
						     /* outputs */
						     float64_t *horiz_ray_intersect_shift_deriv,
						     float64_t *vert_ray_intersect_shift_deriv)
{
  size_t firstidx;
  
  
  firstidx=surfacearray[surfaceid].vertexidx_indices[facetid];
  
  
  // Figure out what happens as we shift to adjacent pixel
  // in image being projected.
  // ray_direc remains constant

  // raysrcloc + horiz_raysrc_shift -> horizshifted_raysrc
  
	
  // Determine intersect position shift in object coordinates with respect to a unit horizontal pixel shift 
  ray_to_plane_raypos_shift(raysrclocobj,rayvecobj,&surfacearray[surfaceid].vertices[3*surfacearray[surfaceid].vertexidx[firstidx + 0]], &surfacearray[surfaceid].facetnormals[facetid*3],horiz_raysrc_shift,horiz_ray_intersect_shift_deriv);

  // Do same for vert, store vert_ray_intersect_shift_deriv 
  // Determine intersect position shift in object coordinates with respect to unit a vertical pixel shift 
  ray_to_plane_raypos_shift(raysrclocobj,rayvecobj,&surfacearray[surfaceid].vertices[3*surfacearray[surfaceid].vertexidx[firstidx + 0]], &surfacearray[surfaceid].facetnormals[facetid*3],vert_raysrc_shift,vert_ray_intersect_shift_deriv);
  
}


void find_intersect_uv(struct projectionsurface *surfacearray,
		       float64_t *focalpointobj,
		       float64_t *rayvecobj,
		       size_t surfaceid,
		       size_t facetid,
		       float64_t dist,
		       float64_t *horiz_ray_intersect_shift_deriv, // only needed if calculating derivatives
		       float64_t *vert_ray_intersect_shift_deriv, // only needed if calculating derivatives
		       /* outputs */
		       float64_t *intersectpoint2uvcanon, // intersection point in texture coordinates (not scaled UV coordinates
		       float64_t *intersectpoint2uvcanon_deriv_horiz, // also flag for whether to calc derivatives
		       float64_t *intersectpoint2uvcanon_deriv_vert)
{
  float64_t intersectpoint[3];
  float64_t intersectpoint3poly[3];
  float64_t intersectpoint2poly[3];
      

  // Each polygon has a refpoint (point centroid),
  // A normal, two in-plane basis vectors
  // (inplanevec1,inplanevec2) as rows
  // of a 2x3 matrix inplanemat, calculated from from SVD of
  // (polygon points - polygon centroid)
  // and a 2x3 matrix that gives texture
  // (intrinsic parameterization) coordinates from
  // (inplanevec1, inplanevec2, 1) 

  //// Evaluate coordinates relative to centroid
  //numvertices=0
  //for CCnt in range(surfacearray[surfaceid].max_vertices):
  //    thisvertexid=surfacearray[surfaceid].vertexids[surfacearray[surfaceid].max_vertices*facetid+CCnt]
  //    if thisvertexid < 0:
  //        // This polygon has no more vertices
  //        break
  //
  //    subvecvec3(&surfacearray[surfaceid].vertices[thisvertexid*3],&surfacearray[surfaceid].refpt,&relcoords[CCnt*3])
  //    
  //    numvertices+=1
  //    pass
  
  // find intersectpoint in object coordinates
  addvecscaledvec3(focalpointobj,dist,rayvecobj,intersectpoint);
	

	
  // Also store the (u,v) coordinates projected so that they can be
  // used to reduce the weightings near edge boundaries.

  
  // Evaluate 3D intersection point relative to polygon
  subvecvec3(intersectpoint,&surfacearray[surfaceid].refpoints[3*facetid],intersectpoint3poly);
  
  // Evaluate 2D intersection point relative to polygon
  multmat23vec(&surfacearray[surfaceid].inplanemats[6*facetid],intersectpoint3poly,intersectpoint2poly);
	
  // be sure to allocate 3 numbers of space
  // for intersectpoint2poly & friends as we will use it
  // in projective form	  
  intersectpoint2poly[2]=1.0; // A point in projective space
	
  // Evaluate 2D polygon to (u,v) (transform in-plane 2D coords -> intrinsic texture coords)
  multmat23vec(&surfacearray[surfaceid].inplane2texcoords[6*facetid],intersectpoint2poly,intersectpoint2uvcanon);  // use ,intersectpoint2tex) once we have a texcoords->uvcanon mapping
  

  if (intersectpoint2uvcanon_deriv_horiz) {
    /* if calculate derivatives... */
    float64_t intersectpoint2poly_deriv_horiz[3]; // needs 3 elements
    float64_t intersectpoint2poly_deriv_vert[3];
    
    // Evaluate derivatives of 2D intersection point

    // be sure to allocate 3 numbers of space
    // for intersectpoint2poly & friends as we will use it
    // in projective form

    multmat23vec(&surfacearray[surfaceid].inplanemats[6*facetid],horiz_ray_intersect_shift_deriv,intersectpoint2poly_deriv_horiz);
    multmat23vec(&surfacearray[surfaceid].inplanemats[6*facetid],vert_ray_intersect_shift_deriv,intersectpoint2poly_deriv_vert);

    intersectpoint2poly_deriv_horiz[2]=0.0; // derivatives are vectors in 2-space, intepreted in a projective space
    intersectpoint2poly_deriv_vert[2]=0.0;

    // Derivatives in (u,v)
    multmat23vec(&surfacearray[surfaceid].inplane2texcoords[6*facetid],intersectpoint2poly_deriv_horiz,intersectpoint2uvcanon_deriv_horiz); // use intersectpoint2tex_deriv_horiz once we have a texcoords->uvcanon coords mapping )
    multmat23vec(&surfacearray[surfaceid].inplane2texcoords[6*facetid],intersectpoint2poly_deriv_vert,intersectpoint2uvcanon_deriv_vert); // use intersectpoint2tex_deriv_vert once we have a texcoords->uvcanon coords mapping )

    
  }
  
  // Now we have uv parameterization coordinates
  // intersectpoint2uvcanon as well as derivatives
  // of those coordinates representing motion of
  // one pixel in the projected image to the right (horiz)
  // and one pixel down (vert).

}


static IPOPS_INLINE void evaluate_zbuffer(struct projectionsurface *surfacearray,
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
					  float32_t *imagedata_horiz_ray_intersect_shift_deriv_cam_z,
					  float32_t *imagedata_vert_ray_intersect_shift_deriv_cam_z,
					  float32_t *imagedata_uvcoords /* in texture coordinates [0,1], not physical scaled coordinates */)
{
  int64_t surfaceycnt;  // must be signed (e.g. not size_t) for MSVC compatibility

#pragma omp parallel default(shared) private(surfaceycnt) // ,src_ny,numsurfaces)
  {
#pragma omp for
    for (surfaceycnt=0;surfaceycnt < src_ny*numsurfaces;surfaceycnt++) {
      
      size_t ycnt,surfacecnt,xcnt;
      float64_t focalpointobj[4];

      float64_t rayvecobj[4]; // rayvec in object coordinates
      float64_t rayvecfactor;


      int trace=FALSE;
      size_t firstidx;
      //float64_t rayvecfactor;
      
      ycnt=surfaceycnt % src_ny; //  numsurfaces;
      surfacecnt=surfaceycnt / src_ny; // numsurfaces;
      
      // calculation of focalpointobj: Coordinates of focal point in object coordinates for this surface.
      // since we are doing this OpenCV style, negate the 2nd and 3rd columns... see imageprojectionmodel.py
      // (except we don't actually do this because we are just accessing the 4th column)
      focalpointobj[0]=surfacearray[surfacecnt].tosurfaceframe[0*4 + 3];
      focalpointobj[1]=surfacearray[surfacecnt].tosurfaceframe[1*4 + 3];
      focalpointobj[2]=surfacearray[surfacecnt].tosurfaceframe[2*4 + 3];
      focalpointobj[3]=surfacearray[surfacecnt].tosurfaceframe[3*4 + 3];
      normalize_wcoord4(focalpointobj);
      
      for (xcnt=0; xcnt < src_nx; xcnt++) {
	// Go through each source image pixel,
	// project it through to surface. If it
	// intersects closer, mark the z-buffer
	// and surface ID. This will give us a map
	// of the source image for each pixel, which surface
	// it maps onto.
	
	// NOTE: These are recalculated later in project_image... formulas should be THE SAME!!!
	
	basic_raycalcs(surfacearray,
		       surfacecnt,
		       cam_mtx_Ut,
		       cam_mtx_Vtsubinv,
		       s1,s2,V31,V32,
		       xcnt,ycnt,
		       /* outputs ... */
		       NULL, /* also a flag for whether we need derivs */
		       NULL, NULL, NULL,
		       rayvecobj, /* 4-vector */
		       &rayvecfactor);
			


	if (FALSE && ycnt==237 && xcnt==352) {
	  int polynum;
	  double dist;
	  int hitcnt=0;
	  struct projectionsurface *surf;

	  /* debugging */
	  fprintf(stderr,"y=%d/%d, x=%d/%d\n",(int)ycnt,(int)src_ny,(int)xcnt,(int)src_nx);
	  trace=TRUE;

	  surf=&surfacearray[surfacecnt];
	  
	  for (polynum = 0; polynum < surf->npolys;polynum++) {
	    firstidx=surf->vertexidx_indices[polynum];
	    dist = ray_to_plane_distance(focalpointobj,rayvecobj, &surf->vertices[3*surf->vertexidx[firstidx + 0]], &surf->facetnormals[polynum*3]);
	    
	    if (ray_intersects_polygon(surf->vertices,&surf->vertexidx[firstidx],surf->numvertices[polynum],surf->facetnormals+3*polynum,surf->refpoints+3*polynum,surf->inplanemats+6*polynum,focalpointobj,rayvecobj,surf->maxradius[polynum],dist,TRUE,trace && FALSE)) {
	      fprintf(stderr,"ray intersects polygon %d; dist=%f\n",(int)polynum,dist);
	      hitcnt++;
	    }
	  }
	  if (!hitcnt) {
	    fprintf(stderr,"ray did not intersect any polygons\n");
	  }
	} else {
	  trace=FALSE;
	}
	
	// find_ray_intersections() should fill out the z-buffer, surface id, and facet id arrays for the closest intersection 
	find_ray_intersections(&surfacearray[surfacecnt],surfacecnt,0,focalpointobj,rayvecobj,imagedata_zbuffer + src_nx*ycnt + xcnt,imagedata_surfaceid + src_nx*ycnt + xcnt,imagedata_facetid + src_nx*ycnt + xcnt,trace);
	
      }
    }
  }

  /* Next phase: evaluate derivatives, etc. */
#pragma omp parallel default(shared) 
  {
    int64_t ycnt; // must be signed (e.g. not size_t) for MSVC compatibility
    static const uint32_t NaNconst=0x7fc00000;

    float32_t NaNval;

    memcpy(&NaNval,&NaNconst,sizeof(NaNval));

#pragma omp for private(ycnt) // ,src_ny)
    for (ycnt=0;ycnt < src_ny;ycnt++) {
      int64_t xcnt; // must be signed (e.g. not size_t) for MSVC compatibility
      float64_t focalpointobj[4];
      double dc1_du0,dc2_du0,dc1_dv0,dc2_dv0;
      
      float64_t rayvecobj[4]; // rayvec in object coordinates
      float64_t rayvecfactor;


      float64_t dist;

      float64_t horiz_ray_intersect_shift_deriv[3];
      float64_t horiz_ray_intersect_shift_deriv_cam[3];
      float64_t vert_ray_intersect_shift_deriv[3];
      float64_t vert_ray_intersect_shift_deriv_cam[3];


      float64_t intersectpoint2uvcanon[3]; // intersection point in texture coordinates (not scaled UV coordinates
      float64_t intersectpoint2uvcanon_deriv_horiz[3];
      float64_t intersectpoint2uvcanon_deriv_vert[3];

      float64_t angle_of_incidence,angle_of_incidence_factor;
      

      uint32_t surfaceid,facetid;



      for (xcnt=0; xcnt < src_nx; xcnt++) {
	// Go through each source image pixel, find derivatives and where it projects in (u,v) space
	
	// Dist was calculated by find_ray_intersections(), above
	dist = imagedata_zbuffer[src_nx*ycnt + xcnt];

	if (!isinf(dist)) { /* if we actually got a ray intersection above */
	  surfaceid=imagedata_surfaceid[ycnt*src_nx + xcnt];
	  facetid=imagedata_facetid[ycnt*src_nx + xcnt];


	  // calculation of focalpointobj: Coordinates of focal point in object coordinates for this surface.
	  // since we are doing this OpenCV style, negate the 2nd and 3rd columns... see imageprojectionmodel.py
	  // (except we don't need to do this because we are just extracting the 4th column) 
	  focalpointobj[0]=surfacearray[surfaceid].tosurfaceframe[0*4 + 3];
	  focalpointobj[1]=surfacearray[surfaceid].tosurfaceframe[1*4 + 3];
	  focalpointobj[2]=surfacearray[surfaceid].tosurfaceframe[2*4 + 3];
	  focalpointobj[3]=surfacearray[surfaceid].tosurfaceframe[3*4 + 3];
	  normalize_wcoord4(focalpointobj);
	  
	  basic_raycalcs(surfacearray,
			 surfaceid,
			 cam_mtx_Ut,
			 cam_mtx_Vtsubinv,
			 s1,s2,V31,V32,
			 xcnt,ycnt,
			 /* outputs ... */
			 &dc1_du0, /* also a flag for whether we need derivs */
			 &dc2_du0,
			 &dc1_dv0,
			 &dc2_dv0,
			 rayvecobj, /* 4-vector */
			 &rayvecfactor);

	  

	  // Evaluate derivatives... We convert them to object coordinates, evaluate the shifts
	  // in object coordinates, then shift back.
	  
	  // (it might be marginally more efficient to evaluate the derivatives in camera coordinates) 
	  
	  
	  deriv_raycalcs(surfacearray,
			 focalpointobj,
			 rayvecobj,rayvecfactor,
			 surfaceid,
			 facetid,
			 dc1_du0,dc2_du0,dc1_dv0,dc2_dv0,
			 /* outputs */
			 horiz_ray_intersect_shift_deriv,
			 vert_ray_intersect_shift_deriv);
	  
	  // Determine intersect position shift in camera coordinates with respect to a unit horizontal pixel shift 
	  multmatvec3(surfacearray[surfaceid].vecfromsurfaceframe,horiz_ray_intersect_shift_deriv,horiz_ray_intersect_shift_deriv_cam);
	  
	  /* evaluate left and right predicted z values... add z coordinate of derivative to this pixel's position */
	  //imagedata_projectzright[src_nx*ycnt + xcnt] = dist + 1.0*horiz_ray_intersect_shift_deriv_cam[2];
	  //imagedata_projectzleft[src_nx*ycnt + xcnt] = dist - 1.0*horiz_ray_intersect_shift_deriv_cam[2];
	  imagedata_horiz_ray_intersect_shift_deriv_cam_z[src_nx*ycnt + xcnt]=1.0*horiz_ray_intersect_shift_deriv_cam[2];
	  
	  
	  
	  // Determine intersect position shift in camera coordinates with respect to a unit vertical pixel shift 
	  multmatvec3(surfacearray[surfaceid].vecfromsurfaceframe,vert_ray_intersect_shift_deriv,vert_ray_intersect_shift_deriv_cam);
	  
	  
	  /* evaluate up and down predicted z values... add z coordinate of derivative to this pixel's position */
	  /* NOTE: We are assuming y ordering is like a raster scan (down is increasing y) */
	  //imagedata_projectzdown[src_nx*ycnt + xcnt] = dist + 1.0*vert_ray_intersect_shift_deriv_cam[2];
	  //imagedata_projectzup[src_nx*ycnt + xcnt] = dist - 1.0*vert_ray_intersect_shift_deriv_cam[2];
	  imagedata_vert_ray_intersect_shift_deriv_cam_z[src_nx*ycnt + xcnt]=vert_ray_intersect_shift_deriv_cam[2];
	  
	  // So horiz_ray_intersect_shift_deriv
	  // and vert_ray_intersect_shift_deriv
	  // are vectors indicating what happens when you move
	  // one pixel to the right and one pixel down
	  // on the image being projected
	  

	
	
	
	/* determine angle of incidence, for weighting use */
	  angle_of_incidence = acos(fabs(dotvecvec3(rayvecobj,&surfacearray[surfaceid].facetnormals[imagedata_facetid[src_nx*ycnt + xcnt]*3]))); // rayvecobj and facetnormals should be unit vectors 
	  
	  // Ignore anything at extreme angles of incidence
	  if (angle_of_incidence > 3*M_PI/8) {
	    angle_of_incidence_factor=0.0;
	  } else {  
	    // Define factor by which we de-emphasize data at larger angles of incidence
	    angle_of_incidence_factor = cos(angle_of_incidence * (M_PI/2)/(3*M_PI/8));	    
	  }
	  imagedata_angleofincidence_weighting[src_nx*ycnt + xcnt] = angle_of_incidence_factor;
	  
	  if (imagedata_angleofincidence) {
	    imagedata_angleofincidence[src_nx*ycnt + xcnt] = angle_of_incidence;
	  }
	  

	  
	  
	  
	  find_intersect_uv(surfacearray,
			    focalpointobj,
			    rayvecobj,
			    surfaceid,
			    imagedata_facetid[src_nx*ycnt + xcnt],
			    dist,
			    horiz_ray_intersect_shift_deriv, // only needed if calculating derivatives
			    vert_ray_intersect_shift_deriv, // only needed if calculating derivatives
			    /* outputs */
			    intersectpoint2uvcanon, // intersection point in texture coordinates (not scaled UV coordinates
			    intersectpoint2uvcanon_deriv_horiz, // also flag for whether to calc derivatives
			    intersectpoint2uvcanon_deriv_vert);
	  
	  
	  imagedata_uvcoords[src_nx*ycnt*2 + xcnt*2]=intersectpoint2uvcanon[0]; // Store u coordinate
	  imagedata_uvcoords[src_nx*ycnt*2 + xcnt*2 + 1]=intersectpoint2uvcanon[1];  // store v coordinate
	  
	  if (surfacearray[surfaceid].angle_of_incidence_factor_buf_uv) {
	    /* if buffer to store angle_of_incidence_factor in uv coordinate frame is provided... */
	    /* use projecttoimgbuf() to map out projection validity region */
	    /* we provide angle_of_incidence_factor as the pixel value
	       so that imgbuf ends up with the sum of projection weighting * angle of incidence factor 
	       whereas validitybuf ends up with the sum of the angle of incidence factors */
	    projecttoimgbuf(angle_of_incidence_factor,1.0,intersectpoint2uvcanon,intersectpoint2uvcanon_deriv_horiz,intersectpoint2uvcanon_deriv_vert,surfacearray[surfaceid].angle_of_incidence_factor_buf_uv,NULL,surfacearray[surfaceid].validitybuf,0,surfacearray[surfaceid].nx,surfacearray[surfaceid].ny,MIN_RADIUS_UV_PIXELS,MIN_RADIUS_SRC_PIXELS,BANDWIDTH_FRACTION);
	  }
	} else {
	  /* No intersection found at all */
	  imagedata_horiz_ray_intersect_shift_deriv_cam_z[src_nx*ycnt + xcnt]=NaNval;
	  imagedata_vert_ray_intersect_shift_deriv_cam_z[src_nx*ycnt + xcnt]=NaNval;
	  imagedata_angleofincidence_weighting[src_nx*ycnt + xcnt] = 0.0;
	  imagedata_angleofincidence[src_nx*ycnt + xcnt]=NaNval;
	  imagedata_uvcoords[src_nx*ycnt*2 + xcnt*2]=NaNval;
	  imagedata_uvcoords[src_nx*ycnt*2 + xcnt*2 + 1]=NaNval;
	}
      }
      
    }
  }
  
  
}

static IPOPS_INLINE void project_image(struct projectionsurface *surfacearray,
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
				       float32_t *imagedata_weighting,int debug)
{

  int64_t ycnt;  // must be signed (e.g. not size_t) for MSVC compatibility

  //omp_lock_t mutex;
  //omp_init_lock(&mutex);
  
  {
    // ***!!!! If you allow this loop to be parallel -- either by
    // enabling OpenMP or by not restricting it to num_threads(1)
    // then (even if you use the commented-out lock to restrain more
    // than one thread from executing in parallel), you get
    // slightly different (!!!) results from the different threads,
    // which are presumably running on different CPU cores.
    // The difference seems to be on the order of machine precision
    // (i.e. O(1e-7) for single precision
    
    // See demos/testregistration_projectionmismatch.py for
    // testcase
#pragma omp parallel for default(shared) private(ycnt)  // num_threads(1)
    for (ycnt=0;ycnt < src_ny;ycnt++) {
      //omp_set_lock(&mutex);
      size_t xcnt;

      for (xcnt=0; xcnt < src_nx; xcnt++) {

	if (!IPOPS_ISINF(imagedata_zbuffer[ycnt*src_nx + xcnt])) {

	  int32_t surfaceid,facetid;
	  float64_t focalpointobj[4];
	  double dc1_du0,dc2_du0,dc1_dv0,dc2_dv0;

	  float64_t rayvecobj[4]; // rayvec in object coordinates
	  float64_t rayvecfactor;
	  
	  
	  float64_t dist;
	  float64_t horiz_ray_intersect_shift_deriv[3];
	  float64_t vert_ray_intersect_shift_deriv[3];
	  float64_t intersectpoint2uvcanon[3]; // intersection point in texture coordinates (not scaled UV coordinates
	  float64_t intersectpoint2uvcanon_deriv_horiz[3];
	  float64_t intersectpoint2uvcanon_deriv_vert[3];
	  //float64_t angle_of_incidence;

	  
	  /* finite z-buffer means this point maps onto the object */
	  surfaceid=imagedata_surfaceid[ycnt*src_nx + xcnt];
	  facetid=imagedata_facetid[ycnt*src_nx + xcnt];

	  /* identify the focal point in object coordinates, as in evaluate_zbuffer */
	  // calculation of focalpointobj: Coordinates of focal point in object coordinates for this surface.
	  // since we are doing this OpenCV style, negate the 2nd and 3rd columns... see imageprojectionmodel.py
	  // (except we don't actually do this because we are just extracting the 4th column)
	  focalpointobj[0]=surfacearray[surfaceid].tosurfaceframe[0*4 + 3];
	  focalpointobj[1]=surfacearray[surfaceid].tosurfaceframe[1*4 + 3];
	  focalpointobj[2]=surfacearray[surfaceid].tosurfaceframe[2*4 + 3];
	  focalpointobj[3]=surfacearray[surfaceid].tosurfaceframe[3*4 + 3];
	  
	  normalize_wcoord4(focalpointobj);

	  // Then we go back through each pixel and
	  // perform the mapping, stamping down sinc-interpolation
	  // (or similar lowpass interpolation) 
	  // patches onto the parameterization output image.

	  // NOTE: These calculations are identical
	  // to the pre-calculations used to do
	  // z-buffering in evaluate_zbuffer, above!
	  basic_raycalcs(surfacearray,
			 surfaceid,
			 cam_mtx_Ut,
			 cam_mtx_Vtsubinv,
			 s1,s2,V31,V32,
			 xcnt,ycnt,
			 /* outputs ... */
			 &dc1_du0, /* also a flag for whether we need derivs */
			 &dc2_du0,
			 &dc1_dv0,
			 &dc2_dv0,
			 rayvecobj, /* 4-vector */
			 &rayvecfactor);

	  
	  // Same for derivatives...

	  deriv_raycalcs(surfacearray,
			 focalpointobj,
			 rayvecobj,rayvecfactor,
			 surfaceid,
			 facetid,
			 dc1_du0,dc2_du0,dc1_dv0,dc2_dv0,
			 /* outputs */
			 horiz_ray_intersect_shift_deriv,
			 vert_ray_intersect_shift_deriv);



	  // So horiz_ray_intersect_shift_deriv
	  // and vert_ray_intersect_shift_deriv
	  // are vectors indicating what happens when you move
	  // one pixel to the right and one pixel down
	  // on the image being projected
	  //
	  // Now we have to figure out where this data
	  // actually goes... i.e. map onto texture coordinates
	  //

	  // // Find intersection point, as in evaluate_zbuffer(), above
	  // firstidx=surfacearray[surfaceid].vertexidx_indices[facetid];
	  //
	  //dist = ray_to_plane_distance(focalpointobj,rayvecobj,&surfacearray[surfaceid].vertices[3*surfacearray[surfaceid].vertexidx[firstidx + 0]], &surfacearray[surfaceid].facetnormals[facetid*3]);

	  //assert(dist==imagedata_zbuffer[src_nx*ycnt + xcnt]); // Should get exactly same result as previous calculation
	  //angle_of_incidence = acos(fabs(dotvecvec3(rayvecobj,&surfacearray[surfaceid].facetnormals[facetid*3]))); // rayvecobj and facetnormals should be unit vectors 
	  dist=imagedata_zbuffer[src_nx*ycnt + xcnt];
	                      
	  // be sure to allocate 3 numbers of space
	  // for intersectpoint2uvcanon & friends as we will use it
	  // in projective form

	  // find intersectpoint
	  
	  find_intersect_uv(surfacearray,
			    focalpointobj,
			    rayvecobj,
			    surfaceid,facetid,
			    dist,
			    horiz_ray_intersect_shift_deriv, // only needed if calculating derivatives
			    vert_ray_intersect_shift_deriv, // only needed if calculating derivatives
			    /* outputs */
			    intersectpoint2uvcanon, // intersection point in texture coordinates (not scaled UV coordinates
			    intersectpoint2uvcanon_deriv_horiz, // also flag for whether to calc derivatives
			    intersectpoint2uvcanon_deriv_vert);
	  

	  // Now we have uv parameterization coordinates
	  // intersectpoint2uvcanon as well as derivatives
	  // of those coordinates representing motion of
	  // one pixel in the projected image to the right (horiz)
	  // and one pixel down (vert).
	  
          
	  // !!!*** Evaluate intrinsic texture coords -> canonical
	  // texture coords  transform and derivatives, here
	  // once implemented ***!!!
          
	  // Evaluate within 1.5px by bandlimited splat 
          
	  projecttoimgbuf(framedata[src_nx*ycnt+xcnt],imagedata_weighting[src_nx*ycnt + xcnt],intersectpoint2uvcanon,intersectpoint2uvcanon_deriv_horiz,intersectpoint2uvcanon_deriv_vert,surfacearray[surfaceid].imgbuf,surfacearray[surfaceid].weightingbuf,surfacearray[surfaceid].validitybuf,framecnt,surfacearray[surfaceid].nx,surfacearray[surfaceid].ny,MIN_RADIUS_UV_PIXELS,MIN_RADIUS_SRC_PIXELS,BANDWIDTH_FRACTION);

	  // ... Check if answer variability is due to floating point mode differences (nope!)
	  //femode_t femode={0};
	  //fegetmode(&femode);
	  //printf("__control_word=0x%4.4x; __glibc_reserved=0x%4.4x; __mxcsr=0x%8.8x\n",(unsigned)femode.__control_word,(unsigned)femode.__glibc_reserved,(unsigned)femode.__mxcsr);
          
	  // Projecttoimgbuf
	  // accumulates not just the patches, but a validity
	  // factor representing how much input there was in
	  // each destination pixel. The accumulated image
	  // can then be normalized by the validity factor
	  // to get a meaningful estimate (and pixels with
	  // validity factor << 1 should be set to NaN) 
	  

	  
	} 
      }
      //omp_unset_lock(&mutex);

    }
    //#pragma omp barrier

    /* Now we can use the validity map to normalize the 
       generated image */
    /* Actually, no we can't because we might have several
       different image sources being collected together, and
       we only want to normalize at the very end. 
       So it wants to be a different function */
    /*    
#pragma omp for
    for (ycnt=0;ycnt < src_ny;ycnt++) {
      size_t xcnt;
      
      for (xcnt=0; xcnt < src_nx; xcnt++) {
	
      }
      }*/
  }
}





static IPOPS_INLINE void evaluate_zbuffer_orthographic(struct projectionsurface *surfacearray,
					  size_t src_nx,
					  size_t src_ny,
					  size_t numsurfaces,
					  float64_t *orthographic_proj_matrix,
					  float32_t *imagedata_zbuffer,
					  uint32_t *imagedata_surfaceid,
					  uint32_t *imagedata_facetid,
					  float32_t *imagedata_angleofincidence,
					  float32_t *imagedata_angleofincidence_weighting,
					  float32_t *imagedata_uvcoords /* in texture coordinates [0,1], not physical scaled coordinates */)
{
  int64_t surfaceycnt;  // must be signed (e.g. not size_t) for MSVC compatibility
  /* orthographic_proj_matrix is defined as follows
   [ u ] = [  fx  0   cu ][ x ]
   [ v ] = [  0   fy  cv ][ y ]
                          [ 1 ]
   where x, y relative to "camera location" (x,y,z) = (0,0,0). Camera looking in +z direction

   where u,v image pixel coordinates 
         fx = pixels/mm in x 
         fy = pixels/mm in y
         cu = center u (in pixels)
         cv = center v (in pixels) 
	 ... z is NOT PASSED to this matrix 


   ... therefore 
   u = fx x + cu
   v = fy y + cv
   u - cu = fx x
   v - cv = fy y 
   x = ( u - cu )/fx
   y = ( v - cv )/fy


   Also, for (u,v) image:
   Step1=1.0/fx
   IniVal1=-cu/fx
   Step2=1.0/fy
   IniVal2=-cv/fy
   or equivalently: 
     fx=1.0/Step1
     cu=-IniVal1*fx
     fy=1.0/Step2
     cv=-IniVal2*fy
*/

  
#pragma omp parallel default(shared) private(surfaceycnt) //,src_ny,numsurfaces)
  {
#pragma omp for
    for (surfaceycnt=0;surfaceycnt < src_ny*numsurfaces;surfaceycnt++) {
      
      size_t ycnt,surfacecnt,xcnt;

      float64_t rayvec[4];
      float64_t rayvecobj[4]; // rayvec in object coordinates
      float64_t raysrcloc[4],raysrclocobj[4];


      int trace=FALSE;
      size_t firstidx;
      //float64_t rayvecfactor;
      
      ycnt=surfaceycnt % src_ny; //  numsurfaces;
      surfacecnt=surfaceycnt / src_ny; // numsurfaces;
      
      
      
      for (xcnt=0; xcnt < src_nx; xcnt++) {
	// Go through each source image pixel,
	// project it through to surface. If it
	// intersects closer, mark the z-buffer
	// and surface ID. This will give us a map
	// of the source image for each pixel, which surface
	// it maps onto.
	
	// NOTE: These are recalculated later in project_image... formulas should be THE SAME!!!


	// rayvec in camera coordinates
	rayvec[0]=0;
	rayvec[1]=0;
	rayvec[2]=1.0;
	rayvec[3]=0.0; // "vector" has 0 w-component

	/* apply correction for camera matrix being in OpenCV coords, but
	   our transforms being in OpenGL coords */
	/* We negate rayvec[1] and rayvec[2] rather than messing with tosurfaceframe. */
	
	rayvec[1]=-rayvec[1];
	rayvec[2]=-rayvec[2];

	// convert to object coordinates
	multmatvec4(surfacearray[surfacecnt].tosurfaceframe,rayvec,rayvecobj);
	// no need to normalize since a unit vector to begin with (rayvecfactor=1)
	

	// calculation of raysrclocobj: Coordinates of ray source in object coordinates: 
	// since we are doing this OpenCV style, negate the 2nd and 3rd columns... see imageprojectionmodel.py
	// (except we don't actually do this (???))
	// ray source location in camera coords 
	raysrcloc[0]=(xcnt - orthographic_proj_matrix[2])/orthographic_proj_matrix[0];  // (u-cu)/fx
	raysrcloc[1]=-(ycnt - orthographic_proj_matrix[5])/orthographic_proj_matrix[4];  // (v-cv)/fy (negate because tosurfaceframe is in opengl coords)
	raysrcloc[2]=-0.0; // negate because tosurfaceframe is in OpenGL coords
	raysrcloc[3]=1.0; /* this is a point */

	multmatvec4(surfacearray[surfacecnt].tosurfaceframe,raysrcloc,raysrclocobj);
	normalize_wcoord4(raysrclocobj);

	if (FALSE && ycnt==237 && xcnt==352) {
	  int polynum;
	  double dist;
	  int hitcnt=0;
	  struct projectionsurface *surf;

	  /* debugging */
	  fprintf(stderr,"y=%d/%d, x=%d/%d\n",(int)ycnt,(int)src_ny,(int)xcnt,(int)src_nx);
	  trace=TRUE;

	  surf=&surfacearray[surfacecnt];
	  
	  for (polynum = 0; polynum < surf->npolys;polynum++) {
	    firstidx=surf->vertexidx_indices[polynum];
	    dist = ray_to_plane_distance(raysrclocobj,rayvecobj, &surf->vertices[3*surf->vertexidx[firstidx + 0]], &surf->facetnormals[polynum*3]);
	    
	    if (ray_intersects_polygon(surf->vertices,&surf->vertexidx[firstidx],surf->numvertices[polynum],surf->facetnormals+3*polynum,surf->refpoints+3*polynum,surf->inplanemats+6*polynum,raysrclocobj,rayvecobj,surf->maxradius[polynum],dist,TRUE,trace && FALSE)) {
	      fprintf(stderr,"ray intersects polygon %d; dist=%f\n",(int)polynum,dist);
	      hitcnt++;
	    }
	  }
	  if (!hitcnt) {
	    fprintf(stderr,"ray did not intersect any polygons\n");
	  }
	} else {
	  trace=FALSE;
	}
	
	// find_ray_intersections() should fill out the z-buffer, surface id, and facet id arrays for the closest intersection 
	find_ray_intersections(&surfacearray[surfacecnt],surfacecnt,0,raysrclocobj,rayvecobj,imagedata_zbuffer + src_nx*ycnt + xcnt,imagedata_surfaceid + src_nx*ycnt + xcnt,imagedata_facetid + src_nx*ycnt + xcnt,trace);
	
      }
    }
  }

  /* Next phase: evaluate derivatives, etc. */
#pragma omp parallel default(shared) // private(src_ny)
  {
    int64_t ycnt; // must be signed (e.g. not size_t) for MSVC compatibility
    static const uint32_t NaNconst=0x7fc00000;

    float32_t NaNval;
    memcpy(&NaNval,&NaNconst,sizeof(NaNval));

#pragma omp for private(ycnt) //src_ny
    for (ycnt=0;ycnt < src_ny;ycnt++) {
      int64_t xcnt; // must be signed (e.g. not size_t) for MSVC compatibility

      float64_t rayvec[4];
      float64_t rayvecobj[4]; // rayvec in object coordinates
      float64_t raysrcloc[4],raysrclocobj[4];


      float64_t dist;

      float64_t horiz_raysrc_shift[3];
      float64_t vert_raysrc_shift[3];
      float64_t horiz_ray_intersect_shift_deriv[3];
      float64_t vert_ray_intersect_shift_deriv[3];


      float64_t intersectpoint2uvcanon[3]; // intersection point in texture coordinates (not scaled UV coordinates
      float64_t intersectpoint2uvcanon_deriv_horiz[3];
      float64_t intersectpoint2uvcanon_deriv_vert[3];

      float64_t angle_of_incidence,angle_of_incidence_factor;
      

      uint32_t surfaceid,facetid;



      for (xcnt=0; xcnt < src_nx; xcnt++) {
	// Go through each source image pixel, find derivatives and where it projects in (u,v) space
	
	// Dist was calculated by find_ray_intersections(), above
	dist = imagedata_zbuffer[src_nx*ycnt + xcnt];

	if (!isinf(dist)) { /* if we actually got a ray intersection above */
	  surfaceid=imagedata_surfaceid[ycnt*src_nx + xcnt];
	  facetid=imagedata_facetid[ycnt*src_nx + xcnt];


	  // rayvec in camera coordinates
	  rayvec[0]=0;
	  rayvec[1]=0;
	  rayvec[2]=1.0;
	  rayvec[3]=0.0; // "vector" has 0 w-component
	  
	  /* apply correction for camera matrix being in OpenCV coords, but
	     our transforms being in OpenGL coords */
	/* We negate rayvec[1] and rayvec[2] rather than messing with tosurfaceframe. */
	  
	  rayvec[1]=-rayvec[1];
	  rayvec[2]=-rayvec[2];
	  
	  // convert to object coordinates
	  multmatvec4(surfacearray[surfaceid].tosurfaceframe,rayvec,rayvecobj);
	  // no need to normalize since a unit vector to begin with (rayvecfactor=1)
	  
	  
	  // calculation of raysrclocobj: Coordinates of ray source in object coordinates: 
	  // since we are doing this OpenCV style, negate the 2nd and 3rd columns... see imageprojectionmodel.py
	  // (except we don't actually do this (???))
	  // ray source location in camera coords 
	  raysrcloc[0]=(xcnt - orthographic_proj_matrix[2])/orthographic_proj_matrix[0];  // (u-cu)/fx
	  raysrcloc[1]=-(ycnt - orthographic_proj_matrix[5])/orthographic_proj_matrix[4];  // (v-cv)/fy (negate because tosurfaceframe is in opengl coords)
	  raysrcloc[2]=-0.0; // negate because tosurfaceframe is in OpenGL coords
	  raysrcloc[3]=1.0; /* this is a point */
	  
	  multmatvec4(surfacearray[surfaceid].tosurfaceframe,raysrcloc,raysrclocobj);
	  normalize_wcoord4(raysrclocobj);
	  

	  // horiz_ray_intersect_shift_deriv it the change of the intersection coordinates per one-pixel
	  // horizontal shift
	  // The ray itself stays in the same direction but shifts laterally in camera coordinates by [ 1/fx, 0, 0 ]

	  // Evaluate derivatives... We convert them to object coordinates, evaluate the shifts
	  // in object coordinates, then shift back.
	  
	  // (it might be marginally more efficient to evaluate the derivatives in camera coordinates) 

	  //horiz_raysrc_shift = surfacearray[surfaceid].tosurfaceframe * [ 1/fx ; 0 ; 0 ; 0]
	  // so extract 1st column of tosurfaceframe and divide by fx
	  horiz_raysrc_shift[0]=surfacearray[surfaceid].tosurfaceframe[0*4 + 0]/orthographic_proj_matrix[0]; // /fx
	  horiz_raysrc_shift[1]=surfacearray[surfaceid].tosurfaceframe[1*4 + 0]/orthographic_proj_matrix[0]; // /fx
	  horiz_raysrc_shift[2]=surfacearray[surfaceid].tosurfaceframe[2*4 + 0]/orthographic_proj_matrix[0]; // /fx

	  // vertical shift
	  // The ray itself stays in the same direction but shifts laterally in camera coordinates by [ 0, 1/fy, 0 ]
	  // because tosurfaceframe in opengl coordinates negate 2nd and 3rd column (this is 2nd)
	  vert_raysrc_shift[0]=-surfacearray[surfaceid].tosurfaceframe[0*4 + 1]/orthographic_proj_matrix[4]; // /fy
	  vert_raysrc_shift[1]=-surfacearray[surfaceid].tosurfaceframe[1*4 + 1]/orthographic_proj_matrix[4]; // /fy 
	  vert_raysrc_shift[2]=-surfacearray[surfaceid].tosurfaceframe[2*4 + 1]/orthographic_proj_matrix[4]; // /fy
	  
	  
	  orthographic_deriv_raycalcs(surfacearray,
				      raysrclocobj,
				      rayvecobj,
				      surfaceid,
				      facetid,
				      horiz_raysrc_shift,
				      vert_raysrc_shift,
				      /* outputs */
				      horiz_ray_intersect_shift_deriv,
				      vert_ray_intersect_shift_deriv);
	  
	  //// Determine intersect position shift in camera coordinates with respect to a unit horizontal pixel shift 
	  //multmatvec3(surfacearray[surfaceid].vecfromsurfaceframe,horiz_ray_intersect_shift_deriv,horiz_ray_intersect_shift_deriv_cam);
	  
	  /* evaluate left and right predicted z values... add z coordinate of derivative to this pixel's position */
	  //imagedata_projectzright[src_nx*ycnt + xcnt] = dist + 1.0*horiz_ray_intersect_shift_deriv_cam[2];
	  //imagedata_projectzleft[src_nx*ycnt + xcnt] = dist - 1.0*horiz_ray_intersect_shift_deriv_cam[2];
	  //imagedata_horiz_ray_intersect_shift_deriv_cam_z[src_nx*ycnt + xcnt]=1.0*horiz_ray_intersect_shift_deriv_cam[2];
	  
	  
	  
	  // Determine intersect position shift in camera coordinates with respect to a unit vertical pixel shift 
	  //multmatvec3(surfacearray[surfaceid].vecfromsurfaceframe,vert_ray_intersect_shift_deriv,vert_ray_intersect_shift_deriv_cam);
	  
	  
	  /* evaluate up and down predicted z values... add z coordinate of derivative to this pixel's position */
	  /* NOTE: We are assuming y ordering is like a raster scan (down is increasing y) */
	  //imagedata_projectzdown[src_nx*ycnt + xcnt] = dist + 1.0*vert_ray_intersect_shift_deriv_cam[2];
	  //imagedata_projectzup[src_nx*ycnt + xcnt] = dist - 1.0*vert_ray_intersect_shift_deriv_cam[2];
	  //imagedata_vert_ray_intersect_shift_deriv_cam_z[src_nx*ycnt + xcnt]=vert_ray_intersect_shift_deriv_cam[2];
	  
	  // So horiz_ray_intersect_shift_deriv
	  // and vert_ray_intersect_shift_deriv
	  // are vectors indicating what happens when you move
	  // one pixel to the right and one pixel down
	  // on the image being projected
	  

	
	
	
	/* determine angle of incidence, for weighting use */
	  angle_of_incidence = acos(fabs(dotvecvec3(rayvecobj,&surfacearray[surfaceid].facetnormals[imagedata_facetid[src_nx*ycnt + xcnt]*3]))); // rayvecobj and facetnormals should be unit vectors 
	  
	  // Ignore anything at extreme angles of incidence
	  if (angle_of_incidence > 3*M_PI/8) {
	    angle_of_incidence_factor=0.0;
	  } else {  
	    // Define factor by which we de-emphasize data at larger angles of incidence
	    angle_of_incidence_factor = cos(angle_of_incidence * (M_PI/2)/(3*M_PI/8));	    
	  }
	  if (imagedata_angleofincidence_weighting) {
	    imagedata_angleofincidence_weighting[src_nx*ycnt + xcnt] = angle_of_incidence_factor;
	  }
	  
	  if (imagedata_angleofincidence) {
	    imagedata_angleofincidence[src_nx*ycnt + xcnt] = angle_of_incidence;
	  }
	  

	  
	  
	  if (imagedata_uvcoords || surfacearray[surfaceid].angle_of_incidence_factor_buf_uv) {
	    find_intersect_uv(surfacearray,
			      raysrclocobj,
			      rayvecobj,
			      surfaceid,
			      imagedata_facetid[src_nx*ycnt + xcnt],
			      dist,
			      horiz_ray_intersect_shift_deriv, // only needed if calculating derivatives
			      vert_ray_intersect_shift_deriv, // only needed if calculating derivatives
			      /* outputs */
			      intersectpoint2uvcanon, // intersection point in texture coordinates (not scaled UV coordinates
			      intersectpoint2uvcanon_deriv_horiz, // also flag for whether to calc derivatives
			      intersectpoint2uvcanon_deriv_vert);
	    
	  
	    if (imagedata_uvcoords) {
	      imagedata_uvcoords[src_nx*ycnt*2 + xcnt*2]=intersectpoint2uvcanon[0]; // Store u coordinate
	      imagedata_uvcoords[src_nx*ycnt*2 + xcnt*2 + 1]=intersectpoint2uvcanon[1];  // store v coordinate
	    }
	    
	    
	    if (surfacearray[surfaceid].angle_of_incidence_factor_buf_uv) {
	      /* if buffer to store angle_of_incidence_factor in uv coordinate frame is provided... */
	      /* use projecttoimgbuf() to map out projection validity region */
	      /* we provide angle_of_incidence_factor as the pixel value
		 so that imgbuf ends up with the sum of projection weighting * angle of incidence factor 
		 whereas validitybuf ends up with the sum of the angle of incidence factors */
	      projecttoimgbuf(angle_of_incidence_factor,1.0,intersectpoint2uvcanon,intersectpoint2uvcanon_deriv_horiz,intersectpoint2uvcanon_deriv_vert,surfacearray[surfaceid].angle_of_incidence_factor_buf_uv,NULL,surfacearray[surfaceid].validitybuf,0,surfacearray[surfaceid].nx,surfacearray[surfaceid].ny,MIN_RADIUS_UV_PIXELS,MIN_RADIUS_SRC_PIXELS,BANDWIDTH_FRACTION);
	    }
	  }
	} else {
	  /* No intersection found at all */
	  //imagedata_horiz_ray_intersect_shift_deriv_cam_z[src_nx*ycnt + xcnt]=NaNval;
	  //imagedata_vert_ray_intersect_shift_deriv_cam_z[src_nx*ycnt + xcnt]=NaNval;
	  if (imagedata_angleofincidence_weighting) {
	    imagedata_angleofincidence_weighting[src_nx*ycnt + xcnt] = 0.0;
	  }

	  if (imagedata_angleofincidence) {
	    imagedata_angleofincidence[src_nx*ycnt + xcnt]=NaNval;
	  }
	  if (imagedata_uvcoords) {
	    imagedata_uvcoords[src_nx*ycnt*2 + xcnt*2]=NaNval;
	    imagedata_uvcoords[src_nx*ycnt*2 + xcnt*2 + 1]=NaNval;
	  }
	}
      }
      
    }
  }
  
  
}    
static IPOPS_INLINE void project_orthographic(struct projectionsurface *surfacearray,
						    float32_t *framedata,
						    size_t framecnt,
						    size_t src_nx,
						    size_t src_ny,
					            float64_t *orthographic_proj_matrix,
						    float32_t *imagedata_zbuffer,
						    uint32_t *imagedata_surfaceid,
						    uint32_t *imagedata_facetid,
						    float32_t *imagedata_weighting)
{

  int64_t ycnt;  // must be signed (e.g. not size_t) for MSVC compatibility

#pragma omp parallel default(shared) private(ycnt) //,src_ny)
  {
#pragma omp for
    for (ycnt=0;ycnt < src_ny;ycnt++) {
      size_t xcnt;
      
      for (xcnt=0; xcnt < src_nx; xcnt++) {

	if (!IPOPS_ISINF(imagedata_zbuffer[ycnt*src_nx + xcnt])) {

	  int32_t surfaceid,facetid;

	  float64_t rayvec[4],rayvecobj[4]; // rayvec in object coordinates
	  float64_t raysrcloc[4],raysrclocobj[4];
	  
	  
	  float64_t horiz_raysrc_shift[3];
	  float64_t vert_raysrc_shift[3];
	  
	  float64_t dist;
	  float64_t horiz_ray_intersect_shift_deriv[3];
	  float64_t vert_ray_intersect_shift_deriv[3];
	  float64_t intersectpoint2uvcanon[3]; // intersection point in texture coordinates (not scaled UV coordinates
	  float64_t intersectpoint2uvcanon_deriv_horiz[3];
	  float64_t intersectpoint2uvcanon_deriv_vert[3];
	  //float64_t angle_of_incidence;

	  
	  /* finite z-buffer means this point maps onto the object */
	  surfaceid=imagedata_surfaceid[ycnt*src_nx + xcnt];
	  facetid=imagedata_facetid[ycnt*src_nx + xcnt];

	  // rayvec in camera coordinates
	  rayvec[0]=0;
	  rayvec[1]=0;
	  rayvec[2]=1.0;
	  rayvec[3]=0.0; // "vector" has 0 w-component
	  
	  /* apply correction for camera matrix being in OpenCV coords, but
	     our transforms being in OpenGL coords */
	/* We negate rayvec[1] and rayvec[2] rather than messing with tosurfaceframe. */
	  
	  rayvec[1]=-rayvec[1];
	  rayvec[2]=-rayvec[2];
	  
	  // convert to object coordinates
	  multmatvec4(surfacearray[surfaceid].tosurfaceframe,rayvec,rayvecobj);
	  // no need to normalize since a unit vector to begin with (rayvecfactor=1)
	  
	  
	  // calculation of raysrclocobj: Coordinates of ray source in object coordinates: 
	  // since we are doing this OpenCV style, negate the 2nd and 3rd columns... see imageprojectionmodel.py
	  // (except we don't actually do this (???))
	  // ray source location in camera coords 
	  raysrcloc[0]=(xcnt - orthographic_proj_matrix[2])/orthographic_proj_matrix[0];  // (u-cu)/fx
	  raysrcloc[1]=-(ycnt - orthographic_proj_matrix[5])/orthographic_proj_matrix[4];  // (v-cv)/fy (negate because tosurfaceframe is in opengl coords)
	  raysrcloc[2]=-0.0; // negate because tosurfaceframe is in OpenGL coords
	  raysrcloc[3]=1.0; /* this is a point */
	  
	  multmatvec4(surfacearray[surfaceid].tosurfaceframe,raysrcloc,raysrclocobj);
	  normalize_wcoord4(raysrclocobj);


	  // horiz_ray_intersect_shift_deriv it the change of the intersection coordinates per one-pixel
	  // horizontal shift
	  // The ray itself stays in the same direction but shifts laterally in camera coordinates by [ 1/fx, 0, 0 ]

	  // Evaluate derivatives... We convert them to object coordinates, evaluate the shifts
	  // in object coordinates, then shift back.
	  
	  // (it might be marginally more efficient to evaluate the derivatives in camera coordinates) 

	  //horiz_raysrc_shift = surfacearray[surfaceid].tosurfaceframe * [ 1/fx ; 0 ; 0 ; 0]
	  // so extract 1st column of tosurfaceframe and divide by fx
	  horiz_raysrc_shift[0]=surfacearray[surfaceid].tosurfaceframe[0*4 + 0]/orthographic_proj_matrix[0]; // /fx
	  horiz_raysrc_shift[1]=surfacearray[surfaceid].tosurfaceframe[1*4 + 0]/orthographic_proj_matrix[0]; // /fx
	  horiz_raysrc_shift[2]=surfacearray[surfaceid].tosurfaceframe[2*4 + 0]/orthographic_proj_matrix[0]; // /fx

	  // vertical shift
	  // The ray itself stays in the same direction but shifts laterally in camera coordinates by [ 0, 1/fy, 0 ]
	  // because tosurfaceframe in opengl coordinates negate 2nd and 3rd column (this is 2nd)
	  vert_raysrc_shift[0]=-surfacearray[surfaceid].tosurfaceframe[0*4 + 1]/orthographic_proj_matrix[4]; // /fy
	  vert_raysrc_shift[1]=-surfacearray[surfaceid].tosurfaceframe[1*4 + 1]/orthographic_proj_matrix[4]; // /fy 
	  vert_raysrc_shift[2]=-surfacearray[surfaceid].tosurfaceframe[2*4 + 1]/orthographic_proj_matrix[4]; // /fy

	  
	  // Then we go back through each pixel and
	  // perform the mapping, stamping down sinc-interpolation
	  // (or similar lowpass interpolation) 
	  // patches onto the parameterization output image.

	  // NOTE: These calculations are identical
	  // to the pre-calculations used to do
	  // z-buffering in evaluate_zbuffer, above!

	  
	  // Same for derivatives...

	  orthographic_deriv_raycalcs(surfacearray,
				      raysrclocobj,
				      rayvecobj,
				      surfaceid,
				      facetid,
				      horiz_raysrc_shift,
				      vert_raysrc_shift,
				      /* outputs */
				      horiz_ray_intersect_shift_deriv,
				      vert_ray_intersect_shift_deriv);


	  // So horiz_ray_intersect_shift_deriv
	  // and vert_ray_intersect_shift_deriv
	  // are vectors indicating what happens when you move
	  // one pixel to the right and one pixel down
	  // on the image being projected
	  //
	  // Now we have to figure out where this data
	  // actually goes... i.e. map onto texture coordinates
	  //

	  // // Find intersection point, as in evaluate_zbuffer(), above
	  // firstidx=surfacearray[surfaceid].vertexidx_indices[facetid];
	  //
	  //dist = ray_to_plane_distance(focalpointobj,rayvecobj,&surfacearray[surfaceid].vertices[3*surfacearray[surfaceid].vertexidx[firstidx + 0]], &surfacearray[surfaceid].facetnormals[facetid*3]);

	  //assert(dist==imagedata_zbuffer[src_nx*ycnt + xcnt]); // Should get exactly same result as previous calculation
	  //angle_of_incidence = acos(fabs(dotvecvec3(rayvecobj,&surfacearray[surfaceid].facetnormals[facetid*3]))); // rayvecobj and facetnormals should be unit vectors 
	  dist=imagedata_zbuffer[src_nx*ycnt + xcnt];
	                      
	  // be sure to allocate 3 numbers of space
	  // for intersectpoint2uvcanon & friends as we will use it
	  // in projective form

	  // find intersectpoint
	  
	  find_intersect_uv(surfacearray,
			    raysrclocobj,
			    rayvecobj,
			    surfaceid,facetid,
			    dist,
			    horiz_ray_intersect_shift_deriv, // only needed if calculating derivatives
			    vert_ray_intersect_shift_deriv, // only needed if calculating derivatives
			    /* outputs */
			    intersectpoint2uvcanon, // intersection point in texture coordinates (not scaled UV coordinates
			    intersectpoint2uvcanon_deriv_horiz, // also flag for whether to calc derivatives
			    intersectpoint2uvcanon_deriv_vert);
	  

	  // Now we have uv parameterization coordinates
	  // intersectpoint2uvcanon as well as derivatives
	  // of those coordinates representing motion of
	  // one pixel in the projected image to the right (horiz)
	  // and one pixel down (vert).
	  
          
	  // !!!*** Evaluate intrinsic texture coords -> canonical
	  // texture coords  transform and derivatives, here
	  // once implemented ***!!!
          
	  // Evaluate within 1.5px by bandlimited splat 
          
	  projecttoimgbuf(framedata[src_nx*ycnt+xcnt],imagedata_weighting[src_nx*ycnt + xcnt],intersectpoint2uvcanon,intersectpoint2uvcanon_deriv_horiz,intersectpoint2uvcanon_deriv_vert,surfacearray[surfaceid].imgbuf,surfacearray[surfaceid].weightingbuf,surfacearray[surfaceid].validitybuf,framecnt,surfacearray[surfaceid].nx,surfacearray[surfaceid].ny,MIN_RADIUS_UV_PIXELS,MIN_RADIUS_SRC_PIXELS,BANDWIDTH_FRACTION);
	  
          
	  // Projecttoimgbuf
	  // accumulates not just the patches, but a validity
	  // factor representing how much input there was in
	  // each destination pixel. The accumulated image
	  // can then be normalized by the validity factor
	  // to get a meaningful estimate (and pixels with
	  // validity factor << 1 should be set to NaN) 
	  

	  
	} 
      }
    }
    //#pragma omp barrier

    /* Now we can use the validity map to normalize the 
       generated image */
    /* Actually, no we can't because we might have several
       different image sources being collected together, and
       we only want to normalize at the very end. 
       So it wants to be a different function */
    /*    
#pragma omp for
    for (ycnt=0;ycnt < src_ny;ycnt++) {
      size_t xcnt;
      
      for (xcnt=0; xcnt < src_nx; xcnt++) {
	
      }
      }*/
  }
}
