#ifdef _MSC_VER
#define GAOPS_INLINE  __inline
#define GAOPS_ALLOCA _alloca
#else
#define GAOPS_INLINE  inline
#define GAOPS_ALLOCA alloca
#endif

static GAOPS_INLINE int point_in_polygon_2d_c(float64_t *vertices_rel_point,size_t numvertices)
{
      
  //# Apply winding number algorithm.
  //# This algorithm is selected -- in its most simple form --
  //# because it is so  simple and robust in the case of the
  //# intersect point being on or near the edge. It may well
  //# be much slower than optimal. It tries to return True
  //# in the edge case. 
    
  //# Should probably implement a faster algorithm then drop
  //# down to this for the special cases.

  //# See Hormann and Agathos, The point in polygon problem
  //# for arbitrary polygons, Computational Geometry 20(3) 131-144 (2001)
  //# http://dx.doi.org/10.1016/S0925-7721(01)00012-8
  //# https://pdfs.semanticscholar.org/e90b/d8865ddb7c7af2b159d413115050d8e5d297.pdf
    
  //# Winding number is sum over segments of
  //# acos((point_to_vertex1 dot point_to_vertex2)/(magn(point_to_vertex1)*magn(point_to_vertex_2))) * sign(det([ point_to_vertex1  point_to_vertex2 ]))
  //# where sign(det) is really: What is the sign of the z
  //# component of (point_to_vertex1 cross point_to_vertex2)
        
  //# Special cases: magn(point_to_vertex1)==0 or
  //#  magn_point_to_vertex2   -> point is on edge
  //# det([ point_to_vertex1  point_to_vertex2 ]) = 0 -> point may be on edge
    
  double windingnum=0.0;
  size_t VertexCnt;
  size_t NextVertex;
  float64_t magn1,magn2;
  float64_t vec1[2],vec2[2];
  float64_t det;
  float64_t cosparam;
  
  for (VertexCnt=0;VertexCnt < numvertices;VertexCnt++) {
    NextVertex=(VertexCnt+1) % numvertices;
    
    // calculate (thisvertex - ourpoint) -> vec1
    //    vec1=vertices_rel_point[VertexCnt,:]
    magn1=normvec2(&vertices_rel_point[2*VertexCnt]);
    
    
        
    // calculate (nextvertex - ourpoint) -> vec2
    //    vec2=vertices_rel_point[NextVertex,:]
    magn2=normvec2(&vertices_rel_point[2*NextVertex]);
    
    if (magn1==0.0 || magn2==0.0){
      // Got it!!!
      return TRUE;
    }
    scalevec2(1.0/magn1,&vertices_rel_point[2*VertexCnt],vec1);
    scalevec2(1.0/magn2,&vertices_rel_point[2*NextVertex],vec2);

    det=vec1[0]*vec2[1]-vec2[0]*vec1[1]; // matrix determinant
        
    cosparam=(vec1[0]*vec2[0]+vec1[1]*vec2[1]); //  /(magn1*magn2);
        
    if (cosparam < -1.0) {
      // Shouldn't be possible...just in case of weird roundoff
      cosparam=-1.0;
    }
        
    if (cosparam > 1.0) {
      // Shouldn't be possible...just in case of weird roundoff
      cosparam=1.0;
    }
    
    if (det > 0) {
      windingnum += acos(cosparam);
    } else if (det < 0) {
      windingnum -= acos(cosparam);
    } else {
      // det==0.0 
      
      // Vectors parallel or anti-parallel 
      
      if (cosparam > 0.9) {
	// Vectors parallel. We are OUTSIDE. Do Nothing
      }
      else if (cosparam < -0.9) {
	// Vectors anti-parallel. We are ON EDGE */
	return TRUE;
      }
      else { 
	assert(0); //# Should only be able to get cosparam = +/- 1.0 if abs(det) > 0.0 */
      }
    }
  }
  
  
  windingnum=fabs(windingnum)*(1.0/(2.0*M_PI)); // divide out radians to number of winds; don't care about clockwise vs. ccw
  if (windingnum > .999 && windingnum < 1.001) {
    // Almost exactly one loop... got it! 
    return TRUE;
  } else if (windingnum >= .001) {
    fprintf(stderr,"spatialnde.geometry.point_in_polygon_2d() Got weird winding number of %le; assuming inaccurate calculation on polygon edge\n",(double)windingnum);
    // Could also be self intersecting polygon 
    // got it !!! 
    return TRUE;
  }  
  // If we got this far, the search failed 
  
  return FALSE;
}

static GAOPS_INLINE int point_in_polygon_3d_c(float64_t *vertices,int32_t *optional_vertexidx, size_t nvertices,float64_t *point, float64_t *inplanemat)
{
  float64_t vert3d_rel_point[3];
  float64_t *vert2d_rel_point;
  size_t cnt,vertexidx;
  
  //vert3d_rel_point=GAOPS_ALLOCA(nvertices*3*sizeof(*vert3d_rel_point));
  vert2d_rel_point=GAOPS_ALLOCA(nvertices*2*sizeof(*vert2d_rel_point));
  
  for (cnt=0;cnt < nvertices;cnt++) {
    if (!optional_vertexidx) {
      vertexidx=cnt;
    } else {
      vertexidx=optional_vertexidx[cnt];
    }

    subvecvec3(&vertices[3*vertexidx],point,vert3d_rel_point);
    vert2d_rel_point[2*cnt+0]=dotvecvec3(vert3d_rel_point,&inplanemat[3*0]);
    vert2d_rel_point[2*cnt+1]=dotvecvec3(vert3d_rel_point,&inplanemat[3*1]);
    
  }
  return point_in_polygon_2d_c(vert2d_rel_point,nvertices);
}


static GAOPS_INLINE int segment_intersects_box_c(float64_t *box_v0,float64_t *box_v1,float64_t *seg_v0, float64_t *seg_v1)
{
  float64_t original_center[3];
  float64_t segvec[3];
  float64_t box_width[3];
  float64_t seg_axisdirections[3];
  int cnt;
  int axis;
  float64_t axisvec[3];
  float64_t surf_normal[3];
  float64_t sn_sign[3];
  float64_t vert0_minus_center[3];
  float64_t directed_box_width[3];

  mean2vec3(box_v0,box_v1,original_center);

  subvecvec3(seg_v1,seg_v0,segvec);
  subvecvec3(box_v1,box_v0,box_width);

  sign_nonzero3(segvec,seg_axisdirections);


  for (cnt=0;cnt < 3;cnt++) {
    //Surfaces at v0 end of the slide
    if (seg_v0[cnt]*seg_axisdirections[cnt]-box_width[cnt]/2.0 > original_center[cnt]*seg_axisdirections[cnt]) {
      return FALSE;
    }

    // Surfaces at v1 end of the slide
    if (seg_v1[cnt]*seg_axisdirections[cnt] + box_width[cnt]/2 < original_center[cnt]*seg_axisdirections[cnt]) {
      return FALSE;
    }

  }

  // Remaining six faces connect the two ends

  for (axis=0;axis < 3; axis++) {
    
    //surf_normal should be normal to axis
    // and normal to segvec
    for (cnt=0;cnt < 3;cnt++) {
      if (cnt==axis) axisvec[cnt]=1;
      else axisvec[cnt]=0;
    }

    crossvecvec3(axisvec,segvec,surf_normal);
    sign_nonzero3(surf_normal,sn_sign);

    subvecvec3(seg_v0,original_center,vert0_minus_center);
    multvecvec3(box_width,sn_sign,directed_box_width);
    if (fabs(dotvecvec3(vert0_minus_center,surf_normal)) > 0.5*dotvecvec3(directed_box_width,surf_normal)) {
      return FALSE;
    }
    
  }
  return TRUE;

}

static GAOPS_INLINE int polygon_intersects_box_3d_c(float64_t *box_v0, float64_t *box_v1, float64_t *vertices, int32_t *optional_vertexidx, size_t nvertices, float64_t *inplanemat, float64_t *facetnormal)
{
  size_t startvertex,endvertex;
  float64_t diagonalvec[3];
  float64_t firstdiagonal[3];
  float64_t normalsigns[3];
  float64_t firstvertex_rel_corner[3];
  float64_t t;
  float64_t intersectioncoords[3];
  float64_t starting_corner[3];
  size_t first_vertexidx,start_vertexidx,end_vertexidx;
  int cnt;
  
  if (!optional_vertexidx) {
    first_vertexidx=0;
  } else {
    first_vertexidx=optional_vertexidx[0];
  }
  
  for (startvertex=0;startvertex < nvertices;startvertex++) {
    endvertex = (startvertex+1) % nvertices;

    if (!optional_vertexidx) {
      start_vertexidx=startvertex;
      end_vertexidx=endvertex;
    } else {
      start_vertexidx=optional_vertexidx[startvertex];
      end_vertexidx=optional_vertexidx[endvertex];
    }
    
    if (segment_intersects_box_c(box_v0,box_v1,&vertices[start_vertexidx*3],&vertices[end_vertexidx*3])) {
      return TRUE;
    }
    
  }

  subvecvec3(box_v1,box_v0,firstdiagonal);
  
  sign_nonzero3(facetnormal,normalsigns);
  multvecvec3(normalsigns,firstdiagonal,diagonalvec);

  for (cnt=0;cnt < 3;cnt++) {
    if (normalsigns[cnt] >= 0) {
      starting_corner[cnt]=box_v0[cnt];
    } else {
      starting_corner[cnt]=box_v1[cnt];

    }
  }
  

  subvecvec3(&vertices[first_vertexidx*3],starting_corner,firstvertex_rel_corner);

  t=dotvecvec3(firstvertex_rel_corner,facetnormal)/dotvecvec3(diagonalvec,facetnormal);

  if (t > 1.0 || t < 0.0) {
    return FALSE;
  }

  addvecscaledvec3(starting_corner,t,diagonalvec,intersectioncoords);
 
  

  return point_in_polygon_3d_c(vertices,optional_vertexidx,nvertices,intersectioncoords,inplanemat);
}

static GAOPS_INLINE int CCW(float64_t *a,float64_t *b,float64_t *c)
{
  // 2D in-plane: are a,b,c in CCW order?
  float64_t b_minus_a[2],c_minus_a[2];
  subvecvec2(b,a,b_minus_a);
  subvecvec2(c,a,c_minus_a);
  
  return (crossvecvec2(b_minus_a,c_minus_a) > 0.0);
}

static GAOPS_INLINE int check_lineseg_intersection(float64_t *a,float64_t *b,float64_t *c,float64_t *d)
{
  // Check if the planar (2D) line segments ab and cd intersect
  // Per http://jeffe.cs.illinois.edu/teaching/373/notes/x06-sweepline.pdf
  // The segments intersect if and only if
  //     endpoints a and b are on opposite sides of cd
  //  and endpoints c and d are on opposite sides of ab
  //
  // a and b are on opposite sides of c and d if and only if 
  // exactly one of (a,c,d) and (b,c,d) are CCW
  //
  //  This is true if the counter clock wise ordering test passes
  //  (CCW(a,c,d) != CCW(b,c,d)) and (CCW(a,b,c) != CCW(a,b,d))

  // ... What if (a,b) and (c,d) are colinear? 
  // ... Then the cross product is 0 and CCW returns False
  // and they count as not intersecting
  // ... That is probably reasonable in this application
  
  return (CCW(a,c,d) != CCW(b,c,d)) && (CCW(a,b,c) != CCW(a,b,d));
}

static GAOPS_INLINE int polygon_intersects_box_2d_c(float64_t *box_v0,float64_t *box_v1,float64_t *vertices,int32_t *optional_vertexidx,size_t numvertices)
{
  float64_t box_v0b[2];
  float64_t box_v1b[2];
  size_t vertexcnt,nextvertex;
  uint32_t cur_vertexidx;
  uint32_t next_vertexidx;

  box_v0b[0]=box_v0[0];
  box_v0b[1]=box_v1[1];

  box_v1b[0]=box_v1[0];
  box_v1b[1]=box_v0[1];

  for (vertexcnt=0;vertexcnt < numvertices;vertexcnt++) {
    nextvertex=(vertexcnt+1) % numvertices;

    if (!optional_vertexidx) {
      cur_vertexidx=vertexcnt;
      next_vertexidx=nextvertex;
    } else {
      cur_vertexidx=optional_vertexidx[vertexcnt];
      next_vertexidx=optional_vertexidx[nextvertex];
    }

    
    // Compare with all of the edges of the box
    if (check_lineseg_intersection(box_v0,box_v1b,&vertices[cur_vertexidx*2],&vertices[next_vertexidx*2]) ||
	check_lineseg_intersection(box_v1b,box_v1,&vertices[cur_vertexidx*2],&vertices[next_vertexidx*2]) ||
	check_lineseg_intersection(box_v1,box_v0b,&vertices[cur_vertexidx*2],&vertices[next_vertexidx*2]) ||
	check_lineseg_intersection(box_v0b,box_v0,&vertices[cur_vertexidx*2],&vertices[next_vertexidx*2])) {
      // Found intersection of this polygon with this box
      return TRUE;
    }
  }   
  return FALSE;
  
}

static GAOPS_INLINE int box_inside_polygon_2d_c(float64_t *box_v0,float64_t *box_v1,float64_t *vertices,int32_t *optional_vertexidx,size_t num_vertices)
{
  size_t cnt;
  int outside_left_x=TRUE;
  int outside_right_x=TRUE;
  int outside_left_y=TRUE;
  int outside_right_y=TRUE;
  float64_t box_v0b[2];
  float64_t box_v1b[2];
  
  float64_t *vertices_rel_v0;
  float64_t *vertices_rel_v1;
  float64_t *vertices_rel_v0b;
  float64_t *vertices_rel_v1b;

  uint32_t cur_vertexidx;

  vertices_rel_v0=GAOPS_ALLOCA(num_vertices*2*sizeof(*vertices_rel_v0));
  vertices_rel_v1=GAOPS_ALLOCA(num_vertices*2*sizeof(*vertices_rel_v1));
  vertices_rel_v0b=GAOPS_ALLOCA(num_vertices*2*sizeof(*vertices_rel_v0b));
  vertices_rel_v1b=GAOPS_ALLOCA(num_vertices*2*sizeof(*vertices_rel_v1b));
  
  
  // Short circuit in case all of the vertices are on one side of a box plane
  for (cnt=0;cnt < num_vertices;cnt++) {
    if (!optional_vertexidx) {
      cur_vertexidx=cnt;
    } else {
      cur_vertexidx=optional_vertexidx[cnt];
    }
    
    outside_left_x = outside_left_x && (vertices[2*cur_vertexidx+0] < box_v0[0]);
    outside_right_x = outside_right_x && (vertices[2*cur_vertexidx+0] > box_v1[0]);
    outside_left_y = outside_left_y && (vertices[2*cur_vertexidx+1] < box_v0[1]);
    outside_right_y = outside_right_y && (vertices[2*cur_vertexidx+1] > box_v1[1]);
  }
  
  if (outside_left_x || outside_right_x || outside_left_y ||outside_right_y) {
    return FALSE;
  }

  // find other box corners
  box_v0b[0]=box_v0[0];
  box_v0b[1]=box_v1[1];
  
  box_v1b[0]=box_v1[0];
  box_v1b[1]=box_v0[1];


  for (cnt=0;cnt < num_vertices;cnt++) {
    if (!optional_vertexidx) {
      cur_vertexidx=cnt;
    } else {
      cur_vertexidx=optional_vertexidx[cnt];
    }
    subvecvec2(&vertices[cur_vertexidx*2],box_v0,&vertices_rel_v0[cnt*2]);
    subvecvec2(&vertices[cur_vertexidx*2],box_v1,&vertices_rel_v1[cnt*2]);
    subvecvec2(&vertices[cur_vertexidx*2],box_v0b,&vertices_rel_v0b[cnt*2]);
    subvecvec2(&vertices[cur_vertexidx*2],box_v1b,&vertices_rel_v1b[cnt*2]);
  }
  
  return (point_in_polygon_2d_c(vertices_rel_v0,num_vertices) &&
	  point_in_polygon_2d_c(vertices_rel_v1,num_vertices) &&
	  point_in_polygon_2d_c(vertices_rel_v0b,num_vertices) &&
	  point_in_polygon_2d_c(vertices_rel_v1b,num_vertices));
  
}

static GAOPS_INLINE int vertices_in_box_2d(float64_t *vertices,int32_t *optional_vertexidxs,size_t numvertices,float64_t *box_v0,float64_t *box_v1)
/* v0 must have lower coordinates than v1 */
/* returns whether all vertices are inside or on the edge of the specified box */
{
  size_t vertexcnt;
  uint32_t vertexidx;
  
  for (vertexcnt=0;vertexcnt < numvertices;vertexcnt++) {
    if (optional_vertexidxs) {
      vertexidx=optional_vertexidxs[vertexcnt];
    } else {
      vertexidx=vertexcnt;
    }

    /* if this vertex is outside the box... */
    if (vertices[vertexidx*2+0] < box_v0[0] ||
	vertices[vertexidx*2+0] > box_v1[0] ||
	vertices[vertexidx*2+1] < box_v0[1] ||
	vertices[vertexidx*2+1] > box_v1[1]) {
      return FALSE;
    }
    
  }
  return TRUE;
}


static GAOPS_INLINE int vertices_in_box_3d(float64_t *vertices,int32_t *optional_vertexidxs,size_t numvertices,float64_t *box_v0,float64_t *box_v1)
/* v0 must have lower coordinates than v1 */
/* returns whether all vertices are inside or on the edge of the specified box */
{
  size_t vertexcnt;
  uint32_t vertexidx;
  
  for (vertexcnt=0;vertexcnt < numvertices;vertexcnt++) {
    if (optional_vertexidxs) {
      vertexidx=optional_vertexidxs[vertexcnt];
    } else {
      vertexidx=vertexcnt;
    }

    /* if this vertex is outside the box... */
    if (vertices[vertexidx*3+0] < box_v0[0] ||
	vertices[vertexidx*3+0] > box_v1[0] ||
	vertices[vertexidx*3+1] < box_v0[1] ||
	vertices[vertexidx*3+1] > box_v1[1] ||
	vertices[vertexidx*3+2] < box_v0[2] ||
	vertices[vertexidx*3+2] > box_v1[2]) {
      return FALSE;
    }
    
  }
  return TRUE;
}

static GAOPS_INLINE void eigvals_2d(float64_t *mat,float64_t *evals_out)
{
  float64_t trace;
  float64_t det;

  det=mat[0*2+0]*mat[1*2+1]-mat[0*2+1]*mat[1*2+0]; // matrix determinant
  trace=mat[0*2+0]+mat[1*2+1];

  evals_out[0] = (trace + sqrt(pow(trace,2.0) - 4*det))/2.0;
  evals_out[1] = (trace - sqrt(pow(trace,2.0) - 4*det))/2.0;
  
}

static GAOPS_INLINE void eigvals_2d_float(float32_t *mat,float32_t *evals_out)
{
  float32_t trace;
  float32_t det;

  det=mat[0*2+0]*mat[1*2+1]-mat[0*2+1]*mat[1*2+0]; // matrix determinant
  trace=mat[0*2+0]+mat[1*2+1];

  evals_out[0] = (trace + sqrt(pow(trace,2.0) - 4*det))/2.0;
  evals_out[1] = (trace - sqrt(pow(trace,2.0) - 4*det))/2.0;
  
}
