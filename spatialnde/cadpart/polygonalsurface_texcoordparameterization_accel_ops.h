#ifdef _MSC_VER
#define PSTCPAOPS_INLINE  __inline
#include <float.h>
#define PSTCPAOPS_ISNAN _isnan
#define PSTCPAOPS_ALLOCA _alloca
#else
#define PSTCPAOPS_INLINE  inline
#define PSTCPAOPS_ISNAN isnan
#define PSTCPAOPS_ALLOCA alloca
#endif

static PSTCPAOPS_INLINE size_t enclosed_or_intersecting_polygons_2d_c(int32_t *polypool,
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
					      uint32_t *num_fully_enclosed) // Returns # of polygons actually returned in retpolys
{
  // retpolys assumed to be at least as big as polypool
  size_t num_returned_polys=0;
  size_t poolidx;
  int32_t idx;
  size_t firstidx;
  int32_t *vertexidxs;
  size_t numvertexidxs;
  int polygon_fully_enclosed;
  uint32_t polynum;

  float64_t box_v0[2];
  float64_t box_v1[2];
  
  (*num_fully_enclosed)=0;

  box_v0[0]=minu;
  box_v0[1]=minv;
  box_v1[0]=maxu;
  box_v1[1]=maxv;

  for (poolidx=0;poolidx < polypoollen;poolidx++) {
    idx=polypool[poolidx];
    if (idx < 0) {
      // masked out polygon
      continue;
    }

    if (idx < num_3d_polygons) {
      firstidx=vertexidx_indices[idx];
      vertexidxs=&texcoordidx[firstidx];
      numvertexidxs=numvertices[idx];
      polygon_fully_enclosed = vertices_in_box_2d(texcoord,vertexidxs,numvertexidxs,box_v0,box_v1); 
    } else {
      // redundant texcoords
      firstidx=texcoordredundant_polystartindexes[idx-num_3d_polygons];
      polynum=texcoordredundant_polystartpolynum[idx-num_3d_polygons];
      vertexidxs = &texcoordredundant_texcoordidx[firstidx];
      numvertexidxs=numvertices[polynum];
      polygon_fully_enclosed = vertices_in_box_2d(texcoord,vertexidxs,numvertexidxs,box_v0,box_v1); 
      
    }

    if (polygon_fully_enclosed) {
      retpolys[num_returned_polys]=idx;
      num_returned_polys++;
      // if it's fully enclosed, nothing else need look at at, so we filter it here from the broader sibling pool
      polypool[poolidx] = -1; // mask out polygon

      (*num_fully_enclosed)++;
    } else {
      /* not polygon_fully_enclosed */

      // does it intersect?

      if (polygon_intersects_box_2d_c(box_v0,box_v1,texcoord,vertexidxs,numvertexidxs)) {
	retpolys[num_returned_polys]=idx;
	num_returned_polys++;
	//Don't filter it out in this case because it must
	// intersect with a sibiling too 

      }
      else if (box_inside_polygon_2d_c(box_v0,box_v1,texcoord,vertexidxs,numvertexidxs)) {
	// What if the box is entirely inside the polygon?
	// Then we should return it also
	retpolys[num_returned_polys]=idx;
	num_returned_polys++;
	
      }
    }
  }

  return num_returned_polys;
}


static PSTCPAOPS_INLINE int test_polynum_uv_c(uint32_t *vertexidx_indices,
				    size_t num_3d_polys,
				    uint32_t *numvertices,
				    int32_t *texcoordidx,
				    uint32_t *texcoordredundant_polystartindexes,
				    uint32_t *texcoordredundant_polystartpolynum,
				    int32_t *texcoordredundant_texcoordidx,
				    float64_t *texcoord,
				    float64_t u,
				    float64_t v,
				    size_t polynum)
{
  // Find whether our (u,v) coordinate is inside this polygon
  // ... is our (u,v) point inside this polygon?
  size_t firstidx;
  size_t numpoints;
  size_t polysurf_polynum;
  int32_t *vertexidxs;
  size_t cnt;

  
  if (polynum < num_3d_polys) {
    firstidx=vertexidx_indices[polynum];
    numpoints=numvertices[polynum];
    vertexidxs=&texcoordidx[firstidx];  // :(firstidx+numpoints)]
            
  }
  else {
    firstidx=texcoordredundant_polystartindexes[polynum-num_3d_polys];
    polysurf_polynum=texcoordredundant_polystartpolynum[polynum-num_3d_polys];
    numpoints=numvertices[polysurf_polynum];
    vertexidxs = &texcoordredundant_texcoordidx[firstidx]; //:(firstidx+polysurf.numvertices[polysurf_polynum])]
  }

  {
    float64_t *vertices_rel_point;
    vertices_rel_point=PSTCPAOPS_ALLOCA(sizeof(*vertices_rel_point)*2*numpoints);

    for (cnt=0;cnt < numpoints;cnt++) {
      vertices_rel_point[cnt*2+0] = texcoord[vertexidxs[cnt]*2+0]-u;
      vertices_rel_point[cnt*2+1] = texcoord[vertexidxs[cnt]*2+1]-v;
    }
    return point_in_polygon_2d_c(vertices_rel_point,numpoints);
    
  }
  
}


#define BOXES_COLUMNS 6

static PSTCPAOPS_INLINE int32_t identify_polynum_uv_c(uint32_t *vertexidx_indices,
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
					    size_t boxpool_mem_nelem)
					    //int32_t *polypool_mem,
					    //size_t polypool_mem_size)
{
  int point_in_poly;
  // size_t num_candidate_polys; // # used in polypool_mem
  size_t num_boxes_to_test; // # used in boxpool_mem
  // int32_t *candidatepolys;
  int32_t *boxes_to_test;
  int32_t curbox;
  int32_t *children;
  int32_t polys_in_curbox_idx;
  int32_t polys_in_curbox_num;
  
  size_t cnt;
  
  if (candidate_polynum >= 0) {
    point_in_poly=test_polynum_uv_c(vertexidx_indices,
				    num_3d_polys,
				    numvertices,
				    texcoordidx,
				    texcoordredundant_polystartindexes,
				    texcoordredundant_polystartpolynum,
				    texcoordredundant_texcoordidx,
				    texcoord,
				    u,v,candidate_polynum);
    if (point_in_poly) return candidate_polynum;
  }

  
  // Search boxes to find candidate polygons
  //candidatepolys=polypool_mem;
  //num_candidatepolys=0;

  boxes_to_test=boxpool_mem;
  boxes_to_test[0]=0; // Always test grandfather box
  num_boxes_to_test=1;
    
  while (num_boxes_to_test > 0) {
    curbox=boxes_to_test[num_boxes_to_test-1];
    num_boxes_to_test--;

    //fprintf(stderr,"curbox=%d\n",(int)curbox);
    
    if (u >= boxcoords[curbox*4+0] &&
	v >= boxcoords[curbox*4+1] &&
	u <= boxcoords[curbox*4+2] &&
	v <= boxcoords[curbox*4+3]) {
      // we are inside box
      children=&boxes[curbox*BOXES_COLUMNS]; // children is length 4
      if (children[0] >= 0) {
	// got children... add to boxes_to_test
	for (cnt=0;cnt < 4;cnt++) {
	  if (num_boxes_to_test >= boxpool_mem_nelem) {
	    // Storage overflow
	    assert(0);
	    return -1;
  	  }
	  //fprintf(stderr,"add child =%d\n",(int)children[cnt]);
	  boxes_to_test[num_boxes_to_test]=children[cnt];
	  num_boxes_to_test++;
        }
      }
      
      if (boxes[curbox*BOXES_COLUMNS+4] >= 0) {
	polys_in_curbox_idx=boxes[curbox*BOXES_COLUMNS+4];
	polys_in_curbox_num=boxes[curbox*BOXES_COLUMNS+5];

	//// Add these polygons to candidatepolys 
	//for (cnt=0;cnt < polys_in_curbox_num;cnt++) {
	//  candidatepolys[num_candidatepolys]=boxpolys[polys_in_curbox_idx+cnt];
	//  num_candidatepolys++;
	//}

	// Test these polygons
	// (if we went back to storing a list,
	// we could potentially optimize by
	// removing duplicates
	for (cnt=0;cnt < polys_in_curbox_num;cnt++) {
	  point_in_poly=test_polynum_uv_c(vertexidx_indices,
					  num_3d_polys,
					  numvertices,
					  texcoordidx,
					  texcoordredundant_polystartindexes,
					  texcoordredundant_polystartpolynum,
					  texcoordredundant_texcoordidx,
					  texcoord,
					  u,v,boxpolys[polys_in_curbox_idx+cnt]);
	  if (point_in_poly) return boxpolys[polys_in_curbox_idx+cnt];
	}
      }
    }
  }   // If we got this far, the search failed!
  return -1;

}

static PSTCPAOPS_INLINE void evaluate_curvature_c(uint32_t *vertexidx_indices,
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
			  float64_t *curvmatout)
{
  // Evaluate the curvature, within polygon # polynum
  // at (u,v) coordinates...  (u,v) in texture coordinate
  // range [0...1]
  size_t polysurf_polynum,firstidx,numpoints;
  float64_t *To2D,*AijMatInv,*centroid;
  float64_t u_v_one[3],TexUVExt[3],Tex3D[3];
  size_t cnt,rowcnt,colcnt,sumcnt,sumcnt2,vertexcnt,axiscnt,princcurvcnt;
  int axis;
  float64_t dist,distssquared,maxdistsquared,eps,sumval;
  float64_t avgcurvmatrix[2*2];
  float64_t asymmetry,totalweights,normval,sumsq,normfactor;
  
  if (polynum >= num_3d_polys) {
    // This polynum corresponds to a redundant texture
    polysurf_polynum=texcoordredundant_polystartpolynum[polynum];
  } else {
    polysurf_polynum=polynum;
  }
  
  To2D=&inplanemats[polysurf_polynum*2*3]; // To2D is 2x3

  AijMatInv=&texcoord2inplane[polynum*3*3]; //texcoord2inplane is 3x3
        
  // Note Capital UV represent the texture  parameterization
  // of the in-plane 3D space of this facet. 
  //TexUVExt = np.inner(AijMatInv,np.array((u,v,1.0)))
  u_v_one[0]=u;
  u_v_one[1]=v;
  u_v_one[2]=1.0;
  multmatvec3(AijMatInv,u_v_one,TexUVExt);
  
  //TexUVExt /= TexUVExt[2] # normalize inhomogeneous coordinates

  TexUVExt[0] /= TexUVExt[2];
  TexUVExt[1] /= TexUVExt[2];
  
  // These coordinates of this (u,v) of this facet are relative to its centroid,
  // and are in terms of the basis vectors in To2D
  //      TexUV = TexUVExt[:2]
  // TexUV equivalent in C context to TexUVExt
        
  //  Get 3D coordinates relative to centroid
  //Tex3D = np.inner(To2D.T,TexUV)
  multvecmat23(TexUVExt,To2D,Tex3D);

  // Need to evaluate 3D vertex coords, relative to centroid,
  // Use them to weight the vertex curvatures
  // according to distance from our point.
        
  centroid = &refpoints[polysurf_polynum*3]; // Centroied in 3d coords
  firstidx=vertexidx_indices[polysurf_polynum];
  numpoints=numvertices[polysurf_polynum];
        
  // Check to see if we have curvatures at all vertices:
  for (cnt=0; cnt < numpoints; cnt++) {
    if (PSTCPAOPS_ISNAN(principal_curvatures[vertexidx[firstidx+cnt]*2+0])) {
      // Curvature not given... abort and return NaN
      curvmatout[0]=curvmatout[1]=curvmatout[2]=curvmatout[3]=my_infnan(0);
      return;
    }
  }

  {
    // For this facet, the 3D coords of the vertices relative to the centroid are coordvals 
    float64_t *coordvals=PSTCPAOPS_ALLOCA(sizeof(*coordvals)*numpoints*3);
    float64_t *dists=PSTCPAOPS_ALLOCA(sizeof(*dists)*numpoints);
    float64_t *rawweights=PSTCPAOPS_ALLOCA(sizeof(*rawweights)*numpoints);
    float64_t *weights=PSTCPAOPS_ALLOCA(sizeof(*weights)*numpoints);
    //float64_t *coordvals2d=PSTCPAOPS_ALLOCA(sizeof(*coordvals2d)*numpoints*2);
    float64_t *CTA_2D=PSTCPAOPS_ALLOCA(sizeof(*CTA_2D)*2*numpoints*2);
    float64_t *curvmatrices=PSTCPAOPS_ALLOCA(sizeof(*curvmatrices)*numpoints*2*2);
    float64_t *weightedcurvmatrices=PSTCPAOPS_ALLOCA(sizeof(*weightedcurvmatrices)*numpoints*2*2);

    for (cnt=0;cnt < numpoints;cnt++) {
      for (axis=0;axis < 3; axis++) {
	coordvals[cnt*3+axis]=vertices[vertexidx[firstidx+cnt]*3+axis]-centroid[axis];
      }
    }
    // coordvals is the coordinates relative to centroid, numpoints x 3 (transpose of version in Python code
    
    // Now coordvals is numvertices x 3, coordinates of the vertices
    // relative to centroid
    // Tex3D is 3 vector, coordinates of our (u,v) location
    // relative to centroid.
    
    // Perform weighted average
    // dists = vecnorm(Tex3D.reshape(3,1) - coordvals,axis=0)
    
    maxdistsquared=0.0;
    for (cnt=0;cnt < numpoints;cnt++) {
      distssquared=0;
      for (axis=0;axis < 3; axis++) {
        dist=coordvals[cnt*3+axis]-Tex3D[axis];
	distssquared += pow(dist,2.0);
      }
      if (distssquared > maxdistsquared) {
        maxdistsquared=distssquared;
      }
      dists[cnt]=sqrt(distssquared);
    }
    //eps = np.max(dists)/10000.0 # small number, so we don't divide by 0
    eps = sqrt(maxdistsquared)/10000.0;


    // rawweights=1.0/(dists+eps)
    // totalweights=np.sum(rawweights)
    totalweights=0.0;
    for (cnt=0;cnt < numpoints;cnt++) {
      rawweights[cnt]=1.0/(dists[cnt]+eps);
      totalweights += rawweights[cnt];
    }
      
        
    // weights=rawweights/totalweights
    for (cnt=0;cnt < numpoints;cnt++) {
      weights[cnt]=rawweights[cnt]/totalweights;
    }
    // The 2D coords of the vertices are 
    //coordvals2d = np.dot(To2D,coordvals) # 2 rows by numpoints cols... in 2D basis relative to centroid   (not used!) (remember our coordvals transposed compared to Python code) 
      
      
    // Likewise 2D coords of the curvature_tangent_axes
    // CTA_2D = np.inner(To2D,polysurf.curvature_tangent_axes[polysurf.vertexidx[firstidx:(firstidx+numpoints)],:,:]).transpose(1,0,2) # Transpose to keep broadcast axis to the left. Pre-transpose axes lengths are: 2 (2D axes) by # of vertices by 2 (principal curvature)
    // *** C implementation does not do transpose
    // so axes are 2 (2d axes) by # of vertices by 2 (principal curvature)
    for (axiscnt=0;axiscnt < 2;axiscnt++) {
      for (vertexcnt=0; vertexcnt < numpoints; vertexcnt++) {
        for(princcurvcnt=0;princcurvcnt < 2; princcurvcnt++) {
          sumval=0.0;
	  for (sumcnt=0;sumcnt < 3;sumcnt++) {
            sumval+=To2D[axiscnt*3+sumcnt]*curvature_tangent_axes[vertexidx[firstidx+vertexcnt]*2*3 + princcurvcnt*3 + sumcnt];
	  }
	  CTA_2D[axiscnt*2*numpoints + vertexcnt*2 + princcurvcnt]=sumval;
        }

      }
    }
      
    // Normalize curvature_tangent_axes (should be unit length)
    //CTA_2D /= vecnormkeepshape(CTA_2D,0) # Axis is axis 0 because it came from To2D
    for (vertexcnt=0; vertexcnt < numpoints; vertexcnt++) {
      for(princcurvcnt=0;princcurvcnt < 2; princcurvcnt++) {
        sumsq=0.0;
	for (axiscnt=0;axiscnt < 2;axiscnt++) {
          sumsq += pow(CTA_2D[axiscnt*numpoints*2 + vertexcnt*2 + princcurvcnt],2);
        }
	normfactor=1.0/sqrt(sumsq);
	
	for (axiscnt=0;axiscnt < 2;axiscnt++) {
          CTA_2D[axiscnt*numpoints*2 + vertexcnt*2 + princcurvcnt] *= normfactor; 
        }
      }
    }

      
    // Construct curvature matrices ...
    // Need to construct V*K*V', broadcasting over which vertex
    // curvmatrices=np.einsum('...ij,...j,...jk->...ik', CTA_2D,polysurf.principal_curvatures[polysurf.vertexidx[firstidx:(firstidx+numpoints)],:],CTA_2D.transpose(0,2,1)) # result is # of vertices by 2x2 curvature matrix

    for (vertexcnt=0;vertexcnt < numpoints;vertexcnt++) { // vertexcnt represented by ellipsis, above
      for (axiscnt=0;axiscnt < 2; axiscnt++) { // axiscnt represented by i, above
        for(princcurvcnt=0;princcurvcnt < 2; princcurvcnt++) { //princcurvcnt represented by k, above
          sumval=0.0;
	  for (sumcnt=0;sumcnt < 2;sumcnt++) { // sumcnt represented by j, above
            sumval+=CTA_2D[axiscnt*numpoints*2 + vertexcnt*2 + sumcnt]*principal_curvatures[vertexidx[firstidx+vertexcnt]*2+sumcnt]*CTA_2D[princcurvcnt*numpoints*2 + vertexcnt*2 + sumcnt];
          }
	    
	  curvmatrices[vertexcnt*2*2 + axiscnt*2 + princcurvcnt] = sumval;
        }
      }
    }

      
    // Weighting of vertices relative to our point (u,v)
    // weightedcurvmatrices = weights.reshape(numpoints,1,1)*curvmatrices
    for (vertexcnt=0;vertexcnt < numpoints;vertexcnt++) { 
      for (axiscnt=0;axiscnt < 2; axiscnt++) { 
        for(princcurvcnt=0;princcurvcnt < 2; princcurvcnt++) { 
          weightedcurvmatrices[vertexcnt*2*2 + axiscnt*2 + princcurvcnt] = curvmatrices[vertexcnt*2*2 + axiscnt*2 + princcurvcnt] * weights[vertexcnt];
        }
      }
    }

    // avgcurvmatrix (weighted average)
    //  avgcurvmatrix = weightedcurvmatrices.sum(axis=0)
    avgcurvmatrix[0*2+0]=avgcurvmatrix[0*2+1]=avgcurvmatrix[1*2+0]=avgcurvmatrix[1*2+1]=0.0;
    for (vertexcnt=0;vertexcnt < numpoints;vertexcnt++) { 
      for (axiscnt=0;axiscnt < 2; axiscnt++) { 
        for(princcurvcnt=0;princcurvcnt < 2; princcurvcnt++) {
          avgcurvmatrix[axiscnt*2+princcurvcnt]+=weightedcurvmatrices[vertexcnt*2*2 + axiscnt*2 + princcurvcnt] ;
        }
      }
    }
    // avgcurvmatrix is a 2x2 which should be close to symmetric
    asymmetry = avgcurvmatrix[1*2+0]-avgcurvmatrix[0*2+1];

    
    //  if abs(asymmetry) > 0.1*np.linalg.norm(avgcurvmatrix):
    normval=0.0;
    for (cnt=0; cnt < 4;cnt++) {
      // Accumulate frobenius norm
      normval+=pow(avgcurvmatrix[cnt],2);
    }
    normval=sqrt(normval);

    if (fabs(asymmetry) > 0.1*normval) {		  
      fprintf(stderr,"evaluate_curvature_c(): WARNING Large asymmetry in mean curvature matrix at (u,v) = (%lg,%lg). Matrix = [ %lg %lg ; %lg %lg ]\n",u,v,avgcurvmatrix[0],avgcurvmatrix[1],avgcurvmatrix[2],avgcurvmatrix[3]);
    }
      
    // correct asymmetry
    avgcurvmatrix[1*2+0] -= asymmetry/2.0;
    avgcurvmatrix[0*2+1] += asymmetry/2.0;
        
    // avgcurvmatrix is in the polysurf.inplanemats (i.e. To2D)
    // orthonormal basis, with units of meters
        
    // We want to return 2D vectors in the AijMat i.e. self.inplane2texcoords) frame
    // with units of texcoords.
    //i.e. AijMat * avgcurv * inv(AijMat)
    // 
    // NOTE: AijMat is in inhomogeneous coordinates
    // ... avgcurvmatrix represents vectors, not coords
    // so use only first two rows
    //return np.dot(self.inplane2texcoords[polynum,:,:2],np.dot(avgcurvmatrix,self.texcoord2inplane[polynum,:2,:2]))
    // i.e ret_il = AijBjkCkl 
    for (rowcnt=0;rowcnt < 2;rowcnt++) {
      for (colcnt=0;colcnt < 2; colcnt++) {
        sumval=0.0;
 	for (sumcnt=0;sumcnt < 2;sumcnt++) {
          for (sumcnt2=0;sumcnt2 < 2;sumcnt2++) {
            sumval += inplane2texcoords[polynum*2*3+rowcnt*3+sumcnt]*avgcurvmatrix[sumcnt*2+sumcnt2]*texcoord2inplane[polynum*3*3+sumcnt2*3+colcnt];
          }
	  
        }
	curvmatout[rowcnt*2+colcnt]=sumval;
      }
    }
    
  }
}


static PSTCPAOPS_INLINE void linelength_sumcurvature_c_one(
	// This operatex on (u,v) in texcoord [0-1] values
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
        int32_t *boxpolys,
        float64_t *boxcoords,
        float64_t *inplanemats,        
        float64_t *inplane2texcoords,
        float64_t *texcoords2inplane,
        float64_t *principal_curvatures,
        float64_t *curvature_tangent_axes,
	float64_t du, // nominal resolution in u
	float64_t dv, // nominal resolution in v
	float64_t u1,
	float64_t v1,
	float64_t u2,
	float64_t v2,
	int32_t *boxpool_mem,
	size_t boxpool_mem_nelem,
	float64_t *linelengthout,
	float64_t *sumcurvatureout)
{

  // Add up the physical length of all the line segments in the straight line
  // (in u,v space) from u1,v1 to u2,v2. Assumes that the straight line is
  // the shortest path, which isn't true in general. (shortest path
  // would be the geodesic.)
  
  float64_t udist_dus,vdist_dvs,uvdist;
  size_t numsteps,stepnum;
  float64_t linevec_uv[2],lineunitvec_uv[2],stepsize_uv[2];
  float64_t segmidu,segmidv; //,stepmag_uv;
  float64_t dXdstep,dYdstep,dRdstep;
  float64_t dR,curvaccum;
  int32_t polynum=-1;
  float64_t curvmat[2*2];
  float64_t curvtempvec[2];
  
  udist_dus = (u2-u1)/du;
  vdist_dvs = (v2-v1)/dv;


  numsteps = (size_t)(1+sqrt(pow(udist_dus,2)+pow(vdist_dvs,2)));

  linevec_uv[0]=(u2-u1);
  linevec_uv[1]=(v2-v1);

  // Distance in uv space
  uvdist = normvec2(linevec_uv);
  //stepmag_uv = uvdist/(numsteps);

    //// unit vector in uv space
  // (avoid this because can be NaN if dist is zero)
  lineunitvec_uv[0]=linevec_uv[0]/uvdist;
  lineunitvec_uv[1]=linevec_uv[1]/uvdist;
  
  //// step size
  //stepsize_uv[0]=stepmag_uv*lineunitvec_uv[0];
  //stepsize_uv[1]=stepmag_uv*lineunitvec_uv[1];
  stepsize_uv[0]=linevec_uv[0]/numsteps;
  stepsize_uv[1]=linevec_uv[1]/numsteps;



  // we split line into a bunch of segments and for each segment
  // just find the midpoint of the segment
  // and use that value as representative...
  curvaccum=0.0;
  
  dR=0.0;  
  for (stepnum=0,polynum=-1;
       stepnum < numsteps;
       stepnum++) {
    //segstartu=u1 + stepsize_uv[0]*stepnum;
    //segstartv=v1 + stepsize_uv[1]*stepnum;
    //segendu=u1 + stepsize_uv[0]*(stepnum+1);
    //segendv=v1 + stepsize_uv[1]*(stepnum+1);
    segmidu = u1 + stepsize_uv[0]*(stepnum+0.5); /* segment midpoint in u */
    segmidv = v1 + stepsize_uv[1]*(stepnum+0.5); /* segment midpoint in v */
    
    polynum=identify_polynum_uv_c(vertexidx_indices,
				  num_3d_polys,
				  numvertices,
				  texcoordidx,
				  texcoordredundant_polystartindexes,
				  texcoordredundant_polystartpolynum,
				  texcoordredundant_texcoordidx,
				  texcoord,
				  boxes,
				  boxpolys,
				  boxcoords,
				  segmidu,
				  segmidv,
				  polynum,
				  boxpool_mem,
				  boxpool_mem_nelem);
    if (polynum < 0) {
      *linelengthout = my_infnan(0);
      *sumcurvatureout = my_infnan(0);
      return;
    }

    // dRdstep=np.linalg.norm(np.dot(self.texcoords2inplane[polynum,:2,:2],stepsize_uv)

    dXdstep = (texcoords2inplane[polynum*3*3 + 0*3 + 0]*stepsize_uv[0] +
	       texcoords2inplane[polynum*3*3 + 0*3 + 1]*stepsize_uv[1]);
    
    dYdstep = (texcoords2inplane[polynum*3*3 + 1*3 + 0]*stepsize_uv[0] +
	       texcoords2inplane[polynum*3*3 + 1*3 + 1]*stepsize_uv[1]);
    
    dRdstep = sqrt(pow(dXdstep,2)+pow(dYdstep,2));
    

    evaluate_curvature_c(vertexidx_indices,
			 num_3d_polys,
			 numvertices,
			 vertexidx,
			 vertices,
			 refpoints,
			 texcoordredundant_polystartpolynum,
			 inplanemats,
			 inplane2texcoords,
			 texcoords2inplane,
			 principal_curvatures,
			 curvature_tangent_axes,
			 polynum,
			 segmidu,
			 segmidv,
			 curvmat);

    // curvmat is a 2x2 in (u,v)...
    // (u,v) are not necessarily orthogonal in real space


    // Really, curvmat is a shape operator.
    // To get a curvature along a line,
    // multiply on the left and right by a tangent vector that maps
    // to 1 meter physical length. The result will be
    // the curvature along that tangent vector (in 1/m).

    // ... But since it's the same vector on the
    // left and right, the scalings will be reciprocals,
    // so we can just multiply on the left and right by unit vectors

    // This then gives us the curvature along the line segment
    // in m^-1

    // So we accumulate a weighted average of the curvatures
    // of the line segments.
    if (!PSTCPAOPS_ISNAN(lineunitvec_uv[0]) && !PSTCPAOPS_ISNAN(lineunitvec_uv[1])) {
      // (If either of these are NaN then our step length is zero so
      // we really don'nt need to accumulate)
      multmatvec2(curvmat,lineunitvec_uv,curvtempvec);

      // dot product is evaluated curvature. dRdstep is weighting factor
      curvaccum += dotvecvec2(lineunitvec_uv,curvtempvec) * dRdstep;
    }
    
    // Curvature is positive when surface is concave,
    // negative when convex... by outward directed
    // selected in curvature calculation code by normal. 
    
    dR += dRdstep;
          
    
  }

  
  *linelengthout=dR; 
  *sumcurvatureout=curvaccum; // divide out sum of weights 
}

static PSTCPAOPS_INLINE void linelength_avgcurvature_c_array(
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
	float64_t *avgcurvatureout)
{
  size_t linecnt;
  int32_t *boxpool_mem;
  float64_t sumcurvature;

  boxpool_mem=malloc(sizeof(*boxpool_mem)*nboxes);

  // this loop can be freely parallelized, but each thread
  // will need its own boxpool_mem
  for (linecnt=0;linecnt < numlines;linecnt++) {
    linelength_sumcurvature_c_one(
				   vertexidx_indices,
				   num_3d_polys,
				   numvertices,
				   vertexidx,
				   vertices,
				   refpoints,
				   texcoordidx,
				   texcoordredundant_polystartindexes,
				   texcoordredundant_polystartpolynum,
				   texcoordredundant_texcoordidx,
				   texcoord,
				   boxes,
				   boxpolys,
				   boxcoords,
				   inplanemats,        
				   inplane2texcoords,
				   texcoords2inplane,
				   principal_curvatures,
				   curvature_tangent_axes,
				   du,dv,
				   u1[linecnt],
				   v1[linecnt],
				   u2[linecnt],
				   v2[linecnt],
				   boxpool_mem,
				   nboxes,
				   &linelengthout[linecnt],
				   &sumcurvature);
    avgcurvatureout[linecnt]=sumcurvature/linelengthout[linecnt];
    
    /***!!! FIXME ... Should handle linelengthout==0 case 
	like mesh based implementation, gracefully degrading 
	to local curvature ****/

  }
  free(boxpool_mem);
}



static PSTCPAOPS_INLINE void linelength_sumcurvature_meshbased_c_one(
							   // This operatex on (u,v) in texcoord [0-1] values 
							   float32_t *curvmats, // indexed c-style (v,u,2,2)  ... must cover entire range of [0,1] texture coordinates
							   float32_t *stepsizearray, // step sizes interpreted as distance per pixel in stepsizearray
							   size_t nv, // number of v steps in curvmats or stepsizearray
							   size_t nu, // number of u steps in curvmats or stepsizearray
							   float64_t param_lowerleft_meaningfulunits_u,
							   float64_t param_lowerleft_meaningfulunits_v,
							   float64_t param_stepsize_u, /* meaningful unit stepsize for curvmats/stepsizearray (nominal value), used to measure u1,v1,etc. */
							   float64_t param_stepsize_v,
	float64_t du, // nominal desired resolution in u, meaning desired step size
	float64_t dv, // nominal desired resolution in v, desired step size, in texcoord [0-1] units
	float64_t u1,
	float64_t v1,
	float64_t u2,
	float64_t v2,
	float64_t *linelengthout,
	float64_t *sumcurvatureout,
	float64_t *sumcrosscurvatureout)
{

  // Add up the physical length of all the line segments in the straight line
  // (in u,v space) from u1,v1 to u2,v2. Assumes that the straight line is
  // the shortest path, which isn't true in general. (shortest path
  // would be the geodesic.)

  
  float64_t udist_dus,vdist_dvs,uvdist;
  size_t numsteps,stepnum;
  float64_t linevec_uv[2],lineunitvec_uv[2],stepsize_uv[2];
  float64_t segmidu,segmidv;// ,stepmag_uv;
  float64_t dRdstep;
  float64_t dR,curvaccum,crosscurvaccum;
  uint32_t uidx,vidx;
  float64_t curvtempvec[2];
  float64_t curvmat[4];
  float64_t crossunitvec[2];
  
  udist_dus = (u2-u1)/du;
  vdist_dvs = (v2-v1)/dv;


  numsteps = (size_t)(1+sqrt(pow(udist_dus,2)+pow(vdist_dvs,2)));

  linevec_uv[0]=(u2-u1);
  linevec_uv[1]=(v2-v1);

  // Distance in uv space
  uvdist = normvec2(linevec_uv);
  //stepmag_uv = uvdist/(numsteps);
  

  if (uvdist==0.0) {
    *linelengthout=0.0;
    *sumcurvatureout=0.0; // end-user needs to divide out sum of weights 
    *sumcrosscurvatureout=0.0;
    return;
  }

  //// unit vector in uv space
  // (avoid using this because can be NaN if dist is zero)
  lineunitvec_uv[0]=linevec_uv[0]/uvdist;
  lineunitvec_uv[1]=linevec_uv[1]/uvdist;

  // step size
  //stepsize_uv[0]=stepmag_uv*lineunitvec_uv[0];
  //stepsize_uv[1]=stepmag_uv*lineunitvec_uv[1];
  stepsize_uv[0]=linevec_uv[0]/numsteps;
  stepsize_uv[1]=linevec_uv[1]/numsteps;
  


  // we split line into a bunch of segments and for each segment
  // just find the midpoint of the segment
  // and use that value as representative...
  curvaccum=0.0;
  crosscurvaccum=0.0; // curvature cross-ways (perpendicular) to the path
  
  dR=0.0;  
  for (stepnum=0;
       stepnum < numsteps;
       stepnum++) {
    //segstartu=u1 + stepsize_uv[0]*stepnum;
    //segstartv=v1 + stepsize_uv[1]*stepnum;
    //segendu=u1 + stepsize_uv[0]*(stepnum+1);
    //segendv=v1 + stepsize_uv[1]*(stepnum+1);
    segmidu = u1 + stepsize_uv[0]*(stepnum+0.5); /* segment midpoint in u */
    segmidv = v1 + stepsize_uv[1]*(stepnum+0.5); /* segment midpoint in v */

    //curvmatpoints are centered at (param_lowerleft_meaningfulunits_u + param_stepsize_u/2 + uidx*param_stepsize_u, ... same in v 

    // solve segmidu = lowerleft + stepsize/2 +uidx*stepsize 
    // uidx = (segmidu-lowerleft)/stepsize -0.5
    // uidx = round((segmidu-lowerleft)/stepsize -0.5)
    // uidx = floor((segmidu-lowerleft)/stepsize - 0.5 + 0.5)
    uidx = (uint32_t)((segmidu-param_lowerleft_meaningfulunits_u)/param_stepsize_u);
    vidx = (uint32_t)((segmidv-param_lowerleft_meaningfulunits_v)/param_stepsize_v);
    if (uidx < 0 || uidx >= nu || vidx < 0 || vidx >= nv || segmidu < param_lowerleft_meaningfulunits_u || segmidv < param_lowerleft_meaningfulunits_v) {
      /* off the world. Assume no curvature, but measure reasonable dR */
      dRdstep = sqrt(stepsize_uv[0]*stepsize_uv[0] + stepsize_uv[1]*stepsize_uv[1]);
  //*sumcurvatureout=0.0;
  //    *sumcrosscurvatureout=0.0;
      //// NaN
      //*linelengthout = my_infnan(0);
      //*sumcurvatureout = my_infnan(0);
      //*sumcrosscurvatureout = my_infnan(0);
      //return;
    } else {
      dRdstep = sqrt(pow(stepsizearray[ vidx*nu*2 + uidx*2 + 0]*stepsize_uv[0]/param_stepsize_u,2.0) + pow(stepsizearray[ vidx*nu*2 + uidx*2 + 1]*stepsize_uv[1]/param_stepsize_v,2.0));

      // &curvmats[vidx*nu*2*2 + uidx*2*2] is a 2x2 curvature matrix in (u,v)...
      // (u,v) are not necessarily orthogonal in real space
      
      curvmat[0]=curvmats[vidx*nu*2*2+uidx*2*2+0];
      curvmat[1]=curvmats[vidx*nu*2*2+uidx*2*2+1];
      curvmat[2]=curvmats[vidx*nu*2*2+uidx*2*2+2];
      curvmat[3]=curvmats[vidx*nu*2*2+uidx*2*2+3];
      
      
      // Really, curvmat is a shape operator.
      // To get a curvature along a line,
      // multiply on the left and right by a tangent vector that maps
      // to 1 meter physical length. The result will be
      // the curvature along that tangent vector (in 1/m).
      
      // ... But since it's the same vector on the
      // left and right, the scalings will be reciprocals,
      // so we can just multiply on the left and right by unit vectors
      
      // This then gives us the curvature along the line segment
      // in m^-1
      
      // So we accumulate a weighted average of the curvatures
      // of the line segments.
      
      if (!PSTCPAOPS_ISNAN(lineunitvec_uv[0]) && !PSTCPAOPS_ISNAN(lineunitvec_uv[1])) {
	// (If either of these are NaN then our step length is zero so
	// we really don'nt need to accumulate)
	multmatvec2(curvmat,lineunitvec_uv,curvtempvec);
	
	// dot product is evaluated curvature. dRdstep is weighting factor
	curvaccum += dotvecvec2(lineunitvec_uv,curvtempvec) * dRdstep;
	
	// cross-curvature: 
	// perpendicular vector in 2D can be found by interchanging
	// the elements and negating one fothem 
	crossunitvec[0]=lineunitvec_uv[1];
	crossunitvec[1]=-lineunitvec_uv[0]
;
	multmatvec2(curvmat,crossunitvec,curvtempvec);
	crosscurvaccum += dotvecvec2(crossunitvec,curvtempvec) * dRdstep;
      }
      // Curvature is positive when surface is concave,
      // negative when convex... by outward directed
      // selected in curvature calculation code by normal. 
    }

    dR += dRdstep;
          
    
  }

  
  *linelengthout=dR; 
  *sumcurvatureout=curvaccum; // end-user needs to divide out sum of weights 
  *sumcrosscurvatureout=crosscurvaccum;
}



// Linelength for a mirrored box in (u,v) space. Line is 
// presumed to start from (u1,v1) inside the box, and 
// connects (u2,v2) which 
// may be outside the box. 
// Locations outside the box are interpreted as mirrored 
// replicas of what is in the box 
static PSTCPAOPS_INLINE void linelength_avgcurvature_mirroredbox_meshbased_c_one(
								      
							   // This operatex on (u,v) in texcoord [0-1] values 
							   float32_t *curvmats, // indexed c-style (v,u,2,2)  ... must cover entire range of [0,1] texture coordinates
							   float32_t *stepsizearray,
							   size_t nv, // number of v steps in curvmats or stepsizearray
							   size_t nu, // number of u steps in curvmats or stepsizearray
							   float64_t param_lowerleft_meaningfulunits_u,
							   float64_t param_lowerleft_meaningfulunits_v,
							   float64_t param_stepsize_u, /* meaningful unit stepsize for curvmats/stepsizearray (nominal value), used to measure u1,v1,etc. */
							   float64_t param_stepsize_v,
							   float64_t boxu1,
							   float64_t boxv1,
							   float64_t boxu2,
							   float64_t boxv2,
							   float64_t du, // nominal resolution in u, meaning desired step size
							   float64_t dv, // nominal resolution in v, desired step size, in texcoord [0-1] units

							   float64_t u1,
							   float64_t v1,
							   float64_t u2,
							   float64_t v2,
							   float64_t *linelengthout,
							   float64_t *avgcurvatureout,
							   float64_t *avgcrosscurvatureout,
							   float64_t *thetaout) // theta=0 means the line is parallel to u at point #2, theta = pi/2 means parallel to v at point #2
{

  // Add up the physical length of all the line segments in the straight line
  // (in u,v space) from u1,v1 to u2,v2. Assumes that the straight line is
  // the shortest path, which isn't true in general. (shortest path
  // would be the geodesic.)

  
  float64_t segstart_u,segstart_v;
  float64_t segstart_mirrored_u,segstart_mirrored_v;
  float64_t segmid_u,segmid_v;
  float64_t segend_u,segend_v;
  float64_t segend_mirrored_u,segend_mirrored_v;

  float64_t segment_linelength;
  float64_t segment_sumcurvature;
  float64_t segment_sumcrosscurvature;
  float64_t total_linelength;
  float64_t sumcurvature;
  float64_t sumcrosscurvature;
  
  float64_t t_bu1,t_bu2,t_bv1,t_bv2;

  uint32_t uidx,vidx;

  // Find intersections between the line from (u1,v1) to (u2,v2)
  // and each of the box boundaries

  // parameterize line segment from (u1,v1) to (u2,v2)
  // as   u = u1 + t*(u2-u1)
  // and  v = v1 + t*(v2-v1) ... segment goes from t=0..1

  // intersections are at u= boxu1, u=boxu2, v=boxv1 and v=boxv2
  //
  // e.g boxu1 = u1 + t*(u2-u1)
  // solve for t
  //  boxu1-u1 = t*(u2-u1)
  // t_bu1 = (boxu1-u1)/(u2-u1)  intersectnum==0
  // t_bu2 = (boxu2-u1)/(u2-u1)  intersectnum==1
  // t_bv1 = (boxv1-u1)/(v2-v1)  intersect
  // t_bv2 = (boxv2-v1)/(v2-v1)

  // up to 4 intersections so 6 times, including 0 and 1
  float64_t times[6];
  float64_t temp;
  int interchangeflag=1;
  int numtimes=1;
  int cnt,segnum;
  float32_t curv_evals[2];

  times[0]=0.0;

  // (u1,v1) is supposed to be inside the box (don't think we actually rely on this
  assert(u1 >= boxu1 && u1 <= boxu2 && v1 >= boxv1 && v1 <= boxv2);
  // ... but won't work if an end is farther than one width from the box edge

  if (u2 != u1) {
    t_bu1 = (boxu1-u1)/(u2-u1);
    if (t_bu1 >= 0.0 && t_bu1 <= 1.0) {
      times[numtimes]=t_bu1;
      numtimes++;
    }

    t_bu2 = (boxu2-u1)/(u2-u1);
    if (t_bu2 >= 0.0 && t_bu2 <= 1.0) {
      times[numtimes]=t_bu2;
      numtimes++;
    }
  }
  if (v2 != v1) {
    t_bv1 = (boxv1-v1)/(v2-v1);
    if (t_bv1 >= 0.0 && t_bv1 <= 1.0) {
      times[numtimes]=t_bv1;
      numtimes++;
    }

    t_bv2 = (boxv2-v1)/(v2-v1);
    if (t_bv2 >= 0.0 && t_bv2 <= 1.0) {
      times[numtimes]=t_bv2;
      numtimes++;
    }

  }
  times[numtimes]=1.0; /* omit final increment */

  /* sort times... bubble sort */
  do {
    interchangeflag=0;
    for (cnt=1; cnt < numtimes-1; cnt++) {
      if (times[cnt+1] < times[cnt]) {
	/* out of order... interchange */
	temp=times[cnt+1];
	times[cnt+1]=times[cnt];
	times[cnt]=temp;
	
	interchangeflag=1;
      } 
    } 
  } while (interchangeflag);

  /* Ok... now we have the segments */
  total_linelength=0.0;
  sumcurvature=0.0;
  sumcrosscurvature=0.0;

  for (segnum=0,segstart_u = u1,segstart_v = v1;
       segnum < numtimes;
       segnum++,segstart_u=segend_u,segstart_v=segend_v) {
    segend_u = u1 + (u2-u1)*times[segnum+1];
    segend_v = v1 + (v2-v1)*times[segnum+1];

    segmid_u = (segstart_u + segend_u)/2.0;
    segmid_v = (segstart_v + segend_v)/2.0;

    segstart_mirrored_u=segstart_u;
    segstart_mirrored_v=segstart_v;

    segend_mirrored_u=segend_u;
    segend_mirrored_v=segend_v;

    // segment midpoint tells us where we are... 
    if (segmid_u > boxu2) {
      // mirrored around boxu2
      segstart_mirrored_u = boxu2 - (segstart_mirrored_u-boxu2);
      segend_mirrored_u = boxu2 - (segend_mirrored_u-boxu2);
    }

    if (segmid_u < boxu1) {
      // mirrored around boxu1
      segstart_mirrored_u = boxu1 + (boxu1-segstart_mirrored_u);
      segend_mirrored_u = boxu1 + (boxu1-segend_mirrored_u);
    }

    if (segmid_v > boxv2) {
      // mirrored around boxv2
      segstart_mirrored_v = boxv2 - (segstart_mirrored_v-boxv2);
      segend_mirrored_v = boxv2 - (segend_mirrored_v-boxv2);
    }

    if (segmid_v < boxv1) {
      // mirrored around boxv1
      segstart_mirrored_v = boxv1 + (boxv1-segstart_mirrored_v);
      segend_mirrored_v = boxv1 + (boxv1-segend_mirrored_v);
    }

    /* calculate distance based on mirrored copy */
    linelength_sumcurvature_meshbased_c_one(curvmats,
					    stepsizearray,
					    nv, nu, 
					    param_lowerleft_meaningfulunits_u,
					    param_lowerleft_meaningfulunits_v,
					    param_stepsize_u, /* meaningful unit stepsize for curvmats/stepsizearray (nominal value), used to measure u1,v1,etc. */
					    param_stepsize_v,
					    du, dv,
					    segstart_mirrored_u,
					    segstart_mirrored_v, 
					    segend_mirrored_u, 
					    segend_mirrored_v,
					    &segment_linelength,
					    &segment_sumcurvature,
					    &segment_sumcrosscurvature);
    total_linelength += segment_linelength;
    sumcurvature += segment_sumcurvature;
    sumcrosscurvature += segment_sumcrosscurvature;
  }

  // simple approximation for theta assuming coordinate axes are all straight
  *thetaout=atan2(v2-v1,u2-u1);
  
  *linelengthout=total_linelength; 
  *avgcurvatureout=sumcurvature/total_linelength; // divide out sum of weights  
  *avgcrosscurvatureout=sumcrosscurvature/total_linelength; // divide out sum of weights  
  if ((PSTCPAOPS_ISNAN(*avgcurvatureout) && !PSTCPAOPS_ISNAN(sumcurvature) && !PSTCPAOPS_ISNAN(total_linelength)) ||
      (PSTCPAOPS_ISNAN(*avgcrosscurvatureout) && !PSTCPAOPS_ISNAN(sumcrosscurvature) && !PSTCPAOPS_ISNAN(total_linelength))) {
    // This would be the case if total linelength is about zero. 
    // i.e. we have a point, not a line
    // ... But what direction? 

    // Evaluate the principal curvatures at this point.
    // assign (arbitrarily) one as avgcurvature and the other as avgcrosscurvature
    //curvmatpoints are centered at (param_lowerleft_meaningfulunits_u + param_stepsize_u/2 + uidx*param_stepsize_u, ... same in v 

    // solve segmidu = lowerleft + stepsize/2 +uidx*stepsize 
    // uidx = (segmidu-lowerleft)/stepsize -0.5
    // uidx = round((segmidu-lowerleft)/stepsize -0.5)
    // uidx = floor((segmidu-lowerleft)/stepsize - 0.5 + 0.5)
    uidx = (uint32_t)((u1-param_lowerleft_meaningfulunits_u)/param_stepsize_u);
    vidx = (uint32_t)((v1-param_lowerleft_meaningfulunits_v)/param_stepsize_v);
    if (uidx < 0 || uidx >= nu || vidx < 0 || vidx >= nv || u1 < param_lowerleft_meaningfulunits_u || v1 < param_lowerleft_meaningfulunits_v) {
      // NaN
      *avgcurvatureout = my_infnan(0);
      *avgcrosscurvatureout = my_infnan(0);
    } else {
      eigvals_2d_float(&curvmats[vidx*nu*2*2 + uidx*2*2],curv_evals);
      *avgcurvatureout= curv_evals[0];
      *avgcrosscurvatureout = curv_evals[1];
    }
  }
  
    
}






static PSTCPAOPS_INLINE void linelength_avgcurvature_meshbased_c_array(
							   // This operatex on (u,v) in texcoord [0-1] values 
							   float32_t *curvmats, // indexed c-style (v,u,2,2)  ... must cover entire range of [0,1] texture coordinates
							   float32_t *stepsizearray,
							   size_t nv, // number of v steps in curvmats or stepsizearray
							   size_t nu, // number of u steps in curvmats or stepsizearray
							   float64_t param_lowerleft_meaningfulunits_u,
							   float64_t param_lowerleft_meaningfulunits_v,
							   float64_t param_stepsize_u, /* meaningful unit stepsize for curvmats/stepsizearray (nominal value), used to measure u1,v1,etc. */
							   float64_t param_stepsize_v,
	float64_t du, // nominal resolution in u, meaning desired step size
	float64_t dv, // nominal resolution in v, desired step size, in texcoord [0-1] units
	float64_t *u1,
	float64_t *v1,
	float64_t *u2,
	float64_t *v2,
	size_t numlines,
	float64_t *linelengthout,
	float64_t *avgcurvatureout,
	float64_t *avgcrosscurvatureout)
{
  int64_t linecnt;  // must be signed (e.g. not size_t) for MSVC compatibility
  float64_t sumcurvature,sumcrosscurvature;

  // this loop can be freely parallelized
#pragma omp parallel default(shared) private(linecnt)
#pragma omp for 
  for (linecnt=0;linecnt < numlines;linecnt++) {
    float32_t curv_evals[2];
    uint32_t uidx,vidx;

    linelength_sumcurvature_meshbased_c_one(
					    curvmats,
					    stepsizearray,
					    nv,nu,
					    param_lowerleft_meaningfulunits_u,
					    param_lowerleft_meaningfulunits_v,
					    param_stepsize_u, /* meaningful unit stepsize for curvmats/stepsizearray (nominal value), used to measure u1,v1,etc. */
					    param_stepsize_v,
					    du,dv,
					    u1[linecnt],
					    v1[linecnt],
					    u2[linecnt],
					    v2[linecnt],
					    &linelengthout[linecnt],
					    &sumcurvature,
					    &sumcrosscurvature);
    
    avgcurvatureout[linecnt]=sumcurvature/linelengthout[linecnt];
    avgcrosscurvatureout[linecnt]=sumcrosscurvature/linelengthout[linecnt];
    if ((PSTCPAOPS_ISNAN(avgcurvatureout[linecnt]) && !PSTCPAOPS_ISNAN(sumcurvature) && !PSTCPAOPS_ISNAN(linelengthout[linecnt])) ||
	(PSTCPAOPS_ISNAN(avgcrosscurvatureout[linecnt]) && !PSTCPAOPS_ISNAN(sumcrosscurvature) && !PSTCPAOPS_ISNAN(linelengthout[linecnt]))) {
      // This would be the case if total linelength is about zero. 
      // i.e. we have a point, not a line
      // ... But what direction? 
      
      // Evaluate the maximum principal curvature at this point.
      // ***!!!! When we start implmenting tangential curvatures
      // we should supply one curvature as max principal 
      // and the other as min principal  ********
      
      //curvmatpoints are centered at (param_lowerleft_meaningfulunits_u + param_stepsize_u/2 + uidx*param_stepsize_u, ... same in v 
      
      // solve segmidu = lowerleft + stepsize/2 +uidx*stepsize 
      // uidx = (segmidu-lowerleft)/stepsize -0.5
      // uidx = round((segmidu-lowerleft)/stepsize -0.5)
      // uidx = floor((segmidu-lowerleft)/stepsize - 0.5 + 0.5)
      uidx = (uint32_t)((u1[linecnt]-param_lowerleft_meaningfulunits_u)/param_stepsize_u);
      vidx = (uint32_t)((v1[linecnt]-param_lowerleft_meaningfulunits_v)/param_stepsize_v);
      if (uidx < 0 || uidx >= nu || vidx < 0 || vidx >= nv || u1[linecnt] < param_lowerleft_meaningfulunits_u || v1[linecnt] < param_lowerleft_meaningfulunits_v) {
	// NaN
	avgcurvatureout[linecnt] = my_infnan(0);
	avgcrosscurvatureout[linecnt] = my_infnan(0);
      } else {
	eigvals_2d_float(&curvmats[vidx*nu*2*2 + uidx*2*2],curv_evals);
	avgcurvatureout[linecnt] = curv_evals[0];
	avgcrosscurvatureout[linecnt] = curv_evals[1];
      }
    }
  }
}



static PSTCPAOPS_INLINE void linelength_avgcurvature_mirroredbox_meshbased_c_array(
							   // This operatex on (u,v) in texcoord [0-1] values 
							   float32_t *curvmats, // indexed c-style (v,u,2,2)  ... must cover entire range of [0,1] texture coordinates
							   float32_t *stepsizearray,
							   size_t nv, // number of v steps in curvmats or stepsizearray
							   size_t nu, // number of u steps in curvmats or stepsizearray
							   float64_t param_lowerleft_meaningfulunits_u,
							   float64_t param_lowerleft_meaningfulunits_v,
							   float64_t param_stepsize_u, /* meaningful unit stepsize for curvmats/stepsizearray (nominal value), used to measure u1,v1,etc. */
							   float64_t param_stepsize_v,
							   float64_t boxu1,
							   float64_t boxv1,
							   float64_t boxu2,
							   float64_t boxv2,
							   float64_t du, // nominal resolution in u, meaning desired step size
							   float64_t dv, // nominal resolution in v, desired step size, in texcoord [0-1] units

							   float64_t *u1,
							   float64_t *v1,
							   float64_t *u2,
							   float64_t *v2,
							   size_t numlines,
							   float64_t *linelengthout,
							   float64_t *avgcurvatureout,
							   float64_t *avgcrosscurvatureout,
							   float64_t *thetaout)
{
  int64_t linecnt; // must be signed (e.g. not size_t) for MSVC compatibility


  // this loop can be freely parallelized
#pragma omp parallel default(shared) private(linecnt)
#pragma omp for 
  for (linecnt=0;linecnt < numlines;linecnt++) {
    linelength_avgcurvature_mirroredbox_meshbased_c_one(
							curvmats,
							stepsizearray,
							nv,nu,
							param_lowerleft_meaningfulunits_u,
							param_lowerleft_meaningfulunits_v,
							param_stepsize_u, /* meaningful unit stepsize for curvmats/stepsizearray (nominal value), used to measure u1,v1,etc. */
							param_stepsize_v,
							boxu1,
							boxv1,
							boxu2,
							boxv2,
							du,dv,
							u1[linecnt],
							v1[linecnt],
							u2[linecnt],
							v2[linecnt],
							&linelengthout[linecnt],
							&avgcurvatureout[linecnt],
							&avgcrosscurvatureout[linecnt],
							&thetaout[linecnt]);
  }
}
