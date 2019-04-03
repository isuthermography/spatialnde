#ifdef _MSC_VER
#define PSAOPS_INLINE  __inline
#else
#define PSAOPS_INLINE  inline
#endif

static PSAOPS_INLINE int enclosed_or_intersecting_polygons_3d_c(int32_t *polypool,
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
					   uint32_t *num_fully_enclosed)
{
  // retpolys assumed to be at least as big as polypool
  size_t num_returned_polys=0;
  size_t poolidx;
  int32_t idx,firstidx;
  int32_t *vertexidxs;
  size_t numvertexidxs;

  int polygon_fully_enclosed;
  (*num_fully_enclosed)=0;

  
  for (poolidx=0;poolidx < polypoollen;poolidx++) {
    idx=polypool[poolidx];
    if (idx < 0) {
      // masked out polygon
      continue;
    }
    
    firstidx=vertexidx_indices[idx];
    vertexidxs=&vertexidx[firstidx];
    numvertexidxs=numvertices[idx];
    //fprintf(stderr,"idx=%d\n",idx);
    polygon_fully_enclosed = vertices_in_box_3d(vertices,vertexidxs,numvertexidxs,box_v0,box_v1);
    //if (idx==266241) {
    //  fprintf(stderr,"266241 v0=%f %f %f v1=%f %f %f fully_enclosed = %d\n",box_v0[0],box_v0[1],box_v0[2],box_v1[0],box_v1[1],box_v1[2],polygon_fully_enclosed);
    //}
    //fprintf(stderr,"idx2=%d\n",idx);
    if (polygon_fully_enclosed) {
      retpolys[num_returned_polys]=idx;
      num_returned_polys++;
      //fprintf(stderr,"fully_enclosed %d\n",idx);
      // if it's fully enclosed, nothing else need look at at, so we filter it here from the broader sibling pool
      polypool[poolidx] = -1; // mask out polygon

      (*num_fully_enclosed)++;

    } else {
      /* not polygon_fully_enclosed */

      // does it intersect?

      if (polygon_intersects_box_3d_c(box_v0,box_v1,vertices,vertexidxs,numvertexidxs,&inplanemats[idx*3*2],&facetnormals[idx*3])) {
	//fprintf(stderr,"returning %d\n",idx);
	retpolys[num_returned_polys]=idx;
	num_returned_polys++;
	//Don't filter it out in this case because it must
	// intersect with a sibiling too 
	//if (idx==266241) {
	//  fprintf(stderr,"266241 intersects_box\n");  
	//}
	
      }
    }
  }
  //fprintf(stderr,"num_returned_polys=%ld\n",num_returned_polys);
  //int cnt;
  //for (cnt=0;cnt < num_returned_polys;cnt++) {
  //  fprintf(stderr,"%d ",retpolys[cnt]);
  //}
  //fprintf(stderr,"\n");
  return num_returned_polys;
  
}
  
