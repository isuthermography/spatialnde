#ifdef _MSC_VER
#define VECOPS_INLINE  __inline
#else
#define VECOPS_INLINE  inline
#endif

static const uint32_t NaNconst=0x7fc00000;
static const uint32_t Infconst=0x7f800000;
static const uint32_t NegInfconst=0xff800000;


double my_infnan(int error) /* be sure to disable SIGFPE */
{
  float32_t retval; 

  assert(sizeof(float32_t)==4);
  if (error==ERANGE) {
    memcpy(&retval,&Infconst,4);
  } else if (error==-ERANGE) {
    memcpy(&retval,&NegInfconst,4);
  } else {
    memcpy(&retval,&NaNconst,4);
  }
  return (double)retval;
}



static VECOPS_INLINE void multmat23vec(float64_t *mat,float64_t *vec,float64_t *out)
/* Multiply 2x3 matrix by 3-vector, giving 2-vector */
{
  int outel,sumidx;

  for (outel=0;outel < 2; outel++) {
    out[outel]=0.0;
    for (sumidx=0;sumidx < 3; sumidx++) {
      out[outel] = out[outel] + mat[ outel*3 + sumidx]*vec[sumidx];
    }
  }
}

static VECOPS_INLINE void multvecmat23(float64_t *vec,float64_t *mat,float64_t *out)
/* Multiply 2-vector by 2x3 matrix giving 3-vector  */
{
  int outel,sumidx;

  for (outel=0;outel < 3; outel++) {
    out[outel]=0.0;
    for (sumidx=0;sumidx < 2; sumidx++) {
      out[outel] = out[outel] + mat[ sumidx*3 + outel]*vec[sumidx];
    }
  }
}

static VECOPS_INLINE void multmatvec4(float64_t *mat,float64_t *vec,float64_t *out)
{
  int outel,sumidx;

  for (outel=0;outel < 4; outel++) {
    out[outel]=0.0;
    for (sumidx=0;sumidx < 4; sumidx++) {
      out[outel] = out[outel] + mat[ outel*4 + sumidx]*vec[sumidx];
    }
  }
}

static VECOPS_INLINE void multmatvec3(float64_t *mat,float64_t *vec,float64_t *out)
{
  int outel,sumidx;

  for (outel=0;outel < 3; outel++) {
    out[outel]=0.0;
    for (sumidx=0;sumidx < 3; sumidx++) {
      out[outel] = out[outel] + mat[ outel*3 + sumidx]*vec[sumidx];
    }
  }
}

static VECOPS_INLINE void multmatvec2(float64_t *mat,float64_t *vec,float64_t *out)
{
  int outel,sumidx;

  for (outel=0;outel < 2; outel++) {
    out[outel]=0.0;
    for (sumidx=0;sumidx < 2; sumidx++) {
      out[outel] = out[outel] + mat[ outel*2 + sumidx]*vec[sumidx];
    }
  }
}


static VECOPS_INLINE float64_t dotvecvec3(float64_t *vec1,float64_t *vec2)
{
  int sumidx;
  float64_t val=0.0;
  for (sumidx=0;sumidx < 3; sumidx++) {
    val = val + vec1[sumidx]*vec2[sumidx];
    
  }
  return val;
}

static VECOPS_INLINE float64_t dotvecvec2(float64_t *vec1,float64_t *vec2)
{
  int sumidx;
  float64_t val=0.0;
  for (sumidx=0;sumidx < 2; sumidx++) {
    val = val + vec1[sumidx]*vec2[sumidx];
    
  }
  return val;
}


static VECOPS_INLINE void scalevec3(float64_t coeff,float64_t *vec1,float64_t *out)
{
  size_t cnt;
  for (cnt=0;cnt < 3; cnt++) {
    out[cnt]=coeff*vec1[cnt];
  }
}


static VECOPS_INLINE void scalevec2(float64_t coeff,float64_t *vec1,float64_t *out)
{
  size_t cnt;
  for (cnt=0;cnt < 2; cnt++) {
    out[cnt]=coeff*vec1[cnt];
  }
}



static VECOPS_INLINE void subvecvec3(float64_t *vec1,float64_t *vec2,float64_t *out)
{
  int outidx;

  for (outidx=0;outidx < 3; outidx++) {
    out[outidx] = vec1[outidx] - vec2[outidx];
    
  }
}

static VECOPS_INLINE void subvecvec2(float64_t *vec1,float64_t *vec2,float64_t *out)
{
  int outidx;

  for (outidx=0;outidx < 2; outidx++) {
    out[outidx] = vec1[outidx] - vec2[outidx];
    
  }
}

static VECOPS_INLINE void addvecscaledvec3(float64_t *vec1,float64_t coeff, float64_t *vec2,float64_t *out)
{
  int outidx;

  for (outidx=0;outidx < 3; outidx++) {
    out[outidx] = vec1[outidx] + coeff* vec2[outidx];
    
  }
}



static VECOPS_INLINE void normalize_wcoord4(float64_t *vec)
/* operates in-place */
{
  vec[0] /= vec[3];
  vec[1] /= vec[3];
  vec[2] /= vec[3];
  vec[3] = 1.0;
  
}



static VECOPS_INLINE double to_unit_vector4(float64_t *vec)
/* operates in-place... returns scaling factor */
{
  float64_t factor;

  factor=1.0/sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  assert(vec[3]==0.0); /* vectors should have no 'w' component */
  vec[0] *= factor;
  vec[1] *= factor;
  vec[2] *= factor;
  //vec[3] *= factor;

  return factor;
  
}


static VECOPS_INLINE float64_t normvec3(float64_t *vec)
/* returns vector norm */
{
  float64_t factor;

  factor=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  return factor;
}

static VECOPS_INLINE float64_t normvec2(float64_t *vec)
/* returns vector norm */
{
  float64_t factor;

  factor=sqrt(vec[0]*vec[0]+vec[1]*vec[1]);
  return factor;
}


static VECOPS_INLINE double to_unit_vector3(float64_t *vec)
/* operates in-place */
{
  float64_t factor;

  factor=1.0/sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  vec[0] *= factor;
  vec[1] *= factor;
  vec[2] *= factor;

  return factor;
  
}


static VECOPS_INLINE void sign_nonzero3(float64_t *input,float64_t *output)
{
  int cnt;
  for (cnt=0;cnt < 3;cnt++) {
    if (input[cnt] < 0.0) output[cnt]=-1.0;
    else output[cnt]=1.0;
  }
}


static VECOPS_INLINE void multvecvec3(float64_t *vec1,float64_t *vec2,float64_t *output)
{
  int cnt;
  for (cnt=0;cnt < 3; cnt++) {
    output[cnt]=vec1[cnt]*vec2[cnt];
  }
}

static VECOPS_INLINE void crossvecvec3(float64_t *vec1,float64_t *vec2,float64_t *output)
{
  /* 
     vec1 cross vec2 

   |   i     j     k    |
   | vec10 vec11  vec12 |
   | vec20 vec21  vec22 |
  */
  output[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1];
  output[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2];
  output[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0];
}

static VECOPS_INLINE float64_t crossvecvec2(float64_t *vec1,float64_t *vec2)
{
  /* 
     vec1 cross vec2 

   |   i     j     k    |
   | vec10 vec11        |
   | vec20 vec21        |
  */
  return vec1[0]*vec2[1]-vec1[1]*vec2[0];
}


static VECOPS_INLINE void mean2vec3(float64_t *vec1,float64_t *vec2,float64_t *out)
{
  int cnt;
  
  for (cnt=0;cnt < 3; cnt++) {
    out[cnt]=(vec1[cnt]+vec2[cnt])/2.0;
  }
}
