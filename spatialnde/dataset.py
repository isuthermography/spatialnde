import sys
import numbers

import numpy as np

# import dg_eval   # now imported at each use
# import dg_metadata as dgm  # now imported at each use

try: 
    from collections.abc import Sequence # python3
    pass
except ImportError:
    from collections import Sequence # python 2.x
    pass

class dataguzzlerpixelarray(object):
    # For use by SurfaceStructuredGridDataSet(object):

    # Handles index reversal, y-reversal, dynamic evaluation, indexing

    # Assign these on creation
    wfminfo=None
    wfmdict=None # Used for procexpr evaluation where we need to refer to other waveforms
    raw=None
    rgba=None

    # Auto-generated:
    ndim=None
    DimLen=None
    IniVal=None
    Step=None
    shape=None


    def __init__(self,**kwargs):
        self.raw=False
        self.rgba=False
        
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
        
            setattr(self,kwarg,kwargs[kwarg])
            pass

        import dg_eval
        (self.ndim,self.DimLen,self.IniVal,self.Step,self.bases)=dg_eval.geom(self.wfminfo,raw=self.raw)
        self.shape=tuple(self.DimLen[::-1])
        
        pass

    def __getitem__(self,slices):
        if not(isinstance(slices,tuple)):
            slices=(slices,)    # make it a tuple
            pass

        IniVal=[]
        Step=[]
        DimLen=[]
        #killaxes=[]

        cnt=0
        for slc in slices:
            if isinstance(slc,slice):
                strt=slc.start
                stop=slc.stop
                stp=slc.step
                if strt is None:
                    strt=0
                    pass
                if strt < 0:
                    strt+=self.DimLen[self.ndim-cnt-1]
                    pass
                if strt < 0:
                    strt=0
                    pass
                if strt > self.DimLen[self.ndim-cnt-1]-1:
                    strt = self.DimLen[self.ndim-cnt-1]-1
                    pass
            
                if stop is None:
                    stop=self.DimLen[self.ndim-cnt-1]
                    pass
                if stop < 0:
                    stop+=self.DimLen[self.ndim-cnt-1]
                    pass
                
                if stop > self.DimLen[self.ndim-cnt-1]:
                    stop=self.DimLen[self.ndim-cnt-1]
                    pass
                if stp is None:
                    stp=1
                    pass
            
                IniVal.insert(0,self.IniVal[self.ndim-cnt-1]+self.Step[self.ndim-cnt-1]*strt)
                Step.insert(0,self.Step[self.ndim-cnt-1]*stp)
                DimLen.insert(0,(stop-strt)//stp)
                cnt+=1
                pass
            elif isinstance(slc,numbers.Integral):
                IniVal.insert(0,self.IniVal[self.ndim-cnt-1]+self.Step[self.ndim-cnt-1]*slc)
                Step.insert(0,self.Step[self.ndim-cnt-1])
                DimLen.insert(0,1)
                
                pass
            else:
                # return self[:][slices]
                raise ValueError("dataguzzler pixel array can only be indexed by integer or slice, got %s. Try evaluating full array with [:], then indexing that." % (slc.__class__.__name__))
            
            pass
        # Fill out remaining axes with full sweep
        while cnt < self.ndim:
            IniVal.insert(0,self.IniVal[self.ndim-cnt-1])
            Step.insert(0,self.Step[self.ndim-cnt-1])
            DimLen.insert(0,self.DimLen[self.ndim-cnt-1])
            cnt+=1
            pass

        import dg_eval
        data=dg_eval.evalv(self.wfminfo,self.wfmdict,self.ndim,IniVal,Step,DimLen,raw=self.raw,rgba=self.rgba).data
        data=data.transpose() # reverse axis interpretation from (x,y,frameno) to (frameno,y,x)

        # Reverse data in y, so now like a raster scan
        dataslices=[ slice(None) ] * (self.ndim - 2) # no change to first indices
        dataslices.append(slice(None,None,-1))  # reverse y
        dataslices.append(slice(None))  # No change to x

        data=data[dataslices]  # apply reversal
        return data
    
    pass
# Register dataguzzlerpixelarray with Sequence abstract base class
Sequence.register(dataguzzlerpixelarray) 



class SurfaceStructuredGridDataSet(object):
    n=None # Number of rows
    m=None # Number of columns
    p=None # number of coordinates per node
    data=None # Numpy data array: ... x n x m  represented as C-stored raster scan, indexed y,x RECOMMEND ALWAYS ACCESS THROUGH SLICING ([])
    mask=None # n x m boolean mask indicating valid locations,
              # or None to indicate 'we don't know'  NOT YET IMPLEMENTED
    coords=None # n x m x p array with coordinates of surface point i,j
    landmarkpixelcoords=None # dictionary by landmark name of p-tuples
                             # indicating coordinates of that landmark, in
                             # pixels measured from the upper-left corner

    xinival=None   # x coordinate of center of upper-left pixel
    yinival=None   # y coordinate of center of upper left pixel
    xstep=None  # x pixel step size
    ystep=None  # y pixel step size

    def __init__(self,**kwargs):
        
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
        
            setattr(self,kwarg,kwargs[kwarg])
            pass

        if self.data is not None:
            (self.n,self.m)=self.data.shape[:2]
            pass

        if self.coords is not None:
            self.p=self.coords.shape[2]
            pass

        pass

    def pixelcoords2physcoords(self,pixelcoords):
        assert(len(pixelcoords)==2)

        return np.array((self.xinival+self.xstep*pixelcoords[0],self.yinival+self.ystep*pixelcoords[1]),dtype='d')
    

    @classmethod
    def fromdataguzzlerpixelimage(cls,wfminfo,wfmdict,raw=False):

        import dg_eval
        import dg_metadata as dgm

        (ndim,DimLen,IniVal,Step,bases)=dg_eval.geom(wfminfo,raw=raw)
        (ndim,Coord,Units,AmplCoord,AmplUnits)=dg_eval.axes(wfminfo,raw=raw)

        #data=wfminfo.data.reshape(wfminfo.data.shape[::-1],order='F')  # reverse axis interpretation from (x,y,frameno) to (frameno,y,x)
        #data=data[:,::-1,:] # Reverse data in y, so now like a raster scan
        data=dataguzzlerpixelarray(wfminfo=wfminfo,wfmdict=wfmdict,raw=raw,rgba=False)   # This effectively reverses axes and flips y so now appears like a raster scan
        n=DimLen[1]
        m=DimLen[0]
        p=2

        #xcoords=np.arange(m,dtype='d')
        #ycoords=np.arange(n,dtype='d')
        #coords=np.concatenate((np.ones((n,1,1),dtype='d')*xcoords.reshape(1,m,1),
        #                      ycoords.reshape(n,1,1)*np.ones((1,m,1),dtype='d')),
        #                      axis=2)
        xcoords = bases[0]
        ycoords = bases[1][::-1]  # reversal because we reversed the data in y
        coords = np.concatenate((ycoords.reshape(n,1,1)*np.ones((1,m,1),dtype='d'),
                                 np.ones((n,1,1),dtype='d')*xcoords.reshape(1,m,1)),axis=2)

        xinival=IniVal[0]
        yinival=IniVal[1]+Step[1]*(DimLen[1]-1)
        xstep=Step[0]
        ystep=-Step[1]

        
        # Treat dataguzzler FIDUCIALs as landmarks
        # Technically, landmarks are probably permanent while
        # fiducials are often temporary
        landmarkpixelcoords={}
        for mdname in wfminfo.MetaData:
            if (mdname.startswith("FIDUCIAL_") or mdname.startswith("LANDMARK_")) and mdname.endswith("_X"):
                fiducialname=mdname[9:-2]
                fiducialmdname=mdname[:-2]

                landmarkx=dgm.GetMetaDatumWIDbl(wfminfo,"%s_X" % (fiducialmdname),np.NaN)
                landmarky=dgm.GetMetaDatumWIDbl(wfminfo,"%s_Y" % (fiducialmdname),np.NaN)

                # Convert landmark coordinates to pixel coords
                landmarkx -= IniVal[0]
                landmarkx /= Step[0]

                landmarky -= IniVal[1]
                landmarky /= Step[1]

                # Flip y axis to evaluate like raster scan (OpenCV) rather than like dataguzzler
                
                landmarky = n-1 - landmarky
                landmarkpixelcoords[fiducialname]=(landmarkx,
                                                   landmarky)
                
                pass
            pass
                

        return cls(data=data,coords=coords,landmarkpixelcoords=landmarkpixelcoords,xinival=xinival,yinival=yinival,xstep=xstep,ystep=ystep)
    pass

