import numpy as np

# Asymmetric patterns have adjacent rows shifted by half a unit
# pattern_size is in (x points, y points) 

def AsymCircPattern(pattern_size,scale):
    pts=[]
    offset=0.0
    hstep=1.0
    vstep=0.5
    for i in range(pattern_size[1]):
        for j in range(pattern_size[0]):
            pts.append([offset+j*hstep,i*vstep,0.0])
            pass
        if offset == 0.0:
            offset = 0.5
            pass
        else:
            offset = 0.0
            pass
        pass
    pattern=np.array(pts,np.float32)
    pattern=scale*pattern
    return(pattern)

# Asym2 starts its first row with offset present
def AsymCircPattern2(pattern_size,scale):
    pts=[]
    offset=0.5
    hstep=1.0
    vstep=0.5

    for i in range(pattern_size[1]):
        for j in range(pattern_size[0]):
            if i%2:
                pts.append([j*hstep+offset, i*vstep, 0.0])
            else:
                pts.append([j*hstep, i*vstep, 0.0])
                
    pattern=np.array(pts,np.float32)
    pattern=scale*pattern
    return(pattern)

def SymCircPattern(pattern_size,scale):
    pts=[]
    hstep=1.0
    vstep=1.0
    for i in range(pattern_size[1]):
        for j in range(pattern_size[0]):
            pts.append([j*hstep,i*vstep,0.0])
            pass
        pass
    pattern=np.array(pts,np.float32)
    pattern=scale*pattern
    return(pattern)

def ChessBoardPattern(pattern_size,scale):
    pts=[]
    hstep=1.0
    vstep=1.0
    for i in range(pattern_size[1]):
        for j in range(pattern_size[0]):
            pts.append([j*hstep,i*vstep,0.0])
            pass
        pass
    pattern=np.array(pts,np.float32)
    pattern=scale*pattern
    return(pattern)

