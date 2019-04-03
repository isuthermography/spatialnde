import numpy as np

def ray_to_plane_distance_OBSOLETE_NOTUSED(raypoints,raydirecs,planepoints,planenormals):
    # raypoints is nx3 array of ray starting points
    # raydirecs is nx3 array of ray directions (unit vectors)

    # planepoints is mx3 array of plane points
    # planenormals is mx3 array of plane normals (unit vectors)
    
    # returns nxm array of distances 
    
    # Vector equation for ray line
    # (x,y,z) = a*direc  + linepoint
    
    # vector equation for surface
    # ((x,y,z) - surfacepoint) dot normal=0
    
    # substitute line eqn into surface eqn
    # (a*direc + linepoint - surfacepoint) dot normal = 0
    # a*direc dot normal + linepoint dot normal - surfacepoint dot normal = 0
    # a*(direc dot normal) = surfacepoint dot normal - linepoint dot normal
    # a = (surfacepoint dot normal - linepoint dot normal)/(direc dot normal)
    # a = ((surfacepoint - linepoint) dot normal)/(direc dot normal)

    if len(raypoints.shape) < 2:
        raypoints=raypoints.reshape(1,raypoints.shape[0])
        pass

    if len(raydirecs.shape) < 2:
        raydirecs=raydirecs.reshape(1,raydirecs.shape[0])
        pass
    
    if len(planepoints.shape) < 2:
        planepoints=planepoints.reshape(1,planepoints.shape[0])
        pass

    if len(planenormals.shape) < 2:
        planenormals=planenormals.reshape(1,planenormals.shape[0])
        pass

    n=raypoints.shape[0]
    m=planepoints.shape[0]
    
    fdists=np.tensordot(planepoints.reshape((1,m,3)) - raypoints.reshape((n,1,3)),planenormals.reshape(1,m,3),axes=((2,),(2,),)) / np.inner(raydirecs,planenormals)

    # Note: Will get NaN in cases where the plane normal is perpendicular to the ray direction

    # (np.inner(planepoints,planenormals) - np.inner(localpoints,frontnormal))/(np.inner(localdirecs,frontnormal))
    return fdists

