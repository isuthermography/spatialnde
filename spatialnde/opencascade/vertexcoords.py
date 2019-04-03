import sys
import numpy as np

from OCC.gp import gp_Pnt2d
from OCC.gp import gp_Vec
from OCC.gp import gp_Dir
from OCC.gp import gp_Pnt
from OCC.TopAbs import TopAbs_IN
from OCC.TopAbs import TopAbs_ON
from OCC.TopAbs import TopAbs_FACE
from OCC.BRepTopAdaptor import BRepTopAdaptor_FClass2d
from OCC.ShapeAnalysis import ShapeAnalysis_Surface


def eval_vertex_caduv_coords(polygonalsurface,Faces,Surfaces,SurfObjs,coarsetol,finetol):
    """Given OpenCascade face info (Faces,Surfaces,Surfobjs from 
    loaders.GetFacesSurfaces()), evaluate the CAD uv coordinates
    and CAD surface index of each vertex from the polygonalsurface.
    
    Note that the CAD coordinate frame and SpatialNDE frame must 
    match exactly except the CAD frame is presumed to be in units of mm,
    whereas the SpatialNDE frame is presumed to be in units of m

    Returns (array of CAD surface numbers, array of (u,v) CAD coordinates)
    """

    # Construct a classifier for each face that will tell us
    # if a given (u,v) point is inside the face boundary
    Classifiers = [ BRepTopAdaptor_FClass2d(Face,coarsetol) for Face in Faces ]
    ShapeAnalyzers = [ ShapeAnalysis_Surface(Surface) for Surface in Surfaces ]

    #vertex_cadsurfnums=np.empty(polygonalsurface.vertices.shape[0],dtype=np.int32)
    #vertex_caduv=np.empty((polygonalsurface.vertices.shape[0],2),dtype=np.float64)

    vertex_cadinfo=[]
    
    for vertexcnt in range(polygonalsurface.vertices.shape[0]):
        CoordsToFind=polygonalsurface.vertices[vertexcnt,:]*1000.0 # 1000 is conversion from meters into mm 
        PointToFind=gp_Pnt(*CoordsToFind)

        #dist_list=[]
        #UVcoords_list=[]
        #surfnum_list=[]

        dist_UVcoords_surfnum_list=[]
        
        dist_fulllist=[]
        surfnum_fulllist=[]
        UVcoords_fulllist=[]

        
        for SurfCnt in range(len(Faces)):
            UVcoords=ShapeAnalyzers[SurfCnt].ValueOfUV(PointToFind,finetol) # when repeating nearby points it may be faster to use NextValueOfUV()...
            reprojected=SurfObjs[SurfCnt].Value(UVcoords.X(),UVcoords.Y())
            reprojected_tup=(reprojected.X(),reprojected.Y(),reprojected.Z())

            dist = np.linalg.norm(CoordsToFind-reprojected_tup)
            dist_fulllist.append(dist)
            UVcoords_fulllist.append((UVcoords.X(),UVcoords.Y()))
            surfnum_fulllist.append(SurfCnt)

            if (dist < coarsetol):
                classification = Classifiers[SurfCnt].Perform(UVcoords)
                inside = (classification == TopAbs_IN or classification == TopAbs_ON)
                if inside:
                    #dist_list.append(dist)
                    #UVcoords_list.append(UVcoords)
                    #surfnum_list.append(SurfCnt)
                    dist_UVcoords_surfnum_list.append((dist,(UVcoords.X(),UVcoords.Y()),SurfCnt))
                    pass
                pass
            pass
        # dist_list should be almost always one but maybe 2 or 3 max when we're very close to an edge or corner
        if len(dist_UVcoords_surfnum_list) > 2:
            sys.stderr.write("eval_vertex_caduv_coords(): Got multiple (%d) (u,v) solutions\n" % (len(dist_UVcoords_surfnum_list)))
            pass
        if len(dist_UVcoords_surfnum_list) < 1:
            raise ValueError("Could not find point in a face or on a face boundary  within coarse tolerance of vertex #%d at (%g,%g,%g) mm. dist_fulllist=%s" % (vertexcnt,CoordsToFind[0],CoordsToFind[1],CoordsToFind[2],str(dist_fulllist)))
        #vertex_cadsurfnums[vertexcnt] = surfnum_list[np.argmin(dist_list)]
        #vertex_caduv[vertexcnt,:]=(UVcoords_list[np.argmin(dist_list)].X(),UVcoords_list[np.argmin(dist_list)].Y())

        vertex_cadinfo.append(dist_UVcoords_surfnum_list)
        
        pass

    #return (vertex_cadsurfnums,vertex_caduv)
    return vertex_cadinfo
        
        
