import sys
import numpy as np
from numpy import linalg
from numpy.linalg import norm


from OCC.GeomLProp import GeomLProp_SLProps
from OCC.gp import gp_Pnt2d
from OCC.gp import gp_Vec
from OCC.gp import gp_Dir
from OCC.gp import gp_Pnt

from .vertexcoords import eval_vertex_caduv_coords
from .loaders import GetFacesSurfaces

def vecnorm(a,axis=-1):
    
    # compatibility for old versions of Numpy
    return np.sqrt(np.sum(a**2.0,axis))

class GeometricInconsistency(Exception):
    pass

def evalcurvature(CADsurface,CADu,CADv,finetol,approxoutwardnormal=None):
    """ Evaluate the curvature of the specified CAD surface
    at the given CAD (u,v) coordinates, with the given tolerance.
    Returns two element vector of principal curvatures (in 1/mm),
    followed by corresponding axes (in 3D space). Returns all NaN
    if curvature is undefined at this point. 

    WARNING: OpenCascade does not always select a normal in a consistent
             (i.e. u cross v pointing outside the object boundary) way.
             Set approxoutwardnormal to force the interpretation of the
             sense of the normal direction
             
"""
    # Create props object with derivatives
    Curvature = GeomLProp_SLProps(CADsurface,
                                  CADu,
                                  CADv,
                                  2, # Calculate 2 derivatives
                                  finetol)

    if not Curvature.IsCurvatureDefined() or not  Curvature.IsTangentUDefined() or not Curvature.IsTangentVDefined():
        # Return all NaN's. 
        return (np.array((np.NaN,np.NaN)),np.ones((3,2),dtype='d')*np.NaN,np.ones(3)*np.NaN,np.ones(3)*np.NaN,np.ones(3)*np.NaN)


    Sigma_u = Curvature.D1U() # vector 1st derivative along U axis
    Sigma_v = Curvature.D1V() # along V axis
    Sigma_uu = Curvature.D2U() # vector 2nd derivative along U axis
    Sigma_vv = Curvature.D2V() # vector 2nd derivative along V axis
    Sigma_uv = Curvature.DUV() # Vector cross-derivative
    Normal = gp_Vec(Curvature.Normal())
    NormalVec = np.array((Normal.X(),
                          Normal.Y(),
                          Normal.Z()),dtype='d')
    NormalVec=NormalVec/np.linalg.norm(NormalVec)

    if approxoutwardnormal is None:
        approxoutwardnormal=NormalVec
        pass
    
    
    
    TangentU_Dir=gp_Dir()
    Curvature.TangentU(TangentU_Dir) # Store tangent u direction in TangentU_Dir
    TangentU=gp_Vec(TangentU_Dir)
    
    TangentV_Dir=gp_Dir()
    Curvature.TangentV(TangentV_Dir) # Store tangent v direction in TangentV_Dir
    TangentV=gp_Vec(TangentV_Dir)


    ## if NormalVec is in the wrong direction we interchange U
    ## and V so as to correctly represent the desired right handedness
    ## of the frame with outward vector facing out.
    ## This also involves flipping NormalVec
    if np.inner(NormalVec,approxoutwardnormal) < 0.0:
        NormalVec=-NormalVec
        Normal=Normal.Reversed()

        temp1=Sigma_u
        Sigma_u=Sigma_v
        Sigma_v=temp1

        temp2=Sigma_uu
        Sigma_uu=Sigma_vv
        Sigma_vv=temp2

        temp3=TangentU
        TangentU=TangentV
        TangentV=temp3
        
        
        pass


    # First fundamental form:
    # See https://en.wikipedia.org/wiki/Differential_geometry_of_surfaces
    # https://en.wikipedia.org/wiki/First_fundamental_form
    # and http://www.cems.uvm.edu/~gswarrin/math295/Chap6.pdf
    
    E=Sigma_u.Dot(Sigma_u)
    F=Sigma_u.Dot(Sigma_v)
    G=Sigma_v.Dot(Sigma_v)
    
    # Second fundamental form:
    L=Sigma_uu.Dot(Normal)
    M=Sigma_uv.Dot(Normal)
    N=Sigma_vv.Dot(Normal)
    
    # Evaluate shape operator
    # https://en.wikipedia.org/wiki/Differential_geometry_of_surfaces#Shape_operator
    #  where e-> L, f->M, g->N
    
    S = (1.0/(E*G-F**2.0)) * np.array(((L*G-M*F,M*G-N*F),
                                       (M*E-L*F,N*E-M*F)),dtype='d') #  (see wikipedia + https://math.stackexchange.com/questions/536563/trouble-computing-the-shape-operator) http://mathworld.wolfram.com/WeingartenEquations.html

    # Shape operator in general will only be symmetric if (u,v) basis
    # is orthonormal, which it may not be
    #print(S)

    # Eigenvalue problem seems to be OK anyway so stuff below commented out
    #
    ## Apply Gram-Schmidt orthonormalization to find an ON basis
    Uvec=np.array((TangentU.X(),
                   TangentU.Y(),
                   TangentU.Z()),dtype='d')
    Uvec=Uvec/np.linalg.norm(Uvec)
    
    Vvec=np.array((TangentV.X(),
                   TangentV.Y(),
                   TangentV.Z()),dtype='d')
    Vvec=Vvec/np.linalg.norm(Vvec)

    # Check that Uvec cross Vvec = + Normal
    # Note that this verifies the handedness of the coordinate frame,
    # i.e. that U cross V gives an outward normal
    # and is thus important
    Wvec = np.cross(Uvec,Vvec)
    Wnormvec = Wvec/np.linalg.norm(Wvec)
    assert(np.abs(np.dot(Wnormvec,np.array((Normal.X(),Normal.Y(),Normal.Z()),dtype='d'))-1.0) < 0.025)
           
    ## Inner product matrix, for converting between covariant
    ## and contravariant components
    #IPMat = np.array( ((np.inner(Uvec,Uvec),np.inner(Uvec,Vvec)),
    #                   (np.inner(Vvec,Uvec),np.inner(Vvec,Vvec))),dtype='d')
    #
    #IPMatInv=np.linalg.inv(IPMat)
    #
    ##def proj(direc,vec): # Projection operator
    ##    return (np.inner(vec,direc)/np.inner(direc,direc))*direc
    #
    ## Vprimevec=Vvec-proj(Uvec,Vvec)
    #
    #Vprimevec=Vvec - (np.inner(Uvec,Vvec)/np.inner(Uvec,Uvec)) * Uvec
    #
    #Vprimenormfactor=1/np.linalg.norm(Vprimevec)
    #Vprimevec=Vprimevec*Vprimenormfactor
    #
    ## Covariant components of a 3D vector in (U,V)
    ## Upos_cov = Uvec dot 3dcoords 
    ## Vpos_cov = Vvec dot 3dcoords
    #
    ## S matrix times contravariant components -> Contravariant components
    ## so S * IPMatInv * (Uvec;Vvec) * 3dcoords -> Contravariant components
    #
    ## Let (Uprimepos,Vprimepos) = (Uvec dot 3dcoords, Vprimevec dot 3dcoords)
    ## ON basis
    #
    ## [ Uprime ] = [           1                                               0 ] [ Uvec  ]
    ## [ Vprime ]   [ -Vprimenormfactor*(Uvec dot Vvec)/(Uvec dot Uvec)    Vprimenormfactor ] [ Vvec  ]
    ##
    ## Therefore
    ## [ Uvec ] =   [          1                                                    0     ]^-1  [ Uprime ]
    ## [ Vvec ] =   [ -Vprimenormfactor*(Uvec dot Vvec)/(Uvec dot Uvec)  Vprimenormfactor ]     [ Vprime ]
    ##
    ## So a vector represented by  some coefficient times Uprime plus
    ## some other coefficient * Vprime, multiplied
    ## on the left by the inverse matrix, gives rise to
    ## some new coefficient * Uvec  plus  some other new coefficient * Vvec 

    ## These new coefficients are contravariant components.
    ##
    ## Let A = [           1                       0 ] 
    ##       = [ -(Uvec dot Vvec)/(Uvec dot Uvec)  1 ]*Vprimenormfactor 
    #
    ## Now A*S*inv(A)  provides a conversion of S to the Uprime,Vprime ON frame.
    ## ... Which should then be symmetric.
    #A = np.array(((1.0,0.0),
    #              (-Vprimenormfactor*np.inner(Uvec,Vvec)/np.inner(Uvec,Uvec),Vprimenormfactor)),dtype='d')
    #
    #Stransformed=np.dot(A,np.dot(S,np.linalg.inv(A)))
    
    ## Confirm that Stransformed is indeed symmetric to reasonable
    ## precision... Use the trace (twice the mean curvature)
    ## as a reference metric
    #assert(abs(Stransformed[0,1]-Stransformed[1,0]) <= 1e-7*np.trace(Stransformed))
    
    ## Force symmetry so our eigenvalue problem comes out real
    #Stransformed[0,1]=Stransformed[1,0]
    
    #(evals,evects) = np.linalg.eig(Stransformed)
    # Should be OK to just take eigenvalues of asymmetric matrix
    # per ShifrinDiffGeo.pdf page 51
    (curvatures,evects) = np.linalg.eig(S)

    
    ## Now evects * diag(evals) * (evects.T) should be Stransformed
    ## Convert to 3D basis
    ##  [ Uprime Vprime ] * evects * diag(evals) * evects.T * [ Uprime ; Vprime ]
    ## So store [ Uprime Vprime ] * evects  as principal curvature tangent directions
    ##   and evals as principal curvatures
    #curvaturetangents = np.dot(np.array((Uvec,Vprimevec),dtype='d').T,evects)
    curvaturetangents=np.dot(np.array((Uvec,Vvec)).T,evects)
    # per ShifrinDiffGeo.pdf page 52 confirm these tangents are orthogonal
    if np.dot(curvaturetangents[:,0],curvaturetangents[:,1]) >= 1e-8:
        raise GeometricInconsistency("Curvature tangents are not orthogonal (inner product = %g)" % (np.dot(curvaturetangents[:,0],curvaturetangents[:,1])))
    # We don't want the eigenframe to be mirrored relative to the (U,V)
    # frame, for consistency in interpreting positive vs. negative curvature.
    # ...  so if the dot/inner product of (UxV) with (TANGENT0xTANGENT1)
    # is negative, that indicates mirroring 
    # Negating one of the eigenvectors will un-mirror it.

    # The exception to this is if our normal and the desired approxoutwardnormal
    # are in opposite directions, then we DO want to mirror it, to
    # correct for OpenCascade's normal being in the wrong direction

    # Note: '^' here operating on booleans acts as a logical XOR operator
    #if not((np.inner(np.cross(Uvec,Vvec),np.cross(curvaturetangents[:,0],curvaturetangents[:,1])) < 0.0) ^ (np.inner(NormalVec,approxoutwardnormal) >= 0.0)):
    if np.inner(np.cross(Uvec,Vvec),np.cross(curvaturetangents[:,0],curvaturetangents[:,1])) < 0.0:
        curvaturetangents[:,0]=-curvaturetangents[:,0]
        pass

    ## Likewise, if NormalVec is in the wrong direction we flip it
    ## and UVec so as to correctly represent the desired right handedness
    ## of the frame with outward vector facing out.
    #if np.inner(NormalVec,approxoutwardnormal) < 0.0:
    #    NormalVec=-NormalVec
    #    Uvec=-Uvec
    #    pass
    

    # transpose curvaturetangents so it is human readable in the serialization
    # axes:  Which Vertex, tangent u or v, coord x,yz
    return (curvatures,curvaturetangents.transpose(1,0),NormalVec,Uvec,Vvec)

def eval_all_curvatures(polygonalsurface,Faces,Surfaces,SurfObjs,coarsetol,finetol,normaltolerance=0.1,curvaturetolerance=np.Inf,usemeshforsharpedges=True):
    # Evaluate curvatures at all vertices, converting CAD units of mm
    # into spatialnde untis of m
    # curvaturetolerance in units of mm^-1

    polygonalsurface.buildconnectivity() # will need polynum_by_vertex field
    
    
    #(vertex_cadsurfnums,vertex_caduv)=eval_vertex_caduv_coords(polygonalsurface,Faces,Surfaces,SurfObjs,coarsetol,finetol)
    dist_UVcoords_surfnum_list=eval_vertex_caduv_coords(polygonalsurface,Faces,Surfaces,SurfObjs,coarsetol,finetol)

    curvatures = np.empty((polygonalsurface.vertices.shape[0],2),dtype='d')
    curvature_tangent_axes = np.empty((polygonalsurface.vertices.shape[0],2,3),dtype='d')
    
    for vertnum in range(polygonalsurface.vertices.shape[0]):
        # For this vertex, iterate over all CAD surfaces in which the
        # vertex is found

        #if vertnum==112570:
        #    import pdb
        #    pdb.set_trace()
        #    pass

        vert_curvatures=[]
        vert_curvaturetangents=[]
        vert_normals=[]
        vert_uvec=[]
        vert_vvec=[]
        for surfacecnt in range(len(dist_UVcoords_surfnum_list[vertnum])):
            (dist,UVcoords,surfnum) = dist_UVcoords_surfnum_list[vertnum][surfacecnt]
            #if vertnum==8304:
            #    import pdb
            #    pdb.set_trace()

            # Evaluate approximate normal... note that this might fail
            # for vertices at a knife edge.... Should not be a problem
            # because we ignore the OpenCascade result for those points
            # anyway (does not give meaningful curvature along edges)
            
            # get roughnormal from 1st facet
            firstpolynum=polygonalsurface.polynum_by_vertex[vertnum][0]
            if polygonalsurface.normalidx is None:
                roughnormal=polygonalsurface.normals[firstpolynum]
                pass
            else:
                firstidx=polygonalsurface.vertexidx_indices[firstpolynum]
                numvertices=polygonalsurface.numvertices[firstpolynum]
                
                polyvertices=polygonalsurface.vertexidx[firstidx:(firstidx+numvertices)]         
                polyvertexidx=np.where(polyvertices==vertnum)[0][0]
            
                roughnormal=polygonalsurface.normals[polygonalsurface.normalidx[firstidx+polyvertexidx],:]
                pass
            

            
            
            # Evaluate the curvature at this vertex on this surface
            try:
                (surf_vert_curvatures,surf_vert_curvaturetangents,NormalVec,UVec,VVec)=evalcurvature(Surfaces[surfnum],
                                                                                                     UVcoords[0],
                                                                                                     UVcoords[1],
                                                                                                     finetol,
                                                                                                     approxoutwardnormal=roughnormal)
                vert_curvatures.append(surf_vert_curvatures)
                vert_curvaturetangents.append(surf_vert_curvaturetangents)
                vert_normals.append(NormalVec)
                vert_uvec.append(UVec)
                vert_vvec.append(VVec)
                pass
            except GeometricInconsistency as e:
                sys.stderr.write("%s\n" % (e.message))
                pass
            pass
        
        if len(vert_curvatures)==1:
            # Single unique curvature for this vertex
            curvatures[vertnum,:]=vert_curvatures[0]*1000.0 # convert 1/mm (CAD) to 1/m (SpatialNDE)
            curvature_tangent_axes[vertnum,:,:] = vert_curvaturetangents[0]
            pass
        elif len(vert_curvatures) > 1:
            # Multiple curvatures and/or surfaces

            # Check tolerance on normal
            normalarray=np.array(vert_normals,dtype='d')
            normalmean=np.mean(normalarray,axis=0)[np.newaxis,:]
            if (vecnorm(normalarray-normalmean) > normaltolerance).any():
                # Too much change in normal... bail and calculate
                # curvature numerically
                curvatures[vertnum,:]=np.NaN
                curvature_tangent_axes[vertnum,:,:]=np.NaN
                continue
                

            # Check tolerance on curvature            
            # Assemble 2x2 curvature matrices

            # ***!!! NOTE: BUG: If the sense of any of the surfaces
            # is reversed relative to another, averaging the
            # curvature matrix this way is not meaningful,
            #
            # The OpenCascade extraction code uses the surface normals
            # to swap (u,v) to get a frame where u cross v points out
            # of the object. But if this is used in environment where
            # those normals are not consistent, this averaging
            # wouldn't make sense.
            #
            # If we wanted to it would probably be possible to
            # flip any surfaces that have a normal in the
            # opposite direction, so as to ensure the
            # curvature matrices can be averaged easily. 
            curvmats=np.empty((len(vert_curvatures),2,2),dtype='d')
            for cnt in range(len(vert_curvatures)):
                uvmat = np.array((vert_uvec[cnt],vert_vvec[cnt]),dtype='d')
                curvmats[cnt,:,:] = np.dot(uvmat,np.dot(vert_curvaturetangents[cnt].T,vert_curvatures[cnt].reshape(2,1)*np.dot(vert_curvaturetangents[cnt],uvmat.T)))
                pass
            curvmatmean=np.mean(curvmats,axis=0)

            if (np.abs(curvmats-curvmatmean[np.newaxis,:,:]) > curvaturetolerance).any():
            
                curvatures[vertnum,:]=np.NaN
                curvature_tangent_axes[vertnum,:,:]=np.NaN
                continue

            # fix any skewness in curvmatmean
            asymmetry = curvmatmean[1,0]-curvmatmean[0,1]
            curvmatmean[1,0]-=asymmetry
            curvmatmean[0,1]+=asymmetry


            # Determine principal curvatures
            (princcurvs,evects) = np.linalg.eig(curvmatmean)

            # use first vert_uvec to convert evects to 3D
            curvtangentaxes=np.dot(np.array((vert_uvec[0],vert_vvec[0]),dtype='d').T,evects)
            # We don't want the eigenframe to be mirrored relative to the (U,V)
            # frame, for consistency in interpreting positive vs. negative curvature.
            # ...  so if the dot/inner product of (UxV) with (TANGENT0xTANGENT1)
            # is negative, that indicates mirroring 
            # Negating one of the eigenvectors will un-mirror it.
            
            if np.inner(np.cross(vert_uvec[0],vert_vvec[0]),np.cross(curvtangentaxes[:,0],curvtangentaxes[:,1])) < 0.0:
                curvtangentaxes[:,0]=-curvtangentaxes[:,0]
                pass

            # Make sure that u, v cross product constraint holds for other surfaces as well
            # (otherwise bail and assign NaN)
            normalalign = [ np.inner(np.cross(vert_uvec[cnt],vert_vvec[cnt]),np.cross(curvtangentaxes[:,0],curvtangentaxes[:,1])) >= 0.0 for cnt in range(len(vert_curvatures)) ]
            if not all(normalalign):
            
                curvatures[vertnum,:]=np.NaN
                curvature_tangent_axes[vertnum,:,:]=np.NaN
                continue
                
            
            
            # Then assign curvature
            curvatures[vertnum,:]=princcurvs*1000.0 # convert 1/mm (CAD) to 1/m (SpatialNDE)
            curvature_tangent_axes[vertnum,:,:] = curvtangentaxes.T
            
            pass
        else:
            # len(vert_curvatures==0).... Curvature evaluation failed.
            # Leave result blank, possibly trigger numerical calculation of 
            # curvature numerically
            curvatures[vertnum,:]=np.NaN
            curvature_tangent_axes[vertnum,:,:]=np.NaN
            pass

        pass

    if usemeshforsharpedges:
        # Fix up sharp (and other failed curvature calculations)
        # 
        for vertnum in range(polygonalsurface.vertices.shape[0]):
            # Find failed curvatures and calculate them numerically
            
            if np.isnan(curvatures[vertnum,0]):
                (curvatures[vertnum,:],curvature_tangent_axes[vertnum,:,:])=polygonalsurface._evalcurvature(vertnum)
            
                pass
            pass
        pass
    
    return (curvatures,curvature_tangent_axes)

def add_curvatures_to_surface(polygonalsurface,BRepShape,coarsetol,finetol,normaltolerance=0.05,curvaturetolerance=0.1,usemeshforsharpedges=True):

    (Faces,Surfaces,SurfObjs) = GetFacesSurfaces(BRepShape)

    (curvatures,curvature_tangent_axes) = eval_all_curvatures(polygonalsurface,Faces,Surfaces,SurfObjs,coarsetol,finetol,normaltolerance,curvaturetolerance,usemeshforsharpedges=usemeshforsharpedges)

    polygonalsurface.assign_curvatures(curvatures,curvature_tangent_axes)
    pass


   
