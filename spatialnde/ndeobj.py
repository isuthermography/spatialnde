import sys
import string
import copy
import traceback
import collections
import json

from .cadpart.cadpart import cadpart
from .landmarks import point_landmarks

try:
    basestring # python 2.x
    pass
except NameError:
    basestring=str # python 3.x
    pass


import numpy as np





#from .uv_parameterization import parameterization

# Need to define parameterization "UV unwrapping" (for faceted representation, break into surfaces 
# that pass vertical line test. Then project surface onto a flat plane and plane
# coordinates become parameterization.
#
# Also need to define bounding region for the parameterization.
# Need to identify adjacent surface edges?
#
# Generate X3D model
#
# NOTE: X3DOM does not in fact support NURBS natively.
# X3DOM NURBS demos are implementing using code from Ayam (Maya backwards)
# that implements and X3DOM plugin, but (unfortunately) does not
# currently implement NURBS texture mapping
# ... for X3D NURBS: 
#	  <!-- # knots should be order + dimension ?  -->
#	  <!-- # weights should be uDimension x vDimension -->

#
# What about accuracy of NURBS surface representations?
#
# Define representation for mapping between parameterizations. (u,v) -> (s,t).... Bijective? 

# Generate x3d geometry gl with gltf externalgeometry nodes

# Notes: "RoadKill" UV unwrapper for Maya, http://www.pullin-shapes.co.uk/page8.htm
# based on blender LSCM and ABF code
# Paper: A. Sheffer, B. Levy, M. Mogilnitsky, A, Bogomyakov, ABF++: Fast and Robust
# Angle Based Flattening, ACM Transactions on Graphics, 24(2), 311-330 2005.
# Blender implementation at https://github.com/dfelinto/blender/tree/master/source/blender/editors/uvedit  UVMap MeshTexturePolyLayer
# http://blender.stackexchange.com/questions/4820/exporting-uv-coordinates

# NOTE: x3d NURBS control points are in homogeneous coordinates
# http://castle-engine.sourceforge.net/x3d_implementation_nurbs.php#section_homogeneous_coordinates

# Nurbs engine: verbnurbs.com

# Brujic, Risti, Ainsworth, Measurement-based modification of NURBS surfaces: http://www3.imperial.ac.uk/portal/pls/portallive/docs/1/44247.PDF

#
# Thoughts on (u,v) parameterizations
#    * Built-in NURBS parameterization from CAD is arbitrary
#    * May want one where (u,v) in meaningful units, not [0,1]
#    * Will need seams in parameterization, perhaps OK if seams
#      are imperfect but ONLY in desired parameterization, not
#      in CAD parameterization.
#    * For rendering, x3d defines (u,v)  parameterization for texture
#      based on same (u,v) as geometry, converting to texture coords.
#      This is conceptually an inverse mapping of texture coords
#      to CAD (u,v)
#    * We probably want a forward mapping instead.
#      Why? Because we would like to be able to support overlapping
#      seams so that defects at a seam are properly analyzed. So
#      from our new (s,t) parameterization to CAD (u,v) would be
#      unique, but CAD (u,v) -> (s,t) would give multiple (s,t)
#      coordinates for the same (u,v) (at the seam overlaps)
#    * Will also need to accommodate where a seam ends or branches.
#      because otherwise we will not be able to represent adequate
#      overlap at the joint. Probably solve by creating a separate
#      patch in a separate portion of (s,t) space to cover the joint.
#    * These new parameterizations can be automatically generated
#      but expect that in most cases it will need to be done
#      manually, as in computer graphics.
#    * May be able to use Blender as a manual tool for parameterization,
#      at least initially.
#
#    * Intentional seams are not a problem -- they should have overlap
#      anyway, and the user chose their location
#    * Intrinsic seams (seams from the intrinsic (u,v) mapping are
#      more problematic because they are an arbitrary implementation
#      detail of the CAD structure. Want to avoid gaps, etc.
#    * Proposed representation: Identify surface location
#      by (patch number,s,t)  (s,t parameterization may overlap across patches)
#    * Each patch's (s,t) space is then carved up by a bunch of seam curves
#      corresponding to intrinsic seams into sub-patches
#    * Each subpatch has its own 2D NURBS representation, which maps
#      from (s,t) to (u,v) within the patch. It is carefully calculated
#      to ensure a precise mapping of its edges and seams to within
#      numerical tolerance, so that there are no gaps in the reassembled
#      structure. 

# Discussion with Rafael, Adarsh, and Norbert, 11/4/16
#   * Consensus: Not sure if (u,v) mapping overlap is really
#     necessary, but supporting it doesn't seem like a bad thing
#   * Question: Use mapping that preserves angle?  YES Believe
#     NURBS mappings always preserve angle.
#   * Question: Can NURBS mapping represent rotation? Maybe...
#     might be limited to < 90 degrees because of
#     horizontal line test (?)
#   * When providing information to NDE interpretation algorithm
#     will need to provide (a) local curvatures, and other local
#     geometry info, and (b) information on the local (u,v)
#     distortion field (symmetric 2x2 tensor) 

#   * Need to be able to parameterize not just a single model, but
#     a __class__ of models
#     * For example define whether coordinates will stretch with
#       a longer object or run further
#       * e.g. a cylinder. Probably want to define (theta,z) param.
#         and distances around theta will stretch with
#         circumference, but a longer cylinder will extend to larger
#         z values
#       * How to deal with thermal expansion/contraction? Want
#         parameterization to identify same physical material
#         within a single specimen. How about bending due to load?
#       * So after initial parameterization, __always__ deal with
#         shape changes by stretching
#         * But exactly where? May need thermal or stress state model!
#         * Probably assume uniform unless we know otherwise,
#           such as from bending strains. 
#     * Need to be able to define procedure for assigning a
#       parameterization
#       * Procedure would include identifying particular points
#     * That procedure, executed, on an accurate model of a
#       reference configuration will define a parameterization
#       of this particular object.
#       * Specifically, parameterization upper boundaries should
#         be assigned by initial parameterization
#       * Following parameterizations may be stretched or deformed
#         compared to initial.
#       * Procedures for identifying landmarks
#     * Distinction between physical shape, which includes various
#       accoutrements and structural shape, which does not (and
#       also may not include e.g. corrosion thinning
#     * Structural shape useful perhaps for analyzing deformation
#       between configurations based on strains


# A data set representing data from or under a surface
#  includes a surface parameterization, but that parameterization
#  is typically unknown. Also, only part of the parameter space
# is valid, and this region of validity may or may not be known. 
#
# Need to identify that parameterization and map it onto
# our reference parameterization, for the object in
# the measured configuration.
#
# Specifically, need a model that relates that
# parameterization to the physical surface in the measured config.
#
# Next, we can solve for unknowns in the model given landmarks
# and/or other measured registration sources.
#
# For image data, this is the n-point pose problem.
# For an ultrasonic scanner, we would need a kinematic model
# for the scan
#
# Then we can relate the data parameterization to the physical
# surface and relate the reference parameterization to the
# physical surface, so we can relate the data parameterization
# to the reference parameterization.
#
#
# Specifically,
# SurfaceStructuredGridDataSet gives a set of coordinates
# for each valid surface element and for each identified
# landmark.
#
# coordinates -> model gives 3D location
# Landmark paramaterization coords -> parameterization gives 3D location
#
# Set these coordinates equal to each other and solve
#
# Problem: If we just minimize residual then our error is 3-space location
# error which may not be appropriate
#   Could request known parameter std dev and add in explicit error term
#   that is solved for alongside everything else, then minimized
#
# WILL NOT WORK FOR PROJECTION MODEL
# instead, let model define metric_from_3D which converts 3-space location
# from transformed parameterization to a vector of "comparables"
# and coord_to_metric, which converts structured data grid coords to a
# comparable metric. Then optimize the sum squared metric error to find
# the transformation.

# !!!*** Ask Rafael: does matrix algebra on affine transforms
# always yield another nicely normalized transform (with
# bottom row (0,0,0,1), if not, what is the normalization procedure

#
#
# Next step: Provide a way to load in landmarks on 3D objects
#
# * DGS File with landmarks identified on texture image
# * Built with special dataguzzler config
# * Need to be able to convert landmarks from texture image
#   coordinates to 3D coordinates in part (and identify non-unique mapping) 
# * Should be able then to have landmarks defined both in
#   part and in structured grid off a live thermal image, 
#   and evaluate the image projection model
# * Use OpenCV projection operation to project thermal image
#   onto UV parameterization model (need to Z-buffer?)
# * Define DG Math function that evaluates the image projection
#   model from defined landmarks and loadable geometry
# * Define DG Math function that performs projection into
#   parameterization space given projection model.
# * Define DG math function that creates VRML waveform for scope viewing





    
# Representation of rays

# vector of 6 doubles and a coordinate frame
# doubles represent
#  Source point (x0,y0,z0) then unit direction vector (x,y,z) components (a,b,c)

# ray points can be represented as
# (x,y,z) = (x0+at,y0+bt,z0+ct) for some parameter t >= 0 


# NOTE see class point_landmark way below!!!!
#class landmark(object):
#    
#    pass





class ndepart(object):
    frame = None   # Applies to landmarks and implpart and everything inside
    landmarks = None
    implpart = None   # implementation of part: class cadpart or similar
    #uvparam = None # canonical re-parameterization of intrinsic cadpart
    #               # parameterization

    def __init__(self,**kwargs):
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
        
            setattr(self,kwarg,kwargs[kwarg])
            pass
        pass

    def ray_intersections(self,frame,rays):
        # rays is a set of n rays, each of which is specified by 6 coordinates, source point then unit direction vector, frame is the coordframe of those rays

        localpoints = frame.transformto(self.frame,rays[:,:3]) # nx3
        localdirecs = frame.transformto(self.frame,rays[:,3:]) # nx3

        (firstintersectpoints,firstintersectnormals,firstintersectdists) = self.implpart.ray_intersections(localpoints,localdirecs)
        
        return (localpoints,localdirecs,firstintersectpoints,firstintersectnormals,firstintersectdists)

    @staticmethod
    def splitparams(params):

        if params is None:
            return (None,None)
        
        # Assume sub-params are all json too...
        #parsed=json.loads(params)

        assert(len(params)==2)
        ## our params are a list or tuple... form a list
        ## and re-encode sub-params
        #components=[ json.dumps(parsedparam) for parsedparam in parsed ]
        #return components
        return (params[0],params[1])

    @staticmethod
    def joinparams(*params):
        # Assume sub-params are all json too...
        parsed=[ json.loads(param) for param in params ]

        return json.dumps(parsed)

    @classmethod
    def fromstl(cls,frame,ndepartparams,filename,implclass=cadpart,tol=1e-4,recalcnormals=False,metersperunit=1.0,defaultappearance=None):

        (landmarks_params,implpartparams)=cls.splitparams(ndepartparams)
        landmarks=point_landmarks.fromparams(landmarks_params)
        
        implpart=implclass.fromstl(filename=filename,implpartparams=implpartparams,tol=tol,recalcnormals=recalcnormals,metersperunit=metersperunit,defaultappearance=defaultappearance)

        #uvparam=None # Reference or canonical uv parameterization
        #if uvparam_params is not None:
        #    uvparam=parameterization.new_from_params(implpart,uvparam_params)
        #    pass
        
        return cls(frame=frame,implpart=implpart,landmarks=landmarks)


    @classmethod
    def fromx3d(cls,frame,ndepartparams,filename,implclass=cadpart,tol=1e-4,recalcnormals=False,defaultappearance=None):

        (landmarks_params,implpartparams)=cls.splitparams(ndepartparams)
        landmarks=point_landmarks.fromparams(landmarks_params)
        
        implpart=implclass.fromx3d(filename=filename,implpartparams=implpartparams,tol=tol,recalcnormals=recalcnormals,defaultappearance=defaultappearance)

        #uvparam=None
        #if uvparam_params is not None:
        #    uvparam=parameterization.new_from_params(implpart,uvparam_params)
        #    pass
        
        return cls(frame=frame,implpart=implpart,landmarks=landmarks)

    
    @classmethod
    def from_implpart_generator(cls,frame,ndepartparams,generator):
        # generator is called with implpartparams
        
        (landmarks_params,implpartparams)=cls.splitparams(ndepartparams)
        landmarks=point_landmarks.fromparams(landmarks_params)
        
        implpart=generator(implpartparams)
        
        #uvparam=None
        #if uvparam_params is not None:
        #    uvparam=parameterization.new_from_params(implpart,uvparam_params)
        #    pass
        
        return cls(frame=frame,implpart=implpart,landmarks=landmarks)
    
    def X3DWrite(self,serializer,destframe,UVparameterization=None):
        self.implpart.X3DWrite(serializer,self.frame,destframe,UVparameterization=UVparameterization)

    def VRMLWrite(self,serializer,destframe,UVparameterization=None):
        self.implpart.VRMLWrite(serializer,self.frame,destframe,UVparameterization=UVparameterization)  

    #def VRMLFile(self,buf,destframe,UVparameterization=None,appearances=None):
    #    Shapes=self.VRMLShapes(destframe,UVparameterization=UVparameterization,appearances=appearances) # Appearance is instance of class VRMLAppearance
    #
    #    if isinstance(buf,basestring):
    #        buf=open(buf,"w")
    #        opened=True
    #        pass
    #    else:
    #        opened=False
    #        pass
    #    
    #    buf.write("#VRML V2.0 utf8\n")
    #    for Shape in Shapes:
    #        buf.write(Shape)
    #        pass
    #
    #    if opened:
    #        buf.close()
    #        pass
    #        
    #    pass

    def assign_appearances(self,appearances):

        self.implpart.assign_appearances(appearances)
        pass

    pass


class ndeassembly(object):
    # Not fully implemented
    parts=None
    
    def __init__(self,**kwargs):
        self.parts=[]
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
        
            setattr(self,kwarg,kwargs[kwarg])
            pass
        pass

    
    pass



