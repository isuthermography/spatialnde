import sys
import copy
import string
import io

from lxml import etree
import numpy as np

if sys.version_info >= (3,0):
    # python3
    translate=str.translate
    x3d_removecomma=str.maketrans(",", " ")
    pass
else:
    translate=string.translate
    x3d_removecomma=string.maketrans(","," ")
    pass

def vecnorm(a,axis=-1):
    
    # compatibility for old versions of Numpy
    return np.sqrt(np.sum(a**2.0,axis))


def read_mfstring(mfstring):
    Strings = [  ]
    StringBuf = None
    pos=0
    in_string = False

    if hasattr(mfstring,"decode"):
        # in py2.7, StringIO() wants unicode strings, not 'str'
        # so since py2.7 str() has a decode(), we call that 
        # to get a unicode object instead
        #... will also make this work on a bytes object.
        mfstring=mfstring.decode('utf-8')
        pass

    last_was_escape = False
    while pos < len(mfstring):
        if not in_string:
            if mfstring[pos]=='"':
                in_string=True
                StringBuf = io.StringIO()
                pass
            else:
                if not mfstring[pos].isspace():
                    raise ValueError("Invalid character %s in between MFString components (\'%s\')" % (mfstring[pos],mfstring))
                pass
            last_was_escape=False
            pass
        else:
            # We are in_string
            if mfstring[pos]=='"' and not last_was_escape:
                # End of the string
                in_string=False
                Strings.append(StringBuf.getvalue())
                StringBuf=None
                pass
            elif mfstring[pos]=='\\' and not last_was_escape:
                # Escape character
                last_was_escape=True
                pass
            elif (mfstring[pos]=='\\' or mfstring[pos]=='"') and last_was_escape:
                # Add escaped character
                StringBuf.write(mfstring[pos])
                last_was_escape=False
                pass
            elif last_was_escape:
                raise ValueError("Invalid escaped character %s in MFString \"%s\"" % (mfstring[pos],mfstring))
            else:
                # not last_was_escape and we have a regular character
                StringBuf.write(mfstring[pos])
                pass
            pass
        pos+=1
        pass

    if in_string:
        raise ValueError("Unterminated string in MFString \"%s\"" % (mfstring))

    return Strings

def ReadX3DIndexedVertexField_OBSOLETE(VertexIndex,numvertices=None,terminatorindices=None,missing_final_terminator=None):

    if numvertices is None:
        # need to determine numvertices, terminatorindices, missing_final_terminator
        polynum=0
        vertexnum=0
        numvertices=[]
        terminatorindices=[] # index into coordIndex of each "-1" polygon terminator
        # find all of the terminator indices, and the #'s of vertices too
        for cnt in range(VertexIndex.shape[0]):
            if VertexIndex[cnt]==-1:
                numvertices.append(vertexnum)
                polynum+=1
                vertexnum=0
                terminatorindices.append(cnt)
                continue
            vertexnum+=1
            pass
        if vertexnum > 0:
            # if coordIndex did not end with -1, need to
            # do cleanup here
            numvertices.append(vertexnum)
            polynum+=1
            missing_final_terminator=True
            pass
        else:
            missing_final_terminator=False
            pass
        numpolys = polynum
        pass
    else:
        numpolys = len(numvertices)
        pass
    
    maxnumvertices=np.max(numvertices)
    terminators_to_add = maxnumvertices - numvertices -1
    
    # add additional terminators to get rectangular matrix
    #vertexids=np.empty((numpolys,maxnumvertices),dtype='i4')
    
    VertexIndexList=list(VertexIndex)
    # do last one manually because of possible missingfinalterminator
    nterms=terminators_to_add[-1]
    if missing_final_terminator:
        nterms+=1
        pass
    if nterms > 0:
        VertexIndexList.extend([-1]*nterms)
        pass
    elif nterms < 0:
        VertexIndexList.pop()
        pass
    
    for cnt in range(numpolys-2,-1,-1):  # Go in reverse so indices are valid
        nterms=terminators_to_add[cnt]
        if nterms > 0:
            for cnt2 in range(nterms):
                VertexIndexList.insert(terminatorindices[cnt],-1)
                pass
            pass
        elif nterms < 0:
            # Remove existing '-1' terminator
            del VertexIndexList[terminatorindices[cnt]]
            pass
        
        pass
    # Now its a rectangular matrix so we can just
    # convert it back into an array and reshape it. 
    vertexids=np.array(VertexIndexList,dtype='i4').reshape(numpolys,maxnumvertices)
    assert(missing_final_terminator is not None)
    return (numvertices,terminatorindices,missing_final_terminator,vertexids)


# ReadX3DIndexedVertexField() converted to:
#   find_vertexidx_indices() in polygonalsurface.py



class x3d_transform(object):
    center=None
    rotation=None
    scale=None
    scaleOrientation=None
    translation=None
    bboxCenter=None
    bboxSize=None
    
    def __init__(self):
        self.center=np.array((0.0,0.0,0.0),dtype='d')
        self.rotation=np.array((0.0,0.0,1.0,0.0),dtype='d')
        self.scale=np.array((1.0,1.0,1.0),dtype='d')
        self.scaleOrientation=np.array((0.0,0.0,1.0,0.0),dtype='d')
        self.translation=np.array((0.0,0.0,0.0),dtype='d')
        self.bboxCenter=np.array((0.0,0.0,0.0),dtype='d')
        self.bboxSize=np.array((-1.0,-1.0,-1.0),dtype='d')
        
        pass

    def eval(self):
        # Evaluate transform as a 4x4 matrix

        
        T = np.matrix(((1.0,0.0,0.0,self.translation[0]),
                      (0.0,1.0,0.0,self.translation[1]),
                      (0.0,0.0,1.0,self.translation[2]),
                      (0.0,0.0,0.0,1.0)),dtype='d')
        C = np.matrix(((1.0,0.0,0.0,self.center[0]),
                      (0.0,1.0,0.0,self.center[1]),
                      (0.0,0.0,1.0,self.center[2]),
                      (0.0,0.0,0.0,1.0)),dtype='d')

        #import pdb
        #pdb.set_trace()


        ## WRONG: X3D rotation is not quaternion, it is (axis, angle)
        #
        ## http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/
        ## http://stackoverflow.com/questions/1556260/convert-quaternion-rotation-to-rotation-matrix
        ## normalize rotation quaternion
        #rotvecnorm=vecnorm(self.rotation)
        #
        #if rotvecnorm != 0.0:
        #    rotnorm = self.rotation/rotvecnorm
        #    pass
        #else:
        #    # assume null quaternion (no rotation)
        #    rotnorm = np.array((0.0,0.0,1.0,0.0),dtype='d')
        #    pass
        #
        ## convert rotation quaternion to matrix
        ##R = np.matrix(((1.0 - 2.0*rotnorm[1]**2.0 - 2.0*rotnorm[2]**2.0,
        #               2.0*rotnorm[0]*rotnorm[1] - 2.0*rotnorm[2]*rotnorm[3],
        #               2.0*rotnorm[0]*rotnorm[2] + 2.0*rotnorm[1]*rotnorm[3],0.0),
        #              (2.0*rotnorm[0]*rotnorm[1] + 2.0*rotnorm[2]*rotnorm[3],
        #               1.0 - 2.0*rotnorm[0]**2.0 - 2.0*rotnorm[2]**2.0,
        #               2.0*rotnorm[1]*rotnorm[2]-2.0*rotnorm[0]*rotnorm[3],0.0),
        #              (2.0*rotnorm[0]*rotnorm[2] - 2.0*rotnorm[1]*rotnorm[3],
        #               2.0*rotnorm[1]*rotnorm[2] + 2.0*rotnorm[0]*rotnorm[3],
        #               1.0-2.0*rotnorm[0]**2.0 - 2.0*rotnorm[1]**2.0,0.0),
        #              (0.0,0.0,0.0,1.0)),dtype='d')
        #
        #scaleorientnorm =self.scaleOrientation/vecnorm(self.scaleOrientation)
        #
        #SR = np.matrix(((1.0 - 2.0*scaleorientnorm[1]**2.0 - 2.0*scaleorientnorm[2]**2.0,
        #               2.0*scaleorientnorm[0]*scaleorientnorm[1] - 2.0*scaleorientnorm[2]*scaleorientnorm[3],
        #               2.0*scaleorientnorm[0]*scaleorientnorm[2] + 2.0*scaleorientnorm[1]*scaleorientnorm[3],0.0),
        #              (2.0*scaleorientnorm[0]*scaleorientnorm[1] + 2.0*scaleorientnorm[2]*scaleorientnorm[3],
        #               1.0 - 2.0*scaleorientnorm[0]**2.0 - 2.0*scaleorientnorm[2]**2.0,
        #               2.0*scaleorientnorm[1]*scaleorientnorm[2]-2.0*scaleorientnorm[0]*scaleorientnorm[3],0.0),
        #              (2.0*scaleorientnorm[0]*scaleorientnorm[2] - 2.0*scaleorientnorm[1]*scaleorientnorm[3],
        #               2.0*scaleorientnorm[1]*scaleorientnorm[2] + 2.0*scaleorientnorm[0]*scaleorientnorm[3],
        #               1.0-2.0*scaleorientnorm[0]**2.0 - 2.0*scaleorientnorm[1]**2.0,0.0),
        #              (0.0,0.0,0.0,1.0)),dtype='d')

        # Apply Rodrigues rotation formula to determine R and SR
        k = self.rotation[:3]
        ang = self.rotation[3]
        kmag = np.linalg.norm(k)
        if kmag == 0.0:
            kmag=1.0  # null rotation
            k=np.array((0.0,0.0,1.0),dtype='d')
            ang=0.0
            pass
        k /= kmag

        # cross product matrix
        RK=np.matrix(((0.0,-k[2],k[1]),
                      (k[2],0.0,-k[0]),
                      (-k[1],k[0],0.0)),dtype='d')
        # R=np.eye(3) + np.sin(ang)*RK + (1.0-np.cos(ang))*np.dot(RK,RK)
        R=np.concatenate((np.concatenate((np.eye(3) + np.sin(ang)*RK + (1.0-np.cos(ang))*np.dot(RK,RK),np.zeros((3,1),dtype='d')),axis=1),
                          np.array(((0.0,0.0,0.0,1.0),),dtype='d')),axis=0)
        


        # Apply Rodrigues rotation formula to determine scale orientation
        SOk = self.scaleOrientation[:3]
        
        SOang = self.scaleOrientation[3]
        SOkmag = np.linalg.norm(SOk)
        if SOkmag == 0.0:
            SOkmag=1.0  # null rotation
            SOk=np.array((0.0,0.0,1.0),dtype='d')
            SOang=0.0
            pass
        SOk /= SOkmag

        # cross product matrix
        SOK=np.matrix(((0.0,-SOk[2],SOk[1]),
                       (SOk[2],0.0,-SOk[0]),
                       (-SOk[1],SOk[0],0.0)),dtype='d')
        SR=np.concatenate((np.concatenate((np.eye(3) + np.sin(SOang)*SOK + (1.0-np.cos(SOang))*np.dot(SOK,SOK),np.zeros((3,1),dtype='d')),axis=1),
                          np.array(((0.0,0.0,0.0,1.0),),dtype='d')),axis=0)

        
        S = np.matrix(((self.scale[0],0.0,0.0,0.0),
                       (0.0,self.scale[1],0.0,0.0),
                       (0.0,0.0,self.scale[2],0.0),
                       (0.0,0.0,0.0,1.0)),dtype='d')

        # Transform components are defined as matrices not arrays,
        # so we can just multiply them
        
        return np.array((T * C * R * SR * S * (-SR) * (-C)),dtype='d')

    @classmethod
    def fromelement(cls,loader,element):
        transform=cls()

        if "center" in element.attrib:
            transform.center=np.fromstring(translate(element.attrib["center"],x3d_removecomma),sep=' \r\n',dtype='d')
            pass
        if "rotation" in element.attrib:
            transform.rotation=np.fromstring(translate(element.attrib["rotation"],x3d_removecomma),sep=' \r\n',dtype='d')
            pass
        if "scale" in element.attrib:
            transform.scale=np.fromstring(translate(element.attrib["scale"],x3d_removecomma),sep=' \r\n',dtype='d')
            pass
        if "scaleOrientation" in element.attrib:
            transform.scaleOrientation=np.fromstring(translate(element.attrib["scaleOrientation"],x3d_removecomma),sep=' \r\n',dtype='d')
            pass
        if "translation" in element.attrib:
            transform.translation=np.fromstring(translate(element.attrib["translation"],x3d_removecomma),sep=' \r\n',dtype='d')
            pass
        if "bboxCenter" in element.attrib:
            transform.bboxCenter=np.fromstring(translate(element.attrib["bboxCenter"],x3d_removecomma),sep=' \r\n',dtype='d')
            pass
        if "bboxSize" in element.attrib:
            transform.bboxSize=np.fromstring(translate(element.attrib["bboxSize"],x3d_removecomma),sep=' \r\n',dtype='d')
            pass
        
        
        
        
        # Handle attributes... our caller handles dispatching to children

        return transform
    pass

    

class x3d_indexedfaceset(object):
    # Reference: http://www.web3d.org/documents/specifications/19775-1/V3.3/Part01/components/geometry3D.html#IndexedFaceSet
    coord=None
    normal=None
    texCoord=None
    texCoordIndex=None
    ccw=None
    solid=None
    coordIndex=None
    normalIndex=None
    normalPerVertex=None
    transform=None
    metadata=None
    principal_curvatures=None
    curvature_tangent_axes=None
    
    def __init__(self):
        self.normalPerVertex=True
        self.ccw=True
        self.solid=True
        pass
    
    @classmethod
    def fromelement(cls,loader,element):
        ifs=cls()
        ifs.transform=copy.copy(loader.transformstack[-1])
        
        for attrname in element.attrib:
            if attrname=="coordIndex":
                ifs.coordIndex=np.fromstring(translate(element.attrib[attrname],x3d_removecomma),sep=' \r\n',dtype='i4')
                pass
            elif attrname=="normalIndex":
                ifs.normalIndex=np.fromstring(translate(element.attrib[attrname],x3d_removecomma),sep=' \r\n',dtype='i4')
                pass
            elif attrname=="texCoordIndex":
                ifs.texCoordIndex=np.fromstring(translate(element.attrib[attrname],x3d_removecomma),sep=' \r\n',dtype='i4')
                pass
            elif attrname=="ccw":
                ifs.ccw=element.attrib[attrname].lower()=="true"
                pass
            elif attrname=="solid":
                ifs.solid=element.attrib[attrname].lower()=="true"
                pass
            elif attrname=="normalPerVertex":
                ifs.normalPerVertex=element.attrib[attrname].lower()=="true"
                pass
            elif attrname=="USE" or attrname=="DEF" or attrname=="class" or attrname=="containerField":
                pass  # handled elsewhere
            else:
                sys.stderr.write("Unknown indexedfaceset attribute: %s\n" % (attrname))
                pass
            pass
                                 
            
        for child in element:
            loader.dispatch_x3d_childnode(child,ifs)
            pass

        return ifs
    pass


class x3d_imagetexture(object):
    metadata=None
    url=None
    repeatS=None
    repeatT=None
    textureProperties=None
    
    def __init__(self):
        self.repeatS=True
        self.repeatT=True
        pass

    @classmethod
    def fromelement(cls,loader,element):
        imagetexture=cls()

        if "repeatS" in element.attrib:
            if element.attrib["repeatS"]=="true":
                imagetexture.repeatS=True
                pass
            else:
                imagetexture.repeatS=False
                pass
            
            pass

        if "repeatT" in element.attrib:
            if element.attrib["repeatT"]=="true":
                imagetexture.repeatT=True
                pass
            else:
                imagetexture.repeatT=False
                pass
            pass

        if "url" in element.attrib:
            if not element.attrib["url"].strip().startswith('"'):
                # Not an MFString: Interpret it as URL directly
                imagetexture.url=element.attrib["url"]
                pass
            else:
                imagetexture.url=read_mfstring(element.attrib["url"])[0]
                pass
            pass
        
        for child in element:
            loader.dispatch_x3d_childnode(child,imagetexture)
            pass
        return imagetexture
    
    pass



class x3d_material(object):
    metadata=None
    ambientIntensity=None
    diffuseColor=None
    shininess=None
    specularColor=None
    transparency=None
    
    def __init__(self):
        self.ambientIntensity=0.2
        self.diffuseColor=(0.8,0.8,0.8)
        self.emissiveColor=(0.0,0.0,0.0)
        self.shininess=0.2
        self.specularColor=(0.0,0.0,0.0)
        self.transparency=0.0
        pass

    @classmethod
    def fromelement(cls,loader,element):
        material=cls()

        if "ambientIntensity" in element.attrib:
            material.ambientIntensity=float(element.attrib["ambientIntensity"])
            pass

        if "diffuseColor" in element.attrib:
            material.diffuseColor=[ float(num) for num in element.attrib["diffuseColor"].split()]
            pass

        if "emissiveColor" in element.attrib:
            material.emissiveColor=[ float(num) for num in element.attrib["emissiveColor"].split()]
            pass
        
        if "shininess" in element.attrib:
            material.shininess=float(element.attrib["shininess"])
            pass

        if "specularColor" in element.attrib:
            material.specularColor=[ float(num) for num in element.attrib["specularColor"].split()]
            pass

        if "transparency" in element.attrib:
            material.transparency=float(element.attrib["transparency"])
            pass
        
        for child in element:
            loader.dispatch_x3d_childnode(child,material)
            pass
        return material
    
    pass




class x3d_appearance(object):
    fillProperties=None
    lineProperties=None   # x3d_indexedfaceset
    material=None
    metadata=None
    shaders=None
    texture=None
    textureTransform=None
    
    def __init__(self):
        pass
    
    @classmethod
    def fromelement(cls,loader,element):
        appearance=cls()
        
        for child in element:
            loader.dispatch_x3d_childnode(child,appearance)
            pass
        
        return appearance
    pass


class x3d_shape(object):
    appearance=None
    geometry=None   # x3d_indexedfaceset
    metadata=None
    
    def __init__(self):
        pass
    
    @classmethod
    def fromelement(cls,loader,element):
        shape=cls()

        for child in element:
            loader.dispatch_x3d_childnode(child,shape)
            pass
        
        return shape
    pass

class x3d_indexedfaceset_loader(object):
    surfaces=None  # list of surfaces (x3d_shape)
    transformstack=None # stack (list) of accumulated transforms
    defindex=None # dictionary, by DEF value of referenceable XML elements
    NSPRE=None
    spatialnde_NSPRE=None
    metersperunit=None
    
    def __init__(self,**kwargs):
        self.surfaces=[]
        self.transformstack=[ np.eye(4,dtype='d') ]
        self.defindex={}
        self.NSPRE=None
        self.spatialnde_NSPRE="{http://spatialnde.org/x3d}"

        self.metersperunit=1.0
        
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
        
            setattr(self,kwarg,kwargs[kwarg])
            pass

        pass

    def parse_x3d_headchild(self,node):
        if node.tag==self.NSPRE+"unit":
            if node.attrib["category"]=="length":
                self.metersperunit=node.attrib["conversionFactor"]
                pass
            # We don't care about other kinds of units
            pass
        elif node.tag==self.NSPRE+"meta":
            pass # ignore
        else:
            sys.stderr.write("x3d_indexedfaceset_loader: unknown header tag: %s\n" % (node.tag))
            pass
        pass

    def parse_indexedfaceset(self,node,parentstruct,containerField):
        if containerField is None:
            containerField="geometry"
            pass
        
        if not hasattr(parentstruct,containerField):
            raise ValueError("Invalid containerField: %s" % (containerField))

        setattr(parentstruct,containerField,x3d_indexedfaceset.fromelement(self,node))
        pass

    def parse_coordinate(self,node,parentstruct,containerField):
        if containerField is None:
            containerField="coord"
            pass

        if not hasattr(parentstruct,containerField):
            raise ValueError("Invalid containerField: %s" % (containerField))
        
        point=None
        if "point" in node.attrib:
            coord=np.fromstring(translate(node.attrib["point"],x3d_removecomma),sep=' \r\n',dtype='d')*self.metersperunit

            coord=coord.reshape(coord.shape[0]//3,3)
            
            setattr(parentstruct,containerField,coord)
            #import pdb
            #pdb.set_trace()
            pass
        pass
    
    def parse_normal(self,node,parentstruct,containerField):
        if containerField is None:
            containerField="normal"
            pass

        if not hasattr(parentstruct,containerField):
            raise ValueError("Invalid containerField: %s" % (containerField))
        
        vector=None
        if "vector" in node.attrib:
            normal=np.fromstring(translate(node.attrib["vector"],x3d_removecomma),sep=' \r\n',dtype='d')
            normal=normal.reshape(normal.shape[0]//3,3)
            setattr(parentstruct,containerField,normal)
            pass
        
        pass

    def parse_texturecoordinate(self,node,parentstruct,containerField):
        if containerField is None:
            containerField="texCoord"
            pass

        if not hasattr(parentstruct,containerField):
            raise ValueError("Invalid containerField: %s" % (containerField))
        
        point=None
        if "point" in node.attrib:
            texcoord=np.fromstring(translate(node.attrib["point"],x3d_removecomma),sep=' \r\n',dtype='d')

            texcoord=texcoord.reshape(texcoord.shape[0]//2,2)
            
            setattr(parentstruct,containerField,texcoord)
            pass
        pass

    def parse_transform(self,node,parentstruct,containerField):
        # Don't forget to multiply by self.metersperunit
        # when appropriate (transforms...)
        transform=x3d_transform.fromelement(self,node).eval()
        self.transformstack.append(np.dot(self.transformstack[-1],transform))

        for child in node:
            self.dispatch_x3d_childnode(child,None)
            
            pass
        self.transformstack.pop()
        
        pass

    def parse_appearance(self,node,parentstruct,containerField):
        if containerField is None:
            containerField="appearance"
            pass
        
        if not hasattr(parentstruct,containerField):
            raise ValueError("Invalid containerField: %s" % (containerField))
        
        setattr(parentstruct,containerField,x3d_appearance.fromelement(self,node))
        
        pass
    
    def parse_imagetexture(self,node,parentstruct,containerField):
        if containerField is None:
            containerField="texture"
            pass
        
        if not hasattr(parentstruct,containerField):
            raise ValueError("Invalid containerField: %s" % (containerField))
        
        setattr(parentstruct,containerField,x3d_imagetexture.fromelement(self,node))
        
        pass


    # WARNING curvatures won't currently transform properly if
    # the indexedfaceset is under a transform node
    def parse_principalcurvatures(self,node,parentstruct,containerField):
        if containerField is None:
            containerField="principal_curvatures"
            pass

        if not hasattr(parentstruct,containerField):
            raise ValueError("Invalid containerField: %s" % (containerField))

        if "curvatures" in node.attrib:
            principalcurvatures = np.fromstring(translate(node.attrib["curvatures"],x3d_removecomma),sep=' \r\n',dtype='d')/self.metersperunit
            principalcurvatures=principalcurvatures.reshape(principalcurvatures.shape[0]//2,2)
            setattr(parentstruct,containerField,principalcurvatures)
            pass
        
        
        pass

    def parse_curvaturetangentaxes(self,node,parentstruct,containerField):
        if containerField is None:
            containerField="curvature_tangent_axes"
            pass
        
        if not hasattr(parentstruct,containerField):
            raise ValueError("Invalid containerField: %s" % (containerField))


        if "axes" in node.attrib:
            curvature_tangent_axes = np.fromstring(translate(node.attrib["axes"],x3d_removecomma),sep=' \r\n',dtype='d')
            curvature_tangent_axes=curvature_tangent_axes.reshape(curvature_tangent_axes.shape[0]//6,2,3)
            setattr(parentstruct,containerField,curvature_tangent_axes)
            pass
        
        pass

    def parse_material(self,node,parentstruct,containerField):
        if containerField is None:
            containerField="material"
            pass
        
        if not hasattr(parentstruct,containerField):
            raise ValueError("Invalid containerField: %s" % (containerField))
        
        setattr(parentstruct,containerField,x3d_material.fromelement(self,node))
        
        pass

    
    def parse_shape(self,node):
        self.surfaces.append(x3d_shape.fromelement(self,node))
        pass
    
    def dispatch_x3d_childnode(self,child,parentstruct):
        containerField=None

        if "containerField" in child.attrib:
            containerField=child.attrib["containerField"]
            pass

        # NOTE BUG: Since USE/DEF are defined at the XML node
        # level rather than at our output level,
        # we end up with multiple copies of e.g. multiply
        # referenced appearance or texture nodes, rather
        # than a single multiply-referenced copy
        
        if "USE" in child.attrib:
            child = self.defindex[child.attrib["USE"]]
            pass

        if "DEF" in child.attrib:
            self.defindex[child.attrib["DEF"]]=child
            pass
        
        if child.tag==self.NSPRE+'x3d' or child.tag==self.NSPRE+'X3D':
            self.parse_x3d_tree(child)
            pass
        elif child.tag==self.NSPRE+'head':
            for grandchild in child:
                self.parse_x3d_headchild(grandchild)
                pass
            pass
        elif child.tag==self.NSPRE+"Scene" or child.tag==self.NSPRE+"scene":
            for grandchild in child:
                self.dispatch_x3d_childnode(grandchild,None)
                pass
            pass
        elif child.tag==self.NSPRE+"Shape" or child.tag==self.NSPRE+"shape":
            self.parse_shape(child)
            pass
        elif child.tag==self.NSPRE+"IndexedFaceSet" or child.tag==self.NSPRE+"indexedfaceset":
            self.parse_indexedfaceset(child,parentstruct,containerField)
            pass
        elif child.tag==self.NSPRE+"Coordinate" or child.tag==self.NSPRE+"coordinate":
            self.parse_coordinate(child,parentstruct,containerField)
            pass
        elif child.tag==self.NSPRE+"Normal" or child.tag==self.NSPRE+"normal":
            self.parse_normal(child,parentstruct,containerField)
            pass
        elif child.tag==self.NSPRE+"TextureCoordinate" or child.tag==self.NSPRE+"texturecoordinate":
            self.parse_texturecoordinate(child,parentstruct,containerField)
            pass
        elif child.tag==self.NSPRE+"Appearance" or child.tag==self.NSPRE+"appearance":
            self.parse_appearance(child,parentstruct,containerField)
            pass
        elif child.tag==self.NSPRE+"ImageTexture" or child.tag==self.NSPRE+"imagetexture":
            self.parse_imagetexture(child,parentstruct,containerField)
            pass
        elif child.tag==self.NSPRE+"Material" or child.tag==self.NSPRE+"material":
            self.parse_material(child,parentstruct,containerField)
            pass
        elif child.tag==self.NSPRE+"Transform" or child.tag==self.NSPRE+"transform":
            self.parse_transform(child,parentstruct,containerField)
            pass
        elif child.tag==self.NSPRE+"Group" or child.tag==self.NSPRE+"group":
            for grandchild in child:
                self.dispatch_x3d_childnode(grandchild,None)
                pass
            pass
        elif child.tag==self.spatialnde_NSPRE+"PrincipalCurvatures"  or child.tag==self.spatialnde_NSPRE+"principalcurvatures":
            self.parse_principalcurvatures(child,parentstruct,containerField)
            pass
        elif child.tag==self.spatialnde_NSPRE+"CurvatureTangentAxes"  or child.tag==self.spatialnde_NSPRE+"curvaturetangentaxes":
            self.parse_curvaturetangentaxes(child,parentstruct,containerField)
            pass
        else:
            sys.stderr.write("x3d_indexedfaceset_loader:  Non-implemented tag %s\n" % (child.tag))
            pass
        pass
    
    def parse_x3d_tree(self,node):

        if self.NSPRE is None:
            # use X3D element tag to determine what namespace to look for
            tag=node.tag
            namespace=None
            self.NSPRE=""
            if tag.startswith('{'):  # has namespace prefix
                namespace=tag[1:tag.index('}')]
                self.NSPRE=tag[:(tag.index('}')+1)]
                tag=tag[(tag.index('}')+1):]
                pass
            pass
        
        
        if tag != 'x3d' and tag != 'X3D':
            raise KeyError("Unknown root element %s in X3D tree" % (node.tag))
        
        for child in node:
            self.dispatch_x3d_childnode(child,None)
            
            pass
        pass

        
        
    @classmethod
    def parsefile(cls,filename):

        # Must create custom parser so that it won't choke
        # on huge files
        hugefileparser=etree.XMLParser(huge_tree=True)
        
        tree=etree.parse(filename,parser=hugefileparser)
        rootnode=tree.getroot()
        
        loader=cls()
        loader.parse_x3d_tree(rootnode)

        return loader
    pass



