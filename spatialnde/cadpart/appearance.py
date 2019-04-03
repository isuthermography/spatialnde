import sys
import copy
from lxml import etree

from ..exporters import escape_mfstring

class appearance(object):
    # Abstract class
    
    pass

class vrml_x3d_appearance(appearance):
    DefName=None
    CrossReference=None

    def generate_vrml_element(self,serialization):
        raise ValueError("vrml_x3d_appearance class is partially abstract. You must instantiate a subclass")

    def generate_x3d_element(self,serialization):
        raise ValueError("vrml_x3d_appearance class is partially abstract. You must instantiate a subclass")

    def generate_vrml(self,serialization):
        if self.DefName is not None:
            if self.DefName not in serialization.UseCnt:
                serialization.UseCnt[self.DefName]=0
                pass
            serialization.UseCnt[self.DefName]+=1
            
            if self.CrossReference:
                if serialization.UseCnt[self.DefName]==1:
                    # generate_vrml_element is implemented by our
                    # subclass which defines what type (image, material, etc.)
                    # of appearance we will generate
                    ret = self.generate_vrml_element(serialization)
                    return b"appearance DEF %s Appearance {\n%s\n}\n" % (self.DefName,ret)
                else:                    
                    return b"appearance USE %s\n" % (self.DefName)
                pass
            else:
                # generate_vrml_element is implemented by our
                # subclass which defines what type (image, material, etc.)
                # of appearance we will generate
                ret = self.generate_vrml_element(serialization)
                if serialization.UseCnt[self.DefName]==1:
                    return b"appearance DEF %s Appearance {\n%s\n}\n" % (self.DefName,ret)
                else:
                    return b"appearance DEF %s+%d Appearance {\n%s\n}\n" % (self.DefName,serialization.UseCnt[self.DefName]-1,ret)
                pass
            pass
        else:
            # generate_vrml_element is implemented by our
            # subclass which defines what type (image, material, etc.)
            # of appearance we will generate
            ret = self.generate_vrml_element(serialization)
            return b"appearance Appearance {\n%s\n}\n" % (ret)
        
        pass

    def generate_x3d(self,serialization):
        if self.DefName is not None:
            if self.DefName not in serialization.UseCnt:
                serialization.UseCnt[self.DefName]=0
                pass
            serialization.UseCnt[self.DefName]+=1

            if self.CrossReference:
                if serialization.UseCnt[self.DefName]==1:
                    # generate_x3d_element is implemented by our
                    # subclass which defines what type (image, material, etc.)
                    # of appearance we will generate
                    ret = self.generate_x3d_element(serialization)
                    ret.attrib["DEF"]=self.DefName
                    return ret
                else:                    
                    ret=etree.Element(serialization.NSPRE+"Appearance")
                    ret.attrib["USE"]=self.DefName
                    return ret
                pass
            else:
                # generate_x3d_element is implemented by our
                # subclass which defines what type (image, material, etc.)
                # of appearance we will generate
                ret = self.generate_x3d_element(serialization)
                if serialization.UseCnt[self.DefName]==1:
                    ret.attrib["DEF"]=self.DefName
                    return ret
                else:
                    ret.attrib["DEF"]="%s+%d" % (self.DefName,serialization.UseCnt[self.DefName]-1)
                    return ret
                pass
            pass
                        
        else:
            ret = self.generate_x3d_element(serialization)
            return ret
        pass
    
    pass

class simple_material(vrml_x3d_appearance):
    diffuseColor=None

    def __init__(self,**kwargs):
        self.CrossReference=True
        for kwarg in kwargs:
            if hasattr(self,kwarg):
                setattr(self,kwarg,kwargs[kwarg])
                pass
            else:
                raise IndexError(kwarg)
            pass
        pass

    def generate_x3d_element(self,serialization):
        ret=etree.Element(serialization.NSPRE+"Appearance")
        material=etree.Element(serialization.NSPRE+"Material")
        material.attrib["diffuseColor"]="%g %g %g" % (self.diffuseColor[0],self.diffuseColor[1],self.diffuseColor[2])
        ret.append(material)
        return ret
        
    def generate_vrml_element(self,serialization):
        return b"material Material { diffuseColor %g %g %g }\n" % (self.diffuseColor[0],self.diffuseColor[1],self.diffuseColor[2])
    
    
    
    @classmethod
    def from_color(cls,diffuseColor,DefName=None):
        return cls(diffuseColor=diffuseColor,DefName=DefName)
    pass



class texture_url(vrml_x3d_appearance):
    texture_url=None

    def __init__(self,**kwargs):
        self.CrossReference=True
        for kwarg in kwargs:
            if hasattr(self,kwarg):
                setattr(self,kwarg,kwargs[kwarg])
                pass
            else:
                raise IndexError(kwarg)
            pass
        pass


    def generate_x3d_element(self,serialization):
        ret=etree.Element(serialization.NSPRE+"Appearance")
        imagetex=etree.Element(serialization.NSPRE+"ImageTexture")
        if self.texture_url is None:
            imagetex.attrib["url"] = ""
            pass
        else:
            imagetex.attrib["url"] = escape_mfstring(self.texture_url)
            pass
        ret.append(imagetex)
        return ret
    
    def generate_vrml_element(self,serialization):
        if self.texture_url is None:
            return b"texture ImageTexture { url [ ] }\n"
        else:
            return b"texture ImageTexture { url [ %s ] }\n" % (escape_mfstring(self.texture_url))
        pass
    
    
    @classmethod
    def from_url(cls,texture_url,DefName=None):
        return cls(texture_url=texture_url,DefName=DefName)
    pass

