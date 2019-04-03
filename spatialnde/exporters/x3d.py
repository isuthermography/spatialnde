import copy
from lxml import etree

from . import serialization,escape_mfstring


class X3DSerialization(serialization):
    NSPRE=None
    spatialnde_NSPRE=None
    NSMAP=None
    tree=None
    root=None
    scene=None # Scene is what the various IndexedFaceSets get appended to
    buf=None
    
    def __init__(self,**kwargs):
        self.UseCnt={}  # dictionary of use counters for named elements

        for kwarg in kwargs:
            if hasattr(self,kwarg):
                setattr(self,kwarg,kwargs[kwarg])
                pass
            else:
                raise IndexError(kwarg)
            pass
        pass

    def finish(self):
        # Prevent self-closing tags for better compatibility with more .x3d readers  .. by giving each element with no text, blank text
        for element in self.tree.iter():
            if element.text is None:
                element.text=""
                pass
            pass
        
        self.tree.write(self.buf,pretty_print=True,encoding='utf-8',xml_declaration=True)
        pass
    

    @classmethod
    def tofileorbuffer(cls,buf,x3dnamespace="http://www.web3d.org/specifications/x3d"):
        NSPRE=""
        NSMAP = {}
        spatialnde_NSPRE="{http://spatialnde.org/x3d}"
        if x3dnamespace is not None:
            NSPRE="{"+x3dnamespace+"}"
            NSMAP[None]=x3dnamespace
            pass

        NSMAP["sndex3d"] = spatialnde_NSPRE[1:-1]

        
        root=etree.Element(NSPRE+"X3D",nsmap=NSMAP)
        root.attrib["version"]="3.3"
        scene=etree.Element(NSPRE+"Scene",nsmap=NSMAP)
        Background=etree.Element(NSPRE+"Background")
        Background.attrib["skyColor"]=".1 .15 .3"
        
        root.append(scene)
        tree=etree.ElementTree(root)
        
        return cls(buf=buf,NSPRE=NSPRE,spatialnde_NSPRE=spatialnde_NSPRE,NSMAP=NSMAP,root=root,scene=scene,tree=tree)
    pass



class X3DMultiFrameSerialization(serialization):
    """This represents a multi-frame script-enabled .x3d file with a single ImageTexture that 
    can be selected with buttons.


    To use: 
      Define the object appearance to be a texture_url with texture_url=None, 
    e.g.
    obj_appearance = appearance.texture_url(texture_url=None)
    X3DObj.assign_appearances(obj_appearance)

    Then Instatiate this object, giving it a file name to write to 
    x3dwriter=X3DMultiFrameSerialization.tofileorbuffer(x3dout_href.getpath())

    Write your object to the writer
    X3DObj.X3DWrite(x3dwriter,objframe)

    Create and assign urls and layer descriptions (lists)
    layer_texture_urls = [ layer_href.attempt_relative_url(x3dout_href.value()) for layer_href in layer_hrefs ]
    x3dwriter.set_textures(layer_texture_urls,layer_descrs)

    Optionally set the rotation center
    x3dwriter.set_rotation_center(centroidy)

    ... and actually write it out
    x3dwriter.finish()

    """
    NSPRE=None
    spatialnde_NSPRE=None
    NSMAP=None
    tree=None
    root=None
    scene=None # Scene is what the various IndexedFaceSets get appended to
    buf=None
    
    template=r"""<?xml version="1.0"?>
    <X3D xmlns:sndex3d="http://spatialnde.org/x3d" profile="Immersive" version="3.3">
      <Scene>
        <NavigationInfo headlight="false"
                        visibilityLimit="0.0"
                        type='"EXAMINE", "ANY"'
                        avatarSize="0.25, 1.75, 0.75"/>
        <Background skyColor=".1 .15 .3"/>
        <Transform DEF="Object"> <!-- Transform tag will have a translation attribute -->
          <Shape/> <!-- This "Shape" tag will be what we specify as the Scene class attribute and will be what gets stuff added to it -->
        </Transform>

        <!-- we use a "proximity sensor as a hack so that we can catch moves and keep our 
             control panel ("HUD") at constant position -->
        <ProximitySensor DEF="PROXIMITY_SENSOR" size="1000000 1000000 1000000"/>
        <Transform DEF="HUD"> <!-- orientation of the HUD relative to world, javascript-extracted from ProximitySensor -->
          <Transform translation="4.5 0 -10"> <!-- translation of HUD on screen -->
            <Transform translation="0 3 0"> <!-- Translation of "Next layer" text -->
              <Shape> 
                <Text string='"Next layer"'>
                  <FontStyle justify='"MIDDLE"' size="0.5"/>
                </Text>
                <Appearance>
                  <Material diffuseColor="0.5 0.5 0.5"/>
                </Appearance>
              </Shape>
            </Transform>
            <Transform translation="0 2 0"> <!-- Translation of "Next layer" button -->
              <TouchSensor DEF="IncreaseButton"/>  <!-- Add "description" attribute? -->
              <Shape>                  
                <Disk2D outerRadius="0.5" innerRadius="0"/>
                <Appearance>
                  <Material diffuseColor="0.5 0.5 0.5" emissiveColor="0.5 0.5 0.5"/>
                </Appearance>
              </Shape>
            </Transform>
            <Transform translation="0 0.3 0"> <!-- Translation of "current layer" label -->
              <Shape>                  
                <Text string='"Current Layer:"'>
                  <FontStyle justify='"MIDDLE"' size="0.5"/>
                </Text>
                <Appearance>
                  <Material diffuseColor="0.5 0.5 0.5"/>
                </Appearance>
              </Shape>
            </Transform>
            <Transform translation="0 -0.5 0"> <!-- Translation of "current layer" readout -->
              <Shape>                  
                <Text DEF="cur_descr" string='""'>
                  <FontStyle justify='"MIDDLE"' size="0.5"/>
                </Text>
                <Appearance>
                  <Material diffuseColor="0.5 0.5 0.5"/>
                </Appearance>
              </Shape>
            </Transform>
            <Transform translation="0 -2 0"> <!-- Translation of "Previous layer" button -->
              <TouchSensor DEF="DecreaseButton"/>  <!-- Add "description" attribute? -->
              <Shape>                  
                <Disk2D outerRadius="0.5" innerRadius="0"/>
                <Appearance>
                  <Material diffuseColor="0.5 0.5 0.5" emissiveColor="0.5 0.5 0.5"/>
                </Appearance>
              </Shape>
            </Transform>
            <Transform translation="0 -3 0"> <!-- Translation of "Previous layer" text -->
              <Shape> 
                <Text string='"Previous Layer"'>
                  <FontStyle justify='"MIDDLE"' size="0.5"/>
                </Text>
                <Appearance>
                  <Material diffuseColor="0.5 0.5 0.5"/>
                </Appearance>
              </Shape>
            </Transform>
          </Transform>
        </Transform>
        <Script DEF='RotationCenter'>
          <!-- This script is here so that if you want, the interactive view
               rotation can be centered around a different point from the
               origin encoded in the underlying data. The object will 
               be offset only after initialization javascript is run -->
          <!-- The desired center of rotation should be stored in this field: -->
          <field accessType='inputOutput' type='SFVec3f' name='rotation_center' value='0 0 0'/>
                        
          <field accessType='outputOnly' type='SFVec3f' name='minus_rotation_center_out'/>
          <![CDATA[
            ecmascript:
              function initialize() {
                minus_rotation_center_out[0] = -rotation_center[0];
                minus_rotation_center_out[1] = -rotation_center[1];
                minus_rotation_center_out[2] = -rotation_center[2];
              }
          ]]>
        </Script>
        <Script DEF='Switch'>
          <field accessType='inputOutput' type='SFInt32' name='index' value='0'/>
          <field accessType='inputOnly' type='SFBool' name='increase'/>
          <field accessType='inputOnly' type='SFBool' name='decrease'/>
          <field accessType='outputOnly' type='MFString' name='new_url'/>
          <field accessType='outputOnly' type='MFString' name='new_descr'/>
          <field accessType='outputOnly' type='SFString' name='up_descr'/>
          <field accessType='outputOnly' type='SFString' name='down_descr'/>
          <field accessType='inputOutput' type="MFString" name="urls"/> <!-- Set value attribute of descrs to MFString array of urls -->
          <field accessType='inputOutput' type="MFString" name="descrs"/> <!-- Set value attribute of descrs to MFString array of descriptions -->
            <![CDATA[
                ecmascript: 
                        function set_url_descrs() {
                            if (index == 0) {
                                down_descr = "";
                            } else {
                                down_descr = descrs[index];
                            }

                            new_descr = new MFString([ descrs[index] ]);
                            new_url = new MFString([urls[index]]);

                            if (index == descrs.length-1) {
                                up_descr = "";
                            } else {
                                up_descr = descrs[index+1];
                            }


                        }
               
                        function initialize() {
                            set_url_descrs();
                        }

                        function increase (value) { // value is boolean which indicates button pressed 
                            if (value && index < urls.length-1) {
                                index = index + 1;
                                set_url_descrs();
                            }
                        }

                        function decrease (value) { // value is boolean which indicates button pressed 
                            if (value && index > 0) {
                                index = index - 1;
                                set_url_descrs();
                            }
                        }
                                
            ]]>
        </Script>
        <ROUTE fromNode='PROXIMITY_SENSOR' fromField='position_changed' toNode='HUD' toField='translation'/>
        <ROUTE fromNode='PROXIMITY_SENSOR' fromField='orientation_changed' toNode='HUD' toField='rotation'/>
        <ROUTE fromNode='IncreaseButton' fromField='isActive' toNode='Switch' toField='increase'/>
        <ROUTE fromNode='DecreaseButton' fromField='isActive' toNode='Switch' toField='decrease'/>
        <ROUTE fromNode='Switch' fromField='new_url' toNode='Texture' toField='url'/>
        <ROUTE fromNode='Switch' fromField='up_descr' toNode='IncreaseButton' toField='description'/>
        <ROUTE fromNode='Switch' fromField='down_descr' toNode='DecreaseButton' toField='description'/>
        <ROUTE fromNode='Switch' fromField='new_descr' toNode='cur_descr' toField='string'/>
        <ROUTE fromNode='RotationCenter' fromField='minus_rotation_center_out' toNode='Object' toField='translation'/>

      </Scene>
    </X3D>
    """


    def __init__(self,**kwargs):
        self.UseCnt={}  # dictionary of use counters for named elements

        for kwarg in kwargs:
            if hasattr(self,kwarg):
                setattr(self,kwarg,kwargs[kwarg])
                pass
            else:
                raise IndexError(kwarg)
            pass
        pass

    def set_rotation_center(self,rotation_center):
        # Look up field element
        assert(len(rotation_center)==3)
        rotation_center_field = self.tree.xpath("Scene/Script[@DEF='RotationCenter']/field[@name='rotation_center']",namespaces=self.NSMAP)[0]

        # Assign value attribute
        rotation_center_field.attrib["value"]="%.16g %.16g %.16g" % (rotation_center[0],rotation_center[1],rotation_center[2])
        
        pass


    def set_textures(self,texture_url_array, description_array):
        if len(texture_url_array) != len(description_array):
            raise ValueError("texture and description array length mismatch: %d vs. %d" % (len(texture_url_array),len(description_array)))

        # Look up field elements
        urls_field = self.tree.xpath("Scene/Script[@DEF='Switch']/field[@name='urls']",namespaces=self.NSMAP)[0]
        descrs_field = self.tree.xpath("Scene/Script[@DEF='Switch']/field[@name='descrs']",namespaces=self.NSMAP)[0]

        quoted_texture_urls = [ escape_mfstring(texture_url) for texture_url in texture_url_array ]
        quoted_descriptions = [ escape_mfstring(description) for description in description_array ]

        # Set value attributes to MFStrings
        urls_field.attrib["value"] = " ".join(quoted_texture_urls)
        descrs_field.attrib["value"] = " ".join(quoted_descriptions)
        pass

    def finish(self):
        # Mark first ImageTexture node with DEF="Texture" attribute, so that it can be ROUTE'd to 
        ImageTextureNode = self.tree.xpath("//ImageTexture",namespaces=self.NSMAP)[0]
        ImageTextureNode.attrib["DEF"]="Texture"


        # Prevent self-closing tags for better compatibility with more .x3d readers  .. by giving each element with no text, blank text
        for element in self.tree.iter():
            if element.text is None:
                element.text=""
                pass
            pass
        
        self.tree.write(self.buf,pretty_print=True,encoding='utf-8',xml_declaration=True)
        pass
    

    @classmethod
    def tofileorbuffer(cls,buf):
        spatialnde_NSPRE="{http://spatialnde.org/x3d}"

        x3dnamespace=None #"http://www.web3d.org/specifications/x3d"
        NSPRE="" # "{"+x3dnamespace+"}"
        
        parser = etree.XMLParser(strip_cdata=False) # use a parser that retains CDATA sections, as at least some X3D viewers treat CDATA in <Script> tags specially
        root=etree.XML(cls.template,parser)
        NSMAP=root.nsmap
        tree=etree.ElementTree(root)

        scene=tree.xpath("Scene/Transform[@DEF='Object']",namespaces=NSMAP)[0]
        
        return cls(buf=buf,NSPRE=NSPRE,spatialnde_NSPRE=spatialnde_NSPRE,NSMAP=NSMAP,root=root,scene=scene,tree=tree)
    pass


