import copy
import numpy as np
import posixpath
import ast
import os
import sys
from io import BytesIO

try:
    # py2.x
    from urllib import pathname2url
    from urllib import url2pathname
    from urllib import quote
    from urllib import unquote
    from urlparse import urlparse
    from urlparse import urlunparse
    from urlparse import urljoin    
    pass
except ImportError:
    # py3.x
    from urllib.request import pathname2url
    from urllib.request import url2pathname
    from urllib.parse import quote
    from urllib.parse import unquote
    from urllib.parse import urlparse
    from urllib.parse import urlunparse
    from urllib.parse import urljoin
    pass

import dataguzzler as dg
import dg_file as dgf
import dg_image
import dg_metadata as dgm
from spatialnde.coordframes import coordframe
from spatialnde.cadpart import appearance
from spatialnde.ndeobj import ndepart,ndeassembly
from spatialnde.exporters.x3d import X3DSerialization,X3DMultiFrameSerialization
from limatix import dc_value

def run(_xmldoc,
        _element,
        _uniquematches,
        _dest_href,
        dc_specimen_str,
        dc_fused_greensinversion_dgsfile_href,
        dc_greensinversion_offset_numericunits,
        dc_greensinversion_scaling_numericunits):
    tikparam=_xmldoc.xpathsinglecontext(_element,"string(dc:fused_greensinversion_dgsfile/@tikparam)")
    
    (junkmd,wfmdict)=dgf.loadsnapshot(dc_fused_greensinversion_dgsfile_href.getpath(),memmapok=True)

    measurements=_xmldoc.xpathcontext(_uniquematches,"..")

    dc_inversion_channel_str=_xmldoc.xpathsinglecontext(measurements[0],"string(dc:inversion_channel)")

    gi_projchan = "Projgreensinversion_" + dc_inversion_channel_str[:-4] # [:-4] strips off trailing _tex
    texchanprefix=""

    while len(dgm.GetMetaDatumWIStr(wfmdict[gi_projchan],"X3DGeomRef","")) > 0:
        texchanprefix+=dgm.GetMetaDatumWIStr(wfmdict[gi_projchan],"TexChanPrefix","")
        gi_projchan = dgm.GetMetaDatumWIStr(wfmdict[gi_projchan],"X3DGeomRef","")
    
        pass

    # dc:inversion_channel_str is a string representing the original dgs channel that was inverted
    gi_texchan="greensinversion_" + dc_inversion_channel_str # Name of the parameterization (texture map) channel
    
    # Read geometry from metadata into an ndepart object
    X3DGeom = dgm.GetMetaDatumWIStr(wfmdict[gi_projchan],"X3DGeom","")
    X3DGeom_fh = BytesIO(X3DGeom.encode('utf-8'))
                         
    objframe=coordframe() 
    X3DObj = ndepart.fromx3d(objframe,None,X3DGeom_fh,tol=1e-6)
    X3DGeom_fh.close()

    # assume single surface
    assert(len(X3DObj.implpart.surfaces)==1)

    # Clear curvature data -- Not needed for simple rendering
    X3DObj.implpart.surfaces[0].clear_curvatures()

    # Check existing appearance tag -- should be a texture_url that matches our gi_texchan
    assert(texchanprefix + X3DObj.implpart.surfaces[0].appearance.texture_url[1:] == gi_texchan)

    # Blank out texture_url (will be rewritten by Javascript
    # in X3DMultiFrameSerialization
    new_obj_appearance = appearance.texture_url(texture_url=None)
    X3DObj.assign_appearances(new_obj_appearance)
    

    dc_inversion_reflectors_str=_xmldoc.xpathsinglecontext(measurements[0],"string(dc:inversion_reflectors)")


    numframes = wfmdict[gi_texchan].dimlen[2]

    approx_dpi_x = 1.0/(dc_value.numericunitsvalue(dgm.GetMetaDatumWIDbl(wfmdict[gi_texchan],"Step1","1.0"),dgm.GetMetaDatumWIStr(wfmdict[gi_texchan],"Units1","m")).value("m")/.0254)
    approx_dpi_y = 1.0/(dc_value.numericunitsvalue(dgm.GetMetaDatumWIDbl(wfmdict[gi_texchan],"Step2","1.0"),dgm.GetMetaDatumWIStr(wfmdict[gi_texchan],"Units2","m")).value("m")/.0254)
    
    
    reflectors=ast.literal_eval(dc_inversion_reflectors_str)


    barename = posixpath.splitext(dc_fused_greensinversion_dgsfile_href.get_bare_unquoted_filename())[0]


    # Should extract .x3d file, purge curvature, and rewrite with different texture URL's for each frame. 

    # output directory for layer images
    layersdir_href = dc_value.hrefvalue(quote(barename+"_layers/"),contexthref=dc_fused_greensinversion_dgsfile_href.leafless())

    # make sure directory exists
    try: 
        os.mkdir(layersdir_href.getpath())
        pass
    except OSError:
        pass




    retval=[ (("dc:greensinversion_layersdir",{"tikparam": str(tikparam)}),layersdir_href) ]

    layer_depths=[]

    layer_descrs=[]
    layer_hrefs=[]

    for frameno in range(numframes):
        if frameno==0:
            layer_depths.append(0.0)
            pass
        else:
            # Note that reflectors array is deepest first, and doesn't include the front surface
            assert(len(reflectors)==numframes-1)
            layer_depths.append(reflectors[len(reflectors)-frameno][0])
            pass

        img = dg_image.toimage(wfmdict[gi_texchan],wfmdict,(frameno,),unitsperintensity=dc_greensinversion_scaling_numericunits.value("J/m^2"),offset=dc_greensinversion_offset_numericunits.value("J/m^2"),colormap="hot")
        
        frameout_href = dc_value.hrefvalue(quote(barename+"_layer%2.2d.png" % (frameno)),contexthref=layersdir_href)
        img.save(frameout_href.getpath(),dpi=(approx_dpi_x,approx_dpi_y))

        #x3dout_href = dc_value.hrefvalue(quote(barename+"_layer%2.2d.x3d" % (frameno)),contexthref=layersdir_href)
        
        #new_obj_appearance = appearance.texture_url(texture_url=frameout_href.attempt_relative_url(x3dout_href.value()))

        #X3DObj.assign_appearances(new_obj_appearance)
        #x3dwriter=X3DSerialization.tofileorbuffer(x3dout_href.getpath(),x3dnamespace=None)
        #X3DObj.X3DWrite(x3dwriter,objframe)
        #x3dwriter.finish()

        layer_hrefs.append(frameout_href)
        layer_descrs.append("Depth = %f mm" % (layer_depths[frameno]*1000.0))
        
        retval.append((("dc:greensinversion_layer",{"tikparam": str(tikparam), "layernum": str(frameno), "layerdepth_meters": str(layer_depths[frameno]) }),frameout_href))
        #retval.append((("dc:greensinversion_layer_3d",{"tikparam": str(tikparam), "layernum": str(frameno), "layerdepth_meters": str(layer_depths[frameno]) }),x3dout_href))


        pass

    # Find some vaguely centroidy location from a weighted average of the
    # centers of the polygons
    centroidy = np.sum(X3DObj.implpart.surfaces[0].refpoints*X3DObj.implpart.surfaces[0].maxradius[:,np.newaxis]**2,axis=0)/np.sum(X3DObj.implpart.surfaces[0].maxradius[:,np.newaxis]**2)

    x3dout_href = dc_value.hrefvalue(quote(barename+".x3d"),contexthref=dc_fused_greensinversion_dgsfile_href.leafless())
    
    x3dwriter=X3DMultiFrameSerialization.tofileorbuffer(x3dout_href.getpath())
    X3DObj.X3DWrite(x3dwriter,objframe)
    layer_texture_urls = [ layer_href.attempt_relative_url(x3dout_href.value()) for layer_href in layer_hrefs ]
    x3dwriter.set_textures(layer_texture_urls,layer_descrs)
    x3dwriter.set_rotation_center(centroidy)
    x3dwriter.finish()

    retval.append((("dc:greensinversion_3d",{"tikparam": str(tikparam)}),x3dout_href))

    return retval
