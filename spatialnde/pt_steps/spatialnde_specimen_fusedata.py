import copy
import numpy as np
import posixpath

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
import dg_metadata as dgm
from limatix import dc_value

from spatialnde.coordframes import coordframe
from spatialnde.dataguzzler.dg_3d import ndepart_from_dataguzzler_wfm


def run(_xmldoc,_element,_uniquematches,_dest_href,dc_specimen_str):
    # _element is the dc:fusion element with a 
    # <dc:specimen> sub-element
    # _uniquematches are the various dc:specimen elements
    # that matched the key (their parents are the 
    # dc:measurement tags)
    # dc_specimen_str is a string representing the 
    # specimen ID
    # dc:inversion_channel_str is a string representing the dgs channel
    # name for the original .dgs file's texture map
    
    measurements=_xmldoc.xpathcontext(_uniquematches,"..")

    dc_inversion_channel_str=_xmldoc.xpathsinglecontext(measurements[0],"string(dc:inversion_channel)")
    dc_cadmodel_channel_str=_xmldoc.xpathsinglecontext(measurements[0],"string(dc:cadmodel_channel)")
    tikparam_nominal=_xmldoc.xpathsinglecontext(measurements[0],"string(dc:greensinversion_dgsfile/@tikparam)")


    dc_prefix_str = "greensinversion_"

    giwfm_name=dc_prefix_str + dc_inversion_channel_str
    giwfm_weights_name=giwfm_name + "_weights"

    

    outwfmdict={}
    firstiter=True

    for measurement in measurements:
        gi_dgsfile_el=_xmldoc.xpathsinglecontext(measurement,"dc:greensinversion_dgsfile")
        greensinversion_dgs_href = dc_value.hrefvalue.fromxml(_xmldoc,gi_dgsfile_el)

        tikparam=_xmldoc.xpathsinglecontext(measurement,"string(dc:greensinversion_dgsfile/@tikparam)")

        assert(tikparam==tikparam_nominal)

        (junkmd,wfmdict)=dgf.loadsnapshot(greensinversion_dgs_href.getpath(),memmapok=True)
        if firstiter:
            outwfmdict[giwfm_name]=dg.wfminfo()
            outwfmdict[giwfm_name].Name=giwfm_name
            outwfmdict[giwfm_name].dimlen=wfmdict[giwfm_name].dimlen
            outwfmdict[giwfm_name].ndim=wfmdict[giwfm_name].dimlen.shape[0]
            outwfmdict[giwfm_name].data=wfmdict[giwfm_name].data
            outwfmdict[giwfm_name].MetaData=copy.deepcopy(wfmdict[giwfm_name].MetaData)
            
            outwfmdict[giwfm_weights_name]=dg.wfminfo()
            outwfmdict[giwfm_weights_name].Name=giwfm_weights_name
            outwfmdict[giwfm_weights_name].dimlen=wfmdict[giwfm_weights_name].dimlen
            outwfmdict[giwfm_weights_name].ndim=wfmdict[giwfm_weights_name].dimlen.shape[0]
            outwfmdict[giwfm_weights_name].data=wfmdict[giwfm_weights_name].data
            outwfmdict[giwfm_weights_name].MetaData=copy.deepcopy(wfmdict[giwfm_weights_name].MetaData)
            firstiter=False
            pass
        else:
            # Data is stored as a scaled sum of weights * data values (sum_i w_i * d_i)/(sum_i w_i) and a separate sum of weights (sum_i w_i)
            # For only a single contribution this maps to just the data value and weight
            
            # In order to accumulate the weighted average, we need to rescale the existing terms,  by multiplying by
            # the old_sum_of_weights divided by the new_sum_of_weights. Then we can add in the the new_data_value*new_weight/new_sum_of_weights
            
            # accumulate data, correcting scaling for previous weights
            # i.e. multiply by old weights, divide by new weights
            weightscaling = np.zeros(wfmdict[giwfm_weights_name].data.shape,dtype='f',order='F')
            
            weightscaling[wfmdict[giwfm_weights_name].data != 0] = outwfmdict[giwfm_weights_name].data[wfmdict[giwfm_weights_name].data != 0]/wfmdict[giwfm_weights_name].data[wfmdict[giwfm_weights_name].data != 0]
            # double check that the scaling factor is finite; Convert any infs or NaNs to 0
            weightscaling[~np.isfinite(weightscaling)]=0.0
            

            
            (nx,ny,nt)=outwfmdict[giwfm_name].data.shape


            accumulation_weights_rs=outwfmdict[giwfm_weights_name].data.reshape(nx*ny,order='F')

            got_accumulation = accumulation_weights_rs > 0
            
            # accumulate weights into new_sum_of_weights
            outwfmdict[giwfm_weights_name].data += wfmdict[giwfm_weights_name].data

            
            weightscaling_newdata = np.zeros((nx,ny),dtype='f',order='F')
            weightscaling_newdata[outwfmdict[giwfm_weights_name].data != 0] = wfmdict[giwfm_weights_name].data[outwfmdict[giwfm_weights_name].data != 0]/outwfmdict[giwfm_weights_name].data[outwfmdict[giwfm_weights_name].data != 0]
            # double check that the scaling factor is finite; Convert any infs or NaNs to 0
            weightscaling_newdata[~np.isfinite(weightscaling_newdata)]=0.0
            
            #paramwfm.data[:,:,:] = accumulation.data[:,:,:]*weightscaling[:,:,np.newaxis] + uv_proj_stack.data[:,:,:]*weightscaling_newdata[:,:,np.newaxis]



            uv_proj_stack_weights_rs=wfmdict[giwfm_weights_name].data.reshape(nx*ny,order='F')
            got_proj = uv_proj_stack_weights_rs > 0

            got_both = got_accumulation & got_proj

            got_only_proj = got_proj & (~got_accumulation)
            

            outwfmdict[giwfm_name].data.reshape(nx*ny,nt,order='F')[got_only_proj,:] = wfmdict[giwfm_name].data.reshape(nx*ny,nt,order='F')[got_only_proj,:]*weightscaling_newdata.reshape(nx*ny,1,order='F')[got_only_proj,:]
            

            outwfmdict[giwfm_name].data.reshape(nx*ny,nt,order='F')[got_both,:] = outwfmdict[giwfm_name].data.reshape(nx*ny,nt,order='F')[got_both,:]*weightscaling.reshape(nx*ny,1,order='F')[got_both,:] + wfmdict[giwfm_name].data.reshape(nx*ny,nt,order='F')[got_both,:]*weightscaling_newdata.reshape(nx*ny,1,order='F')[got_both,:]
            

            pass
        pass

    outwfmdict[dc_cadmodel_channel_str]=copy.deepcopy(wfmdict[dc_cadmodel_channel_str])


    channel3d = "Proj" + giwfm_name[:-4] # Proj + diffstack channel with _tex stripped
    objframe=coordframe()
    (obj, TexChanPrefix) = ndepart_from_dataguzzler_wfm(wfmdict[channel3d],wfmdict,objframe)

    SplitTextureChans=dgm.GetMetaDatumWIStr(wfmdict[dc_cadmodel_channel_str],"TextureChans","").split("|")
    PrefixedTextureChans="|".join([ TexChanPrefix + TexChan for TexChan in SplitTextureChans ])



    gi_3d=dg.wfminfo()
    #gi_3d.Name=dc_prefix_str+dc_cadmodel_channel_str
    gi_3d.Name="Proj"+TexChanPrefix+dc_cadmodel_channel_str
    gi_3d.dimlen=np.array((1,),dtype='i8')
    gi_3d.data=np.array((1,),dtype='f')
    dgm.AddMetaDatumWI(gi_3d,dgm.MetaDatum("VRML97GeomRef",dc_cadmodel_channel_str))
    dgm.AddMetaDatumWI(gi_3d,dgm.MetaDatum("X3DGeomRef",dc_cadmodel_channel_str))
    #texchanprefix=gi_3d.Name[:gi_3d.Name.find(dc_unprefixed_texname_str)]
    dgm.AddMetaDatumWI(gi_3d,dgm.MetaDatum("TexChanPrefix",TexChanPrefix))
    dgm.AddMetaDatumWI(gi_3d,dgm.MetaDatum("TextureChans",PrefixedTextureChans))
    outwfmdict[gi_3d.Name]=gi_3d

    outdgs_fname="%s_fused.dgs" % (posixpath.splitext(greensinversion_dgs_href.get_bare_unquoted_filename())[0])        
    outdgs_href=dc_value.hrefvalue(quote(outdgs_fname),_dest_href)
    dgf.savesnapshot(outdgs_href.getpath(),outwfmdict)

    return [ (("dc:fused_greensinversion_dgsfile",{"tikparam": str(tikparam)}),outdgs_href) ]
    

