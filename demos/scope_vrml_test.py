import numpy as np

import dataguzzler as dg
import dg_metadata as dgm
import dg_file as dgf

fin = "/home/linuxadm/usr_local/src/freecad-git032314/build/data/Mod/Robot/Lib/Kuka/kr125_3.wrl"

fout = "/tmp/robot.dgs"

wfmdict={}

fh=open(fin,"r");

wfmdict["robot"]=dg.wfminfo()
wfmdict["robot"].Name="robot"
wfmdict["robot"].dimlen=np.array((),dtype='i8')
dgm.AddMetaDatumWI(wfmdict["robot"],dgm.MetaDatum("VRML97Geom",fh.read()))

dgf.savesnapshot(fout,wfmdict)

