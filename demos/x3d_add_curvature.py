
import sys
import os
import os.path

from spatialnde.ndeobj import ndepart
from spatialnde.coordframes import coordframe
from spatialnde.opencascade.loaders import load_byfilename
from spatialnde.opencascade.evalcurvature import add_curvatures_to_surface
from spatialnde.exporters.x3d import X3DSerialization

def usage():
    print("Usage: %s <X3DFile_With_UV> <CADfile.brep_step_or_iges> [ -o <Output_X3D_With_UV_and_curvature> ]")
    
    print("Output file name defaults to input file with _curvature appended to base name")
    print("If output file exists, will not overwrite")
    sys.exit(0)
    pass

positionals=[]
OutputX3D=None
CADFile=None
argcnt=1
while argcnt < len(sys.argv):
    if sys.argv[argcnt][0]=='-':
        if sys.argv[argcnt]=="-h" or sys.argv[argcnt]=="--help":
            usage()
            pass
        elif sys.argv[argcnt][1]=='o':
            OutputX3D=sys.argv[argcnt][2:]
            if len(OutputX3D) == 0:
                argcnt+=1
                OutputX3D = sys.argv[argcnt]
                pass
            pass
        else:
            raise ValueError("Unknown command line parameter %s" % (sys.argv[argcnt]))
        pass
    else:
        positionals.append(sys.argv[argcnt])
        pass
    argcnt+=1
    pass

if len(positionals) > 2:
    raise ValueError("Too many positional command line parameters (max 2)")

elif len(positionals) == 2:
    (InputX3D,CADFile) = positionals
    pass

elif len(positionals) == 1:
    InputX3D=positionals[0]

    pass
else:
    usage()
    pass


if OutputX3D is None:
    basename=os.path.splitext(InputX3D)[0]
    OutputX3D = basename + "_curvature.x3d"
    pass


if os.path.exists(OutputX3D):
    raise ValueError("File %s exists. Will not overwrite" % (OutputX3D))

# define object frame
objframe = coordframe()

obj=ndepart.fromx3d(objframe,None,InputX3D,tol=1e-6)

# Only one surface in object is supported for now
assert(len(obj.implpart.surfaces)==1)

if CADFile is not None:
    BRepShape=load_byfilename(CADFile)
    add_curvatures_to_surface(obj.implpart.surfaces[0],BRepShape,coarsetol=1e-3,finetol=1e-7)
else:
    # No CAD file given...
    # mesh-only approach
    obj.implpart.surfaces[0].add_curvatures_to_surface_from_mesh()
    pass


x3dwriter=X3DSerialization.tofileorbuffer(OutputX3D)
obj.X3DWrite(x3dwriter,objframe)
x3dwriter.finish()
