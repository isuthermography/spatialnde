from distutils.core import Extension
from distutils.command.build_ext import build_ext
import sys
import subprocess
import numpy as np
import os

from Cython.Build import cythonize
#from numpy.distutils.core import setup as numpy_setup, Extension as numpy_Extension
from distutils.core import setup


openmp_compiler_opts = {
    "unix": [ "-fopenmp", "-std=c11","-fstack-protector-all"],
    "mingw32": [ "-fopenmp","-std=c11","-fstack-protector-all" ],
    "msvc": [ "/openmp" ]
}

openmp_linker_opts = {
    "unix": [ "-fopenmp","-std=c11" ],
    "mingw32": [ "-fopenmp","-std=c11" ],
}

if os.name=="posix":
    # Check for obsolete gcc, try to use Red Hat Sofware Collections
    # if necessary
    gccversioninfo = subprocess.check_output("gcc --version",shell="True")
    gccver = gccversioninfo.split()[2]
    gccsubvers=gccver.decode('utf-8').split(".")
    if int(gccsubvers[0]) <= 4 and int(gccsubvers[1]) < 9:
        # old version with C11 atomics
        if os.path.exists("/opt/rh/devtoolset-7/root/bin/gcc"):
            os.environ["CC"]="/opt/rh/devtoolset-7/root/bin/gcc"
            pass
        pass
    pass

        

class build_ext_openmp(build_ext):
    # Custom build_ext subclass that adds on extra
    # compilation/link parameters
    # on a compiler-specific basis

    def build_extensions(self):
        # sys.stderr.write("compiler_type=%s\n" % (self.compiler.compiler_type))
        for ext in self.extensions:
            if ext.extra_compile_args is None:
                ext.extra_compile_args=[]
                pass
            
            #ext.extra_compile_args.extend(["-O0", "-g"])
            if self.compiler.compiler_type in openmp_compiler_opts:
                ext.extra_compile_args.extend(openmp_compiler_opts[self.compiler.compiler_type])
                pass

            if ext.extra_link_args is None:
                ext.extra_link_args=[]
                pass
            if self.compiler.compiler_type in openmp_linker_opts:
                ext.extra_link_args.extend(openmp_linker_opts[self.compiler.compiler_type])
                pass
            pass
        return build_ext.build_extensions(self)
    
    
    
    pass

# NOTE: MSVC 2008 must remove "m" from libraries
ext_modules=cythonize([ Extension("spatialnde.imageprojection",["spatialnde/imageprojection.pyx"],
                                  libraries=[ ],
                                  include_dirs=[np.get_include()]),
                        
                        Extension("spatialnde.imageorthographic",["spatialnde/imageorthographic.pyx"],
                                  libraries=[ ],
                                  include_dirs=[np.get_include()]),

                        Extension("spatialnde.geometry_accel",["spatialnde/geometry_accel.pyx"],
                                  libraries=[ ],
                                  include_dirs=[np.get_include()]),
                        Extension("spatialnde.cadpart.polygonalsurface_texcoordparameterization_accel",["spatialnde/cadpart/polygonalsurface_texcoordparameterization_accel.pyx"],
                                  libraries=[ ],
                                  include_dirs=[np.get_include()]),

                        Extension("spatialnde.cadpart.polygonalsurface_accel",["spatialnde/cadpart/polygonalsurface_accel.pyx"],
                                  libraries=[ ],
                                  include_dirs=[np.get_include()])])


setup(name="spatialnde",
      description="spatialnde",
      author="Stephen D. Holland",
      url="http://thermal.cnde.iastate.edu/spatialnde",
      ext_modules=ext_modules,
      packages=["spatialnde",
                "spatialnde.cadpart",
                "spatialnde.cadpart.loaders",
                "spatialnde.exporters",
                "spatialnde.dataguzzler",
                "spatialnde.opencascade"],
      scripts=["scripts/fiducialname","scripts/spatialnde_calib_image"],
      cmdclass = {"build_ext": build_ext_openmp},
)
