import sys
import os
import os.path
import subprocess
import re
from setuptools import setup
from setuptools.command.install_lib import install_lib
from setuptools.command.install import install
from setuptools.command.build_ext import build_ext
from setuptools.extension import Extension
import setuptools.command.bdist_egg
import sys
import distutils.spawn
import numpy as np

from Cython.Build import cythonize
#from numpy.distutils.core import setup as numpy_setup, Extension as numpy_Extension



extra_compile_args = {
    "msvc": ["/openmp"],
    #"unix": ["-O0", "-g", "-Wno-uninitialized","-std=c11","-fstack-protector-all"),    # Replace the line below with this line to enable debugging of the compiled extension
    "unix": ["-fopenmp","-O5","-Wno-uninitialized","-std=c11","-fstack-protector-all"],
    "mingw32": [ "-fopenmp","-O5","-std=c11","-fstack-protector-all" ],
    "clang": ["-fopenmp","-O5","-Wno-uninitialized"],
}

extra_include_dirs = {
    "msvc": [".", np.get_include() ],
    "unix": [".", np.get_include() ],
    "clang": [".", np.get_include() ],
}

extra_libraries = {
    "msvc": [],
    "unix": ["gomp",],
    "clang": [],
}

extra_link_args = {
    "msvc": [],
    "clang": ["-fopenmp=libomp"],
    "unix": [ "-fopenmp","-std=c11" ],
    "mingw32": [ "-fopenmp","-std=c11" ],
    "clang": [".", np.get_include(),"-fopenmp=libomp"],
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

        

class build_ext_compile_args(build_ext):
    def build_extensions(self):
        compiler=self.compiler.compiler_type
        for ext in self.extensions:
            if compiler in extra_compile_args:
                ext.extra_compile_args=extra_compile_args[compiler]
                ext.extra_link_args=extra_link_args[compiler]
                ext.include_dirs.extend(list(extra_include_dirs[compiler]))
                ext.libraries.extend(list(extra_libraries[compiler]))
                pass
            else:
                # use unix parameters as default
                ext.extra_compile_args=extra_compile_args["unix"]
                ext.extra_link_args=extra_link_args["unix"]
                ext.include_dirs.extend(list(extra_include_dirs["unix"]))
                ext.libraries.extend(extra_libraries["unix"])
                pass
                
            pass
            
        
        build_ext.build_extensions(self)
        pass
    pass

spatialnde_package_files = [ "pt_steps/*" ]

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
      zip_safe=False,
      packages=["spatialnde",
                "spatialnde.cadpart",
                "spatialnde.cadpart.loaders",
                "spatialnde.exporters",
                "spatialnde.dataguzzler",
                "spatialnde.opencascade"],
      scripts=["scripts/fiducialname","scripts/spatialnde_calib_image"],
      cmdclass = {"build_ext": build_ext_compile_args},
      package_data={"spatialnde": spatialnde_package_files},
      entry_points={"limatix.processtrak.step_url_search_path": [ "limatix.share.pt_steps = spatialnde:getstepurlpath" ]}
)
