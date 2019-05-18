from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy


ext_modules = [
    Extension("code.nw", ["code/nw.pyx"],
              include_dirs=[".", numpy.get_include()]),
]

setup(name="cython_code",
      packages = find_packages(),
      ext_modules=cythonize(ext_modules),
      #script_args=["build_ext"],
      options={'build_ext': {'inplace': True}}
      )

