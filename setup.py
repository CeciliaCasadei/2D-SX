from distutils.core import setup
from Cython.Build import cythonize
import numpy

print numpy.get_include()

setup(
    ext_modules = cythonize("transform_calculateCCs.pyx", ),
    include_dirs=[numpy.get_include(), ]
)


### PRODUCE .so BY COMMAND LINE COMMAND: python setup.py build_ext --inplace ###
