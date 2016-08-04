from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("transform_calculateCCs_model.pyx")
)


### PRODUCE .so BY COMMAND LINE COMMAND: python setup.py build_ext --inplace ###
