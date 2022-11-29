from distutils.core import setup
from Cython.Build import cythonize

setup(name='Readmapper_module', ext_modules=cythonize('Readmapper_Cython_Convert.pyx'))
