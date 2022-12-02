from distutils.core import setup
from Cython.Build import cythonize

setup(name='Approx_matching', ext_modules=cythonize('Approx_Positions_Cython.pyx'))
