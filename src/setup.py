from distutils.core import setup
from Cython.Build import cythonize

setup(name='Approx_matching', ext_modules=cythonize('/home/runner/work/project-5-python-asg/project-5-python-asg/src/Approx_Positions_Cython.pyx'))
