from distutils.core import setup
import numpy as np
import obspy as obs
from Cython.Build import cythonize

setup(ext_modules=cythonize('auto_fchp.pyx', compiler_directives={'language_level' : "3"}),include_dirs=[np.get_include()])