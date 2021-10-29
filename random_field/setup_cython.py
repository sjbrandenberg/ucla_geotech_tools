from setuptools import setup
import numpy as np
from Cython.Build import cythonize
setup(ext_modules=cythonize('random_field.pyx'),include_dirs=[np.get_include()])
