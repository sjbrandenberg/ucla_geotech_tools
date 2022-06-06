from setuptools import Extension, setup, find_packages
import numpy as np

with open("README.md","r") as fh:
    long_description = fh.read()

sourcefiles = ['auto_fchp.c']
extensions = [Extension("ucla_geotech_tools.auto_fchp",sourcefiles)]

setup(
    name="ucla_geotech_tools-auto_fchp",
    version="0.0.6",
    author="Scott J. Brandenberg",
    author_email="sjbrandenberg@ucla.edu",
    description="Python package for automatically selecting high-pass filter corner frequency.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    ext_modules=[Extension('ucla_geotech_tools.auto_fchp',['auto_fchp.c'])],
    namespace_packages=['ucla_geotech_tools'],	
    include_dirs=[np.get_include()]
)