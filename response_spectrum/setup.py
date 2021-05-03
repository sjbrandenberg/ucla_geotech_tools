from setuptools import Extension, setup, find_packages
import numpy as np

with open("README.md","r") as fh:
    long_description = fh.read()

sourcefiles = ['response_spectrum.c']
extensions = [Extension("ucla_geotech_tools.response_spectrum",sourcefiles)]

setup(
    name="ucla_geotech_tools-response_spectrum",
    version="1.0.0",
    author="Scott J. Brandenberg",
    author_email="sjbrandenberg@ucla.edu",
    description="Python package for computing acceleration response spectra.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    ext_modules=[Extension('ucla_geotech_tools.response_spectrum',['response_spectrum.c'])],
    namespace_packages=['ucla_geotech_tools'],	
    include_dirs=[np.get_include()]
)
