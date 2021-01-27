from setuptools import Extension, setup, find_packages

with open("README.md","r") as fh:
    long_description = fh.read()

sourcefiles = ['pyConsol.c']
extensions = [Extension("pyConsol",sourcefiles)]

setup(
    name="pyConsol",
    version="1.0.6",
    author="Scott J. Brandenberg",
    author_email="sjbrandenberg@ucla.edu",
    description="Python implementation of iConsol program for computing nonlinear consolidation of soil.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    ext_modules=[Extension('pyConsol',['pyConsol.c'])]
)