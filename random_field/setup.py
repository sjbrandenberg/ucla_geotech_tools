from setuptools import Extension, setup, find_packages
import numpy as np

with open("README.md","r") as fh:
    long_description = fh.read()

sourcefiles = ['random_field.c']
extensions = [Extension("ucla_geotech_tools.random_field",sourcefiles)]

setup(
    name='ucla_geotech_tools-random_field',
    version='1.0.4',
    description='Randosm field generation',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Yang Yang',
    author_email='dienyoung@outlook.com',
    license='MIT',
    packages = find_packages(),
    ext_modules=[Extension('ucla_geotech_tools.random_field',['random_field.c'])],
    namespace_packages=['ucla_geotech_tools'],
    include_dirs=[np.get_include()]
)
