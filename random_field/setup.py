from setuptools import setup

classifiers = [
               "Development Status :: 5 - Production/Stable",
               "Intended Audience :: Education",
               "Operating System :: OS Independent",
               "License :: OSI Approved :: MIT License",
               "Programming Language :: Python :: 3"
               ]

setup(
    name='adding',
    version='0.0.3',
    description='this is my first test',
    long_description=open('README.txt').read() + '\n\n' + open('CHANGELOG.txt').read(),
    url='',
    author='Yang Yang',
    author_email='dienyoung@outlook.com',
    license='MIT',
    classifiers=classifiers,
    keywords='',
    py_modules=["adding"],
    package_dir={'': 'src'},
    install_requires=['']
)
