import numpy as np
import os
import shutil
from distutils.core import Distribution, Extension
from Cython.Build import build_ext, cythonize

cython_dir = os.path.join("src","ucla_geotech_tools")
extension1 = Extension(
    "ucla_geotech_tools.ipyconsol",
    [
        os.path.join(cython_dir,"ipyconsol.pyx"),
    ],
    include_dirs=[np.get_include()]
)
extension2 = Extension(
    "ucla_geotech_tools.random_field",
    [
        os.path.join(cython_dir,"random_field.pyx")
    ],
    include_dirs=[np.get_include()]
    #extra_compile_args=["-O3", "-std=c++17"],
)

ext_modules = cythonize([extension1, extension2], include_path=[cython_dir])
dist = Distribution({"ext_modules": ext_modules})
cmd = build_ext(dist)
cmd.ensure_finalized()
cmd.run()

for output in cmd.get_outputs():
    relative_extension = os.path.relpath(output, cmd.build_lib)
    shutil.copyfile(output, 'src/' + relative_extension)