from distutils.core import setup
from distutils.extension import Extension
import os
import sys
import platform

openmm_dir = '@OPENMM_DIR@'
buckinghamplugin_header_dir = '@BUCKINGHAMPLUGIN_HEADER_DIR@'
buckinghamplugin_library_dir = '@BUCKINGHAMPLUGIN_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = []
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

extension = Extension(name='_buckingham',
                      sources=['BuckinghamPluginWrapper.cpp'],
                      libraries=['OpenMM', 'BuckinghamPlugin'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), buckinghamplugin_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), buckinghamplugin_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='buckinghamplugin',
      version='1.0',
      py_modules=['buckinghamplugin'],
      ext_modules=[extension],
     )
