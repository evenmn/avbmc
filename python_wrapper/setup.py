from glob import glob
from distutils.core import setup, Extension

source_files = glob('../src/*.cpp')
source_files.extend(glob('../src/**/*.cpp'))
source_files.remove('../src/main.cpp')
source_files.remove('../src/main_parser.cpp')
source_files.remove('../src/parser.cpp')
source_files.append('wrap.cpp')
print(source_files)

ext_modules = [
    Extension(
    'avbmc_core',
        source_files,
        include_dirs=['/home/evenmn/pybind11/include', '../src/', '../src/boundary', '../src/constraint', '../src/integrator', '../src/moves', '../src/sampler', '../src/rng', '../src/forcefield'],
    language='c++',
    extra_compile_args = ['-std=c++14'],
    ),
]

setup(
    name='avbmc',
    version='0.0.1',
    author='Even Marius Nordhagen',
    author_email='evenmn@mn.uio.no',
    description='Example',
    ext_modules=ext_modules,
)
