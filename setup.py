import os
from glob import glob
from setuptools import setup, find_packages
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext

this_dir = os.path.dirname(os.path.abspath(__file__))

with open("README.md", "r") as fh:
    long_description = fh.read()


source_files = glob('src/*.cpp')
source_files.extend(glob('src/**/*.cpp'))
source_files.remove('src/main.cpp')
source_files.remove('src/main_parser.cpp')
#source_files.remove('src/parser.cpp')
source_files.append('python_wrapper/wrap.cpp')
print(source_files)


ext_modules = [
    Extension(
    'avbmc',
        source_files,
        include_dirs=[os.path.join(this_dir, 'include'),
                      os.path.join(this_dir, 'src'),
                ], #'src/', 'src/boundary', 'src/constraint', 'src/integrator', 'src/moves', 'src/sampler', 'src/rng', 'src/forcefield'],
    language='c++',
    extra_compile_args = ['-std=c++14'],
    ),
]


setup(name='avbmc',
      version="0.0.2",
      description='Python library for atomistic Monte Carlo simulations',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='http://github.com/evenmn/avbmc',
      author='Even Marius Nordhagen',
      author_email='evenmn@mn.uio.no',
      license='MIT',
      include_package_data=True,
      zip_safe=False,
      packages=find_packages(),
      cmdclass={'build_ext': build_ext},
      classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: Apache Software License",
      ],
      python_requires='>=3.5',
      #install_requires=["pybind11"],
      ext_modules=ext_modules,
)
