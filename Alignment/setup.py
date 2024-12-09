# setup.py
from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("filter_orthologs.pyx", compiler_directives={'language_level': "3"})
)

#python setup.py build_ext --inplace
