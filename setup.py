import os
from distutils.core import setup, Extension

setup(name='MLTreeParse',
      version='0.1',
      description='Newick tree parsing utilities for MLTreeMap',
      author='Connor Morgan-Lang',
      author_email='c.morganlang@gmail.com',
      # ext_package='_tree_parser',
      ext_modules=[Extension("_tree_parser", ["tree_parsermodule.cpp"])],
      )
