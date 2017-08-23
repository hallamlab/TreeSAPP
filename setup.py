from distutils.core import setup, Extension

setup(name='TreeSAPP_extend',
      version='0.1',
      description='Newick tree parsing utilities for TreeSAPP',
      author='Connor Morgan-Lang',
      author_email='c.morganlang@gmail.com',
      ext_modules=[Extension("_tree_parser", ["tree_parsermodule.cpp"]),
                   Extension("_fasta_reader", ["fasta_reader.cpp"])],
      )
