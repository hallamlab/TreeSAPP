from distutils.core import setup, Extension

setup(name='TreeSAPP_extend',
      version='0.1',
      description='Parsing utilities for TreeSAPP',
      author='Connor Morgan-Lang',
      author_email='c.morganlang@gmail.com',
      ext_modules=[Extension("_tree_parser", ["sub_binaries/TreeSAPP_extensions/tree_parsermodule.cpp"]),
                   Extension("_fasta_reader", ["sub_binaries/TreeSAPP_extensions/fasta_reader.cpp"])],
      requires=['pygtrie', 'ete3', 'numpy']
      )
