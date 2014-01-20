import os
from setuptools import setup
from distutils.extension import Extension

try:
    import numpy as np
except ImportError:
    raise Exception('''\
EDD has compile time dependencies on numpy. So please install numpy first.
e.g.: pip install --upgrade numpy''')

try:
    import pysam
except ImportError:
    raise Exception('''\
EDD has compile time dependencies on pysam. So please install pysam first.
e.g.: pip install --upgrade pysam''')
try:
    from Cython.Distutils import build_ext # Cython should be installed via pysam
    #from Cython.Distutils.extension import Extension
except ImportError:
    raise Exception('please install cython first, e.g.: pip install --upgrade cython')

setup(name='edd',
      version='0.9rc1',
      description='Enriched domain detector for ChIP-seq data',
      url='http://github.com/eivindgl/edd',
      author='Eivind G. Lund',
      author_email='e.g.lund@medisin.uio.no',
      packages=['edd'],
      scripts=[
          'bin/edd',
          ],
      install_requires=[
          'numpy',
          'pandas',
          'matplotlib',
          'Logbook',
          'pybedtools',
          'rpy2'
          ],
    ext_modules=[
        Extension('edd.read_bam',
                 sources=['edd/read_bam.pyx'],
                 include_dirs=pysam.get_include() + [np.get_include()],
                 define_macros=pysam.get_defines(),
                 ),
        Extension('edd.algorithm.chrom_max_segments',
                  ['edd/algorithm/chrom_max_segments.pyx']),
                 ],
    cmdclass = {'build_ext': build_ext},
    )
