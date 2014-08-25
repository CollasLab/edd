import os, sys
from setuptools import setup
from distutils.extension import Extension

if not sys.version_info[:2] == (2, 7):
    raise Exception('''\
    Edd requires python version 2.7.x, but you are using %d.%d.%d''' %
    sys.version_info[:3])

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
      version='1.1.8',
      description='Enriched domain detector for ChIP-seq data',
      url='http://github.com/CollasLab/edd',
      author='Eivind G. Lund',
      author_email='e.g.lund@medisin.uio.no',
      packages=['eddlib',
                'eddlib.algorithm'],
      # installs into root dir (not what i want)
      #data_files=[('eddlib', ['eddlib/default_parameters.conf'])],
      package_data={'': ['*.conf']},
      scripts=[
          'bin/edd',
        #'bin/edd-tools',
          ],
      install_requires=[
          'Logbook',
          'pybedtools',
          'statsmodels',
          'patsy', # statsmodels dependency
          'pandas',
          'python-dateutil', # pandas dependency
          'scipy',
          'numpy',
          ],
    ext_modules=[
        Extension('eddlib.read_bam',
                 sources=['eddlib/read_bam.pyx'],
                 include_dirs=pysam.get_include() + [np.get_include()],
                 define_macros=pysam.get_defines(),
                 ),
        Extension('eddlib.algorithm.chrom_max_segments',
                  sources=['eddlib/algorithm/chrom_max_segments.pyx'],
                  include_dirs=[np.get_include()]),
                 ],
    cmdclass = {'build_ext': build_ext},
    )
