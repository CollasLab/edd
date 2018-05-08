import os, sys
from setuptools import setup
from distutils.extension import Extension

def get_version(m):
    xs = m.__version__.split('.')
    xs += ['0', '0']
    version, major, minor = xs[:3]
    return version, major, minor

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
    version, major, minor = get_version(pysam)
    too_old_pysam = int(version) == 0 and int(major) < 10
    too_recent_pysam = int(version) != 0 or int(major) >= 12
    if too_old_pysam or too_recent_pysam:
        sys.stderr.write('''\

        ###########
        #  ERROR  #
        ###########
EDD is only compatible with pysam versions from 0.10.0 up to and including 0.11.2.2.
The detected version was %s.
Aborting ...
''' % pysam.__version__)
        sys.exit(1)
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
      version='1.1.19',
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
          'pysam>=0.10,<0.12',
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
