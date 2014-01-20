import os
#from distutils.core import setup
from setuptools import setup
from distutils.extension import Extension

try:
    import pysam
    import numpy as np
except ImportError:
    raise Exception('please install pysam first, e.g.: pip install --upgrade pysam')
try:
    from Cython.Distutils import build_ext # Cython should be installed via pysam
    #from Cython.Distutils.extension import Extension
except ImportError:
    raise Exception('please install cython first, e.g.: pip install --upgrade cython')

setup(name='edd',
      version='0.1',
      description='ChIP-seq broad domain peak caller',
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
          ],
    ext_modules=[
        Extension('edd.read_bam',
                 sources=['edd/read_bam.pyx'],
                 include_dirs=pysam.get_include() + [np.get_include()],
                 define_macros=pysam.get_defines(),
                 ),
                 ],
    cmdclass = {'build_ext': build_ext},
    )
