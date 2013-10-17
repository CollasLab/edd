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
except ImportError:
    raise Exception('please install cython first, e.g.: pip install --upgrade cython')

try:
    csamtools_path = next(p for p in (os.path.join(base_dir, 'csamtools.pxd') 
        for base_dir in pysam.get_include()) if os.path.isfile(p))
except StopIteration:
    print 'unable to locate file csamtools.pxd, which is required for compilation'
    import sys
    sys.exit()

setup(name='edd',
      version='0.1',
      description='ChIP-seq broad domain peak caller',
      url='http://github.com/eivindgl/edd',
      author='Eivind G. Lund',
      author_email='e.g.lund@medisin.uio.no',
      packages=['edd'],
      scripts=[
          'bin/edd',
          'bin/edd-debug',
          'bin/edd-score-cutoff'],
      install_requires=[
          'numpy',
          'pandas',
          'matplotlib',
          ],
    ext_modules=[
            Extension('max_segments', ['edd/max_segments.pyx']),
            Extension('read_bam',
                sources=['edd/read_bam_into_bins.pyx', csamtools_path],
                 include_dirs=pysam.get_include() + [np.get_include()],
                 define_macros=pysam.get_defines()),
                 ],
    cmdclass = {'build_ext': build_ext},
    )
