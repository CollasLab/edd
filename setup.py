#from distutils.core import setup
from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


setup(name='edd',
      version='0.1',
      description='ChIP-seq broad domain peak caller',
      url='http://github.com/eivindgl/edd',
      author='Eivind G. Lund',
      author_email='e.g.lund@medisin.uio.no',
      packages=['edd'],
      scripts=['bin/edd',
               'bin/edd-debug',
               'bin/edd-score-cutoff'],
            install_requires=[
          'numpy',
          'pandas'
      ],

    ext_modules=[
            Extension('max_segments', ['edd/max_segments.pyx']),
            ],
    cmdclass = {'build_ext': build_ext},
    )
