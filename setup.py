from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


setup(name='enriched_domain_caller',
      version='0.1',
      description='ChIP-seq broad domain peak caller',
      url='http://github.com/eivindgl/edc',
      author='Eivind G. Lund',
      author_email='e.g.lund@medisin.uio.no',
      package_dir={'enriched_domain_caller' : 'enriched_domain_caller'},
      packages=['enriched_domain_caller'],
      scripts=['bin/edc'],
      # install_requires=[
      #     'Cython'
      # ],
    ext_modules=[
            Extension('enriched_domain_caller/max_segments', ['enriched_domain_caller/max_segments.pyx']),
            ],
    cmdclass = {'build_ext': build_ext},
    )
