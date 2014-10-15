
import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), 'r') as f:
        return f.read()

setup(name='MUSiCC',
      version='1.0',
      classifiers=['License :: OSI Approved :: BSD License'],
      license=['BSD'],
      description='A toolkit for normalization and correction of gene abundance measurements derived from shotgun metagenomic sequencing',
      long_description=(read('README.rst') + '\n\n' + read('HISTORY.rst') + '\n\n' + read('AUTHORS.rst') + '\n\n' + read('CITE.rst')),
      author='Ohad Manor',
      author_email='omanor@gmail.com',
      url='http://elbo.gs.washington.edu/software_musicc.html',
      packages=['MUSiCC'],
      package_dir={'MUSiCC': '.'},
      package_data={'MUSiCC': ['data/*.tab', 'data/*.lst', 'examples/*.tab']},
      install_requires=['NumPy >= 1.6.1', 'SciPy >= 0.9', 'scikit-learn >= 0.15.2', 'pandas >= 0.14'],
      provides=['MUSiCC'],
      )




