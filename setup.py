
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
      version='1.0.4',
      classifiers=['License :: OSI Approved :: BSD License'],
      license=['BSD'],
      description='MUSICC: A marker genes based framework for metagenomic normalization and accurate profiling of gene abundances in the microbiome.',
      long_description=(read('README.rst') + '\n\n' + read('HISTORY.rst') + '\n\n' + read('AUTHORS.rst') + '\n\n'),
      author='Ohad Manor',
      author_email='omanor@gmail.com',
      url='http://elbo.gs.washington.edu/software_musicc.html',
      packages=['musicc'],
      package_data={'musicc': ['data/*.tab', 'data/*.lst', 'examples/*.tab']},
      install_requires=['NumPy >= 1.17.0', 'SciPy >= 1.3.0', 'scikit-learn >= 0.21.3', 'pandas >= 0.25.0'],
      scripts=['scripts/run_musicc.py', 'tests/test_musicc.py'],
      )




