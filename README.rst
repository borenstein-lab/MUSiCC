
====================
MUSiCC Documentation
====================

MUSiCC is a toolkit for correcting biases in gene abundance measurements derived from shotgun metagenomic sequencing, 
and is available as an online tool and as a Python module. MUSiCC is developed by the Borenstein group at the University of Washington and is available online at: 
http://elbo.gs.washington.edu/software_musicc.html.

MUSiCC is a stand-alone Python module that implements the MUSiCC functionality. It is distributed under a BSD license and can be readily incorporated into custom analysis tools.

=========================
Installation Instructions
=========================
Prerequisites for installing:

In order for MUSiCC to run successfully, the following Python modules should be pre-installed on your system:
* Numpy >= 1.6.1 (http://www.numpy.org/)
* Scipy >= 0.9 (http://www.scipy.org/)
* Scikit-learn >= 0.15.2 (http://scikit-learn.org/stable/)
* Pandas >= 0.14 (http://pandas.pydata.org/)

To install MUSiCC, simply download the package from http://depts.washington.edu/elbogs/MUSiCC/MUSiCC_Python.zip. This is a zip archive containing the following files/directories:
* MUSiCC.py: The MUSiCC Python module
* data/: A directory containing several data files MUSiCC requires to run properly.
* examples/: A directory containing examples of input and output files.
* lic/COPYING.txt: A copy of the BSD License. This is required to be distributed with the MUSiCC package.

OR, install MUSiCC using the PyPI framework by running:
`pip install -U numpy scipy scikit-learn pandas` (for dependencies)
`pip install -U MUSiCC`

============================
Testing the software package
============================

After downloading and installing the software, we recommend testing it by running the following command
from the same directory where MUSiCC was installed:

`python test_MUSiCC.py`

This will invoke a series of tests. A correct output should end with:

Ran 3 tests in 12.972s
OK


==========================
Interface via command line
==========================
The MUSiCC module handles all calculations internally.
MUSiCC offers an interface to the MUSiCC functionality via the command line.
**usage:**

MUSiCC.py [-h] [-o OUTPUT_FILE] [-if {tab,csv,biom}]
                 [-of {tab,csv,biom}] [-n] [-c {use_generic,learn_model}]
                 [-perf] [-v]
                 input_file

positional arguments:
  input_file            Input abundance file to correct

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_FILE, --out OUTPUT_FILE
                        Output destination for corrected abundance (default:
                        MUSiCC.tab)
  -if {tab,csv,biom}, --input_format {tab,csv,biom}
                        Option indicating the format of the input file
                        (default: tab)
  -of {tab,csv,biom}, --output_format {tab,csv,biom}
                        Option indicating the format of the output file
                        (default: tab)
  -n, --normalize       Apply MUSiCC normalization (default: false)
  -c {use_generic, learn_model}, --correct {use_generic,learn_model}
                        Correct abundance per-sample using MUSiCC (default:
                        false)
  -perf, --performance  Calculate model performance on various gene sets (may
                        add to running time) (default: false)
  -v, --verbose         Increase verbosity of module (default: false)


===========================
Interface via python script
===========================
MUSiCC can also be used directly inside a python script. Passing variables and flags to the MUSiCC script is done by
creating a dictionary and passing it to the function *MUSiCC.correct*, as shown below.
**usage:**

>> from MUSiCC import MUSiCC
>> musicc_args = {'input_file': 'lib/python3.3/site-packages/MUSiCC/examples/simulated_ko_relative_abundance.tab',
                  'output_file': 'simulated_ko_MUSiCC_Normalized.tab', 'input_format': 'tab', 'output_format': 'tab', 'MUSiCC_inter': True,
                  'MUSiCC_intra': 'None', 'compute_scores': True, 'verbose': True}
>> MUSiCC.correct(musicc_args)


========
Examples
========
In the Examples directory, the file simulated_ko_relative_abundance.tab contains simulated KO abundance measurements of 20 samples described in the
MUSiCC manuscript. Using this file as input for MUSiCC results in the following files:
simulated_ko_MUSiCC_Normalized.tab (only normalization)
simulated_ko_MUSiCC_Normalized_Corrected_use_generic.tab (normalize and correct using the generic model learned from HMP)
simulated_ko_MUSiCC_Normalized_Corrected_learn_model.tab (normalize and correct learning a new model for each sample)

The commands used are the following (via command line):

`python MUSiCC.py examples/simulated_ko_relative_abundance.tab -n -perf -v -o examples/simulated_ko_MUSiCC_Normalized.tab`
`python MUSiCC.py examples/simulated_ko_relative_abundance.tab -n -c use_generic -perf -v -o examples/simulated_ko_MUSiCC_Normalized_Corrected_use_generic.tab`
`python MUSiCC.py examples/simulated_ko_relative_abundance.tab -n -c learn_model -perf -v -o examples/simulated_ko_MUSiCC_Normalized_Corrected_learn_model.tab`

