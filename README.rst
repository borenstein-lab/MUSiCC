
====================
MUSiCC Documentation
====================

MUSiCC is a marker genes based framework for metagenomic normalization and accurate profiling of gene abundances in the microbiome,
developed and maintained by the Borenstein group at the University of Washington.

============
Availability
============

MUSiCC is available through the following sources:

- As a Python module from GitHub or PyPI (see installation instructions below)
- As an online tool at: http://elbo.gs.washington.edu/software_musicc.html.

=======
License
=======

MUSiCC is distributed under a BSD license and can be readily incorporated into custom analysis tools.

=========================
Installation Instructions
=========================

Prerequisites for installing:

In order for MUSiCC to run successfully, the following Python modules should be pre-installed on your system:

- Numpy >= 1.17.0 (http://www.numpy.org/)
- Scipy >= 1.3.0 (http://www.scipy.org/)
- Scikit-learn >= 0.21.3 (http://scikit-learn.org/stable/)
- Pandas >= 0.25.0 (http://pandas.pydata.org/)

If you have *pip* installed, you can install these packages by running the following command:

``pip install -U numpy scipy scikit-learn pandas``

**Installing MUSiCC:**

To install MUSiCC, download the package from https://github.com/omanor/MUSiCC/archive/1.0.3.tar.gz

After downloading MUSiCC, you’ll need to unzip the file. If you’ve downloaded the release version, do this with the following command:

``tar -xzf MUSiCC-1.0.3.tar.gz``

You’ll then change into the new MUSiCC directory as follows:

``cd MUSiCC-1.0.3``

and install using the following command:

``python setup.py install``

ALTERNATIVELY, you can install MUSiCC directly from PyPI by running:

``pip install -U MUSiCC``

Note for windows users: Under some windows installations, Scipy may fail when importing the Stats module. Workarounds may be found online, such
as `here <https://code.google.com/p/pythonxy/issues/detail?id=745>`_.

============================
Testing the software package
============================

After downloading and installing the software, we recommend testing it by running the following command:

``test_musicc.py``

This will invoke a series of tests. A correct output should end with:

``Ran 3 tests in X.XXXXs``

``OK``

===============================
MUSiCC API via the command line
===============================
The MUSiCC module handles all calculations internally.
MUSiCC offers an interface to the MUSiCC functionality via the command line and the run_musicc script.

Usage:
------

``run_musicc.py input_file [options]``

Required arguments:
-------------------

**input_file**
    Input abundance file to correct

Optional arguments:
-------------------

**-h, --help**
    show help message and exit

**-o OUTPUT_FILE, --out OUTPUT_FILE**
    Output destination for corrected abundance (default: MUSiCC.tab)

**-if {tab,csv}, --input_format {tab,csv}**
    Option indicating the format of the input file (default: tab)

**-of {tab,csv}, --output_format {tab,csv}**
    Option indicating the format of the output file (default: tab)

**-n, --normalize**
    Apply MUSiCC normalization (default: false)

**-c {use_generic, learn_model}, --correct {use_generic,learn_model}**
    Correct abundance per-sample using MUSiCC (default: false)

**-perf, --performance**
    Calculate model performance on various gene sets (may add to running time) (default: false)

**-v, --verbose**
    Increase verbosity of module (default: false)


============================
MUSiCC API via python script
============================
MUSiCC can also be used directly inside a python script. Passing variables and flags to the MUSiCC script is done by
creating a dictionary and passing it to the function *correct_and_normalize*, as shown below.

Usage:
------

>>> from musicc.core import correct_and_normalize
>>> musicc_args = {'input_file': 'test_musicc/lib/python3.3/site-packages/musicc/examples/simulated_ko_relative_abundance.tab', 'output_file': 'MUSiCC.tab','input_format': 'tab', 'output_format': 'tab', 'musicc_inter': True, 'musicc_intra': 'learn_model','compute_scores': True, 'verbose': True}
>>> correct_and_normalize(musicc_args)

Required arguments:
-------------------

**input_file**
    Input abundance file to correct

Optional arguments:
-------------------

**output_file**
    Output destination for corrected abundance (default: MUSiCC.tab)

**input_format {'tab','csv'}**
    Option indicating the format of the input file (default: 'tab')

**output_format {'tab','csv'}**
    Option indicating the format of the output file (default: 'tab')

**musicc_inter {True, False}**
    Apply MUSiCC normalization (default: False)

**musicc_intra {'use_generic', 'learn_model', 'None'}**
    Correct abundance per-sample using MUSiCC (default: 'None')

**compute_scores {True, False}**
    Calculate model performance on various gene sets (may add to running time) (default: False)

**verbose {True, False}**
    Increase verbosity of module (default: False)

========
Examples
========
In the *musicc/examples* directory, the file *simulated_ko_relative_abundance.tab* contains simulated KO abundance measurements of 20 samples described in the
MUSiCC manuscript. Using this file as input for MUSiCC results in the following files:

- simulated_ko_MUSiCC_Normalized.tab (only normalization)
- simulated_ko_MUSiCC_Normalized_Corrected_use_generic.tab (normalize and correct using the generic model learned from HMP)
- simulated_ko_MUSiCC_Normalized_Corrected_learn_model.tab (normalize and correct learning a new model for each sample)

The commands used were the following (via command line):

``run_musicc.py musicc/examples/simulated_ko_relative_abundance.tab -n -perf -v -o musicc/examples/simulated_ko_MUSiCC_Normalized.tab``

``run_musicc.py musicc/examples/simulated_ko_relative_abundance.tab -n -c use_generic -perf -v -o musicc/examples/simulated_ko_MUSiCC_Normalized_Corrected_use_generic.tab``

``run_musicc.py musicc/examples/simulated_ko_relative_abundance.tab -n -c learn_model -perf -v -o musicc/examples/simulated_ko_MUSiCC_Normalized_Corrected_learn_model.tab``

==================
Citing Information
==================

If you use the MUSiCC software, please cite the following paper:

MUSiCC: A marker genes based framework for metagenomic normalization and accurate profiling of gene abundances in the microbiome.
**Ohad Manor and Elhanan Borenstein.** *Genome Biology*

==============
Question forum
==============
For MUSiCC announcements and questions, including notification of new releases, you can visit the `MUSiCC users forum <https://groups.google.com/forum/#!forum/musicc-users>`_.
