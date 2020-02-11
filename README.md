[![Build Status](https://travis-ci.org/genomematt/pylazybam.svg?branch=master)](https://travis-ci.org/genomematt/pylazybam)
[![Coverage Status](https://coveralls.io/repos/genomematt/pylazybam/badge.svg)](https://coveralls.io/r/genomematt/pylazybam)
[![DOI](https://zenodo.org/badge/189527217.svg)](https://zenodo.org/badge/latestdoi/189527217)

Pylazybam
=========

Pylazybam is a pure python library for reading minimal amounts of information from a BAM format mapped sequence
alignment file. It is intended for uses such as filtering reads where the information within the single alignment entry 
is not sufficient to make the filtering decision.

For simple filtering you should consider other approaches such as samtools or sambamba. If editing of the data is
required htslib based solutions, such as pysam should also be considered.

Pylazybam is a minimalist architecture consisting of classes for reading and writing BAM files, and a set of functions 
that can be used on alignments in binary bytestring format to extract information. In filtering applications decisions 
on read alignment output are made on the processed data, and the raw unmodified BAM alignment is written to the output 
BAM file. This minimizes the decoding and encoding work done by the code, thus the pylazybam name. 

Installation
============
Pylazybam requires Python 3.6 or higher and is tested on Linux and MacOS with CPython and PyPy3.

<!-- Installing from the Python Package Index with pip is the easiest option:

    pip3 install pylazybam
-->    
To install from the github repository

    pip install git+git://github.com/genomematt/pylazybam.git

or alternatively by cloning the github repository

    git clone https://github.com/genomematt/pylazybam
    pip install pylazybam
	
Although the repository tests by continuous integration with TravisCI its good practice to run the tests locally and 
check your install works correctly.

The tests are run with the following command:

    python3 -m pylazybam.tests.test_all

Using pylazybam
===============

Pylazybam is a library, and each use case will require the user to construct a bespoke script for their application.

In most applications this will involve opening a compressed BAM file with `gzip`, parsing the header with 
`bam.FileReader` and then extracting information such as the alignment score tag with `bam.get_AS`

For example, a simple script to count the number of primary mappings per reference:

    import gzip
    from collections import Counter
    from pylazybam import bam
    
    counts = Counter()
    
    with bam.FileReader(gzip.open('path/to/bam.bam')) as mybam:    
        for align in mybam:
            if bam.is_flag(align, bam.FLAG['primary']):
                ref_index = bam.get_ref_index[align]
                refname = mybam.index_to_ref[ref_index]
                counts.update([refname,])
    
    print(counts)
    
For more information on available functions and documentation

    from pylazybam import bam
    help(bam)
 
Examples of how to use pylazybam can be found in [example_usage.ipynb](example_usage.ipynb) and a brief example of using
from within R with [reticulate](https://rstudio.github.io/reticulate/) in [reticulate_example.ipynb](reticulate_example.ipynb)

Contributing to pylazybam
=========================
Pylazybam is licensed under the BSD three clause license.  You are free to fork this repository under the terms of that
 license.  If you have suggested changes please start by raising an issue in the issue tracker.  Pull requests are 
welcome and will be included at the discretion of the author, but must have 100% test coverage.

Bug reports should be made to the issue tracker.  Difficulty in understanding how to use the software is a documentation
 bug, and should also be raised on the issue tracker and will be tagged `question` so your question and my response are 
easily found by others.

Pylazybam uses numpy style docstrings, python type annotations, Travis CI, coverage and coveralls. All code should be
compatible with python versions >= 3.6 and contain only pure python code.


Citing pylazybam
================

Pylazybam is in early development and does not yet have a publication. Please cite the github repository.
Each release will have a Zenodo DOI identifier that can be cited. The current DOI for v0.1.0 is [![DOI](https://zenodo.org/badge/189527217.svg)](https://zenodo.org/badge/latestdoi/189527217)

Acknowledgements
----------------
Pylazybam utilizes the excellent bgzf implementation from BioPython written by Peter Cock [@peterjc](https://github.com/peterjc). The slightly modified version is included in this package under the BSD variant of the bgzf codes licensing (this is the same license as pylazybam). The original version of the bgzf code can be found [here](https://github.com/biopython/biopython/blob/master/Bio/bgzf.py)

Thanks to Alan Rubin [@afrubin](https://github.com/afrubin) and Tony Papenfuss [@papenfuss](https://github.com/papenfuss) for helpful discussions and code review