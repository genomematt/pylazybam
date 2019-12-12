#!/usr/bin/env python
from setuptools import setup

setup(
    name='pylazybam',
    version='0.1.0',
    author='Matthew Wakefield',
    author_email='matthew.wakefield@unimelb.edu.au',
    install_requires = [
      'setuptools >= 1.1.6',
    ],
    packages=['pylazybam',
              'pylazybam.tests',
              ],
    zip_safe = True,
    include_package_data = True,
    url='https://git@bitbucket.org/genomematt/ambivert.git',
    license='BSD',
    entry_points={
        'console_scripts': ['pylazybam_unique = pylazybam.unique:main',
                           ]
    },
    test_suite = "pylazybam.tests.test_pylazybam",
    description='pylazybam - a pure python bam parser for rapid content based read sorting',
    long_description=open('README.txt').read(),
    classifiers=[
          'Development Status :: 5 - Production/Stable',
          'License :: OSI Approved :: BSD License',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

)
