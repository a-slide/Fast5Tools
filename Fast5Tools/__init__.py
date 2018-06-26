# -*- coding: utf-8 -*-

# Define self package variable
__version__ = "0.1.a2"
__all__ = ["Fast5Parser", "Fast5", "Helper_fun"]

description = 'Fast5 tools is a collection of tools to manipulate Fast5 files'
long_description = """"""

# Collect info in a dictionary for setup.py
setup_dict = {
    "name": __name__,
    "version": __version__,
    "description": description,
    "long_description": long_description,
    "url": "https://github.com/a-slide/Fast5Tools",
    "author": 'Adrien Leger',
    "author_email": 'aleg {at} ebi.ac.uk',
    "license": "GPLv3",
    "python_requires":'>=3.3',
    "classifiers": [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',],
    "install_requires": ['h5py>=2.7.1', 'numpy>=1.8.1', 'matplotlib>=2.0.0'],
    "packages": [__name__]}
