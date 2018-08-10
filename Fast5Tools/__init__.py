# -*- coding: utf-8 -*-

# Define self package variable
__version__ = "0.3.3"
__all__ = ["Fast5Parse", "Fast5Wrapper", "Fast5", "Basecall", "Alignment", "Eventalign"]

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
    "install_requires": ['h5py>=2.7.0', 'numpy==1.14.0', 'matplotlib>=2.0.0', 'pysam>=0.12.0', 'pandas>=0.23.0'],
    "packages": [__name__],
    "entry_points":{'console_scripts': [
        'Fast5Tools= Fast5Tools.Fast5Tools_Main:main']}}
