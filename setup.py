# -*- coding: utf-8 -*-
import sys
from sys import argv, version_info
import re
from subprocess import call
import threading
import webbrowser
import os
from shutil import copy, move, rmtree
from os.path import join, dirname, realpath, exists
from setuptools import setup, find_packages, Command

directory = dirname(realpath(__file__))
sys.path.insert(0, join(directory, 'metoncofit'))
version = __import__('version').__version__
full_version = __import__('version').__version__
package = __import__('version').package

setup(
    name = 'metoncofit',
    version = full_version,
    author = package['author'],
    url = package['homepage'],
    description=package['description'],
    keywords=', '.joint(package['keywords']),
    license=package['license'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Visualization',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Operating System :: OS Independent'
    ],
    packages = find_packages(),
    ubclude_package_data=True,
    data_files=[
        (
            ''
        )
    ]
    setup_requires = setup_requirements,
    install_requires = [
        "numpy"
        "pandas"
        "scipy"
        "scikit-learn"
        "imbalanced-learn"
        "matplotlib"
        "seaborn"
    ],

    author_email = "csriram@umich.edu"
    maintainer = "Scott Edward Campit"
    maintainer_email = "scampit@umich.edu"
    license = "GNU 3.0"
    keywords = ["Metabolic Modeling", "Cancer metabolism", "Data-Driven"]
    platforms = "GNU/Linux", "Mac OS X", "Microsoft Windows"
)
