# -*- coding: utf-8 -*-
from warnings import warn
from sys import argv, version_info

from setuptools import setup, find_packages

if __name__ == '__main__':
  setup(
    name = 'metoncofit',
    version = '0.1',
    packages = find_packages(),
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
    authors = ["Krishna Dev Oruganty", "Scott Edward Campit", "Sainath Mamade", "Sriram Chandrasekaran"]
    author_email = "csriram@umich.edu"
    maintainer = "Scott Edward Campit"
    maintainer_email = "scampit@umich.edu"
    license = "GNU 3.0"
    keywords = ["Metabolic Modeling", "Cancer metabolism", "Data-Driven"]
    platforms = "GNU/Linux", "Mac OS X", "Microsoft Windows"
  )
