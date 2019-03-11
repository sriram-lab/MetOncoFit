# MetOncoFit - Common biochemical and topological attributes of metabolic genes recurrently dysregulated in tumors

## What is MetOncoFit?
Cancer cells rewire normal cellular metabolism to support cellular differentiation and growth. While several pathways have been associated with cancer metabolism and its increased cellular fitness, there have been no methods to date that describe how metabolic rewiring contributes to cancer. MetOncoFit is a hybrid approach that combines constraint-based genome scale modeling and data-driven approaches to identify metabolic pathways and genes associated with specific cancer tissues.

The documentation is available online at readthedocs.

## Installation
We recommend installing MetOncoFit inside a virtual environment. To install the MetOncoFit package, use PyPI:

``` Python
pip install metoncofit
```

There are 5 core modules in the MetOncoFit package:
  * processing.py: handles function for data processing.
  * random_forest.py: contains code for the random forest classifier and to save/load pickled models.
  * validator.py: outputs statistical measures to assess the model's performance.
  * visualization.py: outputs the figures as .svg files - can be modified to support .jpeg or .png.
  * metoncofit.py: the main module that runs all the functions.
  * do_all.sh: the bash script that runs it all.

## Contributing
Contributions are welcome! Please read the contributions guide to get started.
