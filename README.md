# MetOncoFit - Common biochemical and topological attributes of metabolic genes recurrently dysregulated in tumors

 MetOncoFit is a data-driven approach that uses biochemical and metabolic attributes to predict tumor differential expression, copy number variation, and patient survival.

<<<<<<< HEAD
Visit the documentation to start using MetOncoFit and explore the API. Feel free to submit bugs, feature requests, and Pull requests.

To support the MetOncoFit project, you can cite our publication:
Oruganty, K., Campit, S.E., Mamde, S., & Chandrasekaran, S. Common biochemical and topological attributes of metabolic genes recurrently dysregulated in tumors. *in preparation*

## Installation

MetOncoFit is a Python package that takes curated cancer metabolic models to output differential expression, copy number variation, and patient survival classification predictions and ranked metabolic/biochemical feature importances. We recommend installing MetOncoFit inside a virtual environment. To install the package, use PyPI:
=======
The documentation is available online at readthedocs.

## Installation
We recommend installing MetOncoFit inside a virtual environment. To install the MetOncoFit package, use PyPI:
>>>>>>> bd527f7e6212c61886576c1b7a5c0aeb660e0629

``` Python
pip install metoncofit
```
<<<<<<< HEAD
=======

There are 5 core modules in the MetOncoFit package:
  * processing.py: handles function for data processing.
  * random_forest.py: contains code for the random forest classifier and to save/load pickled models.
  * validator.py: outputs statistical measures to assess the model's performance.
  * visualization.py: outputs the figures as .svg files - can be modified to support .jpeg or .png.
  * metoncofit.py: the main module that runs all the functions.
  * do_all.sh: the bash script that runs it all.

## Contributing
Contributions are welcome! Please read the contributions guide to get started.
>>>>>>> bd527f7e6212c61886576c1b7a5c0aeb660e0629
