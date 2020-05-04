# MetOncoFit
MetOncoFit is a random forest algorithm that uses biochemical and metabolic attributes to predict tumor differential expression, copy number variation, and patient survival.

## Introduction
Tumors reprogram normal cellular metabolism to support uncontrolled proliferation. While some of these metabolic reprogramming strategies are common across most tumors, such as the Warburg effect, there must be diverse metabolic objectives that contribute to tumor heterogeneity.

We hypothesized that cancer cells have few key changes in the metabolic network, and examined frequently dysregulated metabolic genes using a multi-scale systems biology approach to determine commone features that contribute to metabolic dysregulation in tumors. 

Our tumor models contain a broad range of metabolic attributes, including:
  * Enzyme catalytic activity
  * Gene expression and copy number variants
  * Metabolic pathway membership
  * Metabolic subnetwork information
  * Topological connectivity to biomass and medium components
  * Metabolic fluxes obtained from *in silico* knockout experiments.
  
Our study demonstrates how biochemical and metabolic network features are predictive of metabolic gene dysregulation across several cancer types. 

## Installation
To install MetOncoFit, you can fork this GitHub repository onto your local machine. All outputs will be contained in this folder.

We recommend creating a [virtual environment](https://virtualenv.pypa.io/en/latest/) specifically for MetOncoFit and install all the requisite packages in the `requirements.txt` file supplied in the repository. From here, you will be able to run the MetOncoFit package, which is located in the `metoncofit` folder. 

## Usage
To run MetOncoFit with the nine tumor models, you can run the `runMetOncoFit.sh` file. which generates the trained tumor models and outputs the figures in the manuscript.

## Contributing
Contributions are welcome! Please read the contributions guide to get started. Also feel free to submit bugs, feature requests, and pull requests.

Otherwise, you can support the MetOncoFit project by citing our publication:
Oruganty K*, Campit S*, Mamde S, Lyssiotis C, Chandrasekaran S, Common biochemical properties of metabolic genes recurrently dysregulated in tumors, Cancer and Metabolism, 2020 (in press)

