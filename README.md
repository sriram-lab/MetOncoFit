---
layout: default
title: MetOncoFit
---

## MetOncoFit: A Machine Learning Algorithm to Identify Common Biochemical and Topological Attributes of Metabolic Genes Recurrently Dysregulated in Tumors
 MetOncoFit is a data-driven approach that uses biochemical and metabolic attributes to predict tumor differential expression, copy number variation, and patient survival.

Visit the documentation to start using MetOncoFit and explore the API. Feel free to submit bugs, feature requests, and Pull requests.

To support the MetOncoFit project, you can cite our publication:
Oruganty, K., Campit, S.E., Mamde, S., & Chandrasekaran, S. Common biochemical and topological attributes of metabolic genes recurrently dysregulated in tumors. *in preparation*

## Installation
MetOncoFit is a Python package that takes curated cancer metabolic models to output differential expression, copy number variation, and patient survival classification predictions and ranked metabolic/biochemical feature importances.

We recommend installing MetOncoFit inside a [virtual environment](https://virtualenv.pypa.io/en/latest/). To install the package, use PyPI:

``` Python
pip install metoncofit
```

## Contributing
Contributions are welcome! Please read the contributions guide to get started.

## MetOncoFit
{% include Breast_TCGAannotation.html %}

<div class="custom-select" style="width:200px;">
  <select>
    <option value="0">Cancer Type:</option>
    <option value="1">Breast</option>
    <option value="2">Glioma</option>
    <option value="3">Colorectal</option>
    <option value="4">B-cell Lymphoma</option>
    <option value="5">Lung</option>
    <option value="6">Skin</option>
    <option value="7">Renal</option>
    <option value="8">Prostate</option>
    <option value="9">Ovarian</option>
    <option value="10">Pan</option>
  </select>
</div>

<div class="custom-select" style="width:200px;">
  <select>
    <option value="0">Prognostic Marker:</option>
    <option value="1">Differential Expression</option>
    <option value="2">Copy Number Variation</option>
    <option value="3">Patient Survival</option>
  </select>
</div>
