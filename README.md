#  Rare Variant Association Toolkit (RVAT)

Welcome to the RVAT package! 

## Introduction

RVAT was developed to provide an easy to use computational framework for rare variant analyses in big DNA sequencing datasets. We provide a command-line interface that is tailored to large scale analyses on high performance computing clusters, and an R interface that is tailored to interactive use on local machines. RVAT provide a single framework that simplifies the integration of human DNA sequencing datasets with complex annotations (e.g. variant effect predictions and sample phenotypes), whilst also facilitating efficient data querying (eg "lookups" of variants causing specific biological disruptions) and providing user friendly interfaces to an extensive range of rare variant single, aggregate and gene set analyses.

To date RVAT been used to conduct gene discovery analyses in large exome sequencing datasets (n=100,000 samples), to identify disease associated non-coding regulatory elements in large genome sequencing datasets (n=10,000 samples) and to conduct a variety of in-depth rare variant quality control optimizations, rare variant association tests, gene set enrichment analyses, cell type enrichment analyses, spatial clustering analyses, modifier analyses and subcohort analyses. Users can begin their analyses from standard vcf files, import external annotations and then explore their results using our collection of supporting R functions, our Shiny graphical user interface or by passing our Bioconductor compatible RVAT objects / exportable tables to their own custom code and preferred third party tools. RVAT was also designed to make it easy to export small subsets of larger datasets whilst still maintaining the integration of sequencing genotypes with imported annotations and downstream results. This makes it easier to share summary statistics along with the relevant parts of upstream datasets, and facilitates easy downloading of prioritized genomic regions for detailed exploratory analyses on local laptops.
  
## Installation

### Install directly from github

```
remotes::install_github("kennalab/rvat")
```

Many of the examples and tutorials rely on example data which can be installed similarly:

```
remotes::install_github("kennalab/rvatData")
```

## Tutorials

Tutorials and documentation: https://kennalab.github.io/rvat/
