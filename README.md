# Robust longitudinal Differential Expression (RolDE)

## Introduction

Robust longitudinal Differential Expression (`RolDE`), is a composite method, consisting of three independent modules with different approaches to detecting longitudinal differential expression. The **RegROTS** module combines individual regression modelling with the power of the established differential expression method Reproduciblity Optimized Test Statistic (ROTS). In the **DiffROTS** module, the expression between all the individuals in the different conditions is directly compared at all timepoints. The **PolyReg** module uses polynomial regression modelling to evaluate longitudinal differential expression. The combination of these modules allows RolDE to robustly detect differences in longitudinal trends and expression levels in diverse data types and experimental settings.

Contrary to many existing approaches, the `RolDE` does not require prior knowledge concerning the types of differences searched, but can easily be applied even by non-experienced users.


## Installation

####  GitHub

The latest version of `RolDE` can be downloaded from GitHub using the devtools R package:

`devtools::install_github(repo = "https://github.com/elolab/RolDE", build_vignettes = T)`

The `build_vignettes` parameter should be set to `TRUE` to enable the building of instructions and examples for RolDE usage.

####  Bioconductor

The latest release of RolDE can also be installed from Bioconductor:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RolDE")
```


## Instructions

`RolDE` comes with a detailed documentation and a detailed set of instuctions. The documentation of `RolDE` can be simply accessed after installation and loading the package by:

`?RolDE`

For detailed instuctions and examples on how to use RolDE, please view the included vignettes:

`browseVignettes("RolDE")`

## Contact information

If you have comments regarding `RolDE`, please contact us [here](https://github.com/elolab/RolDE/issues). 
