deepSNV
=====

Description
---
This package provides provides a quantitative variant callers for
    detecting subclonal mutations in ultra-deep (>=100x coverage) sequencing
    experiments. The deepSNV algorithm is used for a comparative setup with a
    control experiment of the same loci and uses a beta-binomial model and a
    likelihood ratio test to discriminate sequencing errors and subclonal SNVs.
    The new shearwater algorithm (beta) computes a Bayes classifier based on a
    beta- binomial model for variant calling with multiple samples for
    precisely estimating model parameters such as local error rates and
    dispersion and prior knowledge, e.g. from variation data bases such as
    COSMIC.
    

Note
----
This repository contains the current development snapshot of the deepSNV package 
in the folder deepSNV. It is not guaranteed to work all times.

Installation
--------

### The good way
For unexperienced users it is recommended to use package as provided in bioconductor.
Please follow the instructions at:
http://master.bioconductor.org/packages/devel/bioc/html/deepSNV.html

### The bad way
You can use devtools::github_install() to install from this repository. For advanced users.

	> library(devtools); install_github("mg14/mg14")

### The ugly way
To install this development snapshot of deepSNV, check out the repository and run

	$ make install

Note that this will not install the necessary dependencies.

Previous (running) builds can be found in builds. These will also be mirrored to the
bioconductor development branch.