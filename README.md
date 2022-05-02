# Power Priors for Replication Studies

This repository contains R code and other files to the arXiv preprint
[XXXX.XXXXX](https://arxiv.org/abs/XXXX.XXXXX).

## Reproducing the results

1. Install the `ppRep` package by running
```r
remotes::install_github("SamCH93/ppRep", , subdir = "pkg")
```
from an R session. Alternatively, the package can be built and installed by
running 
```shell
make build install
```
in from the repository root folder in a shell.

2. Install all the dependencies from CRAN by running
```r
install.packages(c("ggplot2", "colorspace", "xtable", "dplyr", "hypergeo" "ReplicationSuccess"))
```
from an R session.

3. Run in the root folder of the repository from a shell
```shell
cd paper
make pdf
```
