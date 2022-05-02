# Power Priors for Replication Studies

This repository contains R code and other files related to the arXiv preprint
[XXXX.XXXXX](https://arxiv.org/abs/XXXX.XXXXX).

## Reproducing the results

1. Install the `ppRep` package by running from an R session
```r
remotes::install_github("SamCH93/ppRep", , subdir = "pkg")
```
This requires the `remotes` package which is available on CRAN.
Alternatively, the package can be built and installed by
running the following command in the repository root folder in a shell
```bash
make build install
```


2. Install all the dependencies from CRAN by running from an R session
```r
install.packages(c("ggplot2", "colorspace", "xtable", "dplyr", "hypergeo" "ReplicationSuccess"))
```


3. Run in the root folder of the repository from a shell
```bash
cd paper
make pdf
```
