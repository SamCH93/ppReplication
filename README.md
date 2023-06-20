# Power Priors for Replication Studies

This repository contains code and data related to the preprint 

Pawel, S., Aust, F., Held, L., and Wagenmakers, E.-J. (2022). Power Priors for
Replication Studies.
doi:[10.48550/arXiv.2206.xxxxx](https://doi.org/10.48550/arXiv.2206.xxxxx)

## Reproducing the results

We offer two ways to reproduce the results

### 1. Reproduction with local computational environment (requires R and LaTeX)

First install the required R packages by running in a shell from the root
directory of the repository

``` sh
## packages from CRAN
R -e 'install.packages(read.delim("CRANpackages.txt", header = FALSE)[,1])'
## requires remotes package
R -e 'remotes::install_github(repo = "SamCH93/ppRep")'
```

Then run

``` sh
cd paper
make pdf
```

this should reproduce all analyses and output the file `ppreplication.pdf` in
the paper directory.

Although our analysis depends on only few dependencies, this approach may lead
to different results (or not even compile successfully) in the future if R or an
R package dependency changes. The R and R package versions which were used in
our analysis can be seen in the output of the sessionInfo command at the bottom
of the manuscript in the snapshot of the GitHub repository at the time of
submission.

### 2. Reproduction within Docker container (requires Docker with root rights)

Run in a shell from the root directory of the repository

``` sh
make drunpdf
```

this should output the file `ppreplication.pdf` in the paper directory. The
Docker approach takes a bit longer but reruns our analyses in a Docker container
which encapsulates the computational environment (R and R package versions) that
was used in the original analysis. 
