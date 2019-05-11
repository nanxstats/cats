# cats

[![Travis build status](https://travis-ci.org/nanxstats/cats.svg?branch=master)](https://travis-ci.org/nanxstats/cats)

`cats` is a fork of the R package `CATS` developed by Anders Albrechtsen <<aalbrechtsen@bio.ku.dk>> for power estimation in two-stage genome-wide association designs.

## Installation

```r
# install.packages("devtools")
devtools::install_github("nanxstats/cats")
```

## Method

The method for power analysis implemented in `CATS` is based on Skol et al. 2006, but more generalized so that the ratio between cases and controls may vary between stages. The allele frequencies, disease prevalence, and relative risk can also vary.

## Changes

The changes made to the original R package `CATS` are only for fixing the package build and check issues (under R 3.6.0 as of May 2019), for example:

- Used unified lowercase names for package and functions
- Fixed C code build errors and warnings (now passes clang)
- Fixed other package formatting issues that blocks building and checking
- Used roxygen2 to automate namespace and documentation generation
- Organized the existing code examples as a vignette

## Source

This package is forked and modified from these resources. They are also available as archives under `bin/`:

- [Code examples](http://www.popgen.dk/software/index.php/CATS)
- [Source tarball](http://popgen.dk/software/download/CATS/CATS_1.02.tar.gz)
- [PDF manual](http://popgen.dk/software/download/CATS/CATS-manual.pdf)

## References

- Skol AD, Scott LJ, Abecasis GR, Boehnke M: Joint analysis is more efficient than replication-based analysis for two-stage genome-wide association studies. Nat Genet 38: 209-213, 2006.
