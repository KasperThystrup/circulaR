
<!-- README.md is generated from README.Rmd. Please edit that file -->

# circulaR

<!-- badges: start -->
<!-- badges: end -->

The goal of circulaR is to identify and annotate circular RNA (circRNA)
contents from total RNA sequencing data. The package contains
functionality for importing chimeric read data generated with the STAR
aligner \[link: Alexander Dobin - STAR github\], filtering backsplice
junctions, and overlapping known gene models with the identified
backsplice junctions.

## Installation

### Software

In order to use the tools provided in this guide, and to install the
required dependencies, the following external software must be
installed:

-   libgit2: [Ubuntu 20.04
    LTS](https://packages.ubuntu.com/source/focal/libgit2) [Arch
    Linux](https://archlinux.org/packages/extra/x86_64/libgit2), [Mac
    Homebrew](https://formulae.brew.sh/formula/libgit2)
-   gcc-fortran: [Ubuntu 20.04
    LTS](https://packages.ubuntu.com/focal/gfortran) [Arch
    Linux](https://archlinux.org/packages/core/x86_64/gcc-fortran) [Mac
    installation
    instructions](https://gcc.gnu.org/wiki/GFortranBinariesMacOS)
-   R: [Ubuntu](https://cran.r-project.org/bin/linux/ubuntu/) \[Arch
    linux\] (<https://archlinux.org/packages/extra/x86_64/r/>)
    [Mac](https://cran.r-project.org/)
-   RStudio: [Ubuntu &
    Mac](https://rstudio.com/products/rstudio/download/) [Arch
    linux](https://aur.archlinux.org/packages/rstudio-desktop-bin/)

### R packages

In order to install the `circulaR` R package, you must ensure to install
the following Biocodncutor packages:

-   AnnotationDbi
-   BiocGenerics
-   BSgenome
-   DESeq2
-   ensembldb
-   GenomeInfoDb
-   GenomicFeatures
-   GenomicRanges
-   Gviz
-   IRanges
-   Rsamtools
-   S4Vectors

``` r
install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "BiocGenerics", "BSgenome", "DESeq2", "ensembldb", "GenomeInfoDb", "GenomicFeatures", "GenomicRanges", "Gviz", "IRanges", "Rsamtools", "S4Vectors"))
```

In addition, `devtools` r package must be installed, to enable easy
installation of `circulaR` from github.

``` r
install.packages("devtools", dependencies = TRUE)
devtools::install_github("https://github.com/KasperThystrup/circulaR")
```

#### Trobuleshooting

This should install all dependencies from the official CRAN repository.
However, if something does not work out, try to install the CRAN
dependencies manually before running the install\_github command again:

-   dplyr
-   ggplot2
-   parallel
-   pbmcapply
-   plyr
-   readr
-   RSQLite
-   stringi
-   tibble
-   tidyr

``` r
install.packages(c("dplyr", "ggplot2", "parallel", "pbmcapply", "plyr", "readr", "RSQLite", "stringi", "tibble", "tidyr", "Hmisc"))
devtools::install_github("https://github.com/KasperThystrup/circulaR")
```
