DEVELOPMENT VERSION: abmR: An R Package for Agent-based Model Analysis of Large-scale Movements Across Taxa
======================================================================================

#### Benjamin Gochanour, Javi Fernandez Lopez, Andrea Contina

#### 2021-04-28


Getting Started
==================

Installation
------------

To use `abmR`, you must first install it from Github using `devtools`
and load the library:

``` r
devtools::install_github("bgoch5/abmR")
# If install gives errors, try running the following:
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
library(abmR,quietly=TRUE,warn.conflicts=FALSE)
```
While this package is still in development, it will be updated frequently, so please be sure to re-install frequently. Installing  abmR will also automaticaly install its dependencies, if you don’t already have them installed. These include `raster`, `sp`, `table1`, `googledrive`, `swfscMisc`, `geosphere`, `kableExtra`, and `gtsummary`, and `ggplot`.

