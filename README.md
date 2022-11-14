[![Build Status](https://app.travis-ci.com/Ali-Mahzarnia/fGMD.svg?branch=master)](https://app.travis-ci.com/Ali-Mahzarnia/fGMD)
# fGMD

This package is built exclusively for the [MFSGrp](https://github.com/Ali-Mahzarnia/MFSGrp) package. This package is a heavily modified version of the [gglasso](https://github.com/cran/gglasso) that was initiated by [Yang Y.and Zou, H. (2015)](http://users.stat.umn.edu/~zouxx019/Papers/gglasso-paper.pdf).The features added to the original gglasso package are: the mixing parameter (alpha) and its net search cross-validation, the curvature penalization for functional regression and its net search cross-validation, the optimized Fortran core function to accelerate the curvature penalization updates, and the progress reports with time estimations. This package does not work independently from the [MFSGrp](https://github.com/Ali-Mahzarnia/MFSGrp) package, and it does not interfere with the functions of gglasso package due to the slight name differences.

## Installation
You can install this package, `fGMD` from [GitHub](https://github.com/Ali-Mahzarnia/fGMD) with the R command:
``` R
install.packages("https://github.com/Ali-Mahzarnia/fGMD/archive/master.tar.gz", repos = NULL, type="source")
```
## Development version
Alternatively, the development version of `fGMD` can be installed with the R command:
``` R
install.packages("pacman")
pacman::p_install_gh("Ali-Mahzarnia/fGMD")
```
## Installation from the source file
If the installation fails with the above methods, install the package from the source file:
``` R
install.packages("https://github.com/Ali-Mahzarnia/fGMD/raw/master/fGMD_1.0.tar.gz",  repos = NULL, type="source")
```
