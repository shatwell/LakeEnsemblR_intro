# Introduction to ensemble lake modelling with LakeEnsemblR

<a href="url"><img src="logo.png" align="right" height="180" width="180"/></a>

------------------------------------------------------------------------

:spiral_calendar: October 13, 2022

:alarm_clock: 9:00-13:00 

:busts_in_silhouette: Tom Shatwell, Karsten Rinke 

:computer: [Material](https://github.com/shatwell/LakeEnsemblR_intro)

:octocat: [Homepage](https://github.com/aemon-j/LakeEnsemblR)

:page_facing_up: [Journal article](https://doi.org/10.1016/j.envsoft.2021.105101)


------------------------------------------------------------------------

## Description

A guided walkthrough of the [R package
"LakeEnsemblR"](https://github.com/aemon-j/LakeEnsemblR) that has been
developed within the GLEON modelling group. Learn how to set up multiple
lake models for your lake physical system, how to calibrate your models
and explore model uncertainty. Brainstorming of further ideas for how
this package can be applied to aid further research. Everyone is welcome
to attend, but the workshop is mostly aimed at people with at least some
experience in lake modelling.

## Acknowledgement

Much of this workshop material was developed by the LakeEnsemblR core
developers Jorrit Mesman, Johannes Feldbauer, Tadhg Moore, and Robert
Ladwig. It was presented at GLEON21.5 in the GSA workshop series. We are
very grateful for their help. 

:computer: [Original GSA workshop material](https://github.com/shatwell/LakeEnsemblR_intro)


## What will this workshop cover?

* Introduction to hydrodynamic modelling of lakes

* Introduction to the FLake model

* Introduction to LakeEnsemblR package: 
- Why use ensembles? 
- What is LakeEnsemblR?

* Using LakeEnsemblR: 
- Standardisation of input data 
- Functions 
- Visualising output & calibration 
- Apply it to YOUR lake! (or on OUR examples)

## Prerequisites

### 1. Install the required software

In this workshop, you will need [R](https://www.r-project.org/) (version>=3.5) and a GUI of your choice (preferably
[RStudio](https://www.rstudio.com/products/rstudio/download/) - desktop
free version). You will also need a decent text editor (e.g.
[Notepad++](https://notepad-plus-plus.org/downloads/)) and netcdf viewer
like PyNcView. Windows users can download PyNcView for example
[here](https://getwinpcsoft.com/PyNcView-2257247/). Linux users can
install PyNcView via [PyPi](https://pypi.python.org/pypi), then with the
command ```pip install pyncview```.

### 2. Set up the LakeEnsemblR R-package 

Clone or download files from this Github repository. To do this, first install the remotes, ggplot2, ggpubr and reshape packages in R.
```
install.packages("remotes")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("reshape")
```
Next, install the specific LakeEnsemblR packages from github:
```
remotes::install_github("GLEON/rLakeAnalyzer")
remotes::install_github("aemon-j/GLM3r", ref = "v3.1.1")
remotes::install_github("USGS-R/glmtools", ref = "ggplot_overhaul")
remotes::install_github("aemon-j/FLakeR", ref = "inflow")
remotes::install_github("aemon-j/GOTMr")
remotes::install_github("aemon-j/gotmtools")
remotes::install_github("aemon-j/SimstratR")
remotes::install_github("aemon-j/MyLakeR")
remotes::install_github("aemon-j/LakeEnsemblR")
```


#### Update: If you experience problems on macOS (we tested the model
binaries only for macOS Catalina) with error messages like 'dyld:
Library not loaded', you can try the following approach:

-   install the missing libraries, e.g. by using
    ['brew'](https://brew.sh):
    `/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"`
-   then you will need the missing libraries or you will need to update
    your current versions of gfortran, netcdf and gc:
    `brew install gcc`, `brew install netcdf`, `brew install gc`
-   check if everything is working: `gcc -v` and `gfortran --version`
    which should give you

> "Apple clang version 12.0.0 (clang-1200.0.32.2)"

and

> "GNU Fortran (Homebrew GCC 10.2.0) 10.2.0"

-   if you still experience library problems, you can continue
    installing these missing dependencies using 'brew'

## Files

### Workshop Files 

The HTML and PDF files in the repository both
contain the information needed for the workshop. You can pick which one
you prefer. You can copy the code into an R script and run it yourself,
or use the pre-made R script ("InventWater_LakeEnsemblR.R"), which
contains all the lines of code that are run during the workshop and some
short comments.

### LakeEnsemblR Examples 

There is a selection of lakes of different
areas, depths and climatic zones that have been collated to show you
different applications of ```LakeEnsemblR```. They can be downloaded from
the [LER_examples](https://github.com/aemon-j/LER_examples) repository
on the "aemon-j" GitHub account.

------------------------------------------------------------------------
