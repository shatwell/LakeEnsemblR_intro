# install_packages.R
# R version >3.5 (Tested on 4.1.2)
# 2022-09-08 Tom Shatwell

## install required libraries
install.packages("devtools")
install.packages("remotes")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("reshape")

remotes::install_github("GLEON/rLakeAnalyzer")
remotes::install_github("aemon-j/GLM3r", ref = "v3.1.1")
remotes::install_github("USGS-R/glmtools", ref = "ggplot_overhaul")
remotes::install_github("aemon-j/FLakeR", ref = "inflow")
remotes::install_github("aemon-j/GOTMr")
remotes::install_github("aemon-j/gotmtools")
remotes::install_github("aemon-j/SimstratR")
remotes::install_github("aemon-j/MyLakeR")
remotes::install_github("aemon-j/LakeEnsemblR")

# END
