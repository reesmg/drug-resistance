# packages required
options(repos=structure(c(CRAN="http://cran.r-project.org")))

install.packages(c("ppcor",
                   "reshape2",
                   "ggplot2",
                   "gghighlight",
                   "glmnet",
                   "magrittr",
                   "dplyr",
                   "devtools",
                   "data.table"))

# to install cdsrmodels see the instructions on github: https://github.com/broadinstitute/cdsr_models
library(devtools)
devtools::install_github("broadinstitute/cdsr_models")
