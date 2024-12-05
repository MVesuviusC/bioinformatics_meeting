# Install everything needed

## Load modules on Franklin:
### ml GCC/9.3.0 OpenMPI/4.0.3 R/4.3.0
### ml miniforge3 OpenSSL/1.1
 # If using a different R version, may need a different version of libffi and
 # SQLite to make sure GCC/**** dependency matches
### ml libffi/3.3
### ml SQLite/3.31.1

#### I still had an error for libffi, which I ignored and it worked fine

## Start R

# install.packages("remotes")
remotes::install_github("rstudio/tensorflow")

library(reticulate)
install_python(version = "3.11")

library(tensorflow)
#install_tensorflow(envname = "r-tensorflow", python_version = "3.9")
#
#install.packages("keras")
library(keras)
# We may want to specify conda = "/path/to/miniforge3/bin/conda" - not sure what it does by default
install_keras(
    method = "conda",
    envname = "r-tensorflow",
    python_version = "3.10"
)

# check if it works
tf$constant("Hello TensorFlow!")
