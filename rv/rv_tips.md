# Getting started
Use `rv init` to start a rv project

This will use whatever version of R you have loaded when you initiate the project

This creates files and folders:
- rv/                  Folder to hold your R environment
  - scripts/           Folder to hold scripts to activate the environment
    - activate.R
    - rvr.R
  - .gitignore         File to ignore the R environment folder in git
- .Rprofile            Forces R to use the rv environment when you start R in this folder
- rv.lock              File to lock the versions of packages you have installed
- rproject.toml        File to specify the packages you want to use in your project

Add repositories to the rproject.toml file such as:

    { alias = "CRAN", url = "https://cran.rstudio.com/" },
    { alias = "PPM", url = "https://packagemanager.posit.co/cran/latest" },
    { alias = "bioconductor", url = "https://bioconductor.org/packages/3.22/bioc" },
    { alias = "BioCann", url = "https://bioconductor.org/packages/3.22/data/annotation"},
    { alias = "BioCexp", url = "https://bioconductor.org/packages/3.22/data/experiment" },
    { alias = "BioCworkflows", url = "https://bioconductor.org/packages/3.22/workflows"},

I don't know why these aren't pre-populated.

# Adding packages
You can do this in two ways. The preferred method is to edit the rproject.toml file and add the packages you want to use. This allows the dependencies to be resolved all at once and ensures that you have compatible versions of all packages. 

The other method is to install the packages using `rv add` on the command line. This will install the packages one at a time and may lead to version conflicts if the packages have incompatible dependencies.

# Testing package compatibility
Use `rv plan` to make sure it can find compatible versions. This runs a dry run of the installation and checks for any version conflicts. If there are no conflicts, it will show you the packages that will be installed and their versions. If there are conflicts, it will show you which packages are causing the conflicts and what versions are incompatible.

# Installing packages
Use `rv sync` to install everything specified in the rproject.toml file. This will install all packages and their dependencies, and update the rv.lock file with the versions that were installed. 

# The toml formatting is strict
Failed to load config at `.` likely means a syntax error in rproject.toml

# Making sure your environment is reproducible
make sure to add these files to your github:
- rv/scripts/activate.R
- rv/scripts/rvr.R
- rv/.gitignore
- .Rprofile
- rv.lock
- rproject.toml

When cloning a new repo that used rv:
- clone the repo
- cd into the repo
- rv sync to install all packages

# Getting rid of it
https://tenor.com/bgZks.gif

Delete these files and folders:
- rv/
- .Rprofile
- rv.lock
- rproject.toml














    { name = "Matrix", git = "https://github.com/cran/Matrix", tag = "1.6-5" },
    "tidyverse",
    { name = "ggrepel", git = "https://github.com/slowkow/ggrepel", tag = "0.9.6"},
    "Seurat",
    "SeuratObject",
    "SingleCellExperiment",