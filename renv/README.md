# TLDR

- Use `renv::init()` to create a new renv environment in the project directory.
- Install packages as normal.
- Packages used in your code need to be added to the lock file.
- Use `renv::snapshot()` to update the lock file with the packages you used.
    - Do this frequently to keep the lock file up to date.
    - Track these files in your git repo:
        - renv.lock
        - .Rprofile
        - renv/.gitignore
        - renv/activate.R
        - renv/settings.json
- Use `renv::restore()` to install all the packages in the lock file to your environment.

# Setup
To get started, initiate the {renv} environment by running `renv::init()` in the R console. This will create a new renv environment in the project directory. This creates metadata files and a folder called `renv` in the project directory that will hold your packages.

Be careful. If you have code in your project that uses any non-base R packages, renv will try to install these when you call `renv::init()`. If you're copying code over from another project and you're not going to use all the code in it, this may install a lot of packages you don't need. You can always remove them later, but it's something to keep in mind.

Make sure you're using the proper R version for the project. We should likely make a note of the R version in the main qmd file.

When using VScode, make sure you enable the "Use Renv Lib Path" setting in the R extension settings, to avoid getting errors about jsonlite not being installed even though it totally is. You'll need to restart vscode for this to take affect.

# Installing packages
Keep in mind you may need to load modules before you install some packages. See https://gist.github.com/MVesuviusC/8aa0cce0a01a72caeacbd048f8046d08 for example.

For the example in this repo, we need to `module load GCC/9.3.0 OpenMPI/4.0.3 R/4.3.0 OpenSSL/1.1` on the Franklin cluster.

Install packages as normal. I like pak::pkg_install() for this.

One nice thing about this. The packages are actually installed in ~/.cache/R/renv/cache/..... which means that if you have multiple projects that use the same package, you only need to install it once and you're not wasting disk space. Keep in mind, though, that each version of each package is stored separately.

# Removing packages
If you no longer need a package, remove all references to it in your code, then run `renv::snapshot()` to update the lock file. This will remove the package from the lock file. You can then call `renv::remove("package_name")` to remove the package from the environment.

# Using renv
Once you have the packages you need, renv will add any package (and all dependencies) to the renv lock file, but only for packages that are used in the project. This includes both packages loaded through `library()` and packages used in the code but not loaded with the `ggplot2::geom_point()` syntax.

Once you have the packages installed and code written in your project, run `renv::snapshot()` to update the lock file with the packages you need. You can also run `renv::status()` to see if any updates are needed. It will be a good idea to run `renv::status()` before each git commit to make sure you have all the packages you need noted in your lock file.

# Re-installing an environment from a lock file
If you are before installing all the packages with `renv::restore()`

To make your environment reproducible by others, add renv.lock, .Rprofile, renv/.gitignore, renv/activate.R and renv/settings.json to your git repo. You don't need to add any other files in the renv folder.

# Getting rid of renv
The easiest way is to run `renv::deactivate(clean = TRUE)`. This deletes the renv folder and renv.lock file. You can also just delete the renv folder and renv.lock file manually. You'll need to restart R afterwards.
