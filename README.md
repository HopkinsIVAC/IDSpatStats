# IDSpatStats



## Rebuilding

To rebuild `IDSpatialStats` after modification, follow these steps:


1. Update `src/IDSpatialStats_init.c`:

    - Run `tools::package_native_routine_registration_skeleton(".")`
    - Copy output in console and replace the full contents of `src/IDSpatialStats_init.c` with the output, saving line 1.
    
2. Update documentation & Check the package:

    - Run `devtools::check(document = TRUE)
    - Check the output for errors. Check that warnings and notes are not critical.
    
3. Push to github 




## Installing

    *NOTE* - Sometimes you need to uninstall the old version. To do this run:
    
    `detach("package:IDSpatialStats")`
    `remove.packages("IDSpatialStats")`
        

    Run `devtools::install_github("https://github.com/HopkinsIVAC/IDSpatStats")` to install from GitHub.
    
    Make sure other instances of RStudio are closed and your .dll files are not being synced by OneDrive or Dropbox or anything.
    
    


### Notes:

As with any R package, make sure you are exporting any new R functions. To do this, make sure to include `@export` in the roxygen2 documentation at the start of the function. 