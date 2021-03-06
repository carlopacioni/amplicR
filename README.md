# amplicR 
An R package to process amplicon data

This package has a number of functions to filter, dereplicate and error and 
chimera check NGS reads stored in fastq files. This funcitonality is currently
limted to single-read analysis. The retained reads after the data processing 
described above (or your own reads if you have done this already in another way) 
can then be compared against reference
sequences and the number of mismatch is reported. This may be useful when, for
example, screening samples for particular taxa (e.g. a pathogen).

`amplicR` is mainly a wrapper of several functions provided by other R packages 
in order to automate some common analyses, so if you use it, please make sure to
also cite the relevant packages (these are generally indicated in the 
documentation of `amplicR`'s functions). To find out the correct citation for a 
package, you can use the function: `citation("package_name")` where you have to 
replace package_name with the actual name of the package you are interested 
into.

## Quickstart 
Install the package from version control from within R: 
``` 
library(devtools) 
install_github("carlopacioni/amplicR") 
``` 
If you have not used `devtools` before, then you have to install it with  
```
install.packages("devtools") 
```
If you are on Windows, before loading `devtools`, shut down R, download the 
Rtools executable file from CRAN webpage and run it.  

Note that the Bioconductor packages, which are required dependencies for 
`amplicR`, are not currently installed along with the package. This may change 
in the future, but for now these can be installed (if not installed already) 
running the following:

```
library(amplicR)
setup()
```

**NOTE (August 2016)**
At the time of writing, the current official release of devtools (1.12.0) has an issue with R 3.3.1 on windows, which results in the dependencies of the package not being installed. This problem is resolved in the dev version. If you have problems installing `amplicR` try the following (from devtools manual for the function `build_github_devtools`):

```
# Install devtools from CRAN if you haven't already
install.packages("devtools")

library(devtools)
build_github_devtools()

#### Restart R before continuing ####
install.packages("./devtools.zip", repos = NULL)

# Remove the package after installation
unlink("./devtools.zip")

library(devtools)
install_github("carlopacioni/amplicR")

```

## Disclamer 
All reasonable care has been taken to ensure that `amplicR` functions report the 
correct results. However there is no guarantee that the package is bug-free. Also, 
`amplicR` was developed on a machine running windows 7 and no testing on other OS 
has been conducted so far. I can't see any reason why it would not work on linux 
or mac, but if in doubt, you may want to replicate the results on a windows 
environment.

## Documentation 
Use `help(amplicR)` `?amplicR` or `??amplicR` to see a broad 
description of the package. Use `help(package = "amplicR")` to see the 
documentations available.
Alternatively, a manual is available [here](https://www.researchgate.net/publication/307545309_amplicR_-_Manual)
and a tutorial is available [here](https://www.researchgate.net/publication/308204446_amplicR_-_tutorial?ev=prf_pub).

## Citation 
If you use `amplicR`, please cite: _Pending_
