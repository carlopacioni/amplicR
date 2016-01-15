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

## Disclamer 
All reasonable care has been taken to ensure that `amplicR` functions report the 
correct results. However there is no garantee that the package is bug-free. Also, 
`amplicR` was developed on a machine running windows 7 and no testing on other OS 
has been conducted so far. I can't see any reason why it would not work on linux 
or mac, but if in doubt, you may want to replicate the results on a windows 
environment.

## Documentation 
Use `help(amplicR)` `?amplicR` or `??amplicR` to see a broad 
description of the package. Use `help(package = "amplicR")` to see the 
documentations available.

## Citation 
If you use `amplicR`, please cite: _Pending_
