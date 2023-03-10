## R CMD check results

0 errors | 0 warnings | 5 notes

* This is a new release.

# Comments

Some functions in the package allow users to compile user-specified functions 
with Rcpp. These functions are compiled using `#include <RcppArmadillo.h>` which
I could only get to work with "Imports: RcppArmadillo". The following note results
from this setup:

```
N checking package dependencies (1.1s)
  Package in Depends/Imports which should probably only be in LinkingTo: 'RcppArmadillo'
```

## Additional Notes:

Depending on the configuration used to check the package, there is a note regarding
links used in the documentation:
```
Found the following (possibly) invalid URLs:
    URL: https://doi.org/10.1111/bmsp.12130
      From: inst/doc/Definition-Variables-and-Multi-Group-SEM.html
            inst/doc/Parameter-transformations.html
            inst/doc/lessSEM.html
            README.md
      Status: 403
      Message: Forbidden
    URL: https://doi.org/10.1111/j.1467-9868.2005.00503.x
      From: inst/doc/lessSEM.html
            README.md
      Status: 403
      Message: Forbidden
    URL: https://doi.org/10.1111/j.2044-8317.1984.tb00802.x
      From: inst/doc/The-Structural-Equation-Model.html
      Status: 403
      Message: Forbidden
    URL: https://doi.org/10.1137/080716542
      From: inst/doc/SCAD-and-MCP.html
            inst/doc/The-optimizer-interface.html
            inst/doc/lessSEM.html
            README.md
      Status: 403
      Message: Forbidden
```
I've checked all of these dois and they work on my system.

```
N checking CRAN incoming feasibility
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    SEM (14:16)
    lavaan (15:6)
    lessSEM (17:59, 17:68)
    lslx (17:4)
    regsem (16:4)
```

The words are all correct; they refer to R packages or common abbreviations of
the models used in the package.

```
N checking installed package size ... NOTE
  installed size is  5.6Mb
  sub-directories of 1Mb or more:
    libs   3.4Mb
``` 

```
N checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.
```
The package uses RcppParallel which requires GNU make.
