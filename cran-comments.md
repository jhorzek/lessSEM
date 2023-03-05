## R CMD check results

0 errors | 0 warnings | 7 notes

* This is a new release.

## 7 Notes:

### N  checking CRAN incoming feasibility

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

Four dois of cited papers are flagged as 404, Forbidden. I've checked these dois manually and 
they link to the correct papers.

### checking package dependencies (1.1s)

```
Package in Depends/Imports which should probably only be in LinkingTo: 'RcppArmadillo'
```

The package allows users to define custom C++ functions which are compiled
using Rcpp with `#include <RcppArmadillo.h>`. Without adding RcppArmadillo to
Imports, compiling these functions seems to not work.

### checking installed package size ... 

```
     installed size is  5.6Mb
     sub-directories of 1Mb or more:
       libs   3.3Mb
```

### checking for future file timestamps ...

```
   unable to verify current time
```

### checking dependencies in R code (1.7s)

```
Namespace in Imports field not imported from: 'RcppArmadillo'
     All declared Imports should be used.
```

See above: The package allows users to define custom C++ functions which are compiled
using Rcpp with `#include <RcppArmadillo.h>`. Without adding RcppArmadillo to
Imports, compiling these functions seems to not work.

### checking R code for possible problems ... [12s] NOTE (11.9s)

```
.compileTransformations: no visible binding for global variable
     'getPtr'
   .compileTransformations: no visible binding for global variable
     'transformationFunction'
   Undefined global functions or variables:
     getPtr transformationFunction
```

See above: The package allows users to define custom C++ functions which are compiled
using Rcpp with `#include <RcppArmadillo.h>`. These functions are called "getPtr" and
"transformationFunction" and are used after compilation throughout the package.

### N  checking for GNU extensions in Makefiles

```
   GNU make is a SystemRequirements.
```

GNU make is used because the package allows for multi-core execution using 
RcppParallel.