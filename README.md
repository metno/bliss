# bliss [![DOI](https://zenodo.org/badge/154810577.svg)](https://zenodo.org/badge/latestdoi/154810577)


BLISS is an R-program for spatial interpolation of near-surface meteorological variables through statistical methods.
It includes several statistical interpolation schemes:

* OI\_multiscale. Multiscale Optimal Interpolation based on observations only.
* OI\_firstguess. Optimal Interpolation based on the combination of observations with a background field. It is possible to specify a data transformation as a pre-processing step (i.e., Gaussian anamorphosis).
* OI\_twosteptemperature. Modified Optimal Interpolation based on observations only.
* letkf. Local Ensemble Transform Kalman Filter.
* hyletkf. Hybrid Local Ensemble Transform Kalman Filter.

Installation instructions
-------------------------
The following R-libraries are needed:

* argparser
* sp
* raster
* igraph
* rgdal
* ncdf4
* dotnc (not available for public download, write a mail to Cristian Lussana)

Install the R-package using either [devtools](https://cran.r-project.org/web/packages/devtools/README.html):

```
devtools::install_github("metno/bliss_rr")
```

or as described [here](https://cran.r-project.org/).

Copyright and license
---------------------
Copyright (C) 2018 MET Norway. TITAN is licensed under [GPL
version 3](https://github.com/metno/TITAN/blob/master/LICENSE) or (at
your option) any later version.

Contact
-------
E-mail: cristianl@met.no
