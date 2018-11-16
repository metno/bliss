# bliss

Multiple statistical interpolation schemes have been implemented, such as:

* OI\_multiscale. Multiscale Optimal Interpolation based on observations only.
* OI\_firstguess. Optimal Interpolation based on the combination of observations with a background field. It is possible to specify a data transformation as a pre-processing step (i.e., Gaussian anamorphosis).
* OI\_twosteptemperature. Modified Optimal Interpolation based on observations only.

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
