# bliss 

BLISS is an R-program for spatial interpolation of near-surface meteorological variables through statistical methods.
It includes several statistical interpolation schemes:

* OI\_multiscale. Multiscale Optimal Interpolation based on observations only.
* OI\_firstguess. Optimal Interpolation based on the combination of observations with a background field. It is possible to specify a data transformation as a pre-processing step (i.e., Gaussian anamorphosis).
* OI\_twosteptemperature. Modified Optimal Interpolation based on observations only.
* OI\_Bratseth
* SC\_Barnes
* Ensemble-based Statitical Interpolation (EnSI)
* EnSI with Gaussian Anamorphosis for Precipitation (EnSIGAP)

and also:

* Rasterize.

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

bliss.r expect the envirnment variable BLISS\_PATH to be defined

```
export BLISS_PATH="/home/cristianl/projects/bliss"
```

Copyright and license
---------------------
Copyright (C) 2018 MET Norway. TITAN is licensed under [GPL
version 3](https://github.com/metno/TITAN/blob/master/LICENSE) or (at
your option) any later version.

Contact
-------
E-mail: cristianl@met.no
