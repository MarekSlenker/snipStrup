# R packages

Required [R](https://www.r-project.org/) packages `ape`, `memuse`, `pinfsc50` and `vcfR` and their dependencies for R version used **must** be installed here.

Installation of source packages in Linux requires basic compilation tools and developmental files for GDAL, GEOS and UDUNITS-2.

Within `R` command line started in `snipStrup` directory use e.g. command

```R
install.packages(pkgs=c("ape", "memuse", "pinfsc50", "vcfR"),
lib="rpackages", repos="https://mirrors.nic.cz/R/", dependencies="Imports")
```

to install needed packages.
