# GermVersity
<p align="center">
<img src = "https://raw.githubusercontent.com/GermVersity/GermVersity/main/inst/app/www/Logo.png" alt = "drawing" align = "center" width = "500" height = "500"/>
</p>

## Installation from GitHub

* First install the next libraries from Github

```
devtools::install_github("silkeszy/Pomona")
```

```
## Require R >= 4.3.x
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
```

* Second install from BiocManager

```
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")
```
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("qvalue")
```
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("snpStats")
```

* Then install from github again

```
devtools::install_github("SFUStatgen/LDheatmap")
```

* Install package maptools from https://cran.r-project.org/src/contrib/Archive/maptools/
  * Download maptools 1.1-8.tar.gz
  * Extract the files in a path
  * Then follow the next instructions

```
## in console R execute the next code

install.packages('/path/to/maptools', repos=NULL, type="source")
```


* Third, install GermVersity

```
devtools::install_github('GermVersity/GermVersity')
```

* Fourth, run GermVersity

```
GermVersity::run_GermVersity()
```
