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

* Third, install GermVersity

```
devtools::install_github('GermVersity/GermVersity')
```
