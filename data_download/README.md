## BUGs:

If you encounter the error message from `GDCprepare` function when download the TCGA data, such as

```R
GDCprepare ERROR: Error: $ operator is invalid for atomic vectors
```

Please, you could update the `TCGAbiolinks` package with the github version with the following code 
```R
devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
```
