## BUGs:

If you encounter the error message in GDCprepare function when download the TCGA data

```R
GDCprepare ERROR: Error: $ operator is invalid for atomic vectors
```

Please, could you update with the github version with the following code 
```R
devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
```
