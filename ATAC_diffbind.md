---
title: "ATAC-seq differential open chromatin analysis"
author: "Brook Wassie"
date: "August 13, 2017"
output: html_document
---


# Statistical Analysis of Differentially Open Chromatin Sites

## Introduction

ATAC-seq (Assay for Transposase-Accessible Chromatin using Sequencing) is a genomic assay that reveals open chromatin sites in the genome using the Tn5 transposase. As with any Omics assay, ATAC-seq can be used to elucidate genomic differences between two conditions. In this vignette, we will use the DiffBind R package to statistically model and test differences in open chromatin sites between iPSC derived motor neurons from ALS and healthy patients. 

For more information on the DiffBind package

User guide: https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf

Reference manual: https://bioconductor.org/packages/release/bioc/manuals/DiffBind/man/DiffBind.pdf

## 1. Load peaks and counts









```r
samples_count = dba.count(samples,score=DBA_SCORE_READS)
samples_count
```

```
## 15 Samples, 104040 sites in matrix:
##              ID    Tissue Factor        Condition     Treatment Caller
## 1   25iCTR_un_1 d32_diMNS    CTR 25iCTR_Untreated CTR_Untreated counts
## 2  83iCTR_ASO_1 d32_diMNS    CTR       83iCTR_ASO   CTR_treated counts
## 3    30ALS_un_1 d32_diMNS    ALS  30ALS_Untreated ALS_Untreated counts
## 4    52ALS_un_2 d32_diMNS    ALS  52ALS_Untreated ALS_Untreated counts
## 5  25iCTR_ASO_1 d32_diMNS    CTR       25iCTR_ASO   CTR_treated counts
## 6    52ALS_un_1 d32_diMNS    ALS  52ALS_Untreated ALS_Untreated counts
## 7   30ALS_ASO_1 d32_diMNS    ALS        30ALS_ASO   ALS_treated counts
## 8   83iCTR_un_2 d32_diMNS    CTR 83iCTR_Untreated CTR_Untreated counts
## 9   52ALS_ASO_2 d32_diMNS    ALS        52ALS_ASO   ALS_treated counts
## 10  83iCTR_un_1 d32_diMNS    CTR 83iCTR_Untreated CTR_Untreated counts
## 11  52ALS_ASO_1 d32_diMNS    ALS        52ALS_ASO   ALS_treated counts
## 12  25iCTR_un_2 d32_diMNS    CTR 25iCTR_Untreated CTR_Untreated counts
## 13 83iCTR_ASO_2 d32_diMNS    CTR       83iCTR_ASO   CTR_treated counts
## 14   30ALS_un_2 d32_diMNS    ALS  30ALS_Untreated ALS_Untreated counts
## 15  30ALS_ASO_2 d32_diMNS    ALS        30ALS_ASO   ALS_treated counts
##    Intervals FRiP
## 1     104040 0.29
## 2     104040 0.20
## 3     104040 0.30
## 4     104040 0.30
## 5     104040 0.21
## 6     104040 0.28
## 7     104040 0.36
## 8     104040 0.23
## 9     104040 0.22
## 10    104040 0.24
## 11    104040 0.29
## 12    104040 0.22
## 13    104040 0.16
## 14    104040 0.26
## 15    104040 0.29
```

## 2. Plotting raw count data


```r
dba.plotHeatmap(samples_count)
```

![plot of chunk plot heatmap and PCA](figure/plot heatmap and PCA-1.svg)

```r
dba.plotPCA(sample_count,DBA_FACTOR,label=DBA_ID)
```

```
## Error in pv.check(DBA, bCheckEmpty = TRUE): object 'sample_count' not found
```


## 3. Setting up contrasts and performing DEseq2 analysis


```r
samples_contrast=dba.contrast(samples_count, categories=DBA_FACTOR, minMembers = 2)
samples_contrast
```

```
## 15 Samples, 104040 sites in matrix:
##              ID    Tissue Factor        Condition     Treatment Caller
## 1   25iCTR_un_1 d32_diMNS    CTR 25iCTR_Untreated CTR_Untreated counts
## 2  83iCTR_ASO_1 d32_diMNS    CTR       83iCTR_ASO   CTR_treated counts
## 3    30ALS_un_1 d32_diMNS    ALS  30ALS_Untreated ALS_Untreated counts
## 4    52ALS_un_2 d32_diMNS    ALS  52ALS_Untreated ALS_Untreated counts
## 5  25iCTR_ASO_1 d32_diMNS    CTR       25iCTR_ASO   CTR_treated counts
## 6    52ALS_un_1 d32_diMNS    ALS  52ALS_Untreated ALS_Untreated counts
## 7   30ALS_ASO_1 d32_diMNS    ALS        30ALS_ASO   ALS_treated counts
## 8   83iCTR_un_2 d32_diMNS    CTR 83iCTR_Untreated CTR_Untreated counts
## 9   52ALS_ASO_2 d32_diMNS    ALS        52ALS_ASO   ALS_treated counts
## 10  83iCTR_un_1 d32_diMNS    CTR 83iCTR_Untreated CTR_Untreated counts
## 11  52ALS_ASO_1 d32_diMNS    ALS        52ALS_ASO   ALS_treated counts
## 12  25iCTR_un_2 d32_diMNS    CTR 25iCTR_Untreated CTR_Untreated counts
## 13 83iCTR_ASO_2 d32_diMNS    CTR       83iCTR_ASO   CTR_treated counts
## 14   30ALS_un_2 d32_diMNS    ALS  30ALS_Untreated ALS_Untreated counts
## 15  30ALS_ASO_2 d32_diMNS    ALS        30ALS_ASO   ALS_treated counts
##    Intervals FRiP
## 1     104040 0.29
## 2     104040 0.20
## 3     104040 0.30
## 4     104040 0.30
## 5     104040 0.21
## 6     104040 0.28
## 7     104040 0.36
## 8     104040 0.23
## 9     104040 0.22
## 10    104040 0.24
## 11    104040 0.29
## 12    104040 0.22
## 13    104040 0.16
## 14    104040 0.26
## 15    104040 0.29
## 
## 1 Contrast:
##   Group1 Members1 Group2 Members2
## 1    CTR        7    ALS        8
```


```r
samples_differential = dba.analyze(samples_contrast, method = DBA_DESEQ2)
```

```
## converting counts to integer mode
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```r
samples_differential
```

```
## 15 Samples, 104040 sites in matrix:
##              ID    Tissue Factor        Condition     Treatment Caller
## 1   25iCTR_un_1 d32_diMNS    CTR 25iCTR_Untreated CTR_Untreated counts
## 2  83iCTR_ASO_1 d32_diMNS    CTR       83iCTR_ASO   CTR_treated counts
## 3    30ALS_un_1 d32_diMNS    ALS  30ALS_Untreated ALS_Untreated counts
## 4    52ALS_un_2 d32_diMNS    ALS  52ALS_Untreated ALS_Untreated counts
## 5  25iCTR_ASO_1 d32_diMNS    CTR       25iCTR_ASO   CTR_treated counts
## 6    52ALS_un_1 d32_diMNS    ALS  52ALS_Untreated ALS_Untreated counts
## 7   30ALS_ASO_1 d32_diMNS    ALS        30ALS_ASO   ALS_treated counts
## 8   83iCTR_un_2 d32_diMNS    CTR 83iCTR_Untreated CTR_Untreated counts
## 9   52ALS_ASO_2 d32_diMNS    ALS        52ALS_ASO   ALS_treated counts
## 10  83iCTR_un_1 d32_diMNS    CTR 83iCTR_Untreated CTR_Untreated counts
## 11  52ALS_ASO_1 d32_diMNS    ALS        52ALS_ASO   ALS_treated counts
## 12  25iCTR_un_2 d32_diMNS    CTR 25iCTR_Untreated CTR_Untreated counts
## 13 83iCTR_ASO_2 d32_diMNS    CTR       83iCTR_ASO   CTR_treated counts
## 14   30ALS_un_2 d32_diMNS    ALS  30ALS_Untreated ALS_Untreated counts
## 15  30ALS_ASO_2 d32_diMNS    ALS        30ALS_ASO   ALS_treated counts
##    Intervals FRiP
## 1     104040 0.29
## 2     104040 0.20
## 3     104040 0.30
## 4     104040 0.30
## 5     104040 0.21
## 6     104040 0.28
## 7     104040 0.36
## 8     104040 0.23
## 9     104040 0.22
## 10    104040 0.24
## 11    104040 0.29
## 12    104040 0.22
## 13    104040 0.16
## 14    104040 0.26
## 15    104040 0.29
## 
## 1 Contrast:
##   Group1 Members1 Group2 Members2 DB.DESeq2
## 1    CTR        7    ALS        8     26071
```


```r
dba.plotMA(samples_differential)
```

![plot of chunk perform MA plot](figure/perform MA plot-1.svg)
