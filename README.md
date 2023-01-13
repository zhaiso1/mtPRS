# Integrating multiple traits for improving polygenic risk prediction in disease and pharmacogenomics (PGx) GWAS

## Overview
mtPRS package (Zhai et al., 2023) implements two novel multi-trait polygenic risk score (mtPRS) methods, mtPRS-PCA and mtPRS-O. Specifically, mtPRS-PCA combines individual single-trait PRSs (stPRSs) with weights calculated by performing principal component analysis (PCA) on the genetic correlation matrix among traits. mtPRS-O aggregates p-values from mtPRS-PCA, mtPRS-ML (Machine Learning; Krapohl et al., 2018), and all stPRSs using Cauchy Combination Test (CCT) (Liu et al., 2019) to provide a robust way for the multi-trait PRS association test.

## Installation

```
library(devtools)
devtools::install_github("zhaiso1/mtPRS")
```

## Usage

Details can be found in **doc** folder for mtPRS user manual (**mtPRS_0.1.0.pdf**), and **vignettes** folder for a demo illustrating how to use our softwares (**README.pdf**).

The mtPRS package includes four functions:

- **generate_dis_data**: simulate disease GWAS summary statistics of K traits in base cohort, individual-level disease genetics data in target cohort, underlying genetic correlation matrix among traits, and underlying true effect sizes
- **generate_pgx_data**: simulate disease GWAS summary statistics of K traits in base cohort, individual-level PGx data in target cohort, underlying genetic correlation matrix among traits, and underlying true effect sizes
- **mtPRS_PCA**: Run mtPRS-PCA algorithm, which will output weights calculated by performing principal component analysis on the genetic correlation matrix, standardized stPRSs, and mtPRS-PCA
- **mtPRS_O**: Run mtPRS-O algorithm. If the phenotype is from disease GWAS, then mtPRS-O will output main prognostic p-value. If the phenotype is from PGx GWAS (with both T and C arms), then mtPRS-O will output a vector of main prognostic p-value from 2df test, interaction predictive p-value from 2df test, main prognostic p-value in T arm only, and main prognostic p-value in C arm only.

## References

Krapohl E, Patel H, Newhouse S, et al. Multi-polygenic score approach to trait prediction. Mol. psychiatry 2018; 23: 1368-1374

Liu Y, Chen S, Li Z, et al. ACAT: a fast and powerful p-value combination method for rare-variant analysis in sequencing studies. Am. J. Hum. Genet. 2019; 104: 410-421.

Zhai S, Guo B, Wu B, Mehrotra DV, and Shen J. Integrating multiple traits for improving polygenic risk prediction in disease and pharmacogenomics GWAS. To be submitted to Briefings in Bioinformatics in 2023.
