# Multi-omics integration tutorial

**Adrià Setó Llorens**, Predoctoral Researcher at the Barcelona Institute for Global Health (ISGlobal).

**Augusto Anguita-Ruiz**, Junior Leader Researcher at the Barcelona Institute for Global Health (ISGlobal).

The multi-omics approach aims to integrate diverse layers of biological information—such as genomics, transcriptomics, proteomics, metabolomics, and epigenomics—to achieve a more comprehensive understanding of biological systems and disease mechanisms. Each omic layer captures a distinct yet interconnected level of cellular regulation, and their integration enables the identification of molecular interactions that cannot be detected through single-omic analyses alone. The main advantage of multi-omics integration over traditional single-omic studies lies in its ability to uncover cross-level biological relationships and multi-factorial drivers of phenotypes, improving prediction accuracy and mechanistic insight. This systems-level perspective supports the discovery of key biomarkers, regulatory networks, and potential therapeutic targets.

There are many multi-omics integration algorithms, each suited for different analytical goals, and they can be classified according to whether they are supervised or unsupervised and whether they perform variable selection—in this session, we will focus on the **RGCCA** (Regularized Generalized Canonical Correlation Analysis) approach.

The objective of this session to offer an introduction to a multi-omics integration analysis using **RGCCA**. We will:
* Load the data
* Preprocess the data
* Perform multi-omics integration
* Understand the results of multi-omics integration
* Evaluate the algorithm’s performance

We will integrate multi-omics data — including proteomics, urine and serum metabolomics, gene expression, and DNA methylation — using Regularized Generalized Canonical Correlation Analysis (RGCCA). The outcome variable will be standardized body mass index (zBMI) at 9 years old. The objective of this analysis is to identify multi-omic signatures predictive of BMI in later childhood while gaining a hands-on understanding of the application of RGCCA to multi-omics data integration.

For this practical tutorial, we will use data from the HELIX exposome study. The HELIX study is a collaborative project between six longitudinal, population-based birth cohort studies from six European countries (France, Greece, Lithuania, Norway, Spain and the UK).

<img width="800" height="150" alt="HELIX" src="https://github.com/user-attachments/assets/090ed53b-dda4-4383-9b9a-0966efc3f90d" />

Note: The data provided in this introductory course were simulated from the HELIX sub-cohort data. Details of the HELIX project and the origin of the data collected can be found in the following publication: BMJ Open - HELIX and on the project website. Additional details about the dataset can be found in the official repository at https://github.com/isglobal-exposomeHub/ExposomeDataChallenge2021.

# Repository guide
The repository contains the following documents:

* The multiomics_integration_tutorial.ipynb. It contains the notebook for the practical tutorial with the code needed to perform the multi-omic integration using RGCCA.
* Functions: This directory contains all the functions used in this session. These functions are stored in separate files to keep the notebook clean and easy to follow. For more details, you can consult the files in this directory.

This is the dataset we will use:

- **Exposome data (n=1301)**: [Rdata file](https://github.com/isglobal-brge/brge_data_large/blob/master/data/ExposomeDataChallenge2021/exposome_NA.RData) containing three objects:
     - 1 object for exposures: `exposome`
     - 1 object for covariates: `covariates`
     - 1 object for outcomes: `phenotype`

The three tables can be linked using **ID** variable. See the [codebook](https://github.com/isglobal-brge/brge_data_large/blob/master/data/ExposomeDataChallenge2021/codebook.xlsx) for variable description (variable name, domain, type of variable, transformation, ...)


- **omic data**: Exposome and omic data can be linked using **ID** variable. 
     - [Proteome](https://github.com/isglobal-brge/brge_data_large/blob/master/data/ExposomeDataChallenge2021/proteome.Rdata): ExpressionSet called `metabol_serum` of **1170 individuals** and **39 proteins** (log-transformed) that are annotated in the `ExpressionSet` object (use `fData(proteome)` after loading `Biobase` Bioconductor package).
     - [Serum Metabolome](https://github.com/isglobal-brge/brge_data_large/blob/master/data/ExposomeDataChallenge2021/metabol_serum.Rdata): ExpressionSet called `metabol_serum` of **1198 individuals** and **177 metabolites** (log-transformed) (see [here](https://github.com/isglobal-brge/brge_data_large/blob/master/data/ExposomeDataChallenge2021/HELIX_serum_metabol_report_IC_v4_APS_2017_04_06.pdf) for a descripton).
     - [Urine Metabolome](https://github.com/isglobal-brge/brge_data_large/blob/master/data/ExposomeDataChallenge2021/metabol_urine.Rdata): ExpressionSet called `metabol_urine` of **1192 individuals** and **44 metabolites** (see [here](https://github.com/isglobal-brge/brgedata/blob/master/data/ExposomeDataChallenge2021/HELIX_urine_metabol_report_IC_v3_CHL_2017_01_26.pdf) for a descripton). 
     - [Gene expression](https://figshare.com/s/571c8cff7acf5167f343): ExpressionSet called `genexpr`  (see [here](https://isglobal-brge.github.io/Master_Bioinformatics/bioconductor.html#expressionset) what an ExpressionSet is) of **1007 individuals** and **28,738 transcripts** with annotated gene symbols. 
     - [Methylation](https://figshare.com/s/46e6a1d66ff135bb15c8): GenomicRatioSet called `methy` (see [here](https://www.rdocumentation.org/packages/minfi/versions/1.18.4/topics/GenomicRatioSet-class) what a GenomicRatioSet is) of **918 individuals** and **386,518 CpGs**

The variables that are available in the metadata are:

> 1. **ID**: identification number
> 2. **e3_sex**: gender (male, female)
> 3. **age_sample_years**: age (in years)
> 4. **h_ethnicity_cauc**: caucasic? (yes, no)
> 5. **ethn_PC1**: first PCA to address population stratification
> 6. **ethn_PC2**: second PCA to address population stratification
> 7. **Cell-type estimates** (only for methylation): NK_6, Bcell_6, CD4T_6, CD8T_6, Gran_6, Mono_6

# Reminder: Introduction to NoteBook
This notebook will guide you step by step, from loading a dataset to analyzing it.

Getting Started:
* Open multiomics_integration_tutorial.ipynb and click “Open in Colab” (sign in with your Google account if needed).
* Select “Open in draft mode” at the top left so you can run the code safely.
* If you see "Warning: This notebook was not created by Google.", don’t worry—just click Run anyway.

How to Use the Notebook:
* The notebook mixes text explanations and code cells for hands-on learning.
* Always run cells in order to avoid errors.
* Click the play button next to a cell, or press Ctrl+Enter (Cmd+Enter on Mac).
* Lines starting with # are comments for guidance, they won’t affect the code.
* Outputs appear below each cell, showing results and any printed messages.
