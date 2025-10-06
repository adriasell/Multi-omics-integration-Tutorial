# Multi-omics integration tutorial

<img width="1419" height="442" alt="2_ATHLETE_logo_subtitle_color" src="https://github.com/user-attachments/assets/b7a1faeb-1263-41e8-a5bb-060e913c28d2" />

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

By the end of the session, you’ll have a hands-on understanding of how RGCCA works and how to apply it to multi-omics data.

For this practical tutorial, we will use data from the HELIX exposome study. The HELIX study is a collaborative project between six longitudinal, population-based birth cohort studies from six European countries (France, Greece, Lithuania, Norway, Spain and the UK).

<img width="1024" height="212" alt="HELIX" src="https://github.com/user-attachments/assets/090ed53b-dda4-4383-9b9a-0966efc3f90d" />

Note: The data provided in this introductory course were simulated from the HELIX sub-cohort data. Details of the HELIX project and the origin of the data collected can be found in the following publication: BMJ Open - HELIX and on the project website. Additional details about the dataset can be found in the official repository at https://github.com/isglobal-exposomeHub/ExposomeDataChallenge2021.

# Repository guide
The repository contains the following documents:

* The multiomics_integration_tutorial.ipynb. It contains the notebook for the practical tutorial with the code needed to perform the multi-omic integration using RGCCA.

**data**: This folder contains the codebook and the datasets that will be used during the session.

* The **exposoma data (n = 1301)** that we will use are contained in a **Rdata** file, which includes the following files:
    * `phenotype` (results)
    * `exposome` (exposome)
    * `covariates` (covariates)
    * `codebook`

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
