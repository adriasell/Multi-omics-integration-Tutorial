# Multi-omics integration tutorial

Adrià Setó Llorens, Predoctoral Researcher at the Barcelona Institute for Global Health (ISGlobal).

Augusto Anguita-Ruiz, Junior Leader Researcher at the Barcelona Institute for Global Health (ISGlobal).

The multi-omics approach aims to integrate diverse layers of biological information—such as genomics, transcriptomics, proteomics, metabolomics, and epigenomics—to achieve a more comprehensive understanding of biological systems and disease mechanisms. Each omic layer captures a distinct yet interconnected level of cellular regulation, and their integration enables the identification of molecular interactions that cannot be detected through single-omic analyses alone. The main advantage of multi-omics integration over traditional single-omic studies lies in its ability to uncover cross-level biological relationships and multi-factorial drivers of phenotypes, improving prediction accuracy and mechanistic insight. This systems-level perspective supports the discovery of key biomarkers, regulatory networks, and potential therapeutic targets.

There are many multi-omics integration algorithms, each suited for different analytical goals, and they can be classified according to whether they are supervised or unsupervised and whether they perform variable selection—in this session, we will focus on the RGCCA (Regularized Generalized Canonical Correlation Analysis) approach.

For this practical tutorial, we will use data from the HELIX exposome study. The HELIX study is a collaborative project between six longitudinal, population-based birth cohort studies from six European countries (France, Greece, Lithuania, Norway, Spain and the UK).

<img width="1024" height="212" alt="HELIX" src="https://github.com/user-attachments/assets/090ed53b-dda4-4383-9b9a-0966efc3f90d" />

Note: The data provided in this introductory course were simulated from the HELIX sub-cohort data. Details of the HELIX project and the origin of the data collected can be found in the following publication: BMJ Open - HELIX and on the project website. Additional details about the dataset can be found in the official repository at "https://github.com/isglobal-exposomeHub/ExposomeDataChallenge2021".

# Repository guide

# Reminder: Introduction to NoteBook
This tutorial (multiomics_integration_tutorial.ipynb file) is a NoteBook object. Within this notebook (NoteBook), you will be guided step by step from loading a dataset to performing analysis of its content.

To access the tutorial:

* Open the multiomics_integration_tutorial.ipynb file.
* Click the “Open in Colab” button at the top of the notebook.

You’ll need to sign in with a Google account to run the Colab tutorial.

The Jupyter (Python) notebook is an approach that combines text blocks (like this one) together with code blocks or cells. The great advantage of this type of cell is its interactivity, as they can be executed to check the results directly within them. Very important: the order of instructions is fundamental, so each cell in this notebook must be executed sequentially. If any are omitted, the program may throw an error, so you should start from the beginning if in doubt.

First of all:

Once inside google colab, it is very very important that at the start you select "Open in draft mode" (draft mode), at the top left. Otherwise, it will not allow you to execute any code block, for security reasons. When the first of the blocks is executed, the following message will appear: "Warning: This notebook was not created by Google.". Do not worry, you should trust the content of the notebook (NoteBook) and click "Run anyway".

Let’s go!

Click the "play" button on the left side of each code cell. Lines of code that begin with a hashtag (#) are comments and do not affect the execution of the program.

You can also click on each cell and press "ctrl+enter" (cmd+enter on Mac).

Each time you run a block, you will see the output just below it. The information is usually always related to the last instruction, along with all the print() commands in the code.
