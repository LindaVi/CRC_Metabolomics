# CRC_Metabolomics
This repository contains the analysis code for the article titled "Untargeted plasma metabolomics and risk of colorectal cancer-an analysis nested within a large-scale prospective cohort".

Files

Here are a short description of the files in this repository:

1.	**Logistic regr_CRC.R**: The R script that performs logistic regression analysis for calculation of odds ratio.

2.	**Classification_MUVR_analysis.R**: The R script that performs random forest analysis for prediction of colorectal cancer cases using the MUVR package.

3.	**targeted.rda**: R data containing information about colorectal cancer metabolite biomarker candidates from literature.

4.	**CRC_metabolite_features_OR.R**: The main R script for replication of potential colorectal cancer biomarkers from literature, including the steps:

    a.	Find candidates in peaklist  
    b.	Investigate association to colorectal cancer

5.	**CRC_metabolite_features_OR_sink.txt**: The output after running CRC_metabolite_features_OR.R

6.	**getHits.R**: R function used for finding candidate molecular weights in actual peak table data.

7.	**clrSubFeat_noenergy.R**:  R script for pre-processing and conditional logistic regression calculations on previous biomarkers identified in our data set.   
