
# Inference and diagnosis in M-quantile models with applications to small area estimation

### M. Bugallo and  D. Morales 

### Center of Operations Research, Miguel Hernández University of Elche, Spain
## Objective

This repository contains the R scripts corresponding to Section 4, "Empirical Results", of the research article "Inference and Diagnosis in M-Quantile Models with Applications to Small Area Estimation". 

Section 4 presents the results of several empirical studies, focusing on:

-- The distribution of the optimal robustness parameters in M-quantile models.

-- Their use as a diagnostic tool for detecting area-level outliers.

The scripts reproduce the analyses and figures reported in this section of the paper.


## Installation and requirements

The R code available in this repository was tested on a 3.3 GHz Intel Xeon W Mac computer with 32 GB RAM. For data manipulation and task processing, the following R packages are required:

-   library(MASS)
-   library(dplyr)


To facilitate code manipulation, we recommend using the **Graphical User Interface (GUI)** [RStudio](https://posit.co/downloads/). Running the scripts and reproducing the results is straightforward. Once all the required libraries are installed, simply download the project sources and execute the main scripts ***Distribution_Optimal_RPs.R*** and ***SimuSizeTest_Optimal_RPS.R***, corresponding to Section 4.1, "Distribution of the Optimal Robustness Parameters", and Section 4.2, "Detection of Atypical Areas", respectively.

## Contents

1. **Main script 1 – *Distribution_Optimal_RPs.R***: This script reproduces the analyses of Section 4.1, "Distribution of the Optimal Robustness Parameters", in the paper. It computes and visualizes the distribution of the optimal robustness parameters for the M-quantile models, generating the figures and summary statistics reported in the manuscript.  

2. **Main script 2 – *SimuSizeTest_Optimal_RPS.R***: This script runs the Monte Carlo simulations described in Section 4.2, "Detection of Atypical Areas". It explores how the optimal robustness parameters can highlight atypical areas and units, producing the figures and results reported in the manuscript.

3. **Auxiliary scripts – *MQ2.R* and *AuxFunctions.R***: Contain essential functions for model fitting, prediction, mean squared error estimacion and bootstrap inference.  
