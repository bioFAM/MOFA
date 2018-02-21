# MOFA: Multi-Omics Factor Analysis

MOFA is a factor analysis model that provides a **general framework for the integration of multi-omic data sets** in a completely unsupervised fashion.  
Intuitively, MOFA can be viewed as a versatile and statistically rigorous generalization of principal component analysis (PCA) to multi-omics data. Given several data matrices with measurements of multiple ‘omics data types on the same or on overlapping sets of samples, MOFA infers an **interpretable low-dimensional data representation in terms of (hidden) factors**. These learnt factors represent the driving sources of variation across data modalities, thus facilitating the identification of cellular states or disease subgroups.  

Once trained, the model output can be used for a range of downstream analyses, including the visualisation of samples in factor space, the automatic annotation of factors using (gene set) enrichment analysis, the identification of outliers (e.g. due to sample swaps) and the imputation of missing values.  

For more details you can read our preprint: https://www.biorxiv.org/content/early/2017/11/10/217554
<p align="center"> 
<img src="images/logo.png" style="width: 50%; height: 50%"/>​
</p>



## News
- 10/11/2017 Paper uploaded to bioRxiv
- 10/11/2017 We created a Slack group to provide personalised help on running and interpreting MOFA, [this is the link](https://join.slack.com/t/mofahelp/shared_invite/enQtMjcxNzM3OTE3NjcxLTkyZmE5YzNiMDc4OTkxYWExYWNlZTRhMWI2OWNkNzhmYmNlZjJiMjA4MjNiYjI2YTc4NjExNzU2ZTZiYzQyNjY)
 

## Installation
MOFA is run exclusively from R, but it requires some python dependencies that you need to install. Here is how to do it:

### Python dependencies 
We are on the process of uploading it to PyPI. For now, you can install it using:
```r
pip install git+git://github.com/PMBio/MOFA
```
Or clone the repository and then install it using the setup.py:
```r
git clone https://github.com/PMBio/MOFA
python setup.py install
```

### R package
The easier way to install the R package is via github:
```r
devtools::install_github("PMBio/MOFA", subdir="MOFAtools")
```

Alternatively, you can clone the repository and install locally:
```r
git clone https://github.com/PMBio/MOFA
R CMD build MOFAtools
R CMD install MOFAtools
```


## MOFA workflow

The workflow of MOFA consists of two steps:  
**(1) Fitting step**: train the model with the multi-omics data to disentangle the heterogeneity into a small number of latent factors.  
**(2) Downstream analysis**: once the factors are inferred they need to be characterised as technical or biological sources of variation by looking at the corresponding weights, doing (gene set) enrichment analysis, plotting the factors, correlating factors with known covariates, etc. Also, one can do imputation of missing values and prediction of clinical outcomes using the latent factors.

<p align="center"> 
<img src="images/workflow.png">
</p>

A list of all **relevant methods** with a short description can be found [here](https://github.com/PMBio/MOFA/blob/master/MOFAtools/Documentation.md)  

### Step 1: Fitting the model
First you need to create the MOFA object with your input data, and subsequently you need to train the model. Everything is explained in [the vignette](http://htmlpreview.github.com/?https://github.com/PMBio/MOFA/blob/master/MOFAtools/vignettes/MOFA_example_CLL.html). 
If everything is successful, you should observe an output analogous to the following:
```
  ###########################################################
  ###                 __  __  ___  _____ _                ###
  ###                |  \/  |/ _ \|  ___/ \               ###
  ###                | |\/| | | | | |_ / _ \              ###
  ###                | |  | | |_| |  _/ ___ \             ###
  ###                |_|  |_|\___/|_|/_/   \_\            ###
  ###                                                     ###
  ###########################################################

##################
## Loading data ##
##################

Loaded /Users/ricard/MOFA/MOFA/test/data/500_0.txt with dim (100,500)...
Loaded /Users/ricard/MOFA/MOFA/test/data/500_1.txt with dim (100,500)...
Loaded /Users/ricard/MOFA/MOFA/test/data/500_2.txt with dim (100,500)...
 

#############################################
## Running trial number 1 with seed 642034 ##
#############################################

Trial 1, Iteration 1: time=0.08 ELBO=-345954.96, Factors=10, Covariates=1
Trial 1, Iteration 2: time=0.10 ELBO=-283729.31, deltaELBO=62225.6421, Factors=10, Covariates=1
Trial 1, Iteration 3: time=0.10 ELBO=-257427.42, deltaELBO=26301.8893, Factors=10, Covariates=1
...
Trial 1, Iteration 100: time=0.07 ELBO=-221171.01, deltaELBO=0.0998, Factors=10, Covariates=1

Converged!
```

There are two important quantities to keep track of: 
* **Number of factors**: you start the model with a large enough amount of factors, and the model will automatically remove the factors that do not explain significant amounts of variation. 
* **deltaELBO**: this is the objective function being maximised which is used to assess model convergence. Once the deltaELBO decreases below a threshold, training will end and the model will be saved as an .hdf5 file. Then, you are ready to start the analysis with the R package.

### Step 2: Disentangle the variability
MOFA disentangles the heterogeneity of a high-dimensional multi-omics data set into a reduced set of latent factors that capture global sources of variation. 
Importantly, these factors can have different activity patterns in different omics. For example, a batch effect might be affecting the RNA data but not the Methylation data. 
Decoupling this heterogeneity is a mandatory first step in the analysis of multi-omics data. For example, this is the variance decomposition plot for the Chronic Lymphocytic Leukemia data set analysed in the paper:

<p align="center"> 
<img src="images/varExplained.png" style="width: 50%; height: 50%"/>​
</p>


### Step 3: Annotation of factors
Once the heterogeneity of the data set is reduced into a set of factors, you need to understand what are they, and whether they capture technical or biological sources of variability. 

We have built a semi-automated pipeline based on our experience annotating factors:  
(1) **Visualisation of the samples in the factor space**: similarly to what is done in Principal Component Analysis, it is useful to plot the factors against each other and color the samples using known covariates such as batch, sex, clinical information, etc.  
(2) **Inspection of top weighted features**: for example, if a factor is associated to the sex of the individual, the mRNA data will have very high loadings for genes located in the X and Y chromosomes.  
(3) **Feature set enrichment analysis**: particularly when having large amounts of features, the inspection of loadings is challenging, and doing gene ontology enrichment analysis can be useful.  

Please refer to the paper for details on the different analysis.  

### Step 4: Using the factors to get biological insights in downstream analysis
The latent factors can be used for several purposes, such as:  
(1) **Dimensionality reduction**: similar to PCA, dimensionality reduction plots can be obtained by plotting the Factors against each other.  
(2) **Imputation**: Factors can be used to predict missing values, including entire missing assays.  
(3) **Predicting clinical response**: if the factors capture phenotypical information, they can capture clinical covariates of interest.  
(4) **Regressing out technical effects**: if a factor is capturing an undesired technical effect, its effect can be regressed out from your original data matrix.  

Please refer to the paper for details on the different analysis. 

## Tutorial
An example workflow is provided in [the vignette](http://htmlpreview.github.com/?https://github.com/PMBio/MOFA/blob/master/MOFAtools/vignettes/MOFA_example_CLL.html).

## Frequently asked questions

**(1) I get the following error when running MOFA:**
```
sh: mofa: command not found
```
This occurs if the mofa binary is not in the $PATH of R. This will be fixed in a new update, a simple workaround is to get the full path of the mofa executable by doing on the terminal:
```
which mofa
```
and then copy the path into runMOFA:
```
runMOFA(object, DirOptions, ..., mofaPath="PUT THE PATH HERE")
```

**(2) I get the following error when installing the R package:**
```
ERROR: dependencies 'pcaMethods', 'MultiAssayExperiment' are not available for package 'MOFAtools'
```
These two packages are available from Bioconductor, not CRAN. You can install them from R as follows:
```
source("https://bioconductor.org/biocLite.R")
biocLite(c('pcaMethods', 'MultiAssayExperiment'))
```

## Contact
The package is maintained by Britta Velten (britta.velten@embl.de) and Ricard Argelaguet (ricard@ebi.ac.uk). 
Please, contact us for problems, comments or suggestions.


