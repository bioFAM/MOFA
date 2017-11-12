# MOFA: Multi-Omics Factor Analysis

MOFA is a factor analysis model that provides a **general framework for the integration of multi-omic data sets** in a completely unsupervised fashion.  
Intuitively, MOFA can be viewed as a versatile and statistically rigorous generalization of principal component analysis (PCA) to multi-omics data. Given several data matrices with measurements of multiple ‘omics data types on the same or on overlapping sets of samples, MOFA infers an **interpretable low-dimensional data representation in terms of (hidden) factors**. These learnt factors represent the driving sources of variation across data modalities, thus facilitating the identification of cellular states or disease subgroups.  

Once trained, the model output can be used for a range of downstream analyses, including the visualisation of samples in factor space, the automatic annotation of factors using (gene set) enrichment analysis, the identification of outliers (e.g. due to sample swaps) and the imputation of missing values.  

For more details you can read our preprint: https://www.biorxiv.org/content/early/2017/11/10/217554
<p align="center"> 
<img src="logo.png">
</p>

## News
- 10/11/2017 Paper uploaded to bioRxiv
- 10/11/2017 We created a Slack group to provide personalised help on running and interpreting MOFA, [this is the link](https://join.slack.com/t/mofahelp/shared_invite/enQtMjcxNzM3OTE3NjcxLTkyZmE5YzNiMDc4OTkxYWExYWNlZTRhMWI2OWNkNzhmYmNlZjJiMjA4MjNiYjI2YTc4NjExNzU2ZTZiYzQyNjY)
 

## Installation
The workflow is splitted in two parts: the training of the model, which is done using a Python framework and the downstream analysis which is done using an R package. If you don't like Python, don't worry, the whole process can be done from R!  
In any case, you need to install both packages, as follows:

### Python package 
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
**(1) Fitting step**: train the model with the multi-omics data to disentangle the heterogeneity into latent factors.  
**(2) Characterisation step**: once the factors are inferred they need to be characterised in terms of technical or biological sources of variation.  

### Step 1: Fitting the model
There are two ways of doing this:
* **Using the command line tool**: modify and run the script [run_basic.sh](mofa/run/run_basic.sh). If you are very familar with the model and want to play with more advanced options, you can use instead [run_advanced.sh](mofa/run/run_advanced.sh).
* **Using the R wrapper**: for the ones not comfortable with the command line we built an R wrapper. See [the vignette](MOFAtools/vignettes/MOFA_example_CLL.Rmd).

Important note: the core framework of MOFA is implemented in Python, so no matter which approach you folllow, you need to install the Python package first.  

No matter which option you went for, if everything is successful, you should observe an output analogous to the following:

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

### Step 2: Downstream analysis: annotation of factors
Once the heterogeneity of the data set is reduced into a set of factors, you need to understand what are they and relate them to technical or biological sources of variability.

We have built a semi-automated pipeline based on our experience annotating factors:  
(1) **Disentangling the heterogeneity**: calculation of variance explained by each factor in each view.  
(2) **Inspection of top weighted features**: for example, if a factor is associated to the presence of a chromosomal duplication, the mRNA data will have very high loadings for genes located in that particular chromosome.  
(4) **Feature set enrichment analysis**: using for example gene ontologies.  
(4) **Visualisation of the samples in the factor space**: similarly to what is done in Principal Component Analysis, it is useful to plot the factors against each other and color using known covariates.  

An example workflow is provided in [the vignette](MOFAtools/vignettes/MOFA_example_CLL.Rmd). From R the vignette can be explored using: 
```r
browseVignettes("MOFAtools")
```

## Contact
The package is maintained by Britta Velten (britta.velten@embl.de) and Ricard Argelaguet (ricard@ebi.ac.uk). 
Please, contact ask for problems, comments or suggestions.


