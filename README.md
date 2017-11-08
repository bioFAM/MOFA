# MOFA: Multi-Omics Factor Analysis

MOFA is a factor analysis model that provides a **general framework for the integration of multi-omic data sets** in a completely unsupervised fashion.  
It calculates a low dimensional representation of the multiple views in terms of a small number of latent factors that capture the main sources of variability, which might be masked in the noisy and complex high-dimensional representation. Furthermore, it reveals whether the different factors are unique to a single -omic or shared between multiple -omics.  
The model can take as input multiple data modalities (continuous, binary and count data) and it is flexible to the presence of missing values, including absence of entire assays.

For more details you can read our preprint:
<p align="center"> 
<img src="logo.png">
</p>

## News
- XX/09/2017 Paper uploaded to bioRxiv and submitted for review

## Installation
The workflow is splitted in two parts: the training of the model, which is done using a Python package and the downstream analysis which is done using an R package.
They both can be installed as follows:

### Python package 
The easiest way to install MOFA is to use PyPI:
```r
```
Alternatively, you can directly install from the repository:
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
Once the heterogeneity of the data set is reduced into a set of factors, you need to understand what they are and relate them to technical or biological sources of variability. 

We have built a semi-automated pipeline based on our experience annotating factors:  
(1) **Disentangling the heterogeneity**: calculation of variance explained by each factor in each view.  
(2) **Inspection of top weighted features**: for example, if a factor is associated to the presence of a chromosomal duplication, the mRNA data will have very high loadings for genes located in that particular chromosome.  
(4) **Feature set enrichment analysis**: using for example gene ontologies.  
(4) **Visualisation of the samples in the factor space**: similarly to what is done in Principal Component Analysis, it is useful to plot the factors against each other and color using known covariates.  

An example workflow is provided in [the vignette](MOFAtools/vignettes/MOFA_example_CLL.Rmd). From R the vignette can be explored using: 
```r
browseVignettes("MOFAtools")
```

