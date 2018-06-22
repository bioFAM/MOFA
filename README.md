# MOFA: Multi-Omics Factor Analysis

MOFA is a factor analysis model that provides a **general framework for the integration of multi-omic data sets** in a completely unsupervised fashion.  
Intuitively, MOFA can be viewed as a versatile and statistically rigorous generalization of principal component analysis (PCA) to multi-omics data. Given several data matrices with measurements of multiple ‘omics data types on the same or on overlapping sets of samples, MOFA infers an **interpretable low-dimensional data representation in terms of (hidden) factors**. These learnt factors represent the driving sources of variation across data modalities, thus facilitating the identification of cellular states or disease subgroups.  

Once trained, the model output can be used for a range of downstream analyses, including the visualisation of samples in factor space, the automatic annotation of factors using (gene set) enrichment analysis, the identification of outliers (e.g. due to sample swaps) and the imputation of missing values.  

For more details you can read our paper: http://msb.embopress.org/cgi/doi/10.15252/msb.20178124
<p align="center"> 
<img src="images/logo.png" style="width: 50%; height: 50%"/>​
</p>



## News
- 21/06/2018 Beta version released
- 10/11/2017 We created a Slack group to provide personalised help on running and interpreting MOFA, [this is the link](https://join.slack.com/t/mofahelp/shared_invite/enQtMjcxNzM3OTE3NjcxLTkyZmE5YzNiMDc4OTkxYWExYWNlZTRhMWI2OWNkNzhmYmNlZjJiMjA4MjNiYjI2YTc4NjExNzU2ZTZiYzQyNjY)
 

## Installation
MOFA is run exclusively from R, but it requires some python dependencies that you need to install. Here is how to do it:

### Python dependencies 
We are on the process of uploading it to PyPI. For now, you can install it using:
```r
pip install git+git://github.com/bioFAM/MOFA
```
Or clone the repository and then install it using the setup.py:
```r
git clone https://github.com/bioFAM/MOFA
python setup.py install
```

### R package
The easier way to install the R package is via github:
```r
devtools::install_github("bioFAM/MOFA", subdir="MOFAtools")
```

Alternatively, you can clone the repository and install locally:
```r
git clone https://github.com/bioFAM/MOFA
R CMD build MOFAtools
R CMD install MOFAtools
```

## Tutorials/Vignettes
We currently provide three example workflows:

* [Integration of multi-omics cancer data](http://htmlpreview.github.com/?https://github.com/bioFAM/MOFA/blob/master/MOFAtools/vignettes/MOFA_example_CLL.html).
* [Integration of single-cell multi-omics data](https://cdn.rawgit.com/bioFAM/MOFA/9eee74b7/MOFAtools/vignettes/MOFA_example_scMT.html).
* [Integration of simulated data](http://htmlpreview.github.com/?https://github.com/bioFAM/MOFA/blob/master/MOFAtools/vignettes/MOFA_example_simulation.html): this tutorial is focused on model selection and robustness. For details on the down-stream analyses have a look at one of the two workflows above.

We are preparing the following workflows, to be released soon:
* Imputation.
* Prediction of clinical covariates.

If there is any tutorial that you would like us to do, or if you want to share your analysis with MOFA, please contact us.


## MOFA workflow

The workflow of MOFA consists of two steps:  
**(1) Fitting step**: train the model with the multi-omics data to disentangle the heterogeneity into a small number of latent factors.  
**(2) Downstream analysis**: once the factors are inferred they need to be characterised as technical or biological sources of variation by looking at the corresponding weights, doing (gene set) enrichment analysis, plotting the factors, correlating factors with known covariates, etc. Also, one can do imputation of missing values and prediction of clinical outcomes using the latent factors.

<p align="center"> 
<img src="images/workflow.png">
</p>

A cheatsheet with all **relevant methods**, together with a short description, can be found [here](https://github.com/bioFAM/MOFA/blob/master/MOFAtools/CheatSheet.md)  

### Step 1: Fitting the model
First you need to create the MOFA object with your input data, and subsequently you need to train the model. Everything is explained in the vignettes.  
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
Trial 1, Iteration 2: time=0.10 ELBO=-283729.31, deltaELBO=62225.6421, Factors=10
Trial 1, Iteration 3: time=0.10 ELBO=-257427.42, deltaELBO=26301.8893, Factors=10
...
Trial 1, Iteration 100: time=0.07 ELBO=-221171.01, deltaELBO=0.0998, Factors=10

Converged!
```

There are two important quantities to keep track of: 
* **Number of factors**: you can choose whether to fix the number or factors or let the model automatically learn the dimensionality of the latent space.
* **deltaELBO**: this is the convergence statistic. Once the deltaELBO decreases below a threshold (close to zero), training will end and the model will be saved as an .hdf5 file. Then, you are ready to start the downstream analysis.

### Step 2: Downstream analysis: disentangle the variability between omics
MOFA disentangles the heterogeneity of a high-dimensional multi-omics data set into a set of latent factors that capture global sources of variation.  
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

Please refer to the vignettes or the paper for details on the different analysis.  

### Step 4: Using the factors to get biological insights in downstream analysis
The latent factors can be used for several purposes, such as:  
(1) **Dimensionality reduction**: similar to PCA, dimensionality reduction visualisations can be obtained by plotting the Factors against each other.  
(2) **Imputation**: Factors can be used to predict missing values, including entire missing assays.  
(3) **Predicting clinical response**: if the factors capture phenotypical information, they can capture clinical covariates of interest.  
(4) **Regressing out technical effects**: if a factor is capturing an undesired technical effect, its effect can be regressed out from your original data matrix.  

Please refer to the vignettes or the paper for details on the different analysis.  

## Frequently asked questions

**(Q) I get the following error when installing the R package:**
```
ERROR: dependencies 'pcaMethods', 'MultiAssayExperiment' are not available for package 'MOFAtools'
```
These two packages are available from Bioconductor, not CRAN. You can install them from R as follows:
```
source("https://bioconductor.org/biocLite.R")
biocLite(c('pcaMethods', 'MultiAssayExperiment'))
```

**(Q) I get the following error when running MOFA:**  
```
AttributeError: 'module' object has no attribute 'core.entry_point
```
This means that either:  
(1) you did not install the mofa Python package (follow instructions above)  
(2) you have multiple python installations and R is not detecting the correct one where mofa is installed. You need to find out the right Python interpreter, which usually will be the one you get when running `which python` in the terminal. You can test if the mofa packaged is installed by running INSIDE python: `import mofa`.  
Once everything is figured out, specify the following at the beginning of your R script:
```
library(reticulate)
use_python("YOUR_PYTHON_PATH")
```
You can read more about the [reticulate](https://rstudio.github.io/reticulate/) package and [how it integrates Python and R](https://rstudio.github.io/reticulate/articles/versions.html)


**(Q) I hate R, can I do MOFA only with Python?**  
Yes you can, and we recommend this for training, as it will be slightly faster. See [this template script](https://github.com/bioFAM/MOFA/blob/master/mofa/run/python_template.py). However, we do not provide downstream analysis functions with Python.


**(Q) How many factors should I use?**  
Similar to other Factor Analysis models, this is a hard question to answer. It depends depends on the data set and the aim of the analysis. As a general rule, the bigger the data set, the higher the number of factors that you will likely retrieve, and the less the variance that will be explained per factor.
If you want to get an overview on the major sources of variability then use a small number of factors (K<=15). If you want to capture small sources of variability, for example to improve imputation performance or for eQTL mapping, then go for a large number of factors (K>=50)


**(Q) How does MOFA handle missing values?**  
It simpy ignores them, there is no a priori imputation step required. In fact, matrix factorisation models are known to be very robust to the presence of large amounts of missing values. 

**(Q) Should I do any filtering to the input data?**  
It is not mandatory, but it is highly recommended to filter lowly variable features. It makes the model more robust and speeds up the training.

**(Q) My data sets have different dimensionalities, does this matter?**  
Yes, this is important. Bigger data modalities will tend to be overrepresent in the MOFA model. It is good practice to filter features (based for example on variance) in order to have the different dimensionalities within the same order of magnitudes. If this is unavoidable, take into account that the model has the risk of missing (small) sources of variation unique to the small data set.


**(Q) Can MOFA automatically learn the number of factors?**  
Yes, MOFA can automatically learn the number of factors, but a hyperparameter needs to be provided. The user needs to specify a minimum value of fraction of variance explained that is considered meaningful. Then, MOFA will actively remove factors (during training) that explain less than the specified amount of variance.
If you have no idea on what to expect, it is better to start with a fixed number of factors.


**(Q) What data modalities can MOFA cope with?**  
* Continuous data: should be modelled using a gaussian likelihood. For example, log normalised RNA-seq data or M-values of bulk methylation data
* Binary data: should be modelled using a bernoulli likelihood. For example, somatic mutations or single-cell methylation data.
* Count data: should be modelled using a poisson likelihood. For example, copy number variation or scRNA-seq UMI data.
The use of non-gaussian likelihoods require further approximations and are not as accurate as the gaussian likelihood. Hence, if your data can be safely transformed to match the gaussian likelihood assumptions, this is always recommended. For example log-transform and variance stabilisation of bulk RNA-seq data or M-value computation in DNA methylation data.

**(Q) How do I assess convergence?**  
MOFA is trained using variational bayes, a fast inference framework that consists on optimising a statistica called the Evidence Lower Bound (ELBO). The model uses the change in ELBO (deltaELBO) to assess convergence. A model is defined to be converged when deltaELBO is close to 0. For a quick exploratory analysis, we suggest a convergence threshold between 1 to 10.

**(Q) What input formats are allowed?**  
The data has to be input in two possible formats: 
* Bioconductor: a [MultiAssayExperiment](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html) object
* Base R approach: a list of matrices where features are rows and samples are columns. Examples are shown in the vignettes.

**(Q) Does MOFA always converge to the same solutions?**  
No, as occurs in most complex Bayesian models, they are not guaranteed to always converge to the smae (optimal) solution.
In practice, however, we observed that the solutions are highly consistent, particularly for strong factors. However, one should always assess the robustness and do a proper model selection. We are currently preparing a vignette on this.


## Contact
The package is maintained by Britta Velten (britta.velten@embl.de) and Ricard Argelaguet (ricard@ebi.ac.uk).  
We created a Slack group to provide personalised help on running and analysing MOFA, [this is the link](https://join.slack.com/t/mofahelp/shared_invite/enQtMjcxNzM3OTE3NjcxLTkyZmE5YzNiMDc4OTkxYWExYWNlZTRhMWI2OWNkNzhmYmNlZjJiMjA4MjNiYjI2YTc4NjExNzU2ZTZiYzQyNjY).  
Please, reach us for problems, comments or suggestions.

