# MOFA: Multi-Omics Factor Analysis

## Important notice: MOFA v1 (this repository) is officially depreciated, please switch to [MOFA v2](https://biofam.github.io/MOFA2/)

<br>

MOFA is a factor analysis model that provides a **general framework for the integration of multi-omic data sets** in a completely unsupervised fashion.  
Intuitively, MOFA can be viewed as a versatile and statistically rigorous generalization of principal component analysis (PCA) to multi-omics data. Given several data matrices with measurements of multiple -omics data types on the same or on overlapping sets of samples, MOFA infers an **interpretable low-dimensional data representation in terms of (hidden) factors**. These learnt factors represent the driving sources of variation across data modalities, thus facilitating the identification of cellular states or disease subgroups.  

Once trained, the model output can be used for a range of downstream analyses, including the visualisation of samples in factor space, the automatic annotation of factors using (gene set) enrichment analysis, the identification of outliers (e.g. due to sample swaps) and the imputation of missing values.  

For more details you can read our paper: http://msb.embopress.org/cgi/doi/10.15252/msb.20178124
<p align="center"> 
<img src="images/logo.png" style="width: 50%; height: 50%"/>​
</p>


## News
- 01/01/2020 Beta version of MOFA+ software is [available](https://github.com/bioFAM/MOFA2). We recommend all users to switch to MOFA+.
- 01/12/2019 The new version of MOFA (MOFA+) manuscript is published in [bioRxiv](https://www.biorxiv.org/content/10.1101/837104v1.article-metrics) 
- 03/05/2019 MOFA is [available in Bioconductor](https://bioconductor.org/packages/devel/bioc/html/MOFA.html) (only for R>=3.6).
- 10/01/2019 Python package uploaded to PyPI (https://pypi.org/project/mofapy/)
- 21/06/2018 Beta version released
- 20/06/2018 Paper published: http://msb.embopress.org/content/14/6/e8124
- 10/11/2017 We created a Slack group to provide personalised help on running and interpreting MOFA, [this is the link](https://join.slack.com/t/mofahelp/shared_invite/enQtMjcxNzM3OTE3NjcxLWNhZmM1MDRlMTZjZWRmYWJjMGFmMDkzNDBmMDhjYmJmMzdlYzU4Y2EzYTI1OGExNzM2MmUwMzJkZmVjNDkxNGI)
 

## Installation
MOFA is run exclusively from R, but it requires some python dependencies.

### Python dependencies 
Python dependencies can be installed using pip (from the Unix terminal)
```r
pip install mofapy
```

Alternatively, they can be installed from R itself using the reticulate package:
```r
library(reticulate)
py_install("mofapy", envname = "r-reticulate", method="auto")
```

### MOFAdata R data package
For illustration purposes we provide several data sets that are used in the vignettes of the MOFA package. Can be installed using R:
```r
devtools::install_github("bioFAM/MOFAdata", build_opts = c("--no-resave-data"))
```

### MOFA R package
This is the core software itself. Can be installed using R:
```r
devtools::install_github("bioFAM/MOFA", build_opts = c("--no-resave-data"))
```


### Reticulate configuration

Before running MOFA, you need to make sure that `reticulate` is pointing to the correct python binary (or conda environment).  
This can become tricky when you have multiple conda environments and versions of Python installed:
```r
library(reticulate)

# Using a specific python binary
use_python("/home/user/python", required = TRUE)

# Using a conda enviroment called "r-reticulate"
use_condaenv("r-reticulate", required = TRUE)
```

For more details on how to set up the reticulate connection, see: https://rstudio.github.io/reticulate/

## Tutorials/Vignettes
We currently provide three example workflows:

* **Integration of multi-omics cancer data**: a cohort of 200 chronic lymphocytic leukaemia patients. This is the main data set analysed in the [paper](http://msb.embopress.org/cgi/doi/10.15252/msb.20178124). Load it using `vignette("MOFA_example_CLL")`.
* **Integration of single-cell multi-omics data**: single-cell profiling of DNA methylation and RNA expression in roughly 100 pluripotent stem cells. This is the secondary data set analysed in the [paper](http://msb.embopress.org/cgi/doi/10.15252/msb.20178124). Load it using `vignette("MOFA_example_scMT")`.
* **Model selection and robustness with simulated data**: this tutorial is focused only on how to perform model selection and assess robustness. Load it using `vignette("MOFA_example_simulated")`

If you have problems loading the vignettes, you can find the html files [here](https://bioconductor.org/packages/devel/bioc/html/MOFA.html)

## Cheatsheet
A list with all relevant functions, together with a short description, can be found at the end of the introductory [vignette](https://bioconductor.org/packages/devel/bioc/vignettes/MOFA/inst/doc/MOFA.html) (`vignette("MOFA")`). 

## MOFA workflow

<p align="center"> 
<img src="images/workflow.png">
</p>

### Step 1: Fitting the model
First you need to create the MOFA object with your input data, and subsequently train the model.
If everything is successful, you should observe an output analogous to the following:
```

#############################################
## Running trial number 1 with seed 642034 ##
#############################################

Trial 1, Iteration 1: time=0.08 ELBO=-345954.96, Factors=10
Trial 1, Iteration 2: time=0.10 ELBO=-283729.31, deltaELBO=62225.6421, Factors=10
Trial 1, Iteration 3: time=0.10 ELBO=-257427.42, deltaELBO=26301.8893, Factors=10
...
Trial 1, Iteration 100: time=0.07 ELBO=-221171.01, deltaELBO=0.0998, Factors=10

Converged!
```

### Step 2: Downstream analysis: disentangle the variability between omics
MOFA disentangles the heterogeneity of a high-dimensional multi-omics data set in terms of a small number of latent factors that capture the global sources of variation. Importantly, MOFA quantififes the variance explained of each of the factors in the different omics. An example is shown in the plot below:
<p align="center"> 
<img src="images/varExplained.png" style="width: 50%; height: 50%"/>​
</p>


### Step 3: Annotation of factors
The next step is to try and interpret what the factors are. We have built a semi-automated pipeline to allow the exploration of the latent space:  
(1) **Visualisation of the samples in the factor space**: as in Principal Component Analysis, it is useful to plot the factors against each other and color the samples using known covariates such as batch, sex, clinical information, etc.  
(2) **Correlation of factors with (clinical) covariates**  
(2) **Inspection of the loadings**: loadings provide a measure of feature importance for each factor.  
(3) **Feature set enrichment analysis**: the inspection of loadings can sometimes be challenging, particularly when having large amounts of features. Summarising genes in terms of biological pathways can be useful in such cases.  

Please refer to the vignettes for details on the different analysis.  

### Step 4: Using the factors in downstream analysis
The latent factors can be used for several purposes, such as:  
(1) **Non-linear dimensionality reduction**: the latent factors can be feed into non-linear dimensionality reduction techniques such as UMAP or t-SNE. This is very powerful because you can detect variability or stratifications beyond the RNA expression!  
(2) **Imputation**: factors can be used to predict missing values, including entire missing assays.  
(3) **Predicting clinical response**: factors can be feed into Cox models to predict patient survival.  
(4) **Regressing out technical variability**: if a factor is capturing an undesired technical effect, its effect can be regressed out from your original data matrix.  
(5) **Clustering**: clustering in the latent space is much more robust than in the high-dimensional space.  
(6) **factor-QTL mapping**: factors are a compressed and denoised representation of your samples. This is a much better proxy for the phenotype than the expression of individual genes. Hence, a very promising area is to do eQTL's with the factors themselves! See [this paper]() for an example (https://www.nature.com/articles/ng.3624).

Again, refer to the vignettes for details on the different analysis.

## Frequently asked questions

**(Q) How do I normalise the data?**  
Always try to remove any technical source of variability before fitting the model.  
For example, for count-based data such as RNA-seq or ATAC-seq we recommend size factor normalisation + variance stabilisation. For microarray DNA methylation data, make sure that samples have no differences in the average intensity.

If this is not done correctly, the model will learn a very strong Factor 1 that will capture this variability, and more subtle sources of variation will be harder to identify.  
We have implemented a function called `regressCovariates` that allows the user to regress out a covariate using linear models. See the documentation and the CLL vignette for examples.

**(Q) I get the following error when installing the R package:**  
```
ERROR: dependencies 'pcaMethods', 'MultiAssayExperiment' are not available for package 'MOFA'
```
You probably tried to install them using `install.packages()`. These packages should be installed from Bioconductor.

**(Q) I get one of the following errors when running MOFA:**  
```
AttributeError: 'module' object has no attribute 'core.entry_point

Error in py_module_import(module, convert = convert) :
 ModuleNotFoundError: No module named 'mofapy'
```
First thing: restart R and try again. If the error still holds, this means that either:  
(1) you did not install the mofa Python package (see instructions above).
(2) you have multiple python installations and R is not detecting the correct one where mofa is installed. You need to find out the right Python interpreter, which usually will be the one you get when running `which python` in the terminal. You can test if the mofa packaged is installed by running INSIDE python: `import mofapy`.  
Once everything is figured out, specify the following at the beginning of your R script:
```
library(reticulate)
use_python("YOUR_PYTHON_PATH", required=TRUE)
```
You can also use `use_conda` instead of `use_python` if you work with conda environments. Read more about the [reticulate](https://rstudio.github.io/reticulate/) package and [how it integrates Python and R](https://rstudio.github.io/reticulate/articles/versions.html)

**(Q) I hate R, can I do MOFA only with Python?**  
Nop. You can use Python to train the model, see [this template script](https://github.com/bioFAM/MOFA/blob/master/mofapy/run/python_template.py). However, we currently do not provide downstream analysis functions in Python. We strongly recommend that you use our MOFA R package for this.

**(Q) How many factors should I learn?**  
Similar to Principal Component Analysis and other latent variable models, this is a hard question to answer. It depends on the data set and the aim of the analysis. If you want to get an overview on the major sources of variability then use a small number of factors (K<=10). If you want to capture small sources of variability, for example to do imputation or eQTL mapping, then go for a large number of factors (K>25).

**(Q) How many samples do I need?**  
At least more than 15.

**(Q) Can MOFA automatically learn the number of factors?**  
Yes, but the user needs to specify a minimum value of % variance explained. Then, MOFA will actively remove factors (during training) that explain less than the specified amount of variance.
If you have no idea on what to expect, it is better to start with a fixed number of factors and set the % variance threshold to 0.

**(Q) Can I put known covariates in the model?**  
Combining known covariates with latent factors is technically possible, but we extensively tested this functionality and it was not yielding good results. The reason is that covariates are usually discrete labels that do not reflect the underlying molecular biology. For example, if you introduce age as a covariate, but the actual age is different from the “molecular age”, the model will simply learn a new factor that corresponds to this “latent” molecular age, and it will drop the covariate from the model.  
We recommend that you learn the factors in a completely unsupervised manner and then relate them to the biological covariates via visualisation or via a simple correlation analysis (see our vignettes). If your covariate of interest is indeed an important driver of variability, do not worry, MOFA will find it! 

**(Q) Should I remove undesired sources of variability (i.e. batch effects) before fitting the model?**  
Yes, if you have clear technical factors, we strongly encourage to regress it out a priori using a simple linear model. The reason for this is that the model will "focus" on the huge variability driven by the technical factors, and smaller sources of variability could be missed.
You can regress out known covaraites using the function `regressCovariates`. See the corresponding documentation and the CLL vignette for details.

**(Q) Should I do any filtering to the input data?**  
You must remove features with zero variance and ideally also features with low variance, as they can cause numerical issues in the model. In practice we generally select the top N most variable features per assay

**(Q) My data sets have different dimensionalities, does this matter?**  
Yes, this is important. Bigger data modalities will tend to be overrepresent in the MOFA model. It is good practice to filter features (based for example on variance, as lowly variable features provide little information) in order to have the different dimensionalities within the same order of magnitudes. If this is unavoidable, take into account that the model has the risk of missing (small) sources of variation unique to the small data set.

**(Q) The weights have different values between runs. Is this expected?**  
This is normal and it happens because of to two reasons. The first one is that the model does not always converge to the same exact solution (see below in the FAQ), although different model instances should be pretty similar. The second reason is that factor analysis models are rotation invariant. This means that you can rotate your factors and your weights and still find the same solution. This implies that the signs of the weight or the factors can NOT be compared across trials, only within a trial.

**(Q) What data modalities can MOFA cope with?**  
* Continuous data: modelled using a gaussian likelihood
* Binary data: modelled using a bernoulli likelihood
* Count data: using a poisson likelihood.
Importantly, the use of non-gaussian likelihoods require further approximations and are not as accurate as the gaussian likelihood. Hence, if your data can be safely transformed to match the gaussian likelihood assumptions, this is ALWAYS recommended. For example RNA-seq data is expected to be normalised and modelled with a gaussian distribution, do not input the counts directly.

**(Q) The model does not converge smoothly, and it oscillates between positive and negative deltaELBO values**  
First, check that you are using the right likelihood model (see above). Second, make sure that you have no features or samples that are full of missing values. Third, check that you have no features with zero (or very little) variance. If the problem does not disappear, please contact us via mail or the Slack group, we will provide (quick!) help.

**(Q) What input formats are allowed?**  
The data has to be input in two possible formats: 
* Bioconductor approach: a [MultiAssayExperiment](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html) object
* Base R approach: a list of matrices where features are rows and samples are columns.
Examples of both are shown in the vignettes.

**(Q) Does MOFA always converge to the same solutions?**  
No, as occurs in most complex Bayesian models, they are not guaranteed to always converge to the same (optimal) solution.
In practice, however, we observed that the solutions are highly consistent, particularly for strong factors. However, one should always assess the robustness and do a proper model selection. For this we recommend to train the model multiple times and check the robustness of the factors across the different solutions. For downstream analysis a single model can be chosen based on the best value of the Evidence Lower Bound (ELBO). We provide functions for these two steps, which are explained in the vignette *Integration of simulated data* (`vignette("MOFA_example_simulated")`).

**(Q) How does MOFA handle missing values?**  
It simpy ignores them, there is no hidden imputation step. Matrix factorisation models are known to be very robust to the presence of missing values!

**(Q) How can I do Gene Set Enrichment Analysis?**  
First, you need to create your binary gene set matrix where rows are feature sets and columns are features (genes). We have manually processed some of Reactome and MSigDB gene sets for mouse and human. Contact us if you would like to use the data.  
Then, you will have to choose a local statistic per feature (the loading, by default), a global statistic per pathway (average loading, by default), and a statistical test. The most trustworthy one is a permutation test with a long number of iterations, but this is slow and a fast parametric tests is also available. However, note that it tends to inflate the p-values due to the correlation structure between related genes (see for example [Gatti2010](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-11-574)).

## Publications using MOFA
- [Single cell multi-omics profiling reveals a hierarchical epigenetic landscape during mammalian germ layer specification](https://www.biorxiv.org/content/10.1101/519207v1)
- [Multi-Omics Factor Analysis: a framework for unsupervised integration of multi‐omics data sets](http://msb.embopress.org/content/14/6/e8124)

## Contact
The package is maintained by Ricard Argelaguet (ricard@ebi.ac.uk) and Britta Velten (britta.velten@embl.de). Please, reach us for problems, comments or suggestions. You can also contact us via a Slack group where we provide quick and personalised help, [this is the link](https://join.slack.com/t/mofahelp/shared_invite/enQtMjcxNzM3OTE3NjcxLWNhZmM1MDRlMTZjZWRmYWJjMGFmMDkzNDBmMDhjYmJmMzdlYzU4Y2EzYTI1OGExNzM2MmUwMzJkZmVjNDkxNGI).  


