# MOFA: Multi-Omics Factor Analysis

MOFA is a factor analysis model that provides a **general framework for the integration of multi-omic data sets** in a completely unsupervised fashion.  
It calculates a low dimensional representation of the multiple views in terms of a small number of latent factors that capture the main sources of variability, which might be masked in the noisy and complex high-dimensional representation. Furthermore, it reveals whether the different factors are unique to a single -omic or shared between multiple -omics. This allows a more comprehensive characterisation of the biological and technicals sources of variation underlying the data.  
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

The workflow fo MOFA consists of two steps: first, you need to train the model, which is done via a Python implementation that is extremely easy to use. Once the model is trained, you should move to our R framework, where we implemented a pipeline for a semi-automated annotation of factors coupled with beautiful visualisations.  
The following sections describe a bit more in detail how to use MOFA.

### Training the model
Once you have installed the Python version of MOFA, you simply need to store your input matrices as text files, and modify the following script accordingly:
[run_basic.sh](MOFA/run/run_basic.sh)  
If you are very familar with the model and want to play with more advaned options, you can edit and run the following script: [run_advanced.sh](MOFA/run/run_advanced.sh)

Then, you run the corresponding script and let the model converge (it can take a while if the data set is big)
```r
bash run_[basic/extended].sh
```
If everything is successful, you should observe the following output:
```

##################
## Loading data ##
##################

Loaded /Users/ricard/MOFA/MOFA/test/data/500_0.txt with dim (100,500)...

Loaded /Users/ricard/MOFA/MOFA/test/data/500_1.txt with dim (100,500)...

Loaded /Users/ricard/MOFA/MOFA/test/data/500_2.txt with dim (100,500)...

#############################################
## Running trial number 1 with seed 90545 ##
#############################################

(A few warnings, don't worry!...)
/Users/ricard/anaconda2/lib/python2.7/site-packages/numpy/lib/scimath.py:262: RuntimeWarning: divide by zero encountered in log
  return nx.log(x)
/Users/ricard/anaconda2/lib/python2.7/site-packages/MOFA/core/updates.py:287: RuntimeWarning: invalid value encountered in multiply
  lb_ps = S*theta['lnE'] + (1.-S)*theta['lnEInv']

Trial 1, Iteration 1: time=0.06 ELBO=-344149.64, Factors=10, Covariates=1

Trial 1, Iteration 2: time=0.09 ELBO=-277174.84, deltaELBO=66974.8028, Factors=10, Covariates=1

Trial 1, Iteration 3: time=0.08 ELBO=-263272.74, deltaELBO=13902.1021, Factors=10, Covariates=1

...

Trial 1, Iteration 100: time=0.07 ELBO=-221171.01, deltaELBO=0.0998, Factors=10, Covariates=1

Converged!


```
There are two important quantities to keep track of: first, the number of latent variables, that should decrease during training if you want the model to automatically learn the number of factors. Second, the difference in evidence lower bound (deltaELBO), which is used to assess model convergence. Once the model reaches a default value of 0.1 the model will converge and it will be saved as an .hdf5 file. Then, you are ready to start the analysis with the R package.

### Downstream analysis: annotation of factors
The most important thing to do first is to understand what your latent factors are. We recommend four approaches:   
(1) Correlate the factors with known covariates (sex, age, etc.)  
(2) Inspect the corresponding loadings: for example, a factor associated with variation due to the presence of a chromosomal duplication will have very high loadings in the Gene expression data for genes located in that particular chromosome.  
(3) Do a Feature Set Enrichment Analysis: if you have feature sets for a given view (for example gene ontologies for mRNA data), then you can perform a feature set enrichment analyisis for each factor:
(4) Visualise the samples in the latent space: similarly to PCA, it is useful to plot the factors against each other to detect clusters of samples.  

The following vignette shows the example we used for the publication: 

