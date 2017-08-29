# MOFA: Multi-Omics Factor Analysis

MOFA is a factor analysis model that provides a general framework for the integration of multi-omic data sets in a completely unsupervised fashion.

MOFA calculates a low dimensional representation of the multiple views in terms of a small number of latent factors that capture the main sources of variability, which might be masked in the noisy and complex high-dimensional representation. Furthermore, it reveals whether the different factors are unique to a single -omic or shared between multiple -omics. This allows a more comprehensive characterisation of the biological and technicals sourdes of variation underlying the data.  
The model can take as input multiple data modalities (continuous, binary and count data) and it is flexible to the presence of missing values, including absence of entire assays.

For more details you can read our preprint:
<p align="center"> 
<img src="logo.png">
</p>

## News
- 01/09/2017 Paper uploaded to bioRxiv and submitted for review

## Installation

### Python package (to train the model)
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

### R package (for downstream analysis)
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

It consists of two steps: first, you need to train the model, which is done via the Python implementation ([MOFA](MOFA)) that is extremely easy to use. Once the model is trained, you should move to our R framework ([MOFAtools](MOFAtools)), where we implemented a pipeline for a semi-automated annotation of factors coupled with beautiful visualisations. The following sections describe a bit more in detail how to use MOFA.

### Training the model
Once you have installed the Python version of MOFA, you simply need to store your input matrices as text files, and modify the following script accordingly:
[run_basic.sh](MOFA/run/run_basic.sh)  
If you are very familar with the model and want to play with more advaned options, you can edit and run the following script:  
[run_advanced.sh](MOFA/run/run_advanced.sh)

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

Trial 1, Iteration 853: time=0.07 ELBO=-221171.01, deltaELBO=0.0998, Factors=10, Covariates=1

Converged!
...


```

The change in ELBO (deltaELBO) is the important quantity that is being monitored, and when this reaches a default value of 0.1 the model will converge and it will be saved as an .hdf5 file.

The model is now trained and you are ready to start the analysis with the R package.

### Downstream analysis: annotation of factors
The most important thing you want to do first is to annotate factors.


