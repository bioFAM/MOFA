# MOFA: Multi-Omics Factor Analysis

MOFA is a factor analysis model that integrates multi-omic data sets in an unsupervised fashion. It automatically discovers the main sources of both technical and biological variability and it identifies whether they are unique or shared between several -omic layers. The model allows for the presence of missing values, including absence of entire assays for a given samples. Also, we provide an accessible and user-friendly package for a semi-automated annotation of factors.

For more details you can read our preprint:


## News
- 01/09/2017 Paper uploaded to bioRxiv and submitted for review


## Installation

### Python package (to train the model)
The easiest way to install DeepCpG is to use PyPI:
```r
```
Alternatively, you can directly install from the repository:
```r
pip install git+git://github.com/PMBio/MOFA
```
Or clone the repository and then install MOFA using the setup.py

```r
git clone https://github.com/PMBio/MOFA
python setup.py install
```

### R package (to analyse the model)
The easier way to install the package is via github:
```r
devtools::install_github("PMBio/MOFA", subdir="MOFAtools")
```

Alternatively, you can clone the repository and insall via R CMD:
```r
git clone https://github.com/PMBio/MOFA
R CMD build MOFAtools
R CMD install MOFAtools
```

## MOFA workflow

### Training the model
The training of the model is performed using the python framework. It can also be rin within the R framework (by indirectly calling python), but is not as fast or stable.

To train the model, simply store your input matrices as text files and modify the following script accordingly:
[run_basic.sh](MOFA/run/run_basic.sh)

If you are very familar with the model and want to play with more advaned options, you can edit and run the following script:
run_advanced.sh

Once the model has converged and is saved as an .hdf5 file, then you should move to the R framework, where we implemented a pipeline for a semi-automated factor annotation.

### Downstream analysis: annotation of factors

