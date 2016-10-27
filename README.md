single-cell Group Factor Analysis (scGFA)
======

The `scgfa` package implements the single-cell Group Factor Analysis problem as described in (ref).

# Introduction

## Motivation
Single-cell RNA sequencing is becoming a well-established routine that is revolutionising our understanding of cellular phenotypes. Interestingly, other data modalities are also starting to be assayed at the single-cell level, including epigenetics, proteomics and metabolomics, raising the question of how to jointly analyse these set of complex high-dimensional data sets using  a statistically rigorous framework.

## The multi-view learning problem
scGFA is based on the traditional latent variable model called factor analysis (FA), in which observed data vectors are assumed to be constructed from a set of hidden variables using a linear mapping:
![eq1](xx)
EXPLAINS VARIANCE
However, FA is a single-view model in the sense that it takes as input a single data matrix. However, in some occasions, data is collected in the from diverse domains or layers and so they can exhibit heterogeneous properties. Each of these layers of information (called views) can be analysed separately using single-view learning methods, but the current challenge is the integration of all the different views into a single model that is able to disentangle the variation unique to a single view and the variation shared between two or more views. This is refered to as the multi-view learning problem (ref).
Traditionally, the multi-view problem has been approached by concatenating the different datasets into a single large matrix, which was subsequently analysed using conventional machine learning algorithms. However, this concatenation should be avoided for several reasons. First, the scale of noise is usually differ- ent between views and one view might end up overepresented in the solution. Second, this concatenation leads to a very wide matrix that causes overfitting in the case of a small number of samples. Third, it is hard to distinguish from which view the signal is coming from and whether it is shared by multiple views or it is unique to a single view. For this reason, statistically riguruous multi-view models are required


# Model definition
Here we present single-cell Bayesian Group Factor Analysis (scBGFA), a multi-view generalisation of factor analysis suited to the analysis of noisy multi-view single-cell sequencing data.
scGFA performs a joint dimensionality reduction for all views, yielding a low dimensional representation of the data which hopefully captures an inherent structure that might be masked by the noisy high-dimensional representation. Furthermore, scBGFA disentangles the variation unique to a single view and the variation shared between two or more views, thereby revealing hidden sources of covariation between different data modalities.

## Group Factor Analysis 
The core of the scGFA model is the Group Factor Analysis (GFA) model as published in (ref).
Briefly, GFA is a multi-view extension of factor analysis that can handle an arbitrary number of views. The key part of the model is a group-wise sparsity imposed on the weight matrix, which allows the latent factors to be active in all views, a subset of views or.


## Extensions to GFA 

### Element-wise spike and slab sparsity
### Non-gaussian likelihoods
### Generalisation of the noise model