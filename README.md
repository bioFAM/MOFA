

# MOFA: Multi-Omics Factor Analysis

MOFA is a factor analysis model that integrates multi-omic data sets in an unsupervised fashion. It automatically discovers the main sources of both technical and biological variability. Moreover, it determines whether each factor of variation is unique or shared between several -omic layers. 

MOFA builds upon the Group Factor Analysis (GFA) statistical framework (Virtanen 2012, Klami 2015, Zhao 2016, Khan 2014, Bunte 2016) and introduces a set of key extensions that enable the broad applicability in multi-omic studies (Table S1). In particular, we combine a two-level sparse model providing interpretability with a scalable inference approach. In addition, we enable non-Gaussian views in the model as these are frequently encountered in multi-omic studies (e.g. binary data such as somatic mutations or count data such as copy number variation). Importantly, as biological data is often incomplete and samples might be missing in a subset of views we allow for presence of missing values, including absence of entire views. Finally, we provide accessible and user-friendly software and downstream analysis functions for a semi-automated annotation of factors, facilitating the wide-ranging application of this method to any kind of multi-omic data set.


![alt text](logo.pdf)


## Why did we develop MOFA?
Despite the increasing availability of measurements from more than one molecular layer, the joint analysis of multi-view data remains challenging for several reasons (Chang, 2013). First, data modalities collected using different techniques generally exhibit heterogeneous properties and have to be modelled under different statistical assumptions. Second, as the number of measured features increases, appropriate regularization strategies are required to avoid over-fitting. Third, traditional approaches based on marginal correlations fail due to amounts of noise in single features as well as a massive multiple testing burden. Finally, missing data is commonly encountered, both in terms of single values as well as individual samples that miss a certain assay. Addressing these challenges and developing robust and powerful methods for the integrative analysis of multi-layer datasets is therefore essential to unfold the full potential of multi-omic studies.  

One popular approach for multi-omic data integration is iCluster (Shen 2009) and its extension (Mo 2012), a latent variable model which jointly clusters samples and identifies cluster-relevant features across data sets. However, due to its focus on clustering, the statistical modelling underpinning iCluster is not well-suited to disentangle the sources of variation and discover which factors are unique to specific modalities and which factors are shared between multiple data modalities. Disentangling the axes of variation is critical for obtaining insight into the biological processes that drive variability in the data. Further, iCluster can be limited because of computationally expensive inference procedures and the inability to handle missing values, which are routinely encountered in multi-omic data sets.

### NEWS
- 01/09/2017 Paper uploaded to bioRxiv and submitted for review

### INSTALLATION

### RUNNING THE MODEL

### VIGNETTE
To demonstrate the potential of MOFA we applied it to a recent multi-omic study on chronic lymphocytic leukaemia (CLL) (Dietrich, Oles, Lu et al). In this data, patient samples are characterized by several omic techniques, including genome sequencing, RNA-Seq, DNA-methylation arrays and ex-vivo drug response assays. Using MOFA we are able to identify major drivers of patient heterogeneity and their molecular basis in a completely unsupervised manner. Naturally, we capture drivers which are related to important clinical disease subtypes, e.g. defined by the mutational status of the immunoglobulin heavy-chain variable-region (IGHV) (Fabbri 2016). Furthermore, we also find new drivers which can be linked to clinical outcome or cancer-associated pathway such as stress response and reactive oxygen species.
