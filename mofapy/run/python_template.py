import pandas as pd
from mofapy.core.entry_point import entry_point

# Load the data 
# Be careful to use the right delimiter, and make sure that you use the right arguments from pandas.read_csv to load the row names and column names, if appropriate.
M = 3 # Number of views
data =  [None]*M
data[0] = pd.read_csv("(...)/view_0.txt", delimiter=" ").astype(pd.np.float32)
data[1] = pd.read_csv("(...).txt", delimiter=" ").astype(pd.np.float32)
data[2] = pd.read_csv("(...).txt", delimiter=" ").astype(pd.np.float32)

# Initialise entry point
ep = entry_point()

# Set data
ep.set_data(data)

## Set model options ##
# factors: number of factors. By default, the model does not automatically learn the number of factors. 
# 	If you want the model to do this (based on a minimum variance explained criteria), set `TrainOptions$dropFactorThreshold` to a non-zero value.
# likelihoods: list with the likelihood for each view. Usually we recommend: 
#	- gaussian for continuous data
# 	- bernoulli for binary data
#  	- poisson for count data
# 	If you are using gaussian likelihoods, we recommend centering the data (specified in data_options) and setting learnIntercept to False. 
# 	However, if you have non-gaussian likelihoods, learning an intercept factor is important
# sparsity: boolean indicating whether to use sparsity. 
# 	This is always recommended, as it will make the loadings more interpretable.
ep.set_model_options(factors=25, likelihoods=["gaussian","bernoulli","poisson"], sparsity=True)

## Set data options ##

# view_names: list with view names
# center_features: boolean indicating whether to center the features to zero mean. 
# 	This only works for gaussian data. Default is TRUE.
# scale_views: boolean indicating whether to scale views to have the same unit variance. 
# 	As long as the scale differences between the data sets is not too high, this is not required. Default is False.
# RemoveIncompleteSamples: boolean indicating whether to remove samples that are not profiled in all omics. 
# 	We recommend this only for testing, as the model can cope with samples having missing assays. Default is False.
ep.set_data_options(center_features=True, scale_views=False, RemoveIncompleteSamples=False)

# Parse the data (optionally center or scale, do some QC, etc.)
ep.parse_data()

# Define training options
# iter: numeric value indicating the maximum number of iterations. 
# 	Default is 1000, we recommend setting this to a large value and using the 'tolerance' as convergence criteria.
# tolerance: numeric value indicating the convergence threshold based on the change in Evidence Lower Bound (deltaELBO). 
# 	For quick exploration we recommend this to be around 1.0, and for a thorough training we recommend a value of 0.01. Default is 0.1
# dropR2: numeric hyperparamter to automatically learn the number of factors. 
# 	It indicates the threshold on fraction of variance explained to consider a factor inactive and automatically drop it from the model during training. 
# 	For example, a value of 0.01 implies that factors explaining less than 1% of variance (in each view) will be dropped. 
# 	Default is 0, which implies that only factors that explain no variance at all will be removed
# elbofreq: frequency of computation of the ELBO. It is useful to assess convergence, but it slows down the training.
# verbose: boolean indicating whether to generate a verbose output.
# seed: random seed. If None, it is sampled randomly
ep.set_train_options(iter=1000, tolerance=0.01, dropR2=0.00, elbofreq=1, verbose=False, seed=2018)

# Define prior distributions
ep.define_priors()

# Define initialisations of variational distributions
ep.define_init()

# Parse intercept factor
ep.parse_intercept()

# Train the model
ep.train_model()

# Save the model
outfile="/Users/mofa/test.hdf5"
ep.save_model(outfile)

# This is the end of the python pipeline.
# Now you are ready to switch to the R framework for downstream analysis:
# (in R) model <- loadModel("/Users/mofa/test.hdf5")
