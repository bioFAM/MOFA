#!/bin/bash

# Input files as plain text format
inFolder="/Users/ricard/data/CLL/views/minView=all"
inFiles=( "$inFolder/mut.txt" "$inFolder/viab.txt" "$inFolder/mRNA.txt" "$inFolder/meth.txt" )
delimiter=' ' # Delimiter 
header_rows=1 # Do the files contain row names?
header_cols=1 # Do the files contain column names?

# Output file as .hdf5 format
outFile=( "/Users/ricard/test/test_k20.hdf5" )

# Data options
center_features=0  # Center the features to zero-mean? (not necessary as long as learnMean=1)
scale_views=0 	   # Scale the views to unit variance (recommended, but not necessary as long as there no massive differences in scale)
RemoveIncompleteSamples=0 # Remove samples that have missing views. The model can cope with incomplete data sets

# Define likelihoods 
likelihoods=( bernoulli gaussian gaussian gaussian )

# Define view names
views=( Mutation Drugs mRNA Methylation )

# Define file with extra covariates
# covariatesFile="/tmp/covariates.txt"

# Define schedule of updates
schedule=( Y SW Z AlphaW Theta Tau )

# Maximum number of iterations
iter=3000

# Convergence criterion
tolerance=0.01 # Training will stop when the change in ELBO is smaller than 0.01
nostop=0 # if nostop=1 the training will complete all iterations even if the convergence criterion is met

# Define frequency for computation of ELBO
elbofreq=1

# Define initial number of latent factors
factors=20

# Define how to drop inactive latent factors
startDrop=1  # initial iteration to start shutting down factors
freqDrop=1 	 # frequency of checking for shutting down factors 
dropR2=0.00  # threshold on variance explained: dropR2=0.03 will shut down a factor if it explains less than 3% of variance (R2) in all views

# Define hyperparameters for feature-wise spike and slab
learnTheta=( 1 1 1 1 ) 	# Each element is a view, 1=active , 0=inactive
initTheta=( 1 1 1 1 ) 	# Initialisation for sparsity levels (1=no sparse)
startSparsity=250 		# Initial iteration to activate the spike and slab

# Learn the feature-wise means? We recommend always using this
learnMean=1

# Use a random seed?
seed=0 # If 0, the seed is automatically generated using the current time

scriptdir="/Users/ricard/MOFA/MOFA/run"

cmd='python $scriptdir/template_run.py
	--delimiter "$delimiter"
	--inFiles ${inFiles[@]}
	--outFile $outFile
	--likelihoods ${likelihoods[@]}
	--views ${views[@]}
	--schedule ${schedule[@]}
	--iter $iter
	--tolerance $tolerance
	--learnTheta ${learnTheta[@]}
	--initTheta ${initTheta[@]}
	--startSparsity ${startSparsity[@]}
	--elbofreq $elbofreq
	--factors $factors
	--startDrop $startDrop
	--freqDrop $freqDrop
	--dropR2 $dropR2
	--seed $seed
'

if [[ $header_rows -eq 1 ]]; then cmd="$cmd --header_rows"; fi
if [[ $header_cols -eq 1 ]]; then cmd="$cmd --header_cols"; fi
if [ -n "$covariatesFile" ]; then cmd="$cmd --covariatesFile $covariatesFile"; fi
if [[ $center_features -eq 1 ]]; then cmd="$cmd --center_features"; fi
if [[ $scale_features -eq 1 ]]; then cmd="$cmd --scale_features"; fi
if [[ $scale_views -eq 1 ]]; then cmd="$cmd --scale_views"; fi
if [[ $nostop -eq 1 ]]; then cmd="$cmd --nostop"; fi
if [[ $learnMean -eq 1 ]]; then cmd="$cmd --learnMean"; fi

eval $cmd

