#!/bin/bash

# Basic script to run MOFA. For a more advanced template with more options, see run_advanced.sh

###################
## START EDITING ##
###################

# Input files as plain text format
inFolder="/Users/ricard/MOFA/MOFA/test/data"
inFiles=( "$inFolder/500_0.txt" "$inFolder/500_1.txt" "$inFolder/500_2.txt")

# Options for the input files
delimiter=" " # Delimiter, such as "\t", "" or " "
header_rows=0 # Set to 1 if the files contain row names
header_cols=0 # Set to 0 if the files contain column names

# Output path, the model saved as .hdf5 format
outFile=( "/Users/ricard/MOFA/MOFA/test/test.hdf5" )

# Define likelihoods ('gaussian' for continuous data, 'bernoulli' for binary data or 'poisson' for count data)
likelihoods=( gaussian gaussian gaussian )

# Define view names
views=( A B C )

# Define the initial number of latent factors and how they are dropped during training.
# The model automatically removes inactive factors while training if they explain no variance.
# 'dropR2' sets the threshold of variance explained (R^2) for a model to shut down a factor. 
# For example, if dropR2=0.03the model will shut down factors if it explains less than 5% of variance in all views
# If you only want to get the most strong drivers of variation (generally less than 5) then we recommend dropR2 to be at least 5%,
# but if you want to capture more subtle sources of variation you should decrease it to 1 or 3%
factors=10
dropR2=0.05

####################
## FINISH EDITING ##
####################

# Prepare command
cmd='mofa
	--delimiter "$delimiter"
	--inFiles ${inFiles[@]}
	--outFile $outFile
	--likelihoods ${likelihoods[@]}
	--views ${views[@]}
	--factors $factors
	--dropR2 $dropR2
	--learnMean
'
if [[ $header_rows -eq 1 ]]; then cmd="$cmd --header_rows"; fi
if [[ $header_cols -eq 1 ]]; then cmd="$cmd --header_cols"; fi

# Run!
eval $cmd

