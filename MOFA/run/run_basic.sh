#!/bin/bash

# Input files as plain text format
inFolder="/Users/ricard/MOFA/MOFA/test/data"
inFiles=( "$inFolder/500_0.txt" "$inFolder/500_1.txt" "$inFolder/500_2.txt")
delimiter=" " # Delimiter 
header_rows=0 # Do the files contain row names?
header_cols=0 # Do the files contain column names?

# Output file as .hdf5 format
outFile=( "/Users/ricard/MOFA/MOFA/test/test.hdf5" )

# Define likelihoods 
likelihoods=( gaussian gaussian gaussian )

# Define view names
views=( A B C )

# Define initial number of latent factors
factors=10

# Define threshold to shut down factors based on the variance explained, 
# dropR2=0.03 will shut down a factor if it explains less than 3% of variance (R2) in all views
dropR2=0.03

# indicate the full path of the 'run' directory inside MOFA package
scriptdir="/Users/ricard/MOFA/MOFA/run"

# Prepare python command with arguments
cmd='python $scriptdir/template_run.py
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

