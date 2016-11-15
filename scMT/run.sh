#!/bin/bash

met_inputfiles=()
expr_inputfiles=("/Users/ricard/data/scMT/expr/processed/filt/tmp/e_matrix.txt")
outdir="/Users/ricard/git/scGFA/scMT/tmp"
view_names=("expr")
k=10
# n_trials=1
iter=10
verbosity=2
# cores=1
sparse=False
lik=("gaussian")
dropK=True
dropK_threshold=0.01
forceiter=True
savefreq=$((iter+1))
savefolder=/tmp/scGFA
elbofreq=1

# echo "python main.py  -m $met_inputfiles -e $expr_inputfiles -o $outdir -k $k -s $view_names -n $n_trials -i $iter -v $verbosity -c $cores"
python main.py  -m $met_inputfiles -e $expr_inputfiles -o $outdir -k $k -s $view_names \
	-i $iter -v $verbosity --sparse $sparse --likelihood $lik \
	-dropK $dropK --dropK_threshold $dropK_threshold --forceiter $forceiter --savefreq $savefreq \
	-elbofreq $elbofreq -savefolder $savefolder
