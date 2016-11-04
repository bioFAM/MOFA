"""
Script to run multiple trials of scGFA
"""


# Import required modules
import argparse
import os
# import pandas as pd


n_trials = 1
cores = 1

met_inputfiles = ().join(" ")
expr_inputfiles = ("/Users/ricard/data/scMT/expr/processed/filt/tmp/e_matrix.txt").join(" ")
outdir = "/Users/ricard/git/scGFA/scMT/tmp"
view_names = ("expr")
k = 10
iterations = 100
verbosity = 2
sparse = False
lik = ("gaussian").join(" ")
dropK = True
dropK_threshold = 0.01
forceiter = True
savefreq = iterations+1
savefolder = "/tmp/scGFA"
elbofreq = 1

for n in xrange(n_trials):
	expr_inputfiles
	command = "python main.py -m %s -e %s -o %s -k %i -s %s -i %i -v %i --sparse %r -l %s -dropK %r --dropK_threshold %0.02f \
	--forceiter %r --savefreq %i --elbofreq %i --savefolder %s" % (met_inputfiles, expr_inputfiles, outdir, k \
		view_names, i, verbosity, sparse, lik, dropK, dropK_threshold, forceiter, savefreq, elbofreq, savefolder)
	os.system(command)