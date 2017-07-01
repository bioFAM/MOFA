from __future__ import division
import numpy as np
import pandas as pd
import numpy.ma as ma
import os
import h5py
from glob import glob

"""
Module to define some useful util functions
"""

# Function to load prior annotations for the weight sparsity level (given by the variable theta)
def loadTheta(data_opts):
	d = data_opts['ThetaDir']
	file_names = glob.glob(d+'/*')

	M = len(data_opts['view_names'])
	ThetaPrior = [None] * M

	for f in file_names:
		annotations = pd.read_csv(f, sep=data_opts['delimiter'], header=0, index_col=0)
		view_name = f.split('/')[-1].split('.')[0]
		tmp = [m for m in xrange(M) if view_name == data_opts['view_names'][m]]
		try:
			i = tmp[0]
			ThetaPrior[i] = annotations.values
		except:
			print 'View names dont match with annotation file names'
			print 'ignoring unrecognised file names'

	return ThetaPrior

# Function to mask the data, mainly to test missing values and imputation
def maskData(data, data_opts):
    print "Masking data with the following options:"
    print "at random:"
    print data_opts['maskAtRandom']
    print "full cases:"
    print data_opts['maskNSamples']

    for m in xrange(len(data)):

        # Mask values at random
        D = data[m].shape[1]
        N = data[m].shape[0]
        p2Mask = data_opts['maskAtRandom'][m]
        if p2Mask != 0:
            idxMask = np.zeros(N*D)
            idxMask[:int(round(N*D*p2Mask))]  = 1
            np.random.shuffle(idxMask)
            idxMask=np.reshape(idxMask, [N, D])
            data[m] = data[m].mask(idxMask==1)

        # Mask samples in a complete view
        Nsamples2Mask = data_opts['maskNSamples'][m]
        if Nsamples2Mask != 0:
            idxMask = np.random.choice(N, size=Nsamples2Mask, replace = False)
            tmp = data[m].copy()
            tmp.ix[idxMask, :] = pd.np.nan
            data[m] = tmp

    return data

# Function to load the data
def loadData(data_opts, verbose=True):

    print "\n"
    print "#"*18
    print "## Loading data ##"
    print "#"*18
    print "\n"

    Y = list()
    for m in xrange(len(data_opts['input_files'])):
        file = data_opts['input_files'][m]

        # Read file (with row and column names)
        Y.append( pd.read_csv(file, delimiter=data_opts["delimiter"], header=data_opts["colnames"], index_col=data_opts["rownames"]) )
        print "Loaded %s with dim (%d,%d)..." % (file, Y[m].shape[0], Y[m].shape[1])

        # Center the data
        if data_opts['center'][m]:
			Y[m] = (Y[m] - Y[m].mean(axis=0))

    return Y

def dotd(A, B, out=None):
    """Diagonal of :math:`\mathrm A\mathrm B^\intercal`.
    If ``A`` is :math:`n\times p` and ``B`` is :math:`p\times n`, it is done in :math:`O(pn)`.
    Args:
        A (array_like): Left matrix.
        B (array_like): Right matrix.
        out (:class:`numpy.ndarray`, optional): copy result to.
    Returns:
        :class:`numpy.ndarray`: Resulting diagonal.
    """
    A = ma.asarray(A, float)
    B = ma.asarray(B, float)
    if A.ndim == 1 and B.ndim == 1:
        if out is None:
            return ma.dot(A, B)
        return ma.dot(A, B, out)

    if out is None:
        out = ma.empty((A.shape[0], ), float)

    out[:] = ma.sum(A * B.T, axis=1)
    return out

def nans(shape, dtype=float):
    a = np.empty(shape, dtype)
    a.fill(np.nan)
    return a

def corr(A,B):
	# Rowwise mean of input arrays & subtract from input arrays themeselves
	A_mA = A - A.mean(1)[:,None]
	B_mB = B - B.mean(1)[:,None]

	# Sum of squares across rows
	ssA = (A_mA**2).sum(1);
	ssB = (B_mB**2).sum(1);

	# Finally get corr coeff
	return np.dot(A_mA,B_mB.T)/np.sqrt(np.dot(ssA[:,None],ssB[None]))

def logdet(X):
	return np.log(np.linalg.det(X))
	# UC = np.linalg.cholesky(X)
	# return 2*sum(np.log(np.diag(UC)))

def sigmoid(X):
	return np.divide(1.,1.+np.exp(-X))
	# return 1./(1.+np.exp(-X))

def ddot(d, mtx, left=True):
	"""Multiply a full matrix by a diagonal matrix.
	This function should always be faster than dot.

	Input:
	  d -- 1D (N,) array (contains the diagonal elements)
	  mtx -- 2D (N,N) array
	  left: is the diagonal matrix on the left or on the right of the product?

	Output:
	  ddot(d, mts, left=True) == dot(diag(d), mtx)
	  ddot(d, mts, left=False) == dot(mtx, diag(d))
	"""
	if left:
		return (d*mtx.T).T
	else:
		return d*mtx

def lambdafn(X):
	return np.tanh(X/2.)/(4.*X)

def saveParameters(model, hdf5, view_names=None):

	# Get nodes from the model
	nodes = model.getNodes()

	# Create groups
	param_grp = hdf5.create_group("parameters")

	# Iterate over nodes
	for node in nodes:

		# Collect node parameters
		parameters = nodes[node].getParameters()

		# Create node subgroup
		node_subgrp = param_grp.create_group(node)

		# Multi-view nodes
		if type(parameters) == list:
			# Loop through the views
			for m in xrange(len(parameters)):
				if view_names is not None:
					tmp = view_names[m]
				else:
					tmp = "%d" % m
				# Create subsubgroup for the view
				view_subgrp = node_subgrp.create_group(tmp)
				# Loop through the parameters of the view
				if parameters[m] is not None:
					# Variational nodes
					if type(parameters[m]) == dict:
						for param_name in parameters[m].keys():
							if parameters[m][param_name] is not None:
								view_subgrp.create_dataset(param_name, data=parameters[m][param_name].T)
					# Non-variational nodes (no distributions)
					elif type(parameters[m]) == np.ndarray:
						   view_subgrp.create_dataset("value", data=parameters[m].T)

		# Single-view nodes
		else:
			for param_name in parameters.keys():
				node_subgrp.create_dataset("%s" % (param_name), data=parameters[param_name].T)
	pass

def saveExpectations(model, hdf5, view_names=None, only_first_moments=True):

	# Get nodes from the model
	nodes = model.getNodes()

	exp_grp = hdf5.create_group("expectations")

	# Iterate over nodes
	for node in nodes:

		# Collect node expectations
		expectations = nodes[node].getExpectations()

		# Create subgroup for the node
		node_subgrp = exp_grp.create_group(node)

		# Multi-view nodes
		if type(expectations) == list:

			# Iterate over views
			for m in xrange(len(expectations)):
				if view_names is not None:
					tmp = view_names[m]
				else:
					tmp = "%d" % m

				# Create subsubgroup for the view
				view_subgrp = node_subgrp.create_group(tmp)

				# Loop through the expectations
				if only_first_moments: expectations[m] = {'E':expectations[m]["E"]}
				if expectations[m] is not None:
					for exp_name in expectations[m].keys():
						view_subgrp.create_dataset(exp_name, data=expectations[m][exp_name].T)

		# Single-view nodes
		else:
			if only_first_moments: expectations = {'E':expectations["E"]}
			for exp_name in expectations.keys():
				node_subgrp.create_dataset("%s" % (exp_name), data=expectations[exp_name].T)

def saveTrainingStats(model, hdf5):
	stats = model.getTrainingStats()
	stats_grp = hdf5.create_group("training_stats")
	stats_grp.create_dataset("activeK", data=stats["activeK"])
	stats_grp.create_dataset("elbo", data=stats["elbo"])
	stats_grp.create_dataset("elbo_terms", data=stats["elbo_terms"].T)
	stats_grp['elbo_terms'].attrs['colnames'] = list(stats["elbo_terms"].columns.values)

def saveTrainingOpts(opts, hdf5):
	# Remove dictionaries from the options
	for k,v in opts.copy().iteritems():
		if type(v)==dict:
			for k1,v1 in v.iteritems():
				opts[str(k)+"_"+str(k1)] = v1
			opts.pop(k)
	hdf5.create_dataset("training_opts", data=np.array(opts.values(), dtype=np.float))
	hdf5['training_opts'].attrs['names'] = opts.keys()

def saveModelOpts(opts, hdf5):
	opts_interest = ["learnMean","schedule","likelihood"]
	opts = dict((k, opts[k]) for k in opts_interest)
	grp=hdf5.create_group('model_opts')
	for k,v in opts.items():
		grp.create_dataset(k, data=v)
	grp[k].attrs['names'] = opts.keys()

def saveTrainingData(model, hdf5, view_names=None, sample_names=None, feature_names=None):
	data = model.getTrainingData()
	data_grp = hdf5.create_group("data")
	featuredata_grp = hdf5.create_group("features")
	hdf5.create_dataset("samples", data=sample_names)
	for m in xrange(len(data)):
		view = view_names[m] if view_names is not None else str(m)
		data_grp.create_dataset(view, data=data[m].data.T)
		if feature_names is not None:
			featuredata_grp.create_dataset(view, data=feature_names[m])

def saveModel(model, outfile, train_opts, model_opts, view_names=None, sample_names=None, feature_names=None):
	assert model.trained == True, "Model is not trained yet"
	assert len(np.unique(view_names)) == len(view_names), 'View names must be unique'
	assert len(np.unique(sample_names)) == len(sample_names), 'Sample names must be unique'

	# For some reason h5py orders the datasets alphabetically, so we have to modify the likelihood accordingly
	idx = sorted(range(len(view_names)), key=lambda k: view_names[k])
	tmp = [model_opts["likelihood"][idx[m]] for m in xrange(len(model_opts["likelihood"]))]
	model_opts["likelihood"] = tmp

	hdf5 = h5py.File(outfile,'w')
	saveExpectations(model,hdf5,view_names)
	saveParameters(model,hdf5,view_names)
	saveTrainingStats(model,hdf5)
	saveTrainingOpts(train_opts,hdf5)
	saveModelOpts(model_opts,hdf5)
	saveTrainingData(model, hdf5, view_names, sample_names, feature_names)
	hdf5.close()
