from __future__ import division
import numpy as np
import os
import h5py

"""
Module to define some useful util functions
"""

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


def saveParameters(model, hdf5, view_names=None):

	# Get nodes from the model
	nodes = model.getNodes()

	# Create groups 
	param_grp = hdf5.create_group("parameters")

	# Iterate over nodes
	for node in nodes:

		# Collect node expectation
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
				# Loop through the parameters
				if parameters[m] is not None:
					for param_name in parameters[m].keys():
						view_subgrp.create_dataset(param_name, data=parameters[m][param_name].T)
		# Single-view nodes
		else:
			for param_name in parameters.keys():
				node_subgrp.create_dataset("%s" % (param_name), data=parameters[param_name].T)
	pass

def saveExpectations(model, hdf5, view_names=None):

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
				if expectations[m] is not None:
					for exp_name in expectations[m].keys():
						view_subgrp.create_dataset(exp_name, data=expectations[m][exp_name].T)

		# Single-view nodes
		else:

			for exp_name in expectations.keys():
				node_subgrp.create_dataset("%s" % (exp_name), data=expectations[exp_name].T)

def saveTrainingStats(model, hdf5):
	stats = model.getTrainingStats()
	stats_grp = hdf5.create_group("training_stats")
	stats_grp.create_dataset("activeK", data=stats["activeK"])
	stats_grp.create_dataset("elbo", data=stats["elbo"])
	stats_grp.create_dataset("elbo_terms", data=stats["elbo_terms"].T)
	stats_grp['elbo_terms'].attrs['colnames'] = list(stats["elbo_terms"].columns.values)

def saveTrainingOpts(model, hdf5):
	opts = model.getTrainingOpts()

	# Remove dictionaries from the options
	opts.pop('dropK')
	hdf5.create_dataset("training_opts", data=opts.values())
	hdf5['training_opts'].attrs['names'] = opts.keys()

# def saveTrainingData(model, hdf5, view_names=None, sample_names=None, feature_names=None):
# 	data = model.getTrainingData()
# 	data_grp = hdf5.create_group("training_data")
# 	for m in xrange(len(data)):
# 		view = view_names[m] if view_names is not None else str(m)
# 		data_grp.create_dataset(view, data=data[m].data.T)
# 		if feature_names is not None: 
# 			data_grp[view].attrs['colnames'] = feature_names[m]
# 		if sample_names is not None: 
# 			data_grp[view].attrs['rownames'] = sample_names

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


def saveModel(model, outfile, view_names=None, sample_names=None, feature_names=None):
	assert model.trained == True, "Model is not trained yet"

	hdf5 = h5py.File(outfile,'w')
	saveExpectations(model,hdf5,view_names)
	saveParameters(model,hdf5,view_names)
	saveTrainingStats(model,hdf5)
	saveTrainingOpts(model,hdf5)
	saveTrainingData(model, hdf5, view_names, sample_names, feature_names)
	hdf5.close()
