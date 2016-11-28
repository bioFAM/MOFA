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
	return 1/(1+np.exp(-X))

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

	pass

def saveTrainingStats(model, hdf5):
	stats = model.getTrainingStats()
	stats_grp = hdf5.create_group("training_stats")
	stats_grp.create_dataset("activeK", data=stats["activeK"])
	stats_grp.create_dataset("elbo", data=stats["elbo"])
	stats_grp.create_dataset("elbo_terms", data=stats["elbo_terms"].T)
	stats_grp['elbo_terms'].attrs['colnames'] = list(stats["elbo_terms"].columns.values)
	pass

def saveTrainingOpts(model, hdf5):
	opts = model.getTrainingOpts()
	hdf5.create_dataset("training_opts", data=opts.values())
	hdf5['training_opts'].attrs['names'] = opts.keys()
	pass


def saveTrainingData(model, hdf5, view_names=None, sample_names=None, feature_names=None):
	data = model.getTrainingData()
	data_grp = hdf5.create_group("training_data")
	for m in xrange(len(data)):
		view = view_names[m] if view_names is not None else str(m)
		data_grp.create_dataset(view, data=data[m].data.T)
		if feature_names is not None: 
			data_grp[view].attrs['colnames'] = feature_names[m]
		if sample_names is not None: 
			data_grp[view].attrs['rownames'] = sample_names
	pass

def saveModel(model, outfile, view_names=None, sample_names=None, feature_names=None):
	assert model.trained == True, "Model is not trained yet"

	hdf5 = h5py.File(outfile,'w')
	saveExpectations(model,hdf5,view_names)
	saveParameters(model,hdf5,view_names)
	saveTrainingStats(model,hdf5)
	saveTrainingOpts(model,hdf5)
	saveTrainingData(model, hdf5, view_names, sample_names, feature_names)
	hdf5.close()
	pass

# def saveTrainingOpts(model, outdir):
#     opts = model.getTrainingOpts()
#     file = os.path.join(outdir,"opts.txt")
#     with open(file, "a") as f:
#         for k,v in opts.iteritems(): 
#             f.write(k + ":" + str(v) + "\n")
#     pass

# def saveModel(model, outdir, compress=False):
# 	# Function to save a trained model to be load in R:
# 	# 	Expectations and parameters are stored as .npy objects to be loaded in R using the RcppCNPy package 
# 	#	
# 	# Inputs: 
# 	# - model (BayesNet class): the trained model
# 	# - outdir (string): output directory
# 	# - compress (bool): compress files using gzip?

# 	# Check that the model is trained
# 	assert model.trained == True, "Model is not trained yet"
# 	nodes = model.getNodes()
# 	# vb_nodes = model.getVariationalNodes()

# 	# Create output folder if it does not exist
# 	if not os.path.exists(outdir): os.makedirs(outdir)

# 	#####################
# 	## Save parameters ##
# 	#####################

# 	# to do...


# 	#######################
# 	## Save expectations ##
# 	#######################

# 	# Iterate over nodes
# 	for node in nodes:
# 		# The data will be saved separately, not here...
# 		if node == "Y": continue

# 		# Collect node expectations
# 		expectations = nodes[node].getExpectations()

# 		# Multi-view nodes
# 		if type(expectations) == list:
# 			# Iterate over views
# 			for m in xrange(len(expectations)):
# 				# Iterate over expectations
# 				for key,value in expectations[m].iteritems():
# 					filename = os.path.join(outdir,"%s_%s_%d.npy" % (node,key,m+1))
# 					# print "\tsaving %s..." % filename
# 					np.save(filename,value)

# 		# Single-view nodes
# 		else:
# 			for key,value in expectations.iteritems():
# 				filename = os.path.join(outdir,"%s_%s.npy" % (node,key))
# 				# print "\tsaving %s..." % filename
# 				np.save(filename,value)

# 	##############
# 	## Compress ##
# 	##############

# 	if compress:
# 		os.system("gzip -f %s/*.npy" % outdir)


# if __name__ == "__main__":
# 	# a = np.arange(20).reshape(5,4)
# 	a = np.arange(20).reshape(5,4).astype(int)
# 	# a.tofile("/tmp/data.bin", sep="", format="%s")
# 	f = open("/tmp/data.bin","wb")
# 	f.write(a)
# 	f.close()
