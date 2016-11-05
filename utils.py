import numpy as np
import struct
import os

"""
Module to define some useful util functions
"""

def logdet(X):
	return np.log(np.linalg.det(X))
	# UC = np.linalg.cholesky(X)
	# return 2*sum(np.log(s.diag(UC)))

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

def saveModel(model, outdir, compress=False):
	# Function to save a trained model to be load in R:
	# 	Expectations and parameters are stored as .npy objects to be loaded in R using the RcppCNPy package 
	# 	Training statistics (lower bound, active factors) are stored as npy files
	#	
	# Inputs: 
	# - model (BayesNet class): the trained model
	# - outdir (string): output directory
	# - compress (bool): compress files using gzip?

	# Check that the model is trained
	assert model.trained == True, "Model is not trained yet"
	nodes = model.getAllNodes()

	# Create output folder if it does not exist
	if not os.path.exists(outdir): os.makedirs(outdir)

	#######################
	## Save expectations ##
	#######################

	# Iterate over nodes
	for node in nodes:
		expectations = nodes[node].getExpectations()

		# Multi-view nodes
		if type(expectations) == list:
			# Iterate over views
			for m in xrange(len(expectations)):
				# Iterate over expectations
				for key,value in expectations[m].iteritems():
					filename = os.path.join(outdir,"%s_%s_%d.npy" % (node,key,m+1))
					print "Saving %s..." % filename
					np.save(filename,value)

		# Single-view nodes
		else:
			for key,value in expectations.iteritems():
				filename = os.path.join(outdir,"%s_%s.npy" % (node,key))
				print "Saving %s..." % filename
				np.save(filename,value)

	if compress:
		os.system("gzip %s/*.npy" % outdir)

	##############################
	## Save training statistics ##
	##############################

	stats = model.getTrainingStats()

	np.save(file=os.path.join(outdir,"activeK.npy"), arr=stats["activeK"])
	np.save(file=os.path.join(outdir,"elbo.npy"), arr=stats["elbo"])
	stats["elbo_terms"].to_csv(os.path.join(outdir,"elbo_terms.txt"), sep="\t", header=True, index=False)



# if __name__ == "__main__":
# 	# a = np.arange(20).reshape(5,4)
# 	a = np.arange(20).reshape(5,4).astype(int)
# 	# a.tofile("/tmp/data.bin", sep="", format="%s")
# 	f = open("/tmp/data.bin","wb")
# 	f.write(a)
# 	f.close()
