import numpy as np
import struct
import os

"""
Module to define some useful util functions
"""

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

def saveTrainingStats(model, outdir):
    stats = model.getTrainingStats()
    np.savetxt(X=stats["activeK"], fmt='%d', fname=os.path.join(outdir,"activeK.txt"), delimiter=' ')
    np.savetxt(X=stats["elbo"], fmt='%0.02f', fname=os.path.join(outdir,"elbo.txt"), delimiter=' ')
    stats["elbo_terms"].to_csv(os.path.join(outdir,"elbo_terms.txt"), sep="\t", header=True, index=False)

def saveTrainingOpts(model, outdir):
    opts = model.getTrainingOpts()
    file = os.path.join(outdir,"opts.txt")
    with open(file, "a") as f:
        for k,v in opts.iteritems(): 
            f.write(k + ":" + str(v) + "\n")
    pass

def saveModel(model, outdir, compress=False, only_first_moments=True):
	# Function to save a trained model to be load in R:
	# 	Expectations and parameters are stored as .npy objects to be loaded in R using the RcppCNPy package 
	#	
	# Inputs: 
	# - model (BayesNet class): the trained model
	# - outdir (string): output directory
	# - compress (bool): compress files using gzip?

	# Check that the model is trained
	assert model.trained == True, "Model is not trained yet"
	nodes = model.getNodes()

	# Create output folder if it does not exist
	if not os.path.exists(outdir): os.makedirs(outdir)

	#######################
	## Save expectations ##
	#######################

	# Iterate over nodes
	for node in nodes:
		# The data will be saved separately, not here...
		if node == "Y": continue

		# Collect node expectations
		if only_first_moments:
			expectations = {'E': nodes[node].getExpectation() }
			if node == "Zeta":
				print expectations
				exit()
		else:
			expectations = nodes[node].getExpectations()

		# Multi-view nodes
		if type(expectations) == list:
			# Iterate over views
			for m in xrange(len(expectations)):
				# Iterate over expectations
				for key,value in expectations[m].iteritems():
					filename = os.path.join(outdir,"%s_%s_%d.npy" % (node,key,m+1))
					print "\tsaving %s..." % filename
					np.save(filename,value)

		# Single-view nodes
		else:
			for key,value in expectations.iteritems():
				filename = os.path.join(outdir,"%s_%s.npy" % (node,key))
				print "\tsaving %s..." % filename
				np.save(filename,value)

	if compress:
		os.system("gzip -f %s/*.npy" % outdir)


# if __name__ == "__main__":
# 	# a = np.arange(20).reshape(5,4)
# 	a = np.arange(20).reshape(5,4).astype(int)
# 	# a.tofile("/tmp/data.bin", sep="", format="%s")
# 	f = open("/tmp/data.bin","wb")
# 	f.write(a)
# 	f.close()
