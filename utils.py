import numpy as np
import struct
import os

# from singleview_nodes import VariationalNode

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

def save_npy(outdir, outprefix, data):
    # Function to save parameter and expectations in .npy format and gzip compressed
    # so that they can be loaded with R using the RcppCNPy package 
    # Inputs: 
    # - outdir (string): output directory
    # - outprefix (string): output prefix
    # - node (instance of Variational_Node): currently only works with unobserved ones

    # assert isinstance(node,VariationalNode)
    # assert isinstance(node,Unobserved_Node), "'node' has to be an instance of Variational_Node"

    # Create output folder if it does not exist
	if not os.path.exists(outdir): os.makedirs(outdir)


	outfile = os.path.join(outdir,outprefix)

	if type(data) == list:
		for m in xrange(len(data)):
			filename = "%s_%d.npy" % (outfile,m+1)
			print "Saving %s..." % filename
			np.save(filename,data[m])
	else:
		filename = "%s.npy" % (outfile)
		print "Saving %s..." % filename
		np.save(filename,data)

	# Compress the files
	# os.system("gzip %s/*" % outdir)

	pass


# if __name__ == "__main__":
# 	# a = np.arange(20).reshape(5,4)
# 	a = np.arange(20).reshape(5,4).astype(int)
# 	# a.tofile("/tmp/data.bin", sep="", format="%s")
# 	f = open("/tmp/data.bin","wb")
# 	f.write(a)
# 	f.close()
