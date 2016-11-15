
"""
Module to define general nodes in a Bayesian network

All nodes (variational or not variational) have two main attributes:
- dim: dimensionality of the node
- markov_blanket: the markov blanket of the node

To-do: 
-  removeFactors should not be here... make it more general to remove any element from any dimensions
"""


class Node(object):
    """ 
    General class for a node in a Bayesian network
    """
    def __init__(self, dim):
    	self.dim = dim
    	pass
    def addMarkovBlanket(self, **kwargs):
    	# Function to define the Markov blanket of the node 
        self.markov_blanket = kwargs

    def removeFactors(self,*idx):
        # General function to remove factors
        pass