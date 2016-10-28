
"""
Module to define general nodes in a Bayesian network

All nodes (variational or not variational) have two main attributes:
- dim: dimensionality of the node
- markov_blanket: the markov blanket of the node
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

# class Observed_Node(Node):
#     """ 
#     General class for an observed node in a Bayesian network
#     """
#     def __init__(self, dim, obs):
#     	Node.__init__(self,dim)
#     	self.obs = obs