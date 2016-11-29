
"""
Module to define general nodes in a Bayesian network

All nodes (variational or not variational) have two main attributes:
- dim: dimensionality of the node
- markov_blanket: the markov blanket of the node

"""
import scipy as s

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

    def update(self):
        pass

    def removeFactors(self, *idx):
        # General function to remove factors
        pass

    def updateParameters(self):
        pass

    def updateExpectations(self):
        pass

    def getParameters(self):
        pass

    def getExpectations(self):
        pass


class Constant_Node(Node):
    """
    """
    def __init__(self, dim, value):
        self.dim = dim
        if isinstance(value,(int,float)):
            self.value = value * s.ones(dim)
        else:
            assert value.shape == dim, "dimensionality mismatch"
            self.value = value

    def getValue(self):
        return self.value 

    def getExpectation(self):
        return self.getValue()
        
    def getExpectations(self):
        return { 'E':self.getValue() }
