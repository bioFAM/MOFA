
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

    def addMarkovBlanket(self, **kwargs):
    	# Method to define the Markov blanket of the node 
        self.markov_blanket = kwargs

    def update(self):
        # General method to update both parameters and expectations
        self.updateParameters()
        self.updateExpectations()

    def updateExpectations(self):
        # General method to update the expectations of a node
        pass

    def getDimensions(self):
        return self.dim
        
    def getExpectations(self):
        # General method to get the expectations of a node
        pass

    def getExpectation(self):
        # General method to get the first moment (expectation) of a node
        pass

    def updateParameters(self):
        # General function to update parameters of the Q distribution
        pass

    def getParameters(self):
        # General function to get the parameters of the distributions
        pass

    def updateDim(self, axis, new_dim):
        # Method to update the dimensionality of a node
        # this seems inefficient but self.dim is a tuple and it cannot be modified using
        # something like self.dim[axis] = new_dim
        dim = list(self.dim)
        dim[axis] = new_dim
        self.dim = tuple(dim)

class Constant_Node(Node):
    """
    General class for a constant node in a Bayesian network
    Constant nodes do not have expectations or parameters but just values.
    However, for technical reasons Expectations are defined to be the same as the values
    """
    def __init__(self, dim, value):
        print dim
        print value.shape
        self.dim = dim
        if isinstance(value,(int,float)):
            self.value = value * s.ones(dim)
        else:
            assert value.shape == dim, "dimensionality mismatch"
            self.value = value

    def getValue(self):
        return self.value 

    # def getParameters(self):
        # return self.value

    def getExpectation(self):
        return self.getValue()
        
    def getExpectations(self):
        return { 'E':self.getValue() }

    def removeFactors(self, idx, axis=None):
        if hasattr(self,"factors_axis"): axis = self.factors_axis
        if axis is not None:
            self.value = s.delete(self.value, idx, axis)
            self.updateDim(axis=axis, new_dim=self.dim[axis]-len(idx))

