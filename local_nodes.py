from __future__ import division
import scipy as s
import scipy.stats as stats
from numpy.linalg import linalg

from nodes import *
from utils import *

"""
"""

#####################################
## General classes for local nodes ##
#####################################

class Local_Node(Node):
    """ 
    Abstract class for local variational parameters in a Bayesian probabilistic model.
    Variational Parameters do not have prior and posterior distributions, but they are updated
    in order to optimise a lower bound.

    Two main attributes:
    - dim: dimensionality of the node
    - value: current estimate of the node
    """
    def __init__(self, dim, initial_value):
        Node.__init__(self,dim)

        # Initialise the values
        if initial_value is None: 
            initial_value = s.zeros(self.dim)
        else: 
            assert initial_value.shape == dim
        self.value = initial_value

    def update(self):
        # General function to update the values of the local node
        pass    
    def getValue(self):
        # General function to return the values of the local node
        return self.value
    def getExpectation(self):
        return self.getValue()
    def getExpectations(self):
        return { 'E':self.getValue(), 'lnE':None, 'E2':None }

# class Multiview_Local_Node(Multiview_Node,Local_Node):
#     """ 
#     General class for multiview variational nodes
#     """
#     def __init__(self, M, *nodes):
#         # dim: list of M tuples with the dimensionality of the node in each view, ex. [ (10,5), (20,5), ...]
#         # nodes: list of M 'Node' instances
#         Multiview_Node.__init__(self, M, *nodes)

#     def update(self):
#         # Update parameters using the current expectations
#         for m in xrange(self.M): self.nodes[m].update()
#         pass
#     def getValue(self):
#         return [ self.nodes[m].getValue() for m in xrange(self.M) ]
#     def getExpectation(self):
#         return self.getValue()

###########################################################
## General class for observed and unobserved local nodes ##
###########################################################

class Observed_Local_Node(Local_Node):
    """ 
    General class for an observed node in a Bayesian network
    """
    def __init__(self, dim, value):
        Local_Node.__init__(self, dim, value)

