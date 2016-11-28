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

    Note: we are defining methods to obtain expectations just for convenience during training.
    However, local nodes do not have expectations.

    Two main attributes:
    - dim: dimensionality of the node
    - value: current estimate of the node
    """
    def __init__(self, dim, value):
        Node.__init__(self,dim)
        self.value = None

        # Initialise the values
        if value is None: 
            value = s.zeros(self.dim)
        else: 
            assert value.shape == dim
        self.value = value

    def getValue(self):
        # General function to return the values of the local node
        return self.value
    def getExpectation(self):
        return self.getValue()
    def getParameters(self):
        return None
    def getExpectations(self):
        # return { 'E':self.getValue(), 'lnE':None, 'E2':None }
        return { 'E':self.getValue() }
    def update(self):
        # General function to update the values of the local node
        pass 

class Unobserved_Local_Node(Local_Node):
    """ 
    Abstract class for local variational parameters in a Bayesian probabilistic model.
    Variational Parameters do not have prior and posterior distributions, but they are updated
    in order to optimise a lower bound.

    Note: we are defining methods to obtain expectations just for convenience during training.
    However, local nodes do not have expectations.

    Two main attributes:
    - dim: dimensionality of the node
    - value: current estimate of the node
    """
    def __init__(self, dim, initial_value):
        Local_Node.__init__(self, dim, initial_value)
   

###########################################################
## General class for observed and unobserved local nodes ##
###########################################################

class Observed_Local_Node(Local_Node):
    """ 
    General class for an observed node in a Bayesian network
    """
    def __init__(self, dim, value):
        Local_Node.__init__(self, dim, value)
    def getExpectations(self):
        # return { 'obs':self.getValue(), 'lnE':None, 'E2':None }
        # return { 'obs':self.getValue() }
        return { 'E':self.getValue() } 
