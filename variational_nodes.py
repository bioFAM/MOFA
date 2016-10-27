
from __future__ import division
import scipy as s

from nodes import *
from distributions import *

"""
This module is used to define Variational Nodes 

TO-DO
- Are we doing well the double parent class? Maybe use super class
"""


###########################################
## General classes for variational nodes ##
###########################################

class Variational_Node(Node):
    """ 
    Abstract class for a multi-view variational node in a Bayesian probabilistic model.
    Variational nodes can be observed or unobserved
    """
    def __init__(self, dim):
        Node.__init__(self,dim)
        self.P = None
        self.Q = None
    def calculateELBO(self):
        # General function to calculate the ELBO of the node
        return 0.
    def updateExpectations(self):
        # General function to update expectations of the Q distribution
        pass
    def getExpectation(self):
        # General function to get the expectated value of the Q distribution
        return self.Q.E
    def getExpectations(self):
        # General function to get all relevant moments
        pass
    def updateParameters(self):
        # General function to update parameters of the Q distribution
        pass
    def update(self):
        # General function to update both parameters and expectations
        self.updateParameters()
        self.updateExpectations()
    def getParameters(self):
        # General function to get the parameters of the distributions
        pass
    def removeFactors(self,*idx):
        # General function to remove factors
        pass

###################################################################
## General classes for observed and unobserved variational nodes ##
###################################################################

class Observed_Variational_Node(Variational_Node):
    """ 
    Abstract class for an observed variational node in a Bayesian probabilistic model.
    Observed variational nodes only contain a single distribution Q(X) which will be stored as 
    instances of Distribution() as a .Q attribute.
    """
    def __init__(self, dim, obs):
        Variational_Node.__init__(self, dim)
        self.obs = obs
    def getObservations(self):
        return self.obs
    def getExpectation(self):
        return self.getObservations()
class Unobserved_Variational_Node(Variational_Node):
    """ 
    Abstract class for an unobserved variational node in a Bayesian probabilistic model.
    Unobserved variational nodes contain a prior P(X) and a variational Q(X) distribution, 
    which will be stored as instances of Distribution() attributes .P and .Q, respectively.
    The distributions are in turn composed of parameters and expectations
    """
    def __init__(self, dim):
        Variational_Node.__init__(self, dim)
        self.P = None
        self.Q = None
    def updateExpectations(self):
        self.Q.updateExpectations()

#######################################################
## Specific classes for unobserved variational nodes ##
#######################################################

class UnivariateGaussian_Unobserved_Variational_Node(Unobserved_Variational_Node):
    """ 
    Abstract class for a variational node where P(x) and Q(x)
    are both univariate Gaussian distributions.

    To save memory, we follow common practice and the the prior distribution is 
    assumed to be the same for all elements
    """
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None, qE2=None):
	    # dim (2d tuple): dimensionality of the node 
	    # pmean (nd array): the mean parameter of the P distribution
	    # qmean (nd array): the mean parameter of the Q distribution
	    # pvar (nd array): the variance parameter of the P distribution
	    # qvar (nd array): the variance parameter of the Q distribution
	    # qE (nd array): the initial first moment of the Q distribution
	    # qE2 (nd array): the initial second moment of the Q distribution
        Unobserved_Variational_Node.__init__(self, dim)

        # Initialise the P and Q distributions
        self.P = UnivariateGaussian(dim=(1,), mean=pmean, var=pvar)
        self.Q = UnivariateGaussian(dim=dim, mean=qmean, var=qvar, E=qE, E2=qE2)

    def getParameters(self):
        return dict({'mean': self.Q.mean, 'var': self.Q.var})
    def getExpectations(self):
        return dict({'E':self.Q.E, 'E2':self.Q.E2, 'lnE':None})
class MultivariateGaussian_Unobserved_Variational_Node(Unobserved_Variational_Node):
    """ 
    Abstract class for a variational node where P(x) and Q(x)
    are both multivariate Gaussian distributions.

    Currently the prior of this distribution is not used anywhere so it is ignored to save memory
    """
    def __init__(self, dim, qmean, qcov, qE=None, qE2=None):
        # dim (2d tuple): dimensionality of the node 
        # qmean (nd array): the mean parameter of the Q distribution
        # qcov (nd array): the covariance parameter of the Q distribution
        # qE (nd array): the initial first moment of the Q distribution
        # qE2 (nd array): the initial second moment of the Q distribution
        Unobserved_Variational_Node.__init__(self, dim)

        # Initialise the P and Q distributions
        # self.P = MultivariateGaussian(dim=dim, mean=pmean, cov=pcov)
        self.Q = MultivariateGaussian(dim=dim, mean=qmean, cov=qcov, E=qE, E2=qE2)

	def getParameters(self):
		return dict({'mean':self.Q.mean, 'cov':self.Q.cov})
    def getExpectations(self):
        return dict({'E':self.Q.E, 'E2':self.Q.E2, 'lnE':None})
class Gamma_Unobserved_Variational_Node(Unobserved_Variational_Node):
    """ 
    Abstract class for a variational node where P(x) and Q(x) are both gamma distributions
    """
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
	    # dim (2d tuple): dimensionality of the node 
	    # pa (nd array): the 'a' parameter of the P distribution
	    # qa (nd array): the 'a' parameter of the Q distribution
	    # pa (nd array): the 'b' parameter of the P distribution
	    # pb (nd array): the 'b' parameter of the Q distribution
	    # qE (nd array): the initial expectation of the Q distribution
        Unobserved_Variational_Node.__init__(self,dim)

        # Initialise the distributions
        self.P = Gamma(dim=(1,), a=pa, b=pb) 
        self.Q = Gamma(dim=dim, a=qa, b=qb, E=qE)

    def getParameters(self):
        return dict({'a':self.Q.a, 'b':self.Q.b})
    def getExpectations(self):
        return dict({'E':self.Q.E, 'lnE':self.Q.lnE, 'E2':None})

class Bernoulli_Unobserved_Variational_Node(Unobserved_Variational_Node):
    """ 
    Abstract class for a variational node where P(x) and Q(x)
    are both bernoulli distributions.

    To save memory, we follow common practice and the the prior distribution is 
    assumed to be the same for all elements
    """
    def __init__(self, dim, ptheta, qtheta, qE=None):
	    # dim (2d tuple): dimensionality of the node 
	    # ptheta (nd array): the 'theta' parameter of the P distribution
	    # qtheta (nd array): the 'theta' parameter of the Q distribution
	    # qE (nd array): the current moment of the Q distribution
        Unobserved_Variational_Node.__init__(self,dim)
       
        # Initialise the distributions
        self.P = Bernoulli(dim=(1,), theta=ptheta)
        self.Q = Bernoulli(dim=dim, theta=qtheta, E=qE)

    def getParameters(self):
        return dict({'theta':self.Q.theta})
    def getExpectation(self):
        return self.Q.E
    def getExpectations(self):
        return dict({'E':self.Q.E, 'E2':None, 'lnE':None})
class BernoulliGaussian_Unobserved_Variational_Node(Unobserved_Variational_Node):
    """ 
    Abstract class for a variational node where P(x) and Q(x)
    are joint gaussian-bernoulli distributions (see paper  Spike and Slab Variational Inference for 
    Multi-Task and Multiple Kernel Learning by Titsias and Gredilla)

    This distribution has several expectations that are required for the variational updates, and they depend
    on other parameters (alpha from the ARD prior). For this reason I decided to define the expectations in the 
    class of the corresponding node, instead of doing it here.

    The only element from the prior distribution that is used is S_ptheta, for this reason I ignore the other elements
    """
    def __init__(self, dim, qmean, qvar, ptheta, qtheta):
	    # dim (2d tuple): dimensionality of the node 
	    # qmean (nd array): the mean parameter of the Q distribution
	    # qvar (nd array): the var parameter of the Q distribution
	    # ptheta (nd array): the theta parameter of the P distribution
	    # qtheta (nd array): the theta parameter of the Q distribution
        Unobserved_Variational_Node.__init__(self,dim)

        # Initialise the distributions
        # self.P = BernoulliGaussian(dim=(1,), theta=S_ptheta, mean=W_pmean, var=W_pvar) 
        self.P_theta = ptheta
        self.Q = BernoulliGaussian(dim=dim, theta=qtheta, mean=qmean, var=qvar)

    def getParameters(self):
        return dict({'theta':self.Q.theta, 'mean':self.Q.mean, 'var':self.Q.var})
    def getExpectation(self):
        return self.Q.ESW
