from __future__ import division
import scipy as s
import numpy.ma as ma

from variational_nodes import Unobserved_Variational_Node
from local_nodes import Unobserved_Local_Node
from utils import sigmoid

"""
Module to define updates for non-conjugate matrix factorisation models using the Seeger approach
Reference: 'Fast Variational Bayesian Inference for Non-Conjugate Matrix Factorisation models' by Seeger and Bouchard (2012)
"""


# class Zeta_Node(Unobserved_Local_Node):
#     """ 
#     Abstract class for the local variational variable zeta that represents the location of the
#     Taylor expansion for the upper bound on the log likelihood.
#     It is required for the approximation that leads to closed-form VB updates in non-conjugate matrix factorisation models

#     Notice that Zeta is a local node, not a variational node!
#     """
#     def __init__(self, dim, initial_value=None):
#         # Inputs:
#         #  dim (ndarray): dimensionality of the node
#         #  initial_value (ndarray): initial value of Zeta
#         Unobserved_Local_Node.__init__(self, dim, initial_value)

#     def update(self):
#         Z = self.markov_blanket["Z"].getExpectation()
#         W = self.markov_blanket["W"].getExpectation()
#         self.value = s.dot(Z,W.T)

################
## Pseudodata ##
################

class PseudoY(Unobserved_Variational_Node):
    """
    General class for pseudodata nodes

    Notice that they are defined as Variational Nodes because they have a lower bound associated with it,
    but we are not using the .P and .Q distribution attributes
    """
    def __init__(self, dim, obs, Zeta=None, E=None):
        # Inputs:
        #  dim (2d tuple): dimensionality of each view
        #  obs (ndarray): observed data
        #  Zeta: parameter
        #  E (ndarray): initial expected value of pseudodata
        Unobserved_Variational_Node.__init__(self, dim)

        # Initialise observed data 
        self.obs = obs

        # Create a boolean mask of the data to hidden missing values
        if type(self.obs) != ma.MaskedArray: 
            self.mask()

        # Precompute some terms
        self.precompute()

        # Initialise Zeta
        self.Zeta = Zeta

        # Initialise expectation
        if E is not None:
            assert E.shape == dim.shape, "Problems with the dimensionalities"
        else:
            E = s.zeros(self.dim)
        self.E = E

    def updateParameters(self):
        Z = self.markov_blanket["Z"].getExpectation()
        W = self.markov_blanket["W"].getExpectation()
        self.Zeta = s.dot(Z,W.T)

    def mask(self):
        # Mask the observations if they have missing values
        self.obs = ma.masked_invalid(self.obs)

    def precompute(self):
        # Precompute some terms to speed up the calculations
        # self.N = self.dim[0]
        # self.D = self.dim[1]
        # self.lbconst = -0.5*self.N*self.D*s.log(2*s.pi)
        self.N = self.dim[0] - ma.getmask(self.obs).sum(axis=0)
        self.D = self.dim[1]
        self.lbconst = -0.5*s.sum(self.N)*s.log(2*s.pi)

    def updateExpectations(self):
        pass

    def getExpectation(self):
        return self.E

    def getObservations(self):
        return self.obs

    def getExpectations(self):
        return { 'E':self.getExpectation() }

    def getParameters(self):
        return { 'zeta':self.Zeta }
        
    # def getExpectations(self):
        # return { 'obs':self.getObservations() }

    def calculateELBO(self):
        # Compute Lower Bound using the Gaussian likelihood with pseudodata
        Z = self.markov_blanket["Z"].getExpectation()
        W = self.markov_blanket["W"].getExpectation()
        kappa = self.markov_blanket["kappa"].getExpectation()

        lb = self.lbconst + s.sum(self.N*s.log(kappa))/2 - s.sum( kappa * (self.E-s.dot(Z,W.T))**2 )/2
        return lb

class Poisson_PseudoY_Node(PseudoY):
    """
    Class for a Poisson pseudodata node with the following likelihood:
        p(y|x) \prop gamma(x) * e^{-gamma(x)}  (1)
    where gamma(x) is a rate function that is chosen to be convex and log-concave
    A simple choise for the rate function is e^{x} but this rate function is non-robust
    in the presence of outliers, so in Seeger et al they chose the function:
        gamma(x) = log(1+e^x)

    The data follows a Poisson distribution, but Followung Seeger et al the pseudodata Yhat_ij
    follows a normal distribution with mean E[W_{i,:}]*E[Z_{j,:}] and precision 'kappa_j'
    where 'kappa_j' is an upper bound of the second derivative of the loglikelihood:
        x_ij = sum_k^k w_{i,k}*z_{k,j}
        f_ij''(x_ij) <= kappa_j for all i,j

    For the Poisson likelihood with rate function (1), the upper bound kappa is calculated as follows:
        f_ij''(x_ij) = 0.25 + 0.17*ymax_j   where ymax_j = max(Y_{:,j})

    Pseudodata is updated as follows:
        yhat_ij = zeta_ij - f'(zeta_ij)/kappa_j = ...
    The bound degrades with the presence of entries with large y_ij, so one should consider
    clipping overly large counts

    """
    def __init__(self, dim, obs, Zeta=None, E=None):
        # - dim (2d tuple): dimensionality of each view
        # - obs (ndarray): observed data
        # - E (ndarray): initial expected value of pseudodata
        PseudoY.__init__(self, dim=dim, obs=obs, Zeta=Zeta, E=E)

        # Initialise the observed data
        assert s.all(s.mod(self.obs, 1) == 0), "Data must not contain float numbers, only integers"
        assert self.obs.shape == dim, "Problems with the dimensionalities"
        assert s.all(self.obs >= 0), "Data must not contain negative numbers"

    def ratefn(self, X):
        # Poisson rate function
        return s.log(1+s.exp(X))

    def clip(self, threshold):
        # The local bound degrades with the presence of large values in the observed data, which should be clipped
        pass

    def updateExpectations(self):
        # Update the pseudodata
        kappa = self.markov_blanket["kappa"].getValue()
        self.E = self.Zeta - sigmoid(self.Zeta)*(1-self.obs/self.ratefn(self.Zeta))/kappa

        pass

    def calculateELBO(self):
        # Compute Lower Bound using the Poisson likelihood with observed data
        Z = self.markov_blanket["Z"].getExpectation()
        W = self.markov_blanket["W"].getExpectation()
        tmp = self.ratefn(s.dot(Z,W.T))
        lb = s.sum( self.obs*s.log(tmp) - tmp)
        return lb
class Bernoulli_PseudoY_Node(PseudoY):
    """
    Class for a Bernoulli (0,1 data) pseudodata node with the following likelihood:
        p(y|x) = (e^{yx}) / (1+e^x)  (1)
        f(x) = -log p(y|x) = log(1+e^x) - yx
    
    The second derivative is upper bounded by kappa=0.25

    Folloiwng Seeger et al, the data follows a Bernoulli distribution but the pseudodata follows a
    normal distribution with mean E[W]*E[Z] and precision 'kappa'

    IMPROVE EXPLANATION

    Pseudodata is updated as follows:
        yhat_ij = zeta_ij - f'(zeta_ij)/kappa 
                = zeta_ij - 4*(sigmoid(zeta_ij) - y_ij)


    """
    def __init__(self, dim, obs, Zeta=None, E=None):
        # - dim (2d tuple): dimensionality of each view
        # - obs (ndarray): observed data
        # - E (ndarray): initial expected value of pseudodata
        PseudoY.__init__(self, dim=dim, obs=obs, Zeta=Zeta, E=E)

        # Initialise the observed data
        assert s.all( (self.obs==0) | (self.obs==1) ), "Data must be binary"

    def updateExpectations(self):
        # Update the pseudodata
        self.E = self.Zeta - 4*(sigmoid(self.Zeta) - self.obs)
        pass

    def calculateELBO(self):
        # Compute Lower Bound using the Bernoulli likelihood with observed data
        Z = self.markov_blanket["Z"].getExpectation()
        W = self.markov_blanket["W"].getExpectation()
        tmp = s.dot(Z,W.T)
        lik = s.sum( self.obs*tmp - s.log(1+s.exp(tmp)) )

        return lik
class Binomial_PseudoY_Node(PseudoY):
    """
    Class for a Binomial pseudodata node with the following likelihood:
        p(x|N,theta) = p(x|N,theta) = binom(N,x) * theta**(x) * (1-theta)**(N-x)
        f(x) = -log p(x|N,theta) = -log(binom(N,x)) - x*theta - (N-x)*(1-theta)
    
    The second derivative is upper bounded:
        f''(x_ij) <= 0.25*max(N_{:,j})

    Folloiwng Seeger et al, the pseudodata yhat_{nd} follows a normal distribution with mean 
    E[w_{i,:}]*E[z_{j,:}] and precision 'kappa_j'

    IMPROVE EXPLANATION

    Pseudodata is updated as follows
        yhat_ij = zeta_ij - f'(zeta_ij)/kappa_j
                = zeta_ij - (N_{ij}*sigmoid(zeta_ij) - y_ij)/kappa_d

    """
    def __init__(self, dim, obs, tot, Zeta=None, E=None):
        # - dim (2d tuple): dimensionality of each view
        # - obs (ndarray): observed data
        # - E (ndarray): initial expected value of pseudodata
        PseudoY.__init__(self, dim=dim, obs=None, Zeta=Zeta, E=E)

        # Initialise the observed data
        assert s.all(s.mod(obs, 1) == 0) and s.all(s.mod(tot, 1) == 0), "Data must not contain float numbers, only integers"
        assert s.all(obs >= 0) and s.all(tot >= 0), "Data must not contain negative numbers"
        assert s.all(obs <= tot), "Observed counts have to be equal or smaller than the total counts"
        self.obs = obs
        self.tot = tot


    def updateExpectations(self):
        # Update the pseudodata
        kappa = self.markov_blanket["kappa"].getValue()
        self.E = self.zeta - s.divide(self.tot*sigmoid(self.Zeta)-self.obs, kappa)
        pass

    def calculateELBO(self):
        # Compute Lower Bound using the Bernoulli likelihood with observed data
        Z = self.markov_blanket["Z"].getExpectation()
        W = self.markov_blanket["W"].getExpectation()

        tmp = sigmoid(s.dot(Z,W.T))

        tmp[tmp==0] = 0.00000001
        tmp[tmp==1] = 0.99999999
        lik = s.log(s.special.binom(self.tot,self.obs)).sum() + s.sum(self.obs*s.log(tmp)) + \
            s.sum((self.tot-self.obs)*s.log(1-tmp))
        return lik


