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







#############################################
## Specific classes for Seeger local nodes ##
#############################################



# class PseudoY_Bernoulli_Node(Local_Variational_Node):
#     """
#     Class for a Bernoulli (0,1 data) pseudodata node with the following likelihood:
#         p(y|x) = (e^{yx}) / (1+e^x)  (1)
#         f(x) = -log p(y|x) = log(1+e^x) - yx
    
#     The second derivative is upper bounded by kappa=0.25

#     Folloiwng Seeger et al, the data follows a Bernoulli distribution but the pseudodata follows a
#     normal distribution with mean E[W]*E[Z] and precision 'kappa'

#     IMPROVE EXPLANATION

#     Pseudodata is updated as follows:
#         yhat_ij = zeta_ij - f'(zeta_ij)/kappa 
#                 = zeta_ij - 4*(sigmoid(zeta_ij) - y_ij)


#     """
#     def __init__(self, dim, obs, zeta, E=None):
#         # - dim (list of integer tuples): dimensionality of each view
#         # - obs (list of ndarray): observed data
#         # - zeta (instance of 'Seeger_Zeta' class)
#         # - E (list of ndarray): initial expected value of pseudodata

#         Multiview_Node.__init__(self, dim=dim)

#         # Initialise kappa
#         self.kappa = 0.25

#         # Initialise the zeta parameter
#         for m in xrange(self.M): 
#             assert zeta.E[m].shape == dim[m], "Problems with the dimensionalities"
#         self.zeta = zeta.E

#         # Initialise the observed data
#         for m in xrange(self.M): 
#             assert obs[m].shape == dim[m], "Problems with the dimensionalities"
#             assert s.all( obs[m]==0 or obs[m]==1 ), "Data must not contain negative numbers"
#         self.obs = obs

#         # Initialise expectation
#         if E is not None:
#             assert len(E) == len(dim), "Problems with the dimensionalities"
#             for m in xrange(M): assert E[m].shape == dim[m].shape, "Problems with the dimensionalities"
#         else:
#             E = [ s.zeros((dim[m][0],dim[m][1])) for m in xrange(self.M) ] 
#         self.E = E


#     def updateExpectations(self):
#         # Update the pseudodata
#         for m in xrange(self.M):
#             self.E[m] = self.zeta[m] - 4*(sigmoid(self.zeta[m])*self.obs[m])
#         pass
            
#     def elbo(self, net):
#         # Compute Lower Bound using the Bernoulli likelihood with observed data
#         M = net.dim['M']
#         K = net.dim['K']
#         N = net.dim['N']
#         D = net.dim['D']

#         Z = net.nodes["Z"].Q.E

#         lik = 0
#         for m in xrange(self.M):
#             W = net.nodes["W"].Q[m].E
#             tmp = s.dot(Z,W.T)
#             lik += s.sum( self.obs[m]*tmp - s.log(1+s.exp(tmp)) )
#         return lik

#     def elbo2(self, net):
#         # Compute Lower Bound using the Gaussian likelihood with pseudodata

#         N = net.dim['N']
#         D = net.dim['D']
#         Z = net.nodes["Z"].Q.E

#         lb = 0
#         for m in xrange(self.M):
#             W = net.nodes["W"].Q[m].E
#             lb += -0.5*N*D[m]*s.log(2*s.pi) + N*0.5*s.log(self.kappa[m]).sum() - 0.5*s.sum( self.kappa[m] * (self.E[m] - s.dot(Z,W.T))**2  )
#         return lb
# # TO-FINISH
# class PseudoY_Binomial_Node(Local_Variational_Node):
#     """
#     Class for a Binomial pseudodata node with the following likelihood:
#         p(x|N,theta) = p(x|N,theta) = binom(N,x) * theta**(x) * (1-theta)**(N-x)
#         f(x) = -log p(x|N,theta) = -log(binom(N,x)) - x*theta - (N-x)*(1-theta)
    
#     The second derivative is upper bounded:
#         f''(x_ij) <= 0.25*max(N_{:,j})

#     Folloiwng Seeger et al, the pseudodata yhat_{nd} follows a normal distribution with mean 
#     E[w_{i,:}]*E[z_{j,:}] and precision 'kappa_j'

#     IMPROVE EXPLANATION

#     Pseudodata is updated as follows
#         yhat_ij = zeta_ij - f'(zeta_ij)/kappa_j
#                 = zeta_ij - (N_{ij}*sigmoid(zeta_ij) - y_ij)/kappa_d

#     """
#     def __init__(self, dim, obs, tot, zeta, E=None):
#         # - dim (list of integer tuples): dimensionality of each view
#         # - obs (list of ndarray): observed number of counts
#         # - tot (list of ndarray): total number of counts
#         # - zeta (instance of 'Seeger_Zeta' class)
#         # - E (list of ndarray): initial expected value of pseudodata
#         Multiview_Node.__init__(self, dim=dim)

#         # Initialise the Zeta parameter
#         for m in xrange(self.M): 
#             assert zeta.E[m].shape == dim[m], "Problems with the dimensionalities"
#         self.zeta = zeta.E

#         # Initialise the observed data
#         for m in xrange(self.M): 
#             assert obs[m].shape == dim[m], "Problems with the dimensionalities"
#             assert s.all(s.mod(obs[m], 1) == 0) and s.all(s.mod(tot[m], 1) == 0), "Data must not contain float numbers, only integers"
#             assert s.all(obs[m] >= 0) and s.all(tot[m] >= 0), "Data must not contain negative numbers"
#             assert s.all(obs[m] <= tot[m]), "Observed counts have to be equal or smaller than the total counts"
#         self.obs = obs
#         self.tot = tot

#         # Initialise kappa
#         self.kappa = [ 0.25*s.amax(self.tot[m],axis=0) for m in xrange(self.M) ]  

#         # Initialise expectation
#         if E is not None:
#             assert len(E) == len(dim), "Problems with the dimensionalities"
#             for m in xrange(M): assert E[m].shape == dim[m].shape, "Problems with the dimensionalities"
#         else:
#             E = [ s.zeros((dim[m][0],dim[m][1])) for m in xrange(self.M) ] 
#         self.E = E


#     def updateExpectations(self):
#         # Update the pseudodata
#         for m in xrange(self.M):
#             self.E[m] = self.zeta[m] - s.divide ( self.tot[m]*sigmoid(self.zeta[m]) - self.obs[m], self.kappa[m] )
#         pass
            
#     def elbo(self, net):
#         # Compute Lower Bound using the Binomial likelihood with observed data
#         M = net.dim['M']
#         K = net.dim['K']
#         N = net.dim['N']
#         D = net.dim['D']

#         Z = net.nodes["Z"].Q.E

#         lik = 0
#         for m in xrange(self.M):
#             W = net.nodes["W"].Q[m].E
#             tmp = s.dot(Z,W.T)
#             # is there a way to avoid that?
#             tmp[tmp==0] = 0.00000001
#             tmp[tmp==1] = 0.99999999

#             lik += s.special.binom(self.tot[m],self.obs[m]).sum() + s.sum(self.obs[m]*s.log(tmp)) + s.sum((self.tot[m]-self.obs[m])*s.log(1-tmp))

#         return lik

#     def elbo2(self, net):
#         # Compute Lower Bound using the Gaussian likelihood with pseudodata
#         N = net.dim['N']
#         D = net.dim['D']
#         Z = net.nodes["Z"].Q.E
#         lb = 0
#         for m in xrange(self.M):
#             W = net.nodes["W"].Q[m].E
#             lb += -0.5*N*D[m]*s.log(2*s.pi) + N*0.5*s.log(self.kappa[m]).sum() - 0.5*s.sum( self.kappa[m] * (self.E[m] - s.dot(Z,W.T))**2  )
#         return lb
