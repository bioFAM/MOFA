from __future__ import division
import numpy.linalg  as linalg

from variational_nodes import *
from utils import *

"""
###########################################
## Updates for the Group Factor Analysis ##
###########################################

(Derivation of equations can be found in Ricard's MSc report)

Current nodes: 
    Y_Node: observed data
    W_Node: weights
    Tau_Node: precision of the noise
    Alpha_Node: ARD precision
    Z_Node: latent variables

Each node is a Variational_Node() class with the following main variables:
    Important methods:
    - precompute: precompute some terms to speed up the calculations
    - calculateELBO: calculate evidence lower bound using current estimates of expectations/params
    - getParameters: return current parameters
    - getExpectations: return current expectations
    - updateParameters: update parameters using current estimates of expectations
    - updateExpectations: update expectations using current estimates of parameters
    - removeFactors: remove a set of latent variables from the node

     Important attributes:
    - markov_blanket: dictionary that defines the set of nodes that are in the markov blanket of the current node
    - Q: an instance of Distribution() which contains the specification of the variational distribution 
    - P: an instance of Distribution() which contains the specification of the prior distribution 
    - dim: dimensionality of the node
"""

class Y_Node(Observed_Variational_Node):
    def __init__(self, dim, obs):
        Observed_Variational_Node.__init__(self, dim, obs)
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        self.N = self.dim[0]
        self.D = self.dim[1]
        self.likconst = -0.5*self.N*self.D*s.log(2*s.pi)

    def calculateELBO(self):
        tau_param = self.markov_blanket["tau"].getParameters()
        tau_exp = self.markov_blanket["tau"].getExpectations()
        # We make the assumption that the prior is so broad that is negligible
        lik = self.likconst + self.N*s.sum(tau_exp["lnE"])/2 - s.dot(tau_exp["E"],tau_param["b"])
        return lik
class W_Node(MultivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, qmean, qcov, qE=None, qE2=None):
        MultivariateGaussian_Unobserved_Variational_Node.__init__(self, dim=dim, qmean=qmean, qcov=qcov, qE=qE)
        self.precompute()

    def precompute(self):
        self.D = self.dim[0]
        self.K = self.dim[1]

    def updateParameters(self):
        Z = self.markov_blanket["Z"].getExpectation()
        ZZ = (self.markov_blanket["Z"].getExpectations()["E2"]).sum(axis=0)
        alpha = self.markov_blanket["alpha"].getExpectation()
        tau = (self.markov_blanket["tau"].getExpectation())[:,None,None]
        Y = self.markov_blanket["Y"].getExpectation()

        self.Q.cov = linalg.inv(tau*s.repeat(ZZ[None,:,:],self.D,0) + s.diag(alpha))
        tmp1 = tau*self.Q.cov
        tmp2 = Y.T.dot(Z)
        self.Q.mean = (tmp1[:,:,:]*tmp2[:,None,:]).sum(axis=2)

        pass

    def calculateELBO(self):
        alpha = self.markov_blanket["alpha"].getExpectations()["E"]
        logalpha = self.markov_blanket["alpha"].getExpectations()["lnE"]

        lb_p = self.D*s.sum(logalpha) - s.sum(self.Q.E2 * s.diag(alpha)[None,:,:])
        # lb_p = self.D*s.sum(logalpha)/2 - s.sum(self.Q.E2 * s.diag(alpha)[None,:,:])
        lb_q = -self.D*self.K - logdet(self.Q.cov).sum()

        # return lb_p - lb_q/2
        return (lb_p - lb_q)/2

    def removeFactors(self, *idx):
        # Method to remove a set of (inactive) latent variables from the node
        keep = s.setdiff1d(s.arange(self.K),idx)
        self.Q.mean = self.Q.mean[:,keep]
        self.Q.cov = self.Q.cov[:,:,keep][:,keep,:]
        self.Q.E = self.Q.E[:,keep]
        self.Q.E2 = self.Q.E2[:,:,keep][:,keep,:]
        self.K = len(keep)
        self.dim = (self.D,self.K)
class Tau_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        self.D = self.dim[0]
        self.lbconst = s.sum(self.D*(self.P.a*s.log(self.P.b) - special.gammaln(self.P.a)))

    def updateParameters(self):
        Z = self.markov_blanket["Z"].getExpectation()
        ZZ = (self.markov_blanket["Z"].getExpectations()["E2"]).sum(axis=0)
        W = self.markov_blanket["W"].getExpectation()
        WW = self.markov_blanket["W"].getExpectations()["E2"]
        Y = self.markov_blanket["Y"].getExpectation()

        # self.Q.a[:] = self.P.a + Z.shape[0]/2
        tmp = (Y**2).sum(axis=0) - 2*(Y*s.dot(Z,W.T)).sum(axis=0) + (WW*ZZ[None,:,:]).sum(axis=(1,2))
        self.Q.b = self.P.b + tmp/2

        pass

    def calculateELBO(self):
        # Calculate Variational Evidence Lower Bound
        p = self.P
        q = self.Q
        lb_p = self.lbconst + (p.a-1)*s.sum(q.lnE) - p.b*s.sum(q.E)
        lb_q = s.sum(q.a*s.log(q.b)) + s.sum((q.a-1)*q.lnE) - s.sum(q.b*q.E) - s.sum(special.gammaln(q.a))

        return lb_p - lb_q
class Alpha_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        self.K = self.dim[0]
        self.lbconst = self.K * ( self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )

    def updateParameters(self):
        WW = self.markov_blanket["W"].getExpectations()["E2"]
        # D = WW.shape[0]
        # self.Q.a[:] = self.P.a + D/2
        self.Q.b = self.P.b + s.diag(WW.sum(axis=0))/2

    def calculateELBO(self):
        # Calculate Variational Evidence Lower Bound
        p = self.P
        q = self.Q

        lb_p = self.lbconst + (p.a-1)*s.sum(q.lnE) - p.b*s.sum(q.E)
        lb_q = s.sum(q.a*s.log(q.b)) + s.sum((q.a-1)*q.lnE) - s.sum(q.b*q.E) - s.sum(special.gammaln(q.a))

        return lb_p - lb_q

    def removeFactors(self, *idx):
        # Method to remove a set of (inactive) latent variables from the node
        keep = s.setdiff1d(s.arange(self.K),idx)
        self.Q.a = self.Q.a[keep]
        self.Q.b = self.Q.b[keep]
        self.Q.E = self.Q.E[keep]
        self.Q.lnE = self.Q.lnE[keep]
        self.K = len(keep)
        self.dim = (self.K,1)
class Z_Node(UnivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None, qE2=None):
        UnivariateGaussian_Unobserved_Variational_Node.__init__(self, dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)
        self.precompute()

    def precompute(self):
        self.N = self.dim[0]
        self.K = self.dim[1]

    def updateParameters(self):
        # Method to update the parameters of the Q distribution of the node Z
        Y = self.markov_blanket["Y"].getExpectation()
        tau = self.markov_blanket["tau"].getExpectation()
        tmp = self.markov_blanket["W"].getExpectations()
        W,WW = tmp["E"], tmp["E2"]

        M = len(Y)

        # covariance
        cov = s.eye(self.K)
        for m in xrange(M):
            cov += (tau[m][:,None,None] * WW[m]).sum(axis=0)
        cov = linalg.inv(cov)
        self.Q.cov = s.repeat(cov[None,:,:],self.N,axis=0)

        # mean
        mean = s.zeros(self.dim[::-1])
        for m in xrange(M):
            mean += s.dot( W[m].T, (tau[m]*Y[m]).T )
        self.Q.mean = cov.dot(mean).T

        pass

    def calculateELBO(self):
        lb_p = -s.trace(self.Q.E2.sum(axis=0))/2
        lb_q = -self.N*logdet(self.Q.cov[0,:,:]) - self.N*self.K/2
        return lb_p - lb_q


    def updateParameters(self):

        ## vectorised ##
        Y = self.markov_blanket["Y"].getExpectation()
        tau = self.markov_blanket["tau"].getExpectation()
        tmp = self.markov_blanket["W"].getExpectations()
        W,WW = tmp["E"], tmp["E2"]

        tau = s.concatenate([net.nodes["tau"].Q[m].E for m in xrange(M)],axis=0)
        SWW = s.concatenate([net.nodes["SW"].Q[m].ESWW for m in xrange(M)],axis=0)
        Y = s.concatenate([net.nodes["Y"].obs[m] for m in xrange(M)],axis=1)
        SW = s.concatenate([net.nodes["SW"].Q[m].ESW for m in xrange(M)],axis=0)

        # Variance
        tmp = 1/((tau*SWW.T).sum(axis=1)+1)
        self.Q.var = s.repeat(tmp[None,:],N,0)

        # Mean: factorised over K
        for k in xrange(K):
            tmp1 = SW[:,k]*tau
            tmp2 = Y - s.dot( self.Q.mean[:,s.arange(K)!=k] , SW[:,s.arange(K)!=k].T )
            self.Q.mean[:,k] = self.Q.var[:,k] * s.dot(tmp2,tmp1)

        # Mean: APPROXIMATED Fully factorised approach
        # It is not correct because we have to use the updated latent variables
        # term1 = (SW.T).dot(s.diag(tau)).dot(s.dot(SW,self.Q.mean.T)) # (K,N)
        # term2 = Y.dot(s.diag(tau)).dot(SW) 
        # term1 = (tau*SW.T).dot(s.dot(SW,self.Q.mean.T)) 
        # term2 = (tau*Y).dot(SW)
        # term3 = s.repeat( (tau*(SW**2).T).sum(axis=1)[None,:] ,N,0) * self.Q.mean
        # self.Q.mean = self.Q.var * (-term1.T + term2 + term3)

        pass
    def removeFactors(self, *idx):
        # Method to remove a set of (inactive) latent variables from the node
        keep = s.setdiff1d(s.arange(self.K),idx)
        self.Q.mean = self.Q.mean[:,keep]
        self.Q.cov = self.Q.cov[:,:,keep][:,keep,:]
        self.Q.E = self.Q.E[:,keep]
        self.Q.E2 = self.Q.E2[:,:,keep][:,keep,:]
        self.K = len(keep)
        self.dim = (self.N,self.K)




