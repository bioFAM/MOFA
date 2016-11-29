from __future__ import division
import numpy.linalg  as linalg
import numpy.ma as ma

# Import manually defined functions
from variational_nodes import *
from utils import *
import pdb
import math
from constant_nodes import Constant_Node

import scipy.special as special

"""
###################################################
## Updates for the Sparse Group Factor Analysis ##
###################################################

Extensions with respect to the Gaussian Group Factor Analysis:
- Element-wise spike and slab

(Derivation of equations can be found in file XX)

Current nodes:
    Y_Node: observed data
    SW_Node: spike and slab weights
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

# TODO: in all remove factor functions, the dimensions of the distribution dont
# seem to be updated. Would dbe neater to have a remove factor implemented in each distribution ?
# which take care of this ?

class Y_Node(Observed_Variational_Node):
    def __init__(self, dim, obs):
        Observed_Variational_Node.__init__(self, dim, obs)

        # Create a boolean mask of the data to hidden missing values
        if type(self.obs) != ma.MaskedArray:
            self.mask()
        # Precompute some terms
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        # self.N = self.dim[0]
        self.N = self.dim[0] - ma.getmask(self.obs).sum(axis=0)
        self.D = self.dim[1]
        # self.likconst = -0.5*self.N*self.D*s.log(2*s.pi)
        self.likconst = -0.5*s.sum(self.N)*s.log(2*s.pi)
        pass

    def mask(self):
        # Mask the observations if they have missing values
        self.obs = ma.masked_invalid(self.obs)
        pass

    def calculateELBO(self):
        tau_param = self.markov_blanket["tau"].getParameters()
        tau_exp = self.markov_blanket["tau"].getExpectations()
        lik = self.likconst + s.sum(self.N*(tau_exp["lnE"]))/2 - s.dot(tau_exp["E"],tau_param["b"])
        return lik

class Z_Node(UnivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None, qE2=None):
        UnivariateGaussian_Unobserved_Variational_Node.__init__(self, dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE)
        self.precompute()

    def precompute(self):
        self.N = self.dim[0]
        self.K = self.dim[1]
        self.covariates = np.zeros(self.K, dtype=bool)

    def setCovariates(self, idx):
        # Method to define which factors are unupdated covariates
        # Input:
        #  idx (integer list): index of the columns of Z that are covariates
        self.covariates[idx] = True

    def updateParameters(self):
        Y = self.markov_blanket["Y"].getExpectation()
        tmp = self.markov_blanket["SW"].getExpectations()
        tau = self.markov_blanket["tau"].getExpectation()

        M = len(Y)
        # Y = s.concatenate([Y[m] for m in xrange(M)],axis=1)
        Y = ma.concatenate([Y[m] for m in xrange(M)],axis=1)
        SW = s.concatenate([tmp[m]["ESW"]for m in xrange(M)],axis=0)
        SWW = s.concatenate([tmp[m]["ESWW"] for m in xrange(M)],axis=0)
        tau = s.concatenate([tau[m] for m in xrange(M)],axis=0)

        # Variance
        # pdb.set_trace()
        tmp = (tau*SWW.T).sum(axis=1)
        tmp = s.repeat(tmp[None,:],self.N,0)
        tmp += 1./self.P.var  # adding the prior precision to the updated precision
        self.Q.var = 1./tmp

        # Mean
        if any(self.covariates):
            oldmean = self.Q.mean[:,self.covariates]

        for k in xrange(self.K):
            tmp1 = SW[:,k]*tau
            tmp2 = Y - s.dot( self.Q.mean[:,s.arange(self.K)!=k] , SW[:,s.arange(self.K)!=k].T )
            # self.Q.mean[:,k] = self.Q.var[:,k] * s.dot(tmp2,tmp1)
            tmp3 = ma.dot(tmp2,tmp1)
            tmp3 += 1./self.P.var[:, k] * self.P.mean[:, k]  # adding contribution from the prior
            self.Q.mean[:,k] = self.Q.var[:,k] * tmp3

        # Do not update the latent variables associated with known covariates
        if any(self.covariates):
            self.Q.mean[:,self.covariates] = oldmean

    def calculateELBO(self):
        # term from the exponential term in the Gaussian
        tmp1 = self.Q.E2/2. - self.P.mean * self.Q.E + self.P.mean**2.0/2.
        tmp1 = -(tmp1/self.P.var).sum()

        # term from the precision factor in front of the Gaussian (TODO should be computed only once)
        tmp2 = - (s.log(self.P.var)/2.).sum()

        lb_p = tmp1 + tmp2
        lb_q = - (s.log(self.Q.var).sum() + self.N*self.K)/2.

        return lb_p-lb_q

    def removeFactors(self, *idx):
        keep = s.setdiff1d(s.arange(self.K),idx)
        #remove variational distribution terms
        self.Q.mean = self.Q.mean[:,keep]
        self.Q.var = self.Q.var[:,keep]
        self.Q.E = self.Q.E[:,keep]
        self.Q.E2 = self.Q.E2[:,keep]
        # remove prior terms
        self.P.mean = self.P.mean[:,keep]
        self.P.var = self.P.var[:,keep]
        # others
        self.covariates = self.covariates[keep]
        self.K = len(keep)
        self.dim = (self.N,self.K)

class Tau_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        self.D = self.dim[0]
        self.lbconst = s.sum(self.D*(self.P.a*s.log(self.P.b) - special.gammaln(self.P.a)))

    def updateParameters(self):
        Y = self.markov_blanket["Y"].getExpectation()
        tmp = self.markov_blanket["SW"].getExpectations()
        SW,SWW = tmp["ESW"], tmp["ESWW"]
        tmp = self.markov_blanket["Z"].getExpectations()
        Z,ZZ = tmp["E"],tmp["E2"]
        ## Vectorised ##
        term1 = (Y**2).sum(axis=0).data
        # term2 = 2*(Y*s.dot(Z,SW.T)).sum(axis=0)
        term2 = 2*(Y*s.dot(Z,SW.T)).sum(axis=0).data
        term3 = (ZZ.dot(SWW.T)).sum(axis=0)
        term4 = s.diag(s.dot( SW.dot(Z.T), Z.dot(SW.T) )) - s.dot(Z**2,(SW**2).T).sum(axis=0)
        tmp = term1 - term2 + term3 + term4

        # self.Q.a[:] = self.P.a + N/2
        self.Q.a = self.P.a + (Y.shape[0] - ma.getmask(Y).sum(axis=0))/2
        self.Q.b = self.P.b + tmp/2

    def calculateELBO(self):
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
        self.K = self.dim[0]
        self.lbconst = self.K * ( self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )

    def updateParameters(self):
        tmp = self.markov_blanket["SW"].getExpectations()
        S,EWW,ESWW = tmp["ES"],tmp["EWW"],tmp["ESWW"]

        # ARD prior on What
        # pdb.set_trace()
        self.Q.b = self.P.b + EWW.sum(axis=0)/2.
        self.Q.a = s.repeat(self.P.a + EWW.shape[0]/2., self.K) # Updated in the initialisation

        # ARD prior on W
        # self.Q.a = self.P.a + S.sum(axis=0)/2
        # self.Q.b = self.P.b + ESWW.sum(axis=0)/2

    def calculateELBO(self):
        p = self.P
        q = self.Q
        lb_p = self.lbconst + (p.a-1)*s.sum(q.lnE) - p.b*s.sum(q.E)
        lb_q = s.sum(q.a*s.log(q.b)) + s.sum((q.a-1)*q.lnE) - s.sum(q.b*q.E) - s.sum(special.gammaln(q.a))
        return lb_p - lb_q

    def removeFactors(self, *idx):
        keep = s.setdiff1d(s.arange(self.K),idx)
        self.Q.a = self.Q.a[keep]
        self.Q.b = self.Q.b[keep]
        self.Q.E = self.Q.E[keep]
        self.Q.lnE = self.Q.lnE[keep]
        self.K = len(keep)
        self.dim = (self.K,1)

class SW_Node(BernoulliGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, qmean, qvar, ptheta, qtheta):
        BernoulliGaussian_Unobserved_Variational_Node.__init__(self, dim=dim, qmean=qmean, qvar=qvar, ptheta=ptheta, qtheta=qtheta)
        self.precompute()

    def precompute(self):
        self.D = self.dim[0]
        self.K = self.dim[1]

    def updateParameters(self):
        tmp = self.markov_blanket["Z"].getExpectations()
        Z,ZZ = tmp["E"],tmp["E2"]
        tau = self.markov_blanket["tau"].getExpectation()
        Y = self.markov_blanket["Y"].getExpectation()
        alpha = self.markov_blanket["alpha"].getExpectation()
        SW = self.Q.ESW[:]
        # TODO make general in mixed node
        theta_lnE = self.markov_blanket['Theta'].getExpectations()['lnE']
        theta_lnEInv = self.markov_blanket['Theta'].getExpectations()['lnEInv']

        # check dimensions of theta and expand if necessary
        if theta_lnE.shape != self.Q.mean.shape:
            theta_lnE = s.repeat(theta_lnE[None,:],self.Q.mean.shape[0],0)
        if theta_lnEInv.shape != self.Q.mean.shape:
            theta_lnEInv = s.repeat(theta_lnEInv[None,:],self.Q.mean.shape[0],0)


        all_term1 = theta_lnE - theta_lnEInv
        # all_term1 = s.log(theta/(1.-theta))
        ## Vectorised ##
        for k in xrange(self.K):

            term1 = all_term1[:, k]
            term2 = 0.5*s.log(s.divide(alpha[k],tau))
            term3 = 0.5*s.log(s.sum(ZZ[:,k]) + s.divide(alpha[k],tau))
            # term41 = ma.dot(Y.T,Z[:,k])
            term41 = ma.dot(Y.T,Z[:,k]).data
            term42 = s.dot( SW[:,s.arange(self.K)!=k] , (Z[:,k]*Z[:,s.arange(self.K)!=k].T).sum(axis=1) )
            term43 = s.sum(ZZ[:,k]) + s.divide(alpha[k],tau)
            term4 = 0.5*tau * s.divide((term41-term42)**2,term43)

            # Update S
            self.Q.theta[:,k] = 1/(1+s.exp(-(term1+term2-term3+term4)))

            # Update W
            self.Q.mean[:,k] = s.divide(term41-term42,term43)
            self.Q.var[:,k] = s.divide(1,tau*term43)

            # Update Expectations for the next iteration
            SW[:,k] = self.Q.theta[:,k] * self.Q.mean[:,k]

    def updateExpectations(self):
        alpha = self.markov_blanket["alpha"].getExpectation()
        self.Q.ES = self.Q.theta[:]
        self.Q.EW = self.Q.mean[:]
        self.Q.ESW = self.Q.ES * self.Q.EW
        self.Q.ESWW = self.Q.ES * (self.Q.EW**2 + self.Q.var)
        self.Q.EWW = self.Q.ES * (self.Q.EW**2 + self.Q.var)  + (1-self.Q.ES)*s.repeat(1/alpha[None,:],self.D,0)
        pass

    def getExpectations(self):
        return dict({'ES':self.Q.ES, 'EW':self.Q.EW, 'ESW':self.Q.ESW,'ESWW':self.Q.ESWW, 'EWW':self.Q.EWW})

    def removeFactors(self, *idx):
        # Method to remove a set of (inactive) latent variables from the node
        keep = s.setdiff1d(s.arange(self.K),idx)
        self.Q.mean = self.Q.mean[:,keep]
        self.Q.var = self.Q.var[:,keep]
        self.Q.theta = self.Q.theta[:,keep]
        self.Q.ES = self.Q.ES[:,keep]
        self.Q.EW = self.Q.EW[:,keep]
        self.Q.ESW = self.Q.ESW[:,keep]
        self.Q.ESWW = self.Q.ESWW[:,keep]
        self.Q.EWW = self.Q.EWW[:,keep]
        self.K = len(keep)
        self.dim = (self.D,self.K)

    def calculateELBO(self):
        # Calculate Variational Evidence Lower Bound
        alpha = self.markov_blanket["alpha"].getExpectations()
        exp = self.getExpectations()
        S = exp["ES"]
        WW = exp["EWW"]
        theta_lnE = self.markov_blanket['Theta'].getExpectations()['lnE']
        theta_lnEInv = self.markov_blanket['Theta'].getExpectations()['lnEInv']

        # Calculate ELBO for W
        lb_pw = (self.D*alpha["lnE"].sum() - s.sum(alpha["E"]*WW))/2
        lb_qw = -0.5*self.K*self.D - 0.5*s.log(S*self.Q.var + ((1-S)/alpha["E"])).sum()
        lb_w = lb_pw - lb_qw

        # Calculate ELBO for S
        # pdb.set_trace()
        # Slower = 0.00001
        # Supper = 0.99999
        # S[S<Slower] = Slower
        # S[S>Supper] = Supper

        lb_ps = s.sum( S*theta_lnE + (1-S)*theta_lnEInv)

        lb_qs_tmp = S*s.log(S) + (1-S)*s.log(1-S)
        lb_qs_tmp[s.isnan(lb_qs_tmp)] = 0

        lb_qs = s.sum(lb_qs_tmp)
        lb_s = lb_ps - lb_qs

        # if math.isnan(lb_w + lb_s):
            # pdb.set_trace()
        return lb_w + lb_s

class Theta_Node_No_Annotation(Beta_Unobserved_Variational_Node):
    """
    This class comtain a Theta node associate to factors for which
    we dont have annotations.

    The inference is done per view and factor, so the dimension of the node is the
    number of non-annotated factors

    the updateParameters function needs to know what factors are non-annotated in
    order to choose from the S matrix
    """

    def __init__(self, dim, pa=1., pb=1., qa=1., qb=1., qE=None):
        Beta_Unobserved_Variational_Node.__init__(self, dim, pa, pb, qa, qb, qE)

    def updateParameters(self, factors_selection=None):
        # get needed node from the markov_blanket
        tmp = self.markov_blanket['SW'].getExpectations()
        S = tmp["ES"]  # S is of dimension D*K

        if factors_selection is not None:
            S = S[:, factors_selection]

        tmp1 = S.sum(axis=0)
        self.Q.a = tmp1 + self.P.a
        self.Q.b = self.P.b - tmp1 + S.shape[0]

    def getExpectations(self):
        return {'E':self.Q.E, 'lnE':self.Q.lnE, 'lnEInv':self.Q.lnEInv}

    def removeFactors(self, *idx):
        keep = s.setdiff1d(s.arange(self.dim[0]),idx)
        #remove variational distribution terms
        self.Q.a = self.Q.a[keep]
        self.Q.b = self.Q.b[keep]
        self.Q.E = self.Q.E[keep]
        self.Q.lnE = self.Q.lnE[keep]
        self.Q.lnEInv = self.Q.lnEInv[keep]
        # remove prior terms
        self.P.a = self.P.a[keep]
        self.P.b = self.P.b[keep]
        self.P.E = self.P.E[keep]
        # update dimensionalities
        self.P.dim = (len(keep),)
        self.Q.dim = (len(keep),)
        self.dim = (len(keep),)

    def calculateELBO(self):
        # minus cross entropy of Q and P
        tmp1 = (self.P.a -1.) * self.Q.lnE + (self.P.b -1.) * self.Q.lnEInv
        tmp1 -= special.betaln(self.P.a, self.P.b)
        lbp = tmp1.sum()

        # minus entropy of Q
        tmp2 = (self.Q.a -1.) * self.Q.lnE + (self.Q.b -1.) * self.Q.lnEInv
        tmp2 -= special.betaln(self.Q.a, self.Q.b)
        lbq = tmp2.sum()

        return lbp - lbq

# inheritance to Variational_Node is purely technical (so that there is an
# update_parameters function for instance)
class Theta_Constant_Node(Constant_Node, Variational_Node):
    """docstring for Theta_Constant_Node."""
    def __init__(self, dim, value):
        super(Theta_Constant_Node, self).__init__(dim, value)
        self.precompute()

    def precompute(self):
        self.E = self.value
        self.lnE = s.log(self.value)
        self.lnEInv = s.log(1-self.value)

    def getExpectations(self):
        return {'E': self.E, 'lnE': self.lnE, 'lnEInv': self.lnEInv}

    # constant tehta does not contribute to ELBO
    def calculateELBO(self):
        return 0

    # TODO what would be neater for all removeFcator functions would be to have
    # a list of class mebers depending on factors and to be updated when removing
    # factors -> only one def of the removeFactors function and then you just
    # fill in the list according to the specific node
    def removeFactors(self, *idx):
        if len(self.dim) == 1:
            keep = s.setdiff1d(s.arange(self.dim[0]),idx)
            self.value = self.value[keep]
            self.dim = (len(self.value),)
            self.E = self.E[keep]
            self.lnE = self.lnE[keep]
            self.lnEInv = self.lnEInv[keep]
        else:
            keep = s.setdiff1d(s.arange(self.dim[1]),idx)
            self.value = self.value[:, keep]
            self.dim = self.value.shape
            self.E = self.E[:, keep]
            self.lnE = self.lnE[:, keep]
            self.lnEInv = self.lnEInv[:, keep]
