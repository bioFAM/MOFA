from __future__ import division
import numpy.linalg  as linalg
import numpy.ma as ma

# Import manually defined functions
from variational_nodes import *
from utils import *
from nodes import Constant_Node

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

    Important attributes:
    - markov_blanket: dictionary that defines the set of nodes that are in the markov blanket of the current node
    - Q: an instance of Distribution() which contains the specification of the variational distribution
    - P: an instance of Distribution() which contains the specification of the prior distribution
    - dim: dimensionality of the node

"""

class Y_Node(Constant_Variational_Node):
    def __init__(self, dim, value):
        Constant_Variational_Node.__init__(self, dim, value)

        # Create a boolean mask of the data to hidden missing values
        if type(self.value) != ma.MaskedArray:
            self.mask()

        # Precompute some terms
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        # self.N = self.dim[0]
        self.N = self.dim[0] - ma.getmask(self.value).sum(axis=0)
        self.D = self.dim[1]
        # self.likconst = -0.5*self.N*self.D*s.log(2*s.pi)
        self.likconst = -0.5*s.sum(self.N)*s.log(2*s.pi)

    def mask(self):
        # Mask the observations if they have missing values
        self.value = ma.masked_invalid(self.value)

    def calculateELBO(self):
        tau_param = self.markov_blanket["tau"].getParameters()
        tau_exp = self.markov_blanket["tau"].getExpectations()
        lik = self.likconst + 0.5*s.sum(self.N*(tau_exp["lnE"])) - s.dot(tau_exp["E"],tau_param["b"])
        return lik

class Z_Node(UnivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None):
        UnivariateGaussian_Unobserved_Variational_Node.__init__(self, dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE)
        self.precompute()

    def precompute(self):
        self.N = self.dim[0]
        self.K = self.dim[1]
        self.covariates = np.zeros(self.K, dtype=bool)

    def setCovariates(self, idx):
        # Method to define which factors are unupdated covariates
        #  - idx (int list): index of the columns of Z that are covariates
        self.covariates[idx] = True

    def updateParameters(self):

        # Collect expectations from other nodes
        # TO DO: MAKE THIS FASTER
        Y = self.markov_blanket["Y"].getExpectation()
        SWtmp = self.markov_blanket["SW"].getExpectations()
        tau = self.markov_blanket["tau"].getExpectation()
        M = len(Y)
        Y = ma.concatenate([Y[m] for m in xrange(M)],axis=1)
        SW = s.concatenate([SWtmp[m]["ESW"]for m in xrange(M)],axis=0)
        SWW = s.concatenate([SWtmp[m]["ESWW"] for m in xrange(M)],axis=0)
        tau = s.concatenate([tau[m] for m in xrange(M)],axis=0)

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pmean, Pvar, Qmean, Qvar = P['mean'], P['var'], Q['mean'], Q['var']

        # Variance
        tmp = (tau*SWW.T).sum(axis=1)
        tmp = s.repeat(tmp[None,:],self.N,0)
        tmp += 1./Pvar  # adding the prior precision to the updated precision
        Qvar = 1./tmp

        # Mean
        if any(self.covariates):
            oldmean = Qmean[:,self.covariates]

        for k in xrange(self.K):
            tmp1 = SW[:,k]*tau
            tmp2 = Y - s.dot( Qmean[:,s.arange(self.K)!=k] , SW[:,s.arange(self.K)!=k].T )
            tmp3 = ma.dot(tmp2,tmp1)
            tmp3 += 1./Pvar[:,k] * Pmean[:,k]
            Qmean[:,k] = Qvar[:,k] * tmp3

        # Do not update the latent variables associated with known covariates
        if any(self.covariates):
            Qmean[:,self.covariates] = oldmean

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean=Qmean, var=Qvar)

    def removeFactors(self, idx):
        self.P.removeDimensions(axis=1, idx=idx)
        self.Q.removeDimensions(axis=1, idx=idx)
        self.updateDim(axis=1, new_dim=self.dim[1]-len(idx))
        self.K -= len(idx)

    def calculateELBO(self):
        # Collect parameters and expectations
        Ppar,Qpar,Qexp = self.P.getParameters(), self.Q.getParameters(), self.Q.getExpectations()
        Pmean, Pvar, Qmean, Qvar = Ppar['mean'], Ppar['var'], Qpar['mean'], Qpar['var']
        QE,QE2 = Qexp['E'],Qexp['E2']

        # compute term from the exponential in the Gaussian
        tmp1 = 0.5*QE2 - Pmean * QE + 0.5*Pmean**2.0
        tmp1 = -(tmp1/Pvar).sum()

        # compute term from the precision factor in front of the Gaussian (TODO should be computed only once)
        tmp2 = - (s.log(Pvar)/2.).sum()

        lb_p = tmp1 + tmp2
        lb_q = - (s.log(Qvar).sum() + self.N*self.K)/2.

        return lb_p-lb_q

class Tau_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        self.D = self.dim[0]
        self.lbconst = s.sum(self.D*(self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a'])))

    def updateParameters(self):

        # Collect expectations from other nodes
        Y = self.markov_blanket["Y"].getExpectation()
        tmp = self.markov_blanket["SW"].getExpectations()
        SW,SWW = tmp["ESW"], tmp["ESWW"]
        Ztmp = self.markov_blanket["Z"].getExpectations()
        Z,ZZ = Ztmp["E"],Ztmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']

        # Calculate terms for the update
        term1 = (Y**2).sum(axis=0).data 
        term2 = 2*(Y*s.dot(Z,SW.T)).sum(axis=0).data
        term3 = (ZZ.dot(SWW.T)).sum(axis=0)
        term4 = s.diag(s.dot( SW.dot(Z.T), Z.dot(SW.T) )) - s.dot(Z**2,(SW**2).T).sum(axis=0)
        tmp = term1 - term2 + term3 + term4

        # Perform updates of the Q distribution
        Qa = Pa + (Y.shape[0] - ma.getmask(Y).sum(axis=0))/2
        Qb = Pb + tmp/2

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.expectations['E'], self.Q.expectations['lnE']

        # Do the calculations
        lb_p = self.lbconst + s.sum((Pa-1)*QlnE) - s.sum(Pb*QE)
        lb_q = s.sum(Qa*s.log(Qb)) + s.sum((Qa-1)*QlnE) - s.sum(Qb*QE) - s.sum(special.gammaln(Qa))

        return lb_p - lb_q

class Alpha_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        # self.K = self.dim[0]
        # self.lbconst = self.K * ( self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )
        self.lbconst = s.sum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) )

    def updateParameters(self):

        # Collect expectations from other nodes
        tmp = self.markov_blanket["SW"].getExpectations()
        S,EWW,ESWW = tmp["ES"],tmp["EWW"],tmp["ESWW"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']

        # ARD prior on What
        Qb = Pb + 0.5*EWW.sum(axis=0)
        Qa = Pa + 0.5*EWW.shape[0]

        # ARD prior on W
        # self.Q.a = self.P.a + S.sum(axis=0)/2
        # self.Q.b = self.P.b + ESWW.sum(axis=0)/2

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.expectations['E'], self.Q.expectations['lnE']

        # Do the calculations
        lb_p = self.lbconst + s.sum((Pa-1)*QlnE) - s.sum(Pb*QE)
        lb_q = s.sum(Qa*s.log(Qb)) + s.sum((Qa-1)*QlnE) - s.sum(Qb*QE) - s.sum(special.gammaln(Qa))

        return lb_p - lb_q

    def removeFactors(self, idx):
        self.P.removeDimensions(axis=0, idx=idx)
        self.Q.removeDimensions(axis=0, idx=idx)
        self.updateDim(axis=0, new_dim=self.dim[0]-len(idx))

class SW_Node(BernoulliGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, ptheta, qtheta):
        BernoulliGaussian_Unobserved_Variational_Node.__init__(self, dim=dim, pmean=pmean, pvar=pvar, ptheta=ptheta, qmean=qmean, qvar=qvar, qtheta=qtheta)
        self.precompute()

    def precompute(self):
        self.D = self.dim[0]
        self.K = self.dim[1]

    def updateParameters(self):

        # Collect expectations from other nodes
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
        # I THINK THIS SHOULD NOT BE HERE...
        if theta_lnE.shape != self.Q.mean.shape:
            theta_lnE = s.repeat(theta_lnE[None,:],self.Q.mean.shape[0],0)
        if theta_lnEInv.shape != self.Q.mean.shape:
            theta_lnEInv = s.repeat(theta_lnEInv[None,:],self.Q.mean.shape[0],0)

        # Collect parameters from the P and Q distributions of this node
        # ....

        all_term1 = theta_lnE - theta_lnEInv
        # all_term1 = s.log(theta/(1.-theta))


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

    def getExpectations(self):
        return dict({'ES':self.Q.ES, 'EW':self.Q.EW, 'ESW':self.Q.ESW,'ESWW':self.Q.ESWW, 'EWW':self.Q.EWW})

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

        return lb_w + lb_s

    def removeFactors(self, idx):
        # Method to remove a set of (inactive) latent variables from the node
        keep = s.setdiff1d(s.arange(self.K),idx)
        # variational distribution terms
        self.Q.mean = self.Q.mean[:,keep]
        self.Q.var = self.Q.var[:,keep]
        self.Q.theta = self.Q.theta[:,keep]
        self.Q.ES = self.Q.ES[:,keep]
        self.Q.EW = self.Q.EW[:,keep]
        self.Q.ESW = self.Q.ESW[:,keep]
        self.Q.ESWW = self.Q.ESWW[:,keep]
        self.Q.EWW = self.Q.EWW[:,keep]
        # prior terms (TO-DO: UPDATE EXPECTATIONS)
        self.P.mean = self.P.mean[:,keep]
        self.P.var = self.P.var[:,keep]
        self.P.theta = self.P.theta[:,keep]
        # others
        self.K = len(keep)
        self.dim = (self.D,self.K)

    # def removeFactors(self, idx):
        # self.P.removeDimensions(axis=1, idx=idx)
        # self.Q.removeDimensions(axis=1, idx=idx)
        # self.updateDim(axis=1, new_dim=self.dim[1]-len(idx))

class Theta_Node_No_Annotation(Beta_Unobserved_Variational_Node):
    """
    This class comtain a Theta node associate to factors for which
    we dont have annotations.

    The inference is done per view and factor, so the dimension of the node is the
    number of non-annotated factors

    the updateParameters function needs to know what factors are non-annotated in
    order to choose from the S matrix
    """

    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        Beta_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

    def updateParameters(self, factors_selection=None):

        # Collect expectations from other nodes
        S = self.markov_blanket['SW'].getExpectations()["ES"]

        # Collect parameters from the P and Q distributions of this node
        Pa, Pb = self.P.params['a'], self.P.params['b']
        Qa, Qb = self.Q.params['a'], self.Q.params['b']

        # Precompute terms
        if factors_selection is not None:
            S = S[:,factors_selection]
        tmp1 = S.sum(axis=0)

        # Perform updates
        Qa = tmp1 + Pa
        Qb = Pb - tmp1 + S.shape[0]

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):

        # Collect parameters and expectations
        Ppar,Qpar,Qexp = self.P.getParameters(), self.Q.getParameters(), self.Q.getExpectations()
        Pa, Pb, Qa, Qb = Ppar['a'], Ppar['b'], Qpar['a'], Qpar['b']
        QE, QlnE, QlnEInv = Qexp['E'], Qexp['lnE'], Qexp['lnEInv']

        # minus cross entropy of Q and P
        tmp1 = (Pa-1)*QlnE + (Pb-1)*QlnEInv
        tmp1 -= special.betaln(Pa,Pb)
        lb_p = tmp1.sum()

        # minus entropy of Q
        tmp2 = (Qa-1)*QlnE + (Qb-1)*QlnEInv
        tmp2 -= special.betaln(Qa, Qb)
        lb_q = tmp2.sum()

        return lb_p - lb_q

    def removeFactors(self, idx):
        self.P.removeDimensions(axis=0, idx=idx)
        self.Q.removeDimensions(axis=0, idx=idx)
        self.updateDim(axis=0, new_dim=self.dim[0]-len(idx))

class Theta_Constant_Node(Constant_Variational_Node):
    def __init__(self, dim, value):
        super(Theta_Constant_Node, self).__init__(dim, value)
        self.precompute()

    def precompute(self):
        self.E = self.value
        self.lnE = s.log(self.value)
        self.lnEInv = s.log(1-self.value)


    # TO-DO...
    # def removeFactors(self, idx):
    #     self.P.removeDimensions(self, idx, axis=1)
    #     self.Q.removeDimensions(self, idx, axis=1)

