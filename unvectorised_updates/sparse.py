from __future__ import division
import numpy.linalg  as linalg
import numpy.ma as ma

# Import manually defined functions
from variational_nodes import *
from utils import *

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

class Y_Node(Observed_Variational_Node):
    def __init__(self, dim, obs):
        Observed_Variational_Node.__init__(self, dim, obs)
        # Create a boolean mask of the data to hidden missing values
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

        ## non-vectorised ##
        # old = self.Q.mean.copy()
        # D = Y.shape[1]
        # for n in xrange(self.N):
        #     for k in xrange(self.K):
        #         # Variance
        #         tmp = 0
        #         for d in xrange(D):
        #             tmp += tau[d]*SWW[d,k]
        #         self.Q.var[n,k] = 1/(tmp + 1)
        #         # Mean
        #         tmp = 0
        #         for d in xrange(D):
        #             tmp += tau[d]*SW[d,k] * (Y[n,d] - s.sum(SW[d,s.arange(self.K)!=k]*self.Q.mean[n,s.arange(self.K)!=k]))
        #         self.Q.mean[n,k] = self.Q.var[n,k]*tmp
        # self.Q.mean = old.copy()


        ## vectorised ##

        # Variance
        tmp = 1/((tau*SWW.T).sum(axis=1)+1)
        self.Q.var = s.repeat(tmp[None,:],self.N,0)

        # Mean
        if any(self.covariates): tmp = self.Q.mean[:,self.covariates]
        for k in xrange(self.K):
            tmp1 = SW[:,k]*tau
            tmp2 = Y - s.dot( self.Q.mean[:,s.arange(self.K)!=k] , SW[:,s.arange(self.K)!=k].T )
            # self.Q.mean[:,k] = self.Q.var[:,k] * s.dot(tmp2,tmp1)
            self.Q.mean[:,k] = self.Q.var[:,k] * ma.dot(tmp2,tmp1)

        # Do not update the latent variables associated with known covariates
        if any(self.covariates):
            self.Q.mean[:,self.covariates] = tmp

        pass

    def calculateELBO(self):
        lb_p = -self.Q.E2.sum()
        lb_q = -s.log(self.Q.var).sum() + self.N*self.K
        return (lb_p-lb_q)/2

    def removeFactors(self, *idx):
        keep = s.setdiff1d(s.arange(self.K),idx)
        self.Q.mean = self.Q.mean[:,keep]
        self.Q.var = self.Q.var[:,keep]
        self.Q.E = self.Q.E[:,keep]
        self.Q.E2 = self.Q.E2[:,keep]
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
        term1 = (Y**2).sum(axis=0)
        term2 = 2*(Y*s.dot(Z,SW.T)).sum(axis=0)
        term3 = (ZZ.dot(SWW.T)).sum(axis=0)
        term4 = s.diag(s.dot( SW.dot(Z.T), Z.dot(SW.T) )) - s.dot(Z**2,(SW**2).T).sum(axis=0)
        tmp = term1 - term2 + term3 + term4 
        
        # self.Q.a[:] = self.P.a + N/2
        self.Q.a = self.P.a + (Y.shape[0] - ma.getmask(Y).sum(axis=0))/2
        self.Q.b = self.P.b + tmp/2


        ## Non-vectorised ##
        # N = Y.shape[0]
        # K = Z.shape[1]
        # mask = ma.getmask(Y)
        # a = 0
        # for d in xrange(self.D):
        #     tmp = 0
        #     for n in xrange(N):
        #             tmptmp = 0
        #             for k in xrange(K):
        #                 a += SWW[d,k]*ZZ[n,k]
        #                 tmptmp += SWW[d,k]*ZZ[n,k]
        #                 if k < (K-1):
        #                     for j in xrange(k+1,K):
        #                         tmptmp += 2*(SW[d,k]*Z[n,k])*(SW[d,j]*Z[n,j])
        #             if not mask[n,d]: 
        #                 tmp += Y[n,d]**2 - 2*Y[n,d]*s.sum(SW[d,:]*Z[n,:])
        #             tmp += tmptmp
        #     self.Q.a[d] = self.P.a + s.sum(~mask[:,d])/2
        #     self.Q.b[d] = self.P.b + tmp/2

        pass

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
        S,SWW = tmp["ES"],tmp["ESWW"]

        # self.Q.a = self.P.a + D/2 # Updated in the initialisation
        self.Q.b = self.P.b + SWW.sum(axis=0)/2

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

        ## Non-vectorised ##
        # oldmean = self.Q.mean.copy()
        # oldtheta = self.Q.theta.copy()
        # oldSW = SW.copy()
        # for d in xrange(self.D):
        #     for k in xrange(self.K):
        #         term1 = s.log(self.P_theta/(1-self.P_theta))
        #         term2 = 0.5*s.log(alpha[k]/tau[d])
        #         term3 = 0.5*s.log(s.sum(ZZ[:,k]) + alpha[k]/tau[d])
        #         foo = ma.dot(Y[:,d],Z[:,k]) - s.dot(SW[d,s.arange(self.K)!=k],(Z[:,k]*Z[:,s.arange(self.K)!=k].T).sum(axis=1))
        #         bar = s.sum(ZZ[:,k]) + alpha[k]/tau[d]
        #         term4 = 0.5*tau[d]*foo**2 / bar

        #         # Update S
        #         self.Q.theta[d,k] = 1/(1+s.exp(-(term1+term2-term3+term4)))

        #         # Update W
        #         self.Q.mean[d,k] = foo/bar
        #         self.Q.var[d,k] = 1/(tau[d]*bar)

        #         # Update expectations
        #         SW[d,k] = self.Q.theta[d,k] * self.Q.mean[d,k]

        # self.Q.mean = oldmean.copy()
        # self.Q.theta = oldtheta.copy()
        # SW = oldSW.copy()


        ## Vectorised ##
        for k in xrange(self.K):
            term1 = s.log(self.P_theta/(1-self.P_theta))
            term2 = 0.5*s.log(s.divide(alpha[k],tau))
            term3 = 0.5*s.log(s.sum(ZZ[:,k]) + s.divide(alpha[k],tau))
            # term41 = s.dot(Y.T,Z[:,k]) 
            term41 = ma.dot(Y.T,Z[:,k]) 
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
        
        pass

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

    # def getExpectation(self):
        # return self.Q.ESW

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

        # Calculate ELBO for W
        lb_pw = (self.D*alpha["lnE"].sum() - s.sum(alpha["E"]*WW))/2
        lb_qw = -0.5*self.K*self.D - 0.5*s.log(S*self.Q.var + ((1-S)/alpha["E"])).sum()
        lb_w = lb_pw - lb_qw

        # Calculate ELBO for S
        Slower = 0.00001
        Supper = 0.99999
        S[S<Slower] = Slower
        S[S>Supper] = Supper
        # theta = self.P_theta
        lb_ps = s.sum( S*s.log(self.P_theta) + (1-S)*s.log(1-self.P_theta) )
        lb_qs = s.sum( S*s.log(S) + (1-S)*s.log(1-S) )
        lb_s = lb_ps - lb_qs

        return lb_w + lb_s