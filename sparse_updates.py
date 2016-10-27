from __future__ import division
# from time import time
import numpy.linalg  as linalg

from variational_nodes import *
from utils import *

"""
This module is used to define the variational nodes and the corresponding updates for
the element-wise spike and slab Group Factor Analysis.
Two layers of sparsity:
(1) Element-wise spike and slab
(2) Group-wise ARD
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

class Z_Node(UnivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None, qE2=None):
        UnivariateGaussian_Unobserved_Variational_Node.__init__(self, dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE)
        self.precompute()

    def precompute(self):
        self.N = self.dim[0]
        self.K = self.dim[1]

    def updateParameters(self):
        Y = self.markov_blanket["Y"].getExpectation()
        tmp = self.markov_blanket["SW"].getExpectations()
        tau = self.markov_blanket["tau"].getExpectation()

        M = len(Y)
        Y = s.concatenate([Y[m] for m in xrange(M)],axis=1)
        SW = s.concatenate([tmp[m]["ESW"]for m in xrange(M)],axis=0)
        SWW = s.concatenate([tmp[m]["ESWW"] for m in xrange(M)],axis=0)
        tau = s.concatenate([tau[m] for m in xrange(M)],axis=0)

        # Variance
        tmp = 1/((tau*SWW.T).sum(axis=1)+1)
        self.Q.var = s.repeat(tmp[None,:],self.N,0)

        # Mean: factorised over K
        for k in xrange(self.K):
            tmp1 = SW[:,k]*tau
            tmp2 = Y - s.dot( self.Q.mean[:,s.arange(self.K)!=k] , SW[:,s.arange(self.K)!=k].T )
            self.Q.mean[:,k] = self.Q.var[:,k] * s.dot(tmp2,tmp1)

        pass

    def calculateELBO(self):
        lb_p = -self.Q.E2.sum()
        lb_q = -s.log(self.Q.var).sum() + self.N*self.K
        return (lb_p-lb_q)/2

    def removeFactors(self, *idx):
        # Method to remove a set of (inactive) latent variables from the node
        keep = s.setdiff1d(s.arange(self.K),idx)
        self.Q.mean = self.Q.mean[:,keep]
        self.Q.var = self.Q.var[:,keep]
        self.Q.E = self.Q.E[:,keep]
        self.Q.E2 = self.Q.E2[:,keep]
        self.K = len(keep)
        self.dim = (self.N,self.K)

class Tau_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        self.D = self.dim[0]
        self.lbconst = s.sum(self.D*(self.P.a*s.log(self.P.b) - special.gammaln(self.P.a)))


    def updateParameters(self):
        Y = self.markov_blanket["Y"].getExpectation()
        tmp = self.markov_blanket["SW"].getExpectations()
        SW,SWW = tmp["ESW"], tmp["ESWW"]
        tmp = self.markov_blanket["Z"].getExpectations()
        Z,ZZ = tmp["E"],tmp["E2"]

        term1 = (Y**2).sum(axis=0)
        term2 = 2*(Y*Z.dot(SW.T)).sum(axis=0)
        term3 = (ZZ.dot(SWW.T)).sum(axis=0)
        term4 = s.diag(s.dot( SW.dot(Z.T), Z.dot(SW.T) )) - s.dot(Z**2,(SW**2).T).sum(axis=0)
        tmp = term1 - term2 + term3 + term4 # THIS MIGHT BE WRONG, CHECK IT
        # tmp = s.split( (term1 - term2 + term3 + term4) , s.cumsum(D))

        # self.Q.a = self.P.a + N/2 # this is already updated in the initialisation
        self.Q.b = self.P.b + tmp/2

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
        # Precompute some terms to speed up the calculations
        self.K = self.dim[0]
        self.lbconst = self.K * ( self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )

    def updateParameters(self):
        tmp = self.markov_blanket["SW"].getExpectations()
        S,SWW = tmp["ES"],tmp["ESWW"]

        # self.Q.a = self.P.a + D/2 # Updated in the initialisation
        # self.Q.a = self.P.a + S.sum(axis=0)/2
        self.Q.b = self.P.b + SWW.sum(axis=0)/2

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

        for k in xrange(self.K):
            # term1 = s.log(self.P.theta/(1-self.P.theta))
            term1 = s.log(self.P_theta/(1-self.P_theta))
            term2 = 0.5*s.log(s.divide(alpha[k],tau))
            term3 = 0.5*s.log(s.sum(ZZ[:,k]) + s.divide(alpha[k],tau))
            term41 = s.dot(Y.T,Z[:,k]) 
            term42 = s.dot( SW[:,s.arange(self.K)!=k] , (Z[:,k]*Z[:,s.arange(self.K)!=k].T).sum(axis=1) )                
            term43 = s.sum(ZZ[:,k]) + s.divide(alpha[k],tau)
            term4 = 0.5*tau * s.divide((term41-term42)**2,term43)
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
        ES = self.Q.ES
        EW = self.Q.EW
        ESW = self.Q.ESW
        ESWW = self.Q.ESWW
        EWW = self.Q.EWW
        return dict({'ES':ES, 'EW':EW, 'ESW':ESW, 'ESWW':ESWW, 'EWW':EWW})


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