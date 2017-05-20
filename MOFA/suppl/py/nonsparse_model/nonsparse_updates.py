from __future__ import division
import numpy.linalg  as linalg
import numpy.ma as ma

from variational_nodes import *
from utils import *


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
        tauQ_param = self.markov_blanket["Tau"].getParameters("Q")
        tauP_param = self.markov_blanket["Tau"].getParameters("P")
        tau_exp = self.markov_blanket["Tau"].getExpectations()
        lik = self.likconst + 0.5*s.sum(self.N*(tau_exp["lnE"])) - s.dot(tau_exp["E"],tauQ_param["b"]-tauP_param["b"])
        return lik
class W_Node(MultivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, qmean, qcov, qE=None, qE2=None):
        MultivariateGaussian_Unobserved_Variational_Node.__init__(self, dim=dim, qmean=qmean, qcov=qcov, qE=qE)
        self.precompute()

    def precompute(self):
        self.D = self.dim[0]
        # self.K = self.dim[1]
        self.factors_axis = 1

    def updateParameters(self):
        Z = self.markov_blanket["Z"].getExpectation()
        ZZ = (self.markov_blanket["Z"].getExpectations()["E2"]).sum(axis=0)
        alpha = self.markov_blanket["alpha"].getExpectation()
        tau = (self.markov_blanket["tau"].getExpectation())[:,None,None]
        Y = self.markov_blanket["Y"].getExpectation()

        Qcov = linalg.inv(tau*s.repeat(ZZ[None,:,:],self.D,0) + s.diag(alpha))
        tmp1 = tau*self.Q.cov
        # tmp2 = Y.T.dot(Z)
        tmp2 = ma.dot(Y.T,Z).data
        Qmean = (tmp1[:,:,:]*tmp2[:,None,:]).sum(axis=2)

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean=Qmean, cov=QCov)

    def calculateELBO(self):

        # Collect parameters and expectations
        alpha = self.markov_blanket["alpha"].getExpectations()["E"]
        logalpha = self.markov_blanket["alpha"].getExpectations()["lnE"]
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()

        lb_p = self.D*s.sum(logalpha) - s.sum(Qexp['E2'] * s.diag(alpha)[None,:,:])
        lb_q = -self.D*self.K - logdet(Qpar['cov']).sum()

        return (lb_p - lb_q)/2
class Tau_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        self.D = self.dim[0]
        self.lbconst = s.sum(self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']))

    def updateParameters(self):
        Z = self.markov_blanket["Z"].getExpectation()
        ZZ = (self.markov_blanket["Z"].getExpectations()["E2"]).sum(axis=0)
        W = self.markov_blanket["W"].getExpectation()
        WW = self.markov_blanket["W"].getExpectations()["E2"]
        Y = self.markov_blanket["Y"].getExpectation()

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        Qa = Pa + (Y.shape[0] - ma.getmask(Y).sum(axis=0))/2
        # tmp = (Y**2).sum(axis=0) - 2*(Y*s.dot(Z,W.T)).sum(axis=0) + (WW*ZZ[None,:,:]).sum(axis=(1,2))
        tmp = (Y**2).sum(axis=0).data - 2*(Y*s.dot(Z,W.T)).sum(axis=0).data + (WW*ZZ[None,:,:]).sum(axis=(1,2))
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
        # Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super(Alpha_Node,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        # self.K = self.dim[0]
        self.lbconst = s.sum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) )
        self.factors_axis = 0

    def updateParameters(self):
        # Collect expectations from other nodes
        WW = self.markov_blanket["W"].getExpectations()["E2"]
        D = self.markov_blanket["W"].getExpectations()["E"].shape[0]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # D = WW.shape[0]
        # self.Q.a[:] = self.P.a + D/2

        Qa = Pa + 0.5*D
        Qb = Pb + 0.5*s.diag(WW.sum(axis=0))

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
class Z_Node(MultivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pcov, qmean, qcov, qE=None, idx_covariates=None):
        MultivariateGaussian_Unobserved_Variational_Node.__init__(self, dim=dim, pmean=pmean, pcov=pcov, qmean=qmean, qcov=qcov, qE=qE)
        # super(Z_Node,self).__init__(dim=dim, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)
        self.precompute()

        if idx_covariates is not None:
            self.covariates[idx_covariates] = True

    def precompute(self):
        self.N = self.dim[0]
        self.covariates = np.zeros(self.dim[1], dtype=bool)
        self.factors_axis = 1

    def updateParameters(self):
        # Method to update the parameters of the Q distribution of the node Z
        Y = self.markov_blanket["Y"].getExpectation()
        M = len(Y)
        tau = self.markov_blanket["tau"].getExpectation()
        tmp = self.markov_blanket["W"].getExpectations()
        W = [ tmp[m]["E"]for m in xrange(M) ]
        WW = [ tmp[m]["E2"]for m in xrange(M) ]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pmean, Qmean = P['mean'], Q['mean']

        # covariance
        cov = s.eye(self.K)
        for m in xrange(M):
            cov += (tau[m][:,None,None] * WW[m]).sum(axis=0)
        cov = linalg.inv(cov)
        Qcov = s.repeat(cov[None,:,:],self.N,axis=0)

        # mean
        if any(self.covariates): 
            oldmean = Qmean[:,self.covariates]
            
        tmp = s.zeros(self.dim[::-1])
        for m in xrange(M):
            # tmp += W[m].T.dot(s.diag(tau[m])).dot(Y[m].T)
            # tmp += ma.dot(W[m].T.dot(s.diag(tau[m])), Y[m].T)
            tmp += ma.dot(tau[m]*W[m].T,Y[m].T)
        Qmean = cov.dot(tmp).T

        # Do not update the latent variables associated with known covariates
        if any(self.covariates): 
            Qmean[:,self.covariates] = oldmean

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean=Qmean, cov=Qcov)

    def calculateELBO(self):
        # lb_p = -s.trace(self.Q.E2.sum(axis=0))/2
        # lb_q = -self.N*logdet(self.Q.cov[0,:,:]) - self.N*self.K/2
        # return lb_p - lb_q

        Qpar, Qexp = self.Q.getParameters(), self.Q.getExpectations()
        lb_p = -s.trace(Qexp['E2'].sum(axis=0))
        lb_q = -self.N*logdet(Qpar['cov'][0,:,:]) - self.N*self.K
        return (lb_p - lb_q)/2


