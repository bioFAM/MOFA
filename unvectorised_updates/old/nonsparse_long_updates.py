from __future__ import division
from time import time
import numpy.linalg  as linalg

from singleview_nodes import *
from multiview_nodes import *
from utils import *

"""
Works perfectly
"""

class Y_Node(Multiview_Observed_Node):
    def __init__(self, dim, obs):
        Multiview_Observed_Node.__init__(self,dim)
        self.obs = obs
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        self.Yconst = [ (self.obs[m]**2).sum(axis=0) for m in xrange(self.M) ]
        self.likconst = -0.5*s.sum( [self.dim[0][0]*self.dim[m][1]*s.log(2*s.pi) for m in xrange(self.M) ] )

    def elbo(self, net):
        M = net.dim['M']
        K = net.dim['K']
        N = net.dim['N']
        D = net.dim['D']

        # Z = net.nodes["Z"].Q.E1
        # W = [ net.nodes["W"].Q[m].E1 for m in xrange(M) ]
        # tau = [ net.nodes["tau"].Q[m] for m in xrange(M) ]
        # mu = [ net.nodes["mu"].Q[m].E1 for m in xrange(M) ]
        # Y = [ self.obs[m] for m in xrange(M) ]

        # concatenated
        # tau_lnE = s.concatenate([tau[m].lnE for m in xrange(M)],axis=0)
        # tau_E = s.concatenate([tau[m].E for m in xrange(M)],axis=0)
        # tau_qb = s.concatenate([tau[m].b for m in xrange(M)],axis=0)
        # tau_pb = net.nodes["tau"].P.b
        # lik = -0.5*N*s.sum(D)*s.log(2*s.pi) + 0.5*N*s.sum(tau_lnE) - s.sum(tau_E * (tau_qb-tau_pb))

        # vectorised
        tau = net.nodes["tau"]
        lik = self.likconst.copy()
        for m in xrange(M):
            lik += 0.5*N*s.sum(tau.Q[m].lnE) - s.dot(tau.Q[m].E,(tau.Q[m].b-tau.P.b))
        return lik

class Z_Node(MultivariateGaussian_Unobserved_Node):
    def __init__(self, dim, qmean, qcov, qE1=None, qE2=None):
    # def __init__(self, dim, pmean, pcov, qmean, qcov, qE1=None, qE2=None):
        MultivariateGaussian_Unobserved_Node.__init__(self, dim=dim, qmean=qmean, qcov=qcov, qE1=qE1, qE2=qE2)
        # MultivariateGaussian_Unobserved_Node.__init__(self, dim=dim, pmean=pmean, pcov=pcov, qmean=qmean, qcov=qcov, qE1=qE1)

    def updateParameters(self, net):
        # Method to update the parameters of the Q distribution of the node Z
        # Inputs:
        # - net: instance of BayesNet() containing the entire bayesian betwork

        M = net.dim['M']
        K = net.dim['K']
        N = net.dim['N']
        D = net.dim['D']

        W = [ net.nodes["W"].Q[m].E1 for m in xrange(M) ]
        WW = [ net.nodes["W"].Q[m].E2 for m in xrange(M) ]
        tau = [ net.nodes["tau"].Q[m].E for m in xrange(M) ]
        Y = [ net.nodes["Y"].obs[m] for m in xrange(M) ]

        ## Non-vectorised ##

        # # covariance
        # tmp = s.zeros((K,K))
        # for m in xrange(M):
        #     for d in xrange(D[m]):
        #         tmp += tau[m][d]*WW[m][d,:,:]
        # cov = linalg.inv(s.eye(K) + tmp)
        # self.Q.cov = s.repeat(cov[None,:,:],N,axis=0)

        # # mean
        # self.Q.mean = s.zeros((N,K))
        # for n in xrange(N):
        #     tmp = 0
        #     for m in xrange(M):
        #         tmp += W[m].T.dot(s.diag(tau[m])).dot(Y[m][n,:])
        #     self.Q.mean[n,:] = cov.dot(tmp)

        ## Vectorised ##

        # covariance
        cov = s.eye(K)
        for m in xrange(M):
            cov += (tau[m][:,None,None] * WW[m]).sum(axis=0)
        cov = linalg.inv(cov)
        self.Q.cov = s.repeat(cov[None,:,:],N,axis=0)

        # mean
        mean = s.zeros((K,N))
        for m in xrange(M):
            mean += s.dot( W[m].T, (tau[m]*Y[m]).T )
        self.Q.mean = cov.dot(mean).T

        ## Concatenated ##

        # W = s.concatenate([net.nodes["W"].Q[m].E1 for m in xrange(M)],axis=0)
        # WW = s.concatenate([net.nodes["W"].Q[m].E2 for m in xrange(M)],axis=0)
        # tau = s.concatenate([net.nodes["tau"].Q[m].E for m in xrange(M)],axis=0)
        # Y = s.concatenate([net.nodes["Y"].obs[m] for m in xrange(M)],axis=1)

        # self.cov = linalg.inv(s.eye(K) + s.tensordot(tau,WW,[0,0]))
        # # self.mean = self.cov.dot( W.T ).dot( s.diag(tau) ).dot( (Y-mu).T ).T 
        # self.mean = self.cov.dot(W.T).dot(ddot(tau,Y.T,left=True)).T

        pass
        

    def elbo(self, net):
        K = net.dim['K']
        N = net.dim['N']

        lb_p = -s.trace(self.Q.E2.sum(axis=0))/2
        # lb_q = -logdet(self.Q.cov).sum() - N*K/2
        lb_q = -N*logdet(self.Q.cov[0,:,:]) - N*K/2

        return lb_p - lb_q

class Tau_Node(Gamma_Unobserved_Multiview_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        Gamma_Unobserved_Multiview_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        self.lbconst = s.sum(self.dim*(self.P.a*s.log(self.P.b) - special.gammaln(self.P.a)))

    def updateParameters(self, net):
        D = net.dim["D"]
        M = net.dim["M"]
        K = net.dim["K"]
        N = net.dim["N"]

        # Vectorised
        Z = net.nodes["Z"].Q.E1
        ZZ = (net.nodes["Z"].Q.E2).sum(axis=0)
        for m in xrange(M):
            W = net.nodes["W"].Q[m].E1
            WW = net.nodes["W"].Q[m].E2
            Y = net.nodes["Y"].obs[m]
            self.Q[m].a[:] = self.P.a + N/2
            tmp = (Y**2).sum(axis=0) - 2*(Y*s.dot(Z,W.T)).sum(axis=0) + (WW*ZZ[None,:,:]).sum(axis=(1,2))
            self.Q[m].b = self.P.b + tmp/2

        # Non-vectorised
        # PROBLEM: I THINK THIS IS NOT WORKING
        # Z = net.nodes["Z"].Q.E1
        # for m in xrange(M):
        #     W = net.nodes["W"].Q[m].E1
        #     Y = net.nodes["Y"].obs[m]
        #     for d in xrange(D[m]):
        #         tmp = 0
        #         ww = net.nodes["W"].Q[m].E2[d,:,:]
        #         # ww = s.outer(W[d,:],W[d,:]) + net.nodes["W"].Q[m].cov[d,:,:]
        #         for n in xrange(N):
        #             zz = net.nodes["Z"].Q.E2[n,:,:]
        #             # zz = s.outer(Z[n,:],Z[n,:]) + net.nodes["Z"].Q.cov[n,:,:]
        #             tmp += Y[n,d]**2 - 2*Y[n,d]*W[d,:].dot(Z[n,:]) + s.trace(ww*zz)
        #             # assert tmp > 0 , "This can never be smaller than zero"
        #             # if tmp < 0:
        #             #     print Y[n,d]**2
        #             #     print s.trace(WW[d,:,:]*zz)
        #             #     print -2*Y[n,d]*W[d,:].dot(Z[n,:])
        #             #     print "Terms:"
        #             #     print Y[n,d]
        #             #     print W[d,:].dot(Z[n,:])
        #             #     print W[d,:]
        #             #     print Z[n,:]
        #             #     exit()
        #         self.Q[m].a[d] = self.P.a + N/2
        #         self.Q[m].b[d] = self.P.b + tmp/2

        pass

    def elbo(self, net):
        # Calculate Variational Evidence Lower Bound
        D = net.dim["D"]
        M = net.dim["M"]

        p = self.P
        q = self.Q

        # Fully vectorised
        # lb_p = s.sum(D*(p.a*s.log(p.b) - special.gammaln(p.a))) + (p.a-1)*sum(map(lambda x:x.lnE.sum(), q)) - p.b*sum(map(lambda x: x.E.sum(), q))
        lb_p = self.lbconst + (p.a-1)*sum(map(lambda x:x.lnE.sum(), q)) - p.b*sum(map(lambda x: x.E.sum(), q))
        lb_q = sum(map(lambda q: s.sum(q.a*s.log(q.b)), q)) + sum(map(lambda q: s.sum((q.a-1)*q.lnE), q)) - sum(map(lambda q: s.sum(q.b*q.E), q)) - sum(map(lambda q: s.sum(special.gammaln(q.a)), q)) 

        # Non-vectorised
        # lb_p = 0
        # lb_q = 0
        # for m in xrange(M):
        #     for d in xrange(D[m]): 
        #         lb_p += p.a*s.log(p.b) + (p.a-1)*q[m].lnE[d] - p.b*q[m].E[d] - special.gammaln(p.a)
        #         lb_q += q[m].a[d]*s.log(q[m].b[d]) + (q[m].a[d]-1)*q[m].lnE[d] - q[m].b[d]*q[m].E[d] - special.gammaln(q[m].a[d])


        return lb_p - lb_q

class Alpha_Node(Gamma_Unobserved_Multiview_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        Gamma_Unobserved_Multiview_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        self.lbconst = self.M*self.dim[0][0]*(self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )

    def updateParameters(self, net):
        D = net.dim["D"]
        M = net.dim["M"]
        K = net.dim["K"]

        # Non-vectorised form
        # for m in xrange(M):
        #     WW = net.nodes["W"].Q[m].E2
        #     for k in xrange(K):
        #         self.Q[m].a[k] = self.P.a + D[m]/2
        #         self.Q[m].b[k] = self.P.b + WW.sum(axis=0)[k,k]/2

        # Vectorised form
        for m in xrange(M):
            WW = net.nodes["W"].Q[m].E2
            self.Q[m].a[:] = self.P.a + D[m]/2
            self.Q[m].b = self.P.b + s.diag(WW.sum(axis=0))/2

    def elbo(self, net):
        # Calculate Variational Evidence Lower Bound
        M = net.dim["M"]
        K = net.dim["K"]

        p = self.P
        q = self.Q

        # fully vectorised
        lb_p = self.lbconst + (p.a-1)*sum(map(lambda x:x.lnE.sum(), q)) - p.b*sum(map(lambda x: x.E.sum(), q))
        # lb_p = M*K*(p.a*s.log(p.b) - special.gammaln(p.a)) + (p.a-1)*sum(map(lambda x:x.lnE.sum(), q)) - p.b*sum(map(lambda x: x.E.sum(), q))
        lb_q = sum(map(lambda q: s.sum(q.a*s.log(q.b)), q)) + sum(map(lambda q: s.sum((q.a-1)*q.lnE), q)) - sum(map(lambda q: s.sum(q.b*q.E), q)) - sum(map(lambda q: s.sum(special.gammaln(q.a)), q)) 

        # non-vectorised
        # lb_p = 0
        # lb_q = 0
        # for m in xrange(M):
        #     for k in xrange(K): 
        #         lb_p += p.a*s.log(p.b) + (p.a-1)*q[m].lnE[k] - p.b*q[m].E[k] - special.gammaln(p.a)
        #         lb_q += q[m].a[k]*s.log(q[m].b[k]) + (q[m].a[k]-1)*q[m].lnE[k] - q[m].b[k]*q[m].E[k] - special.gammaln(q[m].a[k]) 

        return lb_p - lb_q

class W_Node(MultivariateGaussian_Unobserved_Multiview_Node):
    # def __init__(self, dim, pmean, pcov, qmean, qcov, qE1=None, qE2=None):
    def __init__(self, dim, qmean, qcov, qE1=None, qE2=None):
        MultivariateGaussian_Unobserved_Multiview_Node.__init__(self, dim=dim, qmean=qmean, qcov=qcov, qE1=qE1)
        # MultivariateGaussian_Unobserved_Multiview_Node.__init__(self, dim=dim, pmean=pmean, pcov=pcov, qmean=qmean, qcov=qcov, qE1=qE1)

    def updateParameters(self, net):
        D = net.dim["D"]
        M = net.dim["M"]
        K = net.dim["K"]
        N = net.dim["N"]

        ## Vectorised ##
        Z = net.nodes["Z"].Q.E1
        ZZ = (net.nodes["Z"].Q.E2).sum(axis=0)
        for m in xrange(M):
            alpha = net.nodes["alpha"].Q[m].E 
            tau = (net.nodes["tau"].Q[m].E)[:,None,None]
            Y = net.nodes["Y"].obs[m]
            self.Q[m].cov = linalg.inv(tau*s.repeat(ZZ[None,:,:],D[m],0) + s.diag(alpha))
            tmp1 = tau*self.Q[m].cov
            tmp2 = Y.T.dot(Z)
            self.Q[m].mean = (tmp1[:,:,:]*tmp2[:,None,:]).sum(axis=2)

        ## Non-Vectorised ##
        # Z = net.nodes["Z"].Q.E1
        # ZZ = (net.nodes["Z"].Q.E2).sum(axis=0)
        # for m in xrange(M):
        #     alpha = net.nodes["alpha"].Q[m].E 
        #     tau = net.nodes["tau"].Q[m].E
        #     Y = net.nodes["Y"].obs[m]
        #     for d in xrange(D[m]):
        #         self.Q[m].cov[d,:,:] = linalg.inv(tau[d]*ZZ + s.diag(alpha))
        #         self.Q[m].mean[d,:] = tau[d]*self.Q[m].cov[d,:,:].dot(s.dot(Y[:,d].T,Z) )

        pass


    def elbo(self, net):
        M = net.dim["M"]
        D = net.dim["D"]
        K = net.dim["K"]

        # Fully vectorised
        lb_p = 0.5*s.sum(D*map(lambda alpha: s.sum(alpha.lnE), net.nodes["alpha"].Q)) - \
            sum(map(lambda alpha,q: s.sum(q.E2 * s.diag(alpha.E)) , net.nodes["alpha"].Q, self.Q))   
        lb_q = -s.sum(D*K) -sum(map(lambda q: logdet(q.cov).sum(),self.Q))

        # Partially Vectorised
        # lb_p = 0
        # lb_q = -s.sum(D*K)
        # for m in xrange(M):
        #     alpha = net.nodes["alpha"].Q[m]
        #     lb_p += D[m]*s.sum(alpha.lnE)/2 - s.sum(self.Q[m].E2 * s.diag(alpha.E))
        #     lb_q -= logdet(self.Q[m].cov).sum()

        # Non-vectorised
        # lb_p = 0.
        # lb_q = 0.
        # for m in xrange(M):
        #     alpha = net.nodes["alpha"].Q
        #     for d in xrange(D[m]):
        #         lb_p += s.sum(alpha[m].lnE)/2 - s.sum(alpha[m].E*s.diag(self.Q[m].E2[d,:,:]))
        #         lb_q -= K + logdet(self.Q[m].cov[d,:,:])
        # print lb_q
        return lb_p - lb_q/2
