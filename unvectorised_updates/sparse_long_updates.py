from __future__ import division
from time import time
import numpy.linalg  as linalg

from multiview_nodes import *
from singleview_nodes import UnivariateGaussian_Unobserved_Node
from utils import *

"""
This module is used to define the variational nodes and the corresponding updates for
the element-wise spike and slab Group Factor Analysis.
Two layers of sparsity:
(1) Element-wise spike and slab
(2) Group-wise ARD

TO-DO:
- Currently in the updates and lower bound I am assuming some prior distributions as fixed, there is no flexibility.
- Define the lower bound in a more general form
- I have to define every time the required dimensionalities before doing calculations, would be nice to define a general class
to store them. 

"""

class Y_Node(Multiview_Observed_Node):
    def __init__(self, dim, obs):
        self.dim = dim
        self.obs = obs

    def elbo(self, net):
        M = net.dim['M']
        K = net.dim['K']
        N = net.dim['N']
        D = net.dim['D']

        # vectorised
        tau = net.nodes["tau"]
        lik = -s.sum(0.5*N*D*s.log(2*s.pi))
        for m in xrange(M):
            lik += 0.5*N*s.sum(tau.Q[m].lnE) - s.sum(tau.Q[m].E * (tau.Q[m].b-tau.P.b))
        return lik

    # def getParams(self):
        # return dict({'obs': self.obs})

class Z_Node(UnivariateGaussian_Unobserved_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE1=None, qE2=None):
        UnivariateGaussian_Unobserved_Node.__init__(self, dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE1=qE1)

    def removeFactors(self, *args):
        # Method to remove a set of (inactive) latent variables from the node
        K = self.Q.mean.shape[1]
        keep = s.setdiff1d(s.arange(K),args)
        self.Q.mean = self.Q.mean[:,keep].copy()
        self.Q.var = self.Q.var[:,keep].copy()
        self.Q.E1 = self.Q.E1[:,keep].copy()
        self.Q.E2 = self.Q.E2[:,keep].copy()
        self.dim = (self.dim[0],len(keep))


    def updateParameters(self, net):
        # Method to update the parameters of the Q distribution of the node Z
        # Inputs:
        # - net: instance of BayesNet() containing the entire bayesian betwork

        M = net.dim['M']
        K = net.dim['K']
        N = net.dim['N']
        D = net.dim['D']

        ## non-vectorised ##
        # SW = [ net.nodes["SW"].Q[m].ESW for m in xrange(M) ]
        # SWW = [ net.nodes["SW"].Q[m].ESWW for m in xrange(M) ]
        # tau = [ net.nodes["tau"].Q[m].E for m in xrange(M) ]
        # Y = [ net.nodes["Y"].obs[m] for m in xrange(M) ]
        # for n in xrange(N):
        #     for k in xrange(K):
        #         # Variance
        #         tmp = 0
        #         for m in xrange(M):
        #             for d in xrange(D[m]):
        #                 tmp += tau[m][d]*SWW[m][d,k]
        #         self.Q.var[n,k] = 1/(tmp + 1)
        #         # Mean
        #         tmp = 0
        #         for m in xrange(M):
        #             for d in xrange(D[m]):
        #                 tmp += tau[m][d]*SW[m][d,k] * (Y[m][n,d] - s.sum(SW[m][d,s.arange(K)!=k]*self.Q.mean[n,s.arange(K)!=k]))
        #         self.Q.mean[n,k] = self.Q.var[n,k]*tmp


        ## vectorised ##
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
            tmp2 = Y - s.dot( self.Q.mean[:,s.arange(self.K)!=k] , SW[:,s.arange(self.K)!=k].T )
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

    def elbo(self, net):
        # TO-DO: INCLUDE PRIOR DISTRIBUTION
        K = net.dim['K']
        N = net.dim['N']

        lb_p = -self.Q.E2.sum()/2
        lb_q = -(s.sum(s.log(self.Q.var)) + N*K)/2

        return lb_p - lb_q

class Tau_Node(Gamma_Unobserved_Multiview_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        Gamma_Unobserved_Multiview_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

    def updateParameters(self, net):
        D = net.dim["D"]
        M = net.dim["M"]
        K = net.dim["K"]
        N = net.dim["N"]

        # Non-vectorised
        # Z = net.nodes["Z"].Q.E1
        # ZZ = net.nodes["Z"].Q.E2
        # a = [s.zeros(D[m]) for m in xrange(M)]
        # for m in xrange(M):
        #     SW = net.nodes["SW"].Q[m].ESW
        #     SWW = net.nodes["SW"].Q[m].ESWW
        #     Y = net.nodes["Y"].obs[m]
        #     for d in xrange(D[m]):
        #         tmp = 0
        #         for n in xrange(N):
        #             tmptmp = 0
        #             for k in xrange(K):
        #                 tmptmp += SWW[d,k]*ZZ[n,k]
        #                 if k < (K-1):
        #                     for j in xrange(k+1,K):
        #                         tmptmp += 2*(SW[d,k]*Z[n,k])*(SW[d,j]*Z[n,j])
        #             tmp += Y[n,d]**2 - 2*Y[n,d]*s.sum(SW[d,:]*Z[n,:]) + tmptmp
        #         self.Q[m].a[d] = self.P.a + N/2
        #         self.Q[m].b[d] = self.P.b + tmp/2


        # vectorised and concatenated
        Y = s.concatenate([net.nodes["Y"].obs[m] for m in xrange(M)],axis=1)
        SW = s.concatenate([net.nodes["SW"].Q[m].ESW for m in xrange(M)],axis=0)
        SWW = s.concatenate([net.nodes["SW"].Q[m].ESWW for m in xrange(M)],axis=0)
        Z = net.nodes["Z"].Q.E1
        ZZ = net.nodes["Z"].Q.E2

        term1 = (Y**2).sum(axis=0)
        term2 = 2*(Y*Z.dot(SW.T)).sum(axis=0)
        term3 = (ZZ.dot(SWW.T)).sum(axis=0)
        term4 = s.diag(s.dot( SW.dot(Z.T), Z.dot(SW.T) )) - s.dot(Z**2,(SW**2).T).sum(axis=0)
        tmp = s.split( (term1 - term2 + term3 + term4) , s.cumsum(D))

        for m in xrange(M):
            self.Q[m].a = self.P.a + N/2
            self.Q[m].b = self.P.b + tmp[m]/2

        pass

    def elbo(self, net):
        # Calculate Variational Evidence Lower Bound
        D = net.dim["D"]
        M = net.dim["M"]

        p = self.P
        q = self.Q

        # Fully vectorised
        lb_p = s.sum(D*(p.a*s.log(p.b) - special.gammaln(p.a))) + (p.a-1)*sum(map(lambda x:x.lnE.sum(), q)) - p.b*sum(map(lambda x: x.E.sum(), q))
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
    def __init__(self, dim=[(1,)], pa=[1E-5], pb=[1E-5], qa=[1.], qb=[1.], qE=None):
    	Gamma_Unobserved_Multiview_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

    # def removeFactors(self, *args):
    #     # Method to remove a set of (inactive) latent variables from the node
    #     K = self.dim[0][0]
    #     keep = s.setdiff1d(s.arange(K),args)
    #     for m in xrange(self.M):
    #         self.Q[m].a = self.Q[m].a[keep]
    #         self.Q[m].b = self.Q[m].b[keep]
    #         self.updateExpectations()
    #         # self.Q[m].E = self.Q[m].E[keep]
    #         # self.Q[m].lnE = self.Q[m].lnE[keep]
    #         self.dim[m] = (len(keep),1)


    def updateParameters(self, net):
        
        D = net.dim["D"]
        M = net.dim["M"]
        K = net.dim["K"]


        # RICARD
        # for m in xrange(M):
        #     WW = net.nodes["SW"].Q[m].EWW
        #     for k in xrange(K):
        #         self.Q[m].a[k] = self.P.a + D[m]/2
        #         self.Q[m].b[k] = self.P.b + s.sum(WW[:,k])/2

        # FLORIAN
        # Non-vectorised form
        # for m in xrange(M):
        #     SWW = net.nodes["SW"].Q[m].ESWW
        #     S = net.nodes["SW"].Q[m].ES
        #     for k in xrange(K):
        #         self.Q[m].a[k] = self.P.a + s.sum(S[:,k])/2
        #         self.Q[m].b[k] = self.P.b + 0.5*s.sum(SWW[:,k])

        # Vectorised form
        for m in xrange(M):
            SWW = net.nodes["SW"].Q[m].ESWW
            S = net.nodes["SW"].Q[m].ES
            self.Q[m].a = self.P.a + D[m]/2
            # self.Q[m].a = self.P.a + S.sum(axis=0)/2
            self.Q[m].b = self.P.b + SWW.sum(axis=0)/2

        
    def elbo(self, net):
        # Calculate Variational Evidence Lower Bound
        M = net.dim["M"]
        K = net.dim["K"]

        p = self.P
        q = self.Q

        # fully vectorised
        lb_p = M*K*(p.a*s.log(p.b) - special.gammaln(p.a)) + (p.a-1)*sum(map(lambda x:x.lnE.sum(), q)) - p.b*sum(map(lambda x: x.E.sum(), q))
        lb_q = sum(map(lambda q: s.sum(q.a*s.log(q.b)), q)) + sum(map(lambda q: s.sum((q.a-1)*q.lnE), q)) - sum(map(lambda q: s.sum(q.b*q.E), q)) - sum(map(lambda q: s.sum(special.gammaln(q.a)), q)) 


        # non-vectorised
        # lb_p = 0
        # lb_q = 0
        # for m in xrange(M):
        #     for k in xrange(K): 
        #         lb_p += p.a*s.log(p.b) + (p.a-1)*q[m].lnE[k] - p.b*q[m].E[k] - special.gammaln(p.a)
        #         lb_q += q[m].a[k]*s.log(q[m].b[k]) + (q[m].a[k]-1)*q[m].lnE[k] - q[m].b[k]*q[m].E[k] - special.gammaln(q[m].a[k]) 

        return lb_p - lb_q

class SW_Node(BernoulliGaussian_Unobserved_Multiview_Node):
    # def __init__(self, dim, W_pmean, W_pvar, W_qmean, W_qvar, S_ptheta, S_qtheta, alpha):
    def __init__(self, dim, W_qmean, W_qvar, S_ptheta, S_qtheta, alpha):
        self.alpha = alpha
        # BernoulliGaussian_Unobserved_Multiview_Node.__init__(self, dim=dim, W_pmean=W_pmean, W_pvar=W_pvar, W_qmean=W_qmean, W_qvar=W_qvar, S_ptheta=S_ptheta, S_qtheta=S_qtheta)
        BernoulliGaussian_Unobserved_Multiview_Node.__init__(self, dim=dim, W_qmean=W_qmean, W_qvar=W_qvar, S_ptheta=S_ptheta, S_qtheta=S_qtheta)

    # def removeFactors(self, *args):
    #     # Method to remove a set of (inactive) latent variables from the node
    #     K = self.dim[0][1]
    #     keep = s.setdiff1d(s.arange(K),args)
    #     for m in xrange(self.M):
    #         self.Q[m].mean = self.Q[m].mean[:,keep]
    #         self.Q[m].var = self.Q[m].var[:,keep]
    #         self.Q[m].theta = self.Q[m].theta[:,keep]
    #         self.updateExpectations()
    #         self.dim[m] = (self.dim[m][0],keep)

    def updateParameters(self, net):
        
        D = net.dim["D"]
        M = net.dim["M"]
        K = net.dim["K"]

        # Non-vectorised form
        # Z = net.nodes["Z"].Q.E1
        # ZZ = net.nodes["Z"].Q.E2
        # for m in xrange(M):
        #     tau = net.nodes["tau"].Q[m].E
        #     Y = net.nodes["Y"].obs[m]
        #     alpha = self.alpha.Q[m].E
        #     SW = self.Q[m].ESW.copy()
        #     for d in xrange(D[m]):
        #         for k in xrange(K):
        #             # Update S
        #             term1 = s.log(self.P_theta/(1-self.P_theta))
        #             term2 = 0.5*s.log(alpha[k]/tau[d])
        #             term3 = 0.5*s.log(s.sum(ZZ[:,k]) + alpha[k]/tau[d])
        #             foo = s.dot(Y[:,d],Z[:,k]) - s.dot(SW[d,s.arange(K)!=k],(Z[:,k]*Z[:,s.arange(K)!=k].T).sum(axis=1))
        #             bar = s.sum(ZZ[:,k]) + alpha[k]/tau[d]
        #             term4 = 0.5*tau[d]*foo**2 / bar
        #             self.Q[m].theta[d,k] = 1/(1+s.exp(-(term1+term2-term3+term4)))

        #             # Update W
        #             self.Q[m].mean[d,k] = foo/bar
        #             self.Q[m].var[d,k] = 1/(tau[d]*bar)

        #             # Update expectations
        #             SW[d,k] = self.Q[m].theta[d,k] * self.Q[m].mean[d,k]


        # Vectorised form over K (I haven't checked that it matchea with the non-vectorised form)
        Z = net.nodes["Z"].Q.E1
        ZZ = net.nodes["Z"].Q.E2
        for m in xrange(M):
            tau = net.nodes["tau"].Q[m].E
            Y = net.nodes["Y"].obs[m]
            alpha = self.alpha.Q[m].E
            SW = self.Q[m].ESW.copy()
            for k in xrange(K):
                # term1 = s.log(self.P.theta/(1-self.P.theta))
                term1 = s.log(self.P_theta/(1-self.P_theta))
                term2 = 0.5*s.log(s.divide(alpha[k],tau))
                term3 = 0.5*s.log(s.sum(ZZ[:,k]) + s.divide(alpha[k],tau))
                term41 = s.dot(Y.T,Z[:,k]) 
                term42 = s.dot( SW[:,s.arange(K)!=k] , (Z[:,k]*Z[:,s.arange(K)!=k].T).sum(axis=1) )                
                term43 = s.sum(ZZ[:,k]) + s.divide(alpha[k],tau)
                term4 = 0.5*tau * s.divide((term41-term42)**2,term43)
                # self.Q[m].theta[:,k]= s.divide(1,1+s.exp(-(term1+term2-term3+term4)))
                self.Q[m].theta[:,k] = 1/(1+s.exp(-(term1+term2-term3+term4)))

                # Update W
                self.Q[m].mean[:,k] = s.divide(term41-term42,term43)
                self.Q[m].var[:,k] = s.divide(1,tau*term43)

                # Update Expectations
                SW[:,k] = self.Q[m].theta[:,k] * self.Q[m].mean[:,k]


        # Approximated Fully vectorised form
        # fINISH THIS
        # Z = net.nodes["Z"].Q.E1
        # ZZ = net.nodes["Z"].Q.E2
        # for m in xrange(M):
        #     tau = net.nodes["tau"].Q[m].E
        #     Y = net.nodes["Y"].obs[m]
        #     alpha = self.alpha.Q[m].E
        #     SW = self.Q[m].ESW
        #     term1 = s.log(self.P.theta/(1-self.P.theta))
        #     term2tmp = s.repeat(alpha[None,:],D[m],0) / s.repeat(tau[:,None],K,1)
        #     term2 = 0.5*s.log(term2tmp)
        #     term3 = 0.5*s.log(ZZ.sum(axis=0) + term2tmp)
        #     term41 = s.dot(Y.T,Z)

        #     term42 = s.dot(SW[:,s.arange(K)!=k],(Z[:,k]*Z[:,s.arange(K)!=k].T).sum(axis=1))

        #     term43 = s.sum(ZZ[:,k]) + alpha[k]/tau
        #     term4 = (term41 + term42)**2 / term43
        #     self.Q[m].theta[:,k] = 1/(1+s.exp(-(term1 + term2 - term3 + 0.5*tau*term4)))

        #     # Update W
        #     self.Q[m].mean[:,k] = (term41+term42)/term43
        #     self.Q[m].var[:,k] = 1/(tau*term43)


        pass

    def updateExpectations(self):
        # Non-vectorised
        # K = self.dim[0][1]
        # D = [ self.dim[m][0] for m in xrange(self.M) ]
        # for m in xrange(self.M):
        #     self.Q[m].ES = s.zeros((D[m],K))
        #     self.Q[m].ESW = s.zeros((D[m],K))
        #     self.Q[m].ESWW = s.zeros((D[m],K))
        #     self.Q[m].EWW = s.zeros((D[m],K))
        #     for d in xrange(D[m]):
        #         for k in xrange(K):
        #             alpha = self.alpha.Q[m].E[k]
        #             self.Q[m].ES[d,k] = self.Q[m].theta[d,k]
        #             self.Q[m].ESW[d,k] = self.Q[m].ES[d,k] * self.Q[m].mean[d,k]
        #             self.Q[m].ESWW[d,k] = self.Q[m].ES[d,k] * (self.Q[m].mean[d,k]**2 + self.Q[m].var[d,k])
        #             self.Q[m].EWW[d,k] = self.Q[m].ES[d,k] * (self.Q[m].mean[d,k]**2 + self.Q[m].var[d,k]) + (1-self.Q[m].ES[d,k])/alpha

        # Vectorised
        for m in xrange(self.M):
            alpha = self.alpha.Q[m].E
            self.Q[m].ES = self.Q[m].theta.copy()
            self.Q[m].EW = self.Q[m].mean.copy()
            self.Q[m].ESW = self.Q[m].ES * self.Q[m].EW
            self.Q[m].ESWW = self.Q[m].ES * (self.Q[m].EW**2 + self.Q[m].var)
            self.Q[m].EWW = self.Q[m].ES * (self.Q[m].EW**2 + self.Q[m].var)  + (1-self.Q[m].ES)*s.repeat(1/alpha[None,:],self.dim[m][0],0)

        pass

    def getExpectations(self):
        ES = [ self.Q[m].ES for m in xrange(self.M) ]
        EW = [ self.Q[m].EW for m in xrange(self.M) ]
        ESW = [ self.Q[m].ESW for m in xrange(self.M) ]
        ESWW = [ self.Q[m].ESWW for m in xrange(self.M) ]
        EWW = [ self.Q[m].EWW for m in xrange(self.M) ]
        return dict({'ES':ES, 'EW':EW, 'ESW':ESW, 'ESWW':ESWW, 'EWW':EWW})

    def elbo(self, net):
        # Calculate Variational Evidence Lower Bound
        M = net.dim["M"]
        K = net.dim["K"]
        D = net.dim["D"]

        # Vectorised
        lb_pw = 0
        lb_qw = 0
        for m in xrange(M):
            alpha = self.alpha.Q[m]
            WW = net.nodes["SW"].Q[m].EWW
            S = self.Q[m].ES
            varW = self.Q[m].var
            lb_pw += 0.5*D[m]*alpha.lnE.sum() - 0.5*s.sum(alpha.E*WW)
            lb_qw += -0.5*K*D[m] - 0.5*s.log(S*varW + ((1-S)/alpha.E)).sum()
        lb_w = lb_pw - lb_qw

        Slower = 0.00001
        Supper = 0.99999
        lb_ps = 0
        lb_qs = 0
        for m in xrange(M):
            S = self.Q[m].ES
            S[S<Slower] = Slower
            S[S>Supper] = Supper
            theta = self.P_theta
            # theta = self.P.theta
            lb_ps += s.sum( S*s.log(theta) + (1-S)*s.log(1-theta) )
            lb_qs += s.sum( S*s.log(S) + (1-S)*s.log(1-S) )
        lb_s = lb_ps - lb_qs

        # Non-vectorised
        # lb_pw = 0
        # lb_qw = 0
        # for m in xrange(M):
        #     alpha = self.alpha.Q[m]
        #     WW = net.nodes["SW"].Q[m].EWW
        #     S = self.Q[m].ES
        #     varW = self.Q[m].var
        #     for d in xrange(D[m]):
        #         for k in xrange(K):
        #             lb_pw +=  -0.5*s.log(2*s.pi) + 0.5*alpha.lnE[k] - 0.5*alpha.E[k]*WW[d,k] 
        #             lb_qw +=  -0.5*(1+s.log(2*s.pi)) - 0.5*s.log(S[d,k]*varW[d,k] + (1-S[d,k])/alpha.E[k])
        # lb_w = lb_pw - lb_qw

        # Slower = 0.00001
        # Supper = 0.99999
        # lb_ps = 0
        # lb_qs = 0
        # for m in xrange(M):
        #     S = self.Q[m].ES
        #     S[S<Slower] = Slower
        #     S[S>Supper] = Supper
        #     theta = self.P.theta
        #     for d in xrange(D[m]):
        #         for k in xrange(K):
        #             lb_ps += S[d,k]*s.log(theta) + (1-S[d,k])*s.log(1-theta)
        #             lb_qs += S[d,k]*s.log(S[d,k]) + (1-S[d,k])*s.log(1-S[d,k])
        # lb_s = lb_ps - lb_qs

        return lb_w + lb_s