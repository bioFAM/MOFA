
from __future__ import division
import scipy as s
from scipy.stats import bernoulli, norm, gamma, uniform, poisson, binom
from random import sample
from utils import sigmoid

"""
Script to simulate data to test the group factor analysis models
"""

class Simulate(object):
    def __init__(self, M, N, D, K):

        # Sanity checks
        assert len(D) == M
        assert K < min(D)
        assert K < N

        self.M = M
        self.N = N
        self.K = K
        self.D = D

    def initAlpha(self, alpha=None):
        if alpha is None:
            alpha = bernoulli.rvs(p=0.5, size=self.M*self.K).reshape((self.M,self.K))
            alpha[alpha==1] = 1.
            alpha[alpha==0] = 1E5
        else:
            assert (alpha.shape[0] == self.M) and (alpha.shape[1] == self.K)
        return alpha

    def initW_ard(self, alpha=None):
        if alpha is None:
            alpha = self.initAlpha()
        W = [ s.zeros((self.D[m],self.K)) for m in xrange(self.M) ]
        for m in xrange(self.M):
            for k in xrange(self.K):
                W[m][:,k] = norm.rvs(loc=0, scale=1/s.sqrt(alpha[m,k]), size=self.D[m])
        return W,alpha

    def initW_spikeslab(self, theta, alpha=None):
        S = [ s.zeros((self.D[m],self.K)) for m in xrange(self.M) ]
        for m in xrange(self.M):
            for k in xrange(self.K):
                S[m][:,k] = bernoulli.rvs(p=theta[m][k], size=self.D[m])

        if alpha is None:
            alpha = self.initAlpha()

        W_hat = [ s.empty((self.D[m],self.K)) for m in xrange(self.M) ]
        W = [ s.empty((self.D[m],self.K)) for m in xrange(self.M) ]
        for m in xrange(self.M):
            for k in xrange(self.K):
                W_hat[m][:,k] = norm.rvs(loc=0, scale=s.sqrt(1/alpha[m,k]), size=self.D[m])
            W[m] = W_hat[m] * S[m]
        return S,W,W_hat,alpha

    def initZ(self, Z=None):
        if Z is None:
            Z = s.empty((self.N,self.K))
            for n in xrange(self.N):
                for k in xrange(self.K):
                    Z[n,k] = norm.rvs(loc=0, scale=1, size=1)
        return Z

    def initTau(self, tau=None):
        if tau is None:
            tau = [ uniform.rvs(loc=1,scale=3,size=self.D[m]) for m in xrange(self.M) ]
        else:
            assert len(tau) == self.M
            for m in range(self.M): assert tau[m].shape[0] == self.D[m]
        return tau

    def initMu(self, mu=None):
        if mu is None:
            mu = [ s.zeros(self.D[m]) for m in xrange(self.M) ]
        else:
            assert len(mu) == self.M
            for m in xrange(self.M): assert mu[m].shape[0] == self.D[m]
        return mu

    def generateData(self, W, Z, Tau, Mu, likelihood, min_trials=None, max_trials=None):
        Y = [ s.zeros((self.N,self.D[m])) for m in xrange(self.M) ]

        # Sample observations using a gaussian likelihood
        if likelihood == "gaussian":
            # Fast way
            # for m in xrange(self.M):
                # Y[m] = s.dot(Z,W[m].T) + Mu[m] + norm.rvs(loc=0, scale=1/s.sqrt(Tau[m])).T
            for m in xrange(self.M):
                for n in xrange(self.N):
                    for d in xrange(self.D[m]):
                        Y[m][n,d] = s.dot(Z[n,:],W[m][d,:].T) + Mu[m][d] + norm.rvs(loc=0,scale=1/s.sqrt(Tau[m][d]))

        # Sample observations using a poisson likelihood
        elif likelihood == "poisson":
            # Slow way 
            for m in xrange(self.M):
                for n in xrange(self.N):
                    for d in xrange(self.D[m]):
                        f = s.dot(Z[n,:],W[m][d,:].T)
                        # f = s.dot(Z[n,:],W[m][d,:].T) + norm.rvs(loc=0,scale=s.sqrt(1/Tau[m][d]))
                        rate = s.log(1+s.exp(f))
                        # Sample from the Poisson distribution
                        # Y[m][n,d] = poisson.rvs(rate)
                        # Use the more likely values
                        Y[m][n,d] = s.special.round(rate)

            # # Fast way
            # for m in xrange(self.M):
            #     F = s.dot(Z,W[m].T)
            #     # F = s.dot(Z,W[m].T) + norm.rvs(loc=0,scale=s.sqrt(1/Tau[m]))
            #     rate = s.log(1+s.exp(F))
            #     # Sample from the Poisson distribution
            #     # MAYBE THIS REQUIRES RESHAPING
            #     # Y[m] = poisson.rvs(rate)
            #     # Use the more likely values
            #     Y[m] = s.special.round(rate)
                
        # Sample observations using a bernoulli likelihood
        elif likelihood == "bernoulli":
            # Slow way 
            for m in xrange(self.M):
                for n in xrange(self.N):
                    for d in xrange(self.D[m]):
                        # Without noise
                        f = sigmoid( s.dot(Z[n,:],W[m][d,:].T) )
                        # With noise, problem: it shifts the sigmoid...
                        # f = sigmoid( s.dot(Z[n,:],W[m][d,:].T) + norm.rvs(loc=0,scale=1/s.sqrt(Tau[m][d])) )

                        # Sample from the Bernoulli distributionn
                        # Y[m][n,d] = bernoulli.rvs(f)
                        # Use the more likely state
                        Y[m][n,d] = s.special.round(f)

        elif likelihood == "binomial":
            Y = dict(tot=[s.zeros((self.N,self.D[m])) for m in xrange(self.M)], 
                     obs=[s.zeros((self.N,self.D[m])) for m in xrange(self.M)] )
            # Slow way 
            for m in xrange(self.M):
                for n in xrange(self.N):
                    for d in xrange(self.D[m]):
                        # Sample the total number of trials
                        Y["tot"][m][n,d] = s.random.random_integers(low=min_trials, high=max_trials, size=1)
                        # Sample the total number of successes
                        f = sigmoid( s.dot(Z[n,:],W[m][d,:].T) )
                        Y["obs"][m][n,d] = binom.rvs(Y["tot"][m][n,d], f) 

        return Y


