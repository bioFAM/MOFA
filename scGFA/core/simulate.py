from __future__ import division
import scipy as s
from scipy.stats import bernoulli, norm, gamma, uniform, poisson, binom
from random import sample
import numpy.ma as ma
import warnings

from utils import sigmoid

"""
Module to simulate test data

To-do:
- Currently non-gaussian likelihoods have no noise
"""

class Simulate(object):
    def __init__(self, M, N, D, K):
        # M (int): number of views
        # N (int): number of samples
        # D (list/tuple of length M): dimensionality of each view
        # K (int): number of latent variables

        # Sanity checks
        assert len(D) == M
        assert K < min(D)
        assert K < N

        self.M = M
        self.N = N
        self.K = K
        self.D = D

    def initAlpha(self):
        # ARD precision is initialised randomly using a Bernoulli distribution with p=0.5
        alpha = [ s.zeros(self.K,) for m in xrange(self.M) ]
        for m in xrange(self.M):
            tmp = bernoulli.rvs(p=0.5, size=self.K)
            tmp[tmp==1] = 1.
            tmp[tmp==0] = 1E5
            alpha[m] = tmp
        return alpha

    def initW_ard(self, alpha=None):
        # ARD weights are initialised
        if alpha is None:
            alpha = self.initAlpha()
        W = [ s.zeros((self.D[m],self.K)) for m in xrange(self.M) ]
        for m in xrange(self.M):
            for k in xrange(self.K):
                W[m][:,k] = norm.rvs(loc=0, scale=1/s.sqrt(alpha[m][k]), size=self.D[m])
        return W,alpha

    # TODO: here is a quick hack for taking into account informative priors
    # need to work out how to do this better
    def initW_spikeslab(self, theta, alpha=None, annotation=False):
        # checking there is no zero in alpha input
        if alpha is not None:
            assert not any([0 in toto for toto in alpha]), 'alpha cannot be zero'

        # Simulate bernoulli variable S
        S = [ s.zeros((self.D[m],self.K)) for m in xrange(self.M) ]
        if annotation:
            # if annotation is True, theta is an informative prior and its dimensions
            # are M * K * D[m], so theta[m][k] is already of length D
            for m in xrange(self.M):
                for k in xrange(self.K):
                    # S[m][:,k] = bernoulli.rvs(p=theta[m][:, k])
                    S[m][:,k] = (theta[m][:, k] > .7) *1.

        else:
            for m in xrange(self.M):
                for k in xrange(self.K):
                    S[m][:,k] = bernoulli.rvs(p=theta[m][k], size=self.D[m])

        # Simualte ARD precision
        if alpha is None:
            alpha = self.initAlpha()

        # Simulate gaussian weights
        W_hat = [ s.empty((self.D[m],self.K)) for m in xrange(self.M) ]
        W = [ s.empty((self.D[m],self.K)) for m in xrange(self.M) ]
        for m in xrange(self.M):
            for k in xrange(self.K):
                W_hat[m][:,k] = norm.rvs(loc=0, scale=s.sqrt(1/alpha[m][k]), size=self.D[m])
            W[m] = W_hat[m] * S[m]
        return S,W,W_hat,alpha

    def initZ(self):
        # Latent variables are initialised by default using a spherical gaussian distribution
        Z = s.empty((self.N,self.K))
        for n in xrange(self.N):
            for k in xrange(self.K):
                Z[n,k] = norm.rvs(loc=0, scale=1, size=1)
        return Z

    def initTau(self):
        # Precision of noise is initialised by default using a uniform distribution
        return [ uniform.rvs(loc=1,scale=3,size=self.D[m]) for m in xrange(self.M) ]

    def initMu(self, mu=None):
        # Means are initialised to zero by default
        return [ s.zeros(self.D[m]) for m in xrange(self.M) ]

    def generateData(self, W, Z, Tau, Mu, likelihood, min_trials=None, max_trials=None, missingness=0.0):
        # W (list of length M where each element is a np array with shape (Dm,K)): weights
        # Z (np array with shape (N,K): latent variables
        # Tau (list of length M where each element is a np array with shape (Dm,)): precision of the normally-distributed noise
        # Mu (list of length M where each element is a np array with shape (Dm,)): feature-wise means
        # likelihood (str): type of likelihood
        # min_trials (int): only for binomial likelihood, minimum number of total trials
        # max_trials (int): only for binomial likelihood, maximum number of total trials
        # missingness (float): percentage of missing values

        Y = [ s.zeros((self.N,self.D[m])) for m in xrange(self.M) ]

        # Sample observations using a gaussian likelihood
        if likelihood == "gaussian":
            for m in xrange(self.M):
                Y[m] = s.dot(Z,W[m].T) + Mu[m] + norm.rvs(loc=0, scale=1/s.sqrt(Tau[m]), size=[self.N, self.D[m]])

            # Non-vectorised, slow
            # for m in xrange(self.M):
            #     for n in xrange(self.N):
            #         for d in xrange(self.D[m]):
            #             Y[m][n,d] = s.dot(Z[n,:],W[m][d,:].T) + Mu[m][d] + norm.rvs(loc=0,scale=1/s.sqrt(Tau[m][d]))

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
                # for n in xrange(self.N):
                    # for d in xrange(self.D[m]):
                        # Without noise
                        # f = sigmoid( s.dot(Z[n,:],W[m][d,:].T) )
                        # With noise, problem: it shifts the sigmoid...
                        # f = sigmoid( s.dot(Z[n,:],W[m][d,:].T) + norm.rvs(loc=0,scale=1/s.sqrt(Tau[m][d])) )

                        # Sample from the Bernoulli distributionn
                        # Y[m][n,d] = bernoulli.rvs(f)
                        # Use the more likely state
                        # Y[m][n,d] = s.special.round(f)
                f = sigmoid( s.dot(Z,W[m].T) )
                Y[m] = s.special.round(f)
        # Sample observations using a binomial likelihood
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

        # Introduce missing values into the data
        # DOESNT WORK FOR BINOMIAL RIGHT NOW
        if missingness > 0.0:
            for m in xrange(self.M):
                nas = s.random.randint(0, self.N*self.D[m], missingness*self.N*self.D[m])
                tmp = Y[m].flatten()
                tmp[nas] = s.nan
                Y[m] = tmp.reshape((self.N,self.D[m]))

        # Create a mask
        for m in xrange(self.M):
            Y[m] = ma.masked_invalid(Y[m])
        return Y
