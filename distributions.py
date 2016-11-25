from __future__ import division
import scipy as s
import logging as L
import numpy.linalg as linalg
import scipy.special as special
import scipy.stats as stats

import pdb

from utils import *

"""
This module is used to define statistical distributions in a general sense

Each Distribution class can store an arbitrary number of distributions of the same type, this is specified in the 'dim' argument

A Distribution class has two main types of attributes: parameters and expectations. Both have to be defined when initialising a class.
Note that in some distributions (Gaussian mainly) a parameter is equal to an expectation. However, they are stored as separate
attributes and are not always necessarily equal (i.e. VBEM algorithm).

There are cases (mainly with the BernoulliGaussian distribution) where some expectations are quite complex and
require expectations from other distributions. In that case, I prefer to define the expectations together with the updates
of the corresponding node, instead of doing it here.

"""

# General class for probability distributions
class Distribution(object):
    """ Abstract class for a distribution """
    def __init__(self, dim):
        self.dim = dim
    def density(self):
        pass
    def loglik(self):
        pass
    def sample(self):
        pass
    def entropy(self):
        pass
    def updateExpectations(self):
        pass

# Specific classes for probability distributions
class MultivariateGaussian(Distribution):
    """
    Class for a multivariate Gaussian distribution. This class can store N multivate
    Gaussian Distributions of dimensionality D each.
    Here we assume that the data factorises over samples (rows),
    so we have N distributions of dimensionality D each.

    Equations:
    p(X|Mu,Sigma) = 1/(2pi)^{D/2} * 1/(|Sigma|^0.5) * exp( -0.5*(X-Mu)^{T} Sigma^{-1} (X-Mu) )
    log p(X|Mu,Sigma) = -D/2*log(2pi) - 0.5log(|Sigma|) - 0.5*(X-Mu)^{T} Sigma^{-1} (X-Mu)
    E[X] = Mu
    E[X^2] = E[X]E[X] + Cov[X] = MuMu + Sigma
    cov[X] = Sigma
    H[X] = 0.5*log(|Sigma|) + D/2*(1+log(2pi))

    Dimensionalities:
    - X: (N,D)
    - Mu: (N,D)
    - Sigma: (N,D,D)
    - E[X]: (N,D)
    - E[X^2]: (N,D,D)
    """
    def __init__(self, dim, mean, cov, E=None, E2=None):
        Distribution.__init__(self, dim)

        # Check dimensions are correct
        assert len(dim) == 2, "You need to specify two dimensions for 'dim': (number of distributions, dimensionality) "
        assert not (dim[0]==1 and dim[1]==1), "A 1x1 multivariate normal reduces to a Univariate normal "

        ## Initialise the parameters ##

        # Initialise the mean
        # If 'mean' is a scalar, broadcast it to all dimensions
        if isinstance(mean,(int,float)): mean = s.ones( (dim[0],dim[1]) ) * mean
        # If 'mean' has dim (D,) and we have N distributions, broadcast it to all N distributions
        if len(mean.shape)==1 and mean.shape[0]==dim[0]: mean = s.repeat(mean,N,0)
        assert sum(mean.shape) > 2, "The mean has to be a matrix with shape (N,D) "

        # Initialise the covariance
        # If 'cov' is a matrix and not a tensor, broadcast it along the zeroth axis
        if len(cov.shape) == 2: cov = s.repeat(cov[None,:,:],dim[0],0)
        assert (cov.shape[1]==cov.shape[2]) and (sum(cov.shape[1:])>1), "The covariance has to be a tensor with shape (N,D,D)"

        # Check that the dimensionalities of 'mean' and 'cov' match
        assert cov.shape[1] == mean.shape[1] == dim[1], "Error in the dimensionalities"
        assert cov.shape[0] == mean.shape[0] == dim[0], "Error in the dimensionalities"

        self.mean = mean
        self.cov = cov

        ## Initialise expectations ##
        if E is not None: self.E = E
        if E2 is not None: self.E2 = E2

    def updateExpectations(self):
        # Update first and second moments using current parameters
        self.E = self.mean.copy()

        # self.E2 = s.empty( (self.dim[0],self.dim[1],self.dim[1]) )
        # for i in xrange(self.dim[0]):
        #     self.E2[i,:,:] = s.outer(self.E[i,:],self.E[i,:]) + self.cov[i,:,:]

        self.E2 = self.cov.copy()
        for i in xrange(self.dim[0]):
            self.E2[i,:,:] += s.outer(self.E[i,:],self.E[i,:])
        pass
    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        return s.sum( stats.multivariate_normal.pdf(x, mean=self.mean[n,:], cov=self.cov[n,:,:]) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        l = 0.
        D = self.dim[1]
        for n in xrange(self.dim[0]):
            qterm = (x[n,:]-self.mean[n,:]).T.dot(linalg.det(self.cov[n,:,:])).dot(x[n,:]-self.mean[n,:])
            l += -0.5*D*s.log(2*s.pi) - 0.5*s.log(linalg.det(self.cov[n,:,:])) -0.5*qterm
        return l
        # return s.sum( s.log(stats.multivariate_normal.pdf(x, mean=self.mean[n,:], cov=self.cov[n,:,:])) )

    # def entropy(self):
        # CHECK THIs Is CORRECT
        # tmp = sum( [ logdet(self.cov[i,:,:]) for i in xrange(self.dim[0]) ] )
        # return ( 0.5*(tmp + (self.dim[0]*self.dim[1])*(1+s.log(2*pi)) ).sum() )
class UnivariateGaussian(Distribution):
    """
    Class for an arbitrary number of univariate Gaussian distributions

    Equations:
    Class for a univariate Gaussian distributed node
    p(x|mu,sigma^2) = 1/sqrt(2*pi*sigma^2) * exp(-0.5*(x-mu)^2/(sigma^2) )
    log p(x|mu,sigma^2) =
    E[x] = mu
    var[x] = sigma^2
    H[x] = 0.5*log(sigma^2) + 0.5*(1+log(2pi))

    """
    def __init__(self, dim, mean, var, E=None, E2=None):
        Distribution.__init__(self, dim)

        ## Initialise parameters ##
        # If 'mean' or 'var' are scalars, broadcast it to all dimensions
        if isinstance(mean,(int,float)): mean = s.ones(dim) * mean
        self.mean = mean
        if isinstance(var,(int,float)): var = s.ones(dim) * var
        self.var = var

        ## Initialise expectations ##
        # DO IT PROPERLY, COPY FROM GAMMA OR POISSON
        if E is not None: self.E = E
        if E2 is not None: self.E2 = E2

        # Check that dimensionality match
        assert self.mean.shape == self.var.shape == self.dim, "Dimensionalities do not match"

    def updateExpectations(self):
        # Update first and second moments using current parameters
        self.E = self.mean.copy()
        self.E2 = self.E**2 + self.var

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        # print stats.norm.pdf(x, loc=self.mean, scale=s.sqrt(self.var))
        return s.sum( (1/s.sqrt(2*s.pi*self.var)) * s.exp(-0.5*(x-self.mean)**2/self.var) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        # return s.log(stats.norm.pdf(x, loc=self.mean, scale=s.sqrt(self.var)))
        return s.sum( -0.5*s.log(2*s.pi) - 0.5*s.log(self.var) -0.5*(x-self.mean)**2/self.var )

    def entropy(self):
        return s.sum( 0.5*s.log(self.var) + 0.5*(1+s.log(2*s.pi)) )


class Gamma(Distribution):
    """
    This class can store an arbitrary number of Gamma distributions

    Equations:
    p(x|a,b) = (1/Gamma(a)) * b^a * x^(a-1) * e^(-b*x)
    log p(x|a,b) = -log(Gamma(a)) + a*log(b) + (a-1)*log(x) - b*x
    E[x] = a/b
    var[x] = a/b^2
    E[ln(x)] = digamma(a) - ln(b)
    H[x] = ln(Gamma(a)) - (a-1)*digamma(a) - ln(b) + a
    """

    def __init__(self, dim=(1,), a=1E-3, b=1E-3, E=None, lnE=None):
        Distribution.__init__(self, dim)

        ## Initialise parameters ##
        # If 'a' or 'b' are scalars, broadcast it to all dimensions
        if isinstance(a,(int,float)): a = s.ones(dim) * a
        if isinstance(b,(float,int)): b = s.ones(dim) * b
        self.a = a
        self.b = b
        # pdb.set_trace()

        ## Initialise expectations ##
        if E is not None:
            if isinstance(E,(float,int)): E = s.ones(dim) * E
        self.E = E

        if lnE is not None:
            if isinstance(lnE,(float,int)): lnE = s.ones(dim) * lnE
        self.lnE = lnE

    def updateExpectations(self):
        self.E = self.a/self.b
        self.lnE = special.digamma(self.a) - s.log(self.b)

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        return s.prod( (1/special.gamma(self.a)) * self.b**self.a * x**(self.a-1) * s.exp(-self.b*x) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        return s.sum( -s.log(special.gamma(self.a)) + self.a*s.log(self.b) * (self.a-1)*s.log(x) -self.b*x )


class Poisson(Distribution):
    """
    Class for a Poisson distribution.
    This class can store an arbitrary number of Poisson distributions

    Equations:
    p(x|theta) = theta**x * exp(-theta) * 1/theta!
    log p(x|a,b) = x*theta - theta - log(x!)
    E[x] = theta
    var[x] = theta
    """

    def __init__(self, dim, theta, E=None):
        Distribution.__init__(self, dim)

        # Initialise parameters
        # If 'theta' is a scalar, broadcast it to all dimensions
        if isinstance(theta,(float,int)): a = s.ones(dim) * theta
        self.theta = theta

        # Initialise expectations
        if E is not None:
            if isinstance(E,(float,int)): E = s.ones(dim) * E
        else:
            E = s.zeros(dim)
        self.E = E

        # Check that dimensionality match
        assert self.theta.shape == self.E.shape == self.dim, "Dimensionalities do not match"

    def updateExpectations(self):
        self.E = self.theta.copy()

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        assert x.dtype == int, "x has to be an integer array"
        theta = self.theta.flatten()
        x = x.flatten()
        # return s.prod (stats.poisson.pmf(x,theta) )
        return s.prod( s.divide(theta**x * s.exp(-theta),s.misc.factorial(x)) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        assert x.dtype == int, "x has to be an integer array"
        theta = self.theta.flatten()
        x = x.flatten()
        # return s.log( s.prod (stats.poisson.pmf(x,theta) ))
        return s.sum( x*s.log(theta) - theta - s.log(s.misc.factorial(x)) )
class Bernoulli(Distribution):
    """
    Class to store an arbitrary number of Bernoulli distributions

    Equations:

    """
    def __init__(self, dim, theta, E=None):
        Distribution.__init__(self, dim)

        ## Initialise parameters ##
        # If 'theta' is a scalar, broadcast it to all dimensions
        if isinstance(theta,(float,int)): theta = s.ones(dim) * theta
        self.theta = theta

        ## Initialise expectations ##
        if E is None:
            E = self.theta
        else:
            if isinstance(E,(float,int)): E = s.ones(dim) * E
        self.E = E

        # Check that dimensionality match
        assert self.theta.shape == self.E.shape == self.dim, "Dimensionalities do not match"

    def updateExpectations(self):
        self.E = self.theta

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        return s.prod( self.theta**x * (1-self.theta)**(1-x) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        return s.sum( x*self.theta + (1-x)*(1-self.theta) )
class BernoulliGaussian(Distribution):
    """
    Class to store an arbitrary number of Bernoulli-Gaussian distributions (see paper Spike and Slab
    Variational Inference for Multi-Task and Multiple Kernel Learning by Titsias and Gredilla)

    Equations:
    p(w,s) = Normal(w|mean,var) * Bernoulli(s|theta)

    This distribution has several expectations that are required for the variational updates, and they depend
    on other parameters (alpha from the ARD prior). For this reason I decided to define the expectations in the
    class of the corresponding node, instead of doing it here.

    """
    def __init__(self, dim, mean, var, theta):
        Distribution.__init__(self,dim)

        ## Initialise parameters ##
        # Broadcast input scalars to all dimensions

        # If 'theta' is a scalar, broadcast it to all dimensions
        if isinstance(theta,(float,int)): theta = s.ones(dim) * theta
        if isinstance(mean,(float,int)): mean = s.ones(dim) * mean
        if isinstance(var,(float,int)): var = s.ones(dim) * var

        self.theta = theta
        self.mean = mean
        self.var = var

    def density(self, x):
        pass

    def loglik(self, x):
        pass
class Binomial(Distribution):
    """
    Class to store an arbitrary number of Binomial distributions

    Equations:
    p(x|N,theta) = binom(N,x) * theta**(x) * (1-theta)**(N-x)
    log p(x|N,theta) = log(binom(N,x)) + x*theta + (N-x)*(1-theta)
    E[x] = N*theta
    var[x] = N*theta*(1-theta)
    """
    def __init__(self, dim, N, theta, E=None):
        Distribution.__init__(self, dim)

        ## Initialise parameters ##
        # Broadcast scalars to all dimensions
        if isinstance(theta,(float,int)): theta = s.ones(dim) * theta
        if isinstance(N,(float,int)): N = s.ones(dim) * N
        # if isinstance(K,(float,int)): K = s.ones(dim) * K
        self.theta = theta
        self.N = N
        # self.K = K

        ## Initialise expectations ##
        if E is None:
            self.updateExpectations()
        else:
            if isinstance(E,(float,int)): E = s.ones(dim) * E
            self.E = E

        # Check that dimensionalities match
        assert self.theta.shape == self.E.shape == self.N.shape == self.dim, "Dimensionalities do not match"

    def updateExpectations(self):
        self.E = self.N * self.theta

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        assert x.dtype == int, "x has to be an integer array"
        # return s.prod( stats.binom.pmf(x, self.N, self.theta) )
        return s.prod( special.binom(self.N,x) * self.theta**x * (1-self.theta)**(self.N-x) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        assert x.dtype == int, "x has to be an integer array"
        # print s.sum (stats.binom.logpmf(x, self.N, self.theta) )
        return s.sum( s.log(special.binom(self.N,x)) + x*s.log(self.theta) + (self.N-x)*s.log(1-self.theta) )
class Beta(Distribution):
    """
    Class to store an arbitrary number of Beta distributions

    Equations:
    p(x|a,b) =
    log p(x|a,b) =
    E[x] = a/(a+b)
    var[x] = a*b / ((a+b)**2 * (a+b+1))
    """
    def __init__(self, dim, a, b, E=None):
        Distribution.__init__(self, dim)

        ## Initialise parameters ##
        # Broadcast scalars to all dimensions
        if isinstance(a,(float,int)): a = s.ones(dim) * a
        if isinstance(b,(float,int)): b = s.ones(dim) * b
        self.a = a
        self.b = b

        ## Initialise expectations ##
        if E is None:
            self.updateExpectations()
        else:
            if isinstance(E,(float,int)): E = s.ones(dim) * E
            self.E = E
        # Check that dimensionalities match
        assert self.a.shape == self.E.shape == self.b.shape == self.dim, "Dimensionalities do not match"

    def updateExpectations(self):
        self.E = s.divide(self.a,self.a+self.b)
        self.lnE = special.digamma(self.a) - special.digamma(self.a + self.b)
        # expectation of ln(1-X)
        self.lnEInv = special.digamma(self.b) - special.digamma(self.a + self.b)

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        assert x.dtype == int, "x has to be an integer array"
        pass
        # return s.prod( stats.binom.pmf(x, self.N, self.theta) )
        # return s.prod( special.binom(self.N,x) * self.theta**x * (1-self.theta)**(self.N-x) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        assert x.dtype == int, "x has to be an integer array"
        pass
        # return s.sum (stats.binom.logpmf(x, self.N, self.theta) )
        # return s.sum( s.log(special.binom(self.N,x)) + x*s.log(self.theta) + (self.N-x)*s.log(1-self.theta) )
