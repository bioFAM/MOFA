from __future__ import division
import scipy as s
import logging as L
import numpy.linalg as linalg
import scipy.special as special
import scipy.stats as stats

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

TO-DO:
- we should pass parametrs and expectations to Distribution() and perform sanity checks there
- Improve initialisation of Multivariate Gaussian
- Sanity checks on setter and getter functions
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

    def getParameters(self):
        # Getter function for parameters
        return self.params

    def setParameters(self,**params):
        # Setter function for parameters
        self.params = params

    def getExpectation(self):
        # Getter function for expectations
        return self.expectations['E']

    def getExpectations(self):
        # Getter function for expectations
        return self.expectations

    def CheckDimensionalities(self):
        # Method to do a sanity check on the dimensionalities
        p_dim = set(map(s.shape, self.params.values()))
        e_dim = set(map(s.shape, self.expectations.values()))
        assert len(p_dim) == 1, "Parameters have different dimensionalities"
        assert len(e_dim) == 1, "Expectations have different dimensionalities"
        assert e_dim == p_dim, "Parameters and Expectations have different dimensionality"

    def removeDimensions(self, axis, idx):
        # Method to remove undesired dimensions
        # - axis (int): axis from where to remove the elements
        # - idx (numpy array): indices of the elements to remove
        assert axis <= len(self.dim)
        assert s.all(idx < self.dim[axis])
        for k in self.params.keys(): self.params[k] = s.delete(self.params[k], idx, axis)
        for k in self.expectations.keys(): self.expectations[k] = s.delete(self.expectations[k], idx, axis)
        self.updateDim(axis=axis, new_dim=self.dim[axis]-len(idx))

    def updateDim(self, axis, new_dim):
        # Function to update the dimensionality of a particular axis
        # this seems inefficient but tuples cannot be modified...
        dim = list(self.dim)
        dim[axis] = new_dim
        self.dim = tuple(dim)

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
    def __init__(self, dim, mean, cov, E=None):
        Distribution.__init__(self, dim)

        # Check dimensions are correct
        assert len(dim) == 2, "You need to specify two dimensions for 'dim': (number of distributions, dimensionality) "
        assert not (dim[0]==1 and dim[1]==1), "A 1x1 multivariate normal reduces to a Univariate normal "

        ## Initialise the parameters ##

        # Initialise the mean
        # If 'mean' is a scalar, broadcast it to all dimensions
        if isinstance(mean,(int,float)): mean = s.ones( (dim[0],dim[1]) ) * mean
        # If 'mean' has dim (D,) and we have N distributions, broadcast it to all N distributions
        if len(mean.shape)==1 and mean.shape[0]==dim[1]: mean = s.repeat(mean,dim[0],0)
        assert sum(mean.shape) > 2, "The mean has to be a matrix with shape (N,D) "

        # Initialise the covariance
        # If 'cov' is a matrix and not a tensor, broadcast it along the zeroth axis
        if len(cov.shape) == 2: cov = s.repeat(cov[None,:,:],dim[0],0)
        assert (cov.shape[1]==cov.shape[2]) and (sum(cov.shape[1:])>1), "The covariance has to be a tensor with shape (N,D,D)"

        # Check that the dimensionalities of 'mean' and 'cov' match
        assert cov.shape[1] == mean.shape[1] == dim[1], "Error in the dimensionalities"
        assert cov.shape[0] == mean.shape[0] == dim[0], "Error in the dimensionalities"

        self.params = {'mean':mean, 'cov':cov }

        # Initialise expectations
        if E is None:
            self.updateExpectations()
        else:
            self.expectations = { 'E':E }

        # TO-DO: SANITY CHECKS ON EXPECTATIONS


    def updateExpectations(self):
        # Update first and second moments using current parameters
        E = self.params['mean']

        # self.E2 = s.empty( (self.dim[0],self.dim[1],self.dim[1]) )
        # for i in xrange(self.dim[0]):
        #     self.E2[i,:,:] = s.outer(self.E[i,:],self.E[i,:]) + self.cov[i,:,:]

        E2 = self.params['cov'].copy()
        for i in xrange(self.dim[0]):
            E2[i,:,:] += s.outer(E[i,:],E[i,:])

        self.expectations = {'E':E, 'E2':E2}

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        return s.sum( stats.multivariate_normal.pdf(x, mean=self.params['mean'][n,:], cov=self.params['cov'][n,:,:]) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        l = 0.
        D = self.dim[1]
        for n in xrange(self.dim[0]):
            qterm = (x[n,:]-self.params['mean'][n,:]).T.dot(linalg.det(self.params['cov'][n,:,:])).dot(x[n,:]-self.params['mean'][n,:])
            l += -0.5*D*s.log(2*s.pi) - 0.5*s.log(linalg.det(self.params['cov'][n,:,:])) -0.5*qterm
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

        # Initialise parameters
        mean = s.ones(dim) * mean
        var = s.ones(dim) * var
        self.params = { 'mean':mean, 'var':var }

        # Initialise expectations 
        self.expectations = {}
        if E is None:
            self.updateExpectations()
        else:
            self.expectations['E'] = s.ones(dim)*E 

        if E2 is not None:
            self.expectations['E2'] = s.ones(dim)*E2 


        # Check that dimensionalities match
        self.CheckDimensionalities()

    def updateExpectations(self):
        # Update first and second moments using current parameters
        E = self.params['mean']
        E2 = E**2 + self.params['var']
        self.expectations = { 'E':E, 'E2':E2 }


    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        # print stats.norm.pdf(x, loc=self.mean, scale=s.sqrt(self.var))
        return s.sum( (1/s.sqrt(2*s.pi*self.params['var'])) * s.exp(-0.5*(x-self.params['mean'])**2/self.params['var']) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        # return s.log(stats.norm.pdf(x, loc=self.mean, scale=s.sqrt(self.var)))
        return s.sum( -0.5*s.log(2*s.pi) - 0.5*s.log(self.params['var']) -0.5*(x-self.params['mean'])**2/self.params['var'] )

    def entropy(self):
        return s.sum( 0.5*s.log(self.params['var']) + 0.5*(1+s.log(2*s.pi)) )
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

    def __init__(self, dim, a, b, E=None):
        Distribution.__init__(self, dim)

        # Initialise parameters
        a = s.ones(dim) * a
        b = s.ones(dim) * b
        self.params = { 'a':a, 'b':b }

        # Initialise expectations
        if E is None:
            self.updateExpectations()
        else:
            self.expectations = { 'E':s.ones(dim)*E }

        # Check that dimensionalities match
        self.CheckDimensionalities()

    def updateExpectations(self):
        E = self.params['a']/self.params['b']
        lnE = special.digamma(self.params['a']) - s.log(self.params['b'])
        self.expectations = { 'E':E, 'lnE':lnE }

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        return s.prod( (1/special.gamma(self.params['a'])) * self.params['b']**self.params['a'] * x**(self.params['a']-1) * s.exp(-self.params['b']*x) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        return s.sum( -s.log(special.gamma(self.params['a'])) + self.params['a']*s.log(self.params['b']) * (self.params['a']-1)*s.log(x) -self.params['b']*x )
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
        theta = s.ones(dim) * theta
        self.params = { 'theta':theta }

        # Initialise expectations
        if E is None:
            self.updateExpectations()
        else:
            self.expectations = { 'E':s.ones(dim)*E }

        # Check that dimensionalities match
        self.CheckDimensionalities()

    def updateExpectations(self):
        E = self.params['theta']
        self.expectations = { 'E':E }

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        assert x.dtype == int, "x has to be an integer array"
        theta = self.params['theta'].flatten()
        x = x.flatten()
        # return s.prod (stats.poisson.pmf(x,theta) )
        return s.prod( s.divide(theta**x * s.exp(-theta),s.misc.factorial(x)) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        assert x.dtype == int, "x has to be an integer array"
        theta = self.params['theta'].flatten()
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

        # Initialise parameters
        theta = s.ones(dim) * theta
        self.params = { 'theta':theta }

        # Initialise expectations
        if E is None:
            self.updateExpectations()
        else:
            self.expectations = { 'E':s.ones(dim)*E }

        # Check that dimensionalities match
        self.CheckDimensionalities()

    def updateExpectations(self):
        E = self.params['theta']
        self.expectations = { 'E':E }

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        return s.prod( self.params['theta']**x * (1-self.params['theta'])**(1-x) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        return s.sum( x*self.params['theta'] + (1-x)*(1-self.params['theta']) )
class BernoulliGaussian(Distribution):
    """
    Class to store an arbitrary number of Bernoulli-Gaussian distributions (see paper Spike and Slab
    Variational Inference for Multi-Task and Multiple Kernel Learning by Titsias and Gredilla)

    The best way to characterise a joint Bernoulli-Gaussian distribution P(w,s) is by considering
    its factorisation p(w|s)p(s) where s is a bernoulli distribution and w|s=0 and w|s=1 are normal distributions

    Equations:
    p(w,s) = Normal(w|mean,var) * Bernoulli(s|theta)
    FINISH EQUATIONS

    ROOM FOR IMPROVEMENT: i think the current code is inefficient because you have to keep track of the params
    and expectations in both the factorised distributions and the joint one. I think
    perhaps is better not to define new Bernoulli and UnivariateGaussian but work directly with the joint model

    """
    def __init__(self, dim, mean_S0, mean_S1, var_S0, var_S1, theta, EW_S0=None, EW_S1=None, ES=None):
        Distribution.__init__(self,dim)
        self.S = Bernoulli(dim=dim, theta=theta, E=ES)
        self.W_S0 = UnivariateGaussian(dim=dim, mean=mean_S0, var=var_S0, E=EW_S0)
        self.W_S1 = UnivariateGaussian(dim=dim, mean=mean_S1, var=var_S1, E=EW_S1)

        # Collect parameters
        self.params = { 'mean_S0':mean_S0, 'mean_S1':mean_S1, 'var_S0':var_S0, 'var_S1':var_S1, 'theta':theta }
        
        # Collect expectations
        self.updateExpectations()

    def setParameters(self,**params):
        # Setter function for parameters
        self.S.setParameters(theta=params['theta'])
        self.W_S0.setParameters(mean=params['mean_S0'], var=params['var_S0'])
        self.W_S1.setParameters(mean=params['mean_S1'], var=params['var_S1'])
        self.params = params

    def updateParameters(self):
        # Method to update the parameters of the joint distribution based on its constituent distributions
        self.params = { 'theta':self.S.params["theta"], 
                        'mean_S0':self.W_S0.params["mean"], 'var_S0':self.W_S0.params["var"],
                        'mean_S1':self.W_S1.params["mean"], 'var_S1':self.W_S1.params["var"] }

    def updateExpectations(self):
        # Method to calculate the expectations based on the current estimates for the parameters

        # Update expectations of the constituent distributions
        self.S.updateExpectations()
        self.W_S0.updateExpectations()
        self.W_S1.updateExpectations()

        # Calculate expectations of the joint distribution
        ES = self.S.getExpectation()
        EW = self.W_S1.getExpectation()
        ESW = ES * EW
        ESWW = ES * (EW**2 + self.params["var_S1"])
        EWW = ES*(EW**2+self.params["var_S1"]) + (1-ES)*self.params["var_S0"]

        # Collect expectations
        self.expectations = {'E':ESW, 'ES':ES, 'EW':EW, 'ESW':ESW, 'ESWW':ESWW, 'EWW':EWW }

    def removeDimensions(self, axis, idx):
        # Method to remove undesired dimensions
        # - axis (int): axis from where to remove the elements
        # - idx (numpy array): indices of the elements to remove
        assert axis <= len(self.dim)
        assert s.all(idx < self.dim[axis])
        self.S.removeDimensions(axis,idx)
        self.W_S0.removeDimensions(axis,idx)
        self.W_S1.removeDimensions(axis,idx)
        self.updateParameters()
        self.updateExpectations()

    def updateDim(self, axis, new_dim):
        # Function to update the dimensionality of a particular axis
        self.S.updateDim(axis,new_dim)
        self.W_S0.updateDim(axis,new_dim)
        self.W_S1.updateDim(axis,new_dim)
        dim = list(self.dim)
        dim[axis] = new_dim
        self.dim = tuple(dim)
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

        # Initialise parameters
        theta = s.ones(dim)*theta
        N = s.ones(dim)*N
        self.params = { 'theta':theta, 'N':N }

        # Initialise expectations 
        if E is None:
            self.updateExpectations()
        else:
            E = s.ones(dim)*E
            self.expectations = { 'E':E }

        # Check that dimensionalities match
        self.CheckDimensionalities()

    def updateExpectations(self):
        E = self.params["N"] * self.params["N"]
        self.expectations = { 'E':E }

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        assert x.dtype == int, "x has to be an integer array"
        # return s.prod( stats.binom.pmf(x, self.params["N"], self.theta) )
        return s.prod( special.binom(self.params["N"],x) * self.params["theta"]**x * (1-self.params["theta"])**(self.params["N"]-x) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        assert x.dtype == int, "x has to be an integer array"
        # print s.sum (stats.binom.logpmf(x, self.params["N"], self.theta) )
        return s.sum( s.log(special.binom(self.params["N"],x)) + x*s.log(self.params["theta"]) + (self.params["N"]-x)*s.log(1-self.params["theta"]) )
class Beta(Distribution):
    """
    Class to store an arbitrary number of Beta distributions

    Equations:
    p(x|a,b) = GammaF(a+b)/(GammaF(a)*GammaF(b)) * x**(a-1) * (1-x)**(b-1)
    log p(x|a,b) = log[GammaF(a+b)/(GammaF(a)*GammaF(b))] + (a-1)*x + (b-1)*(1-x)
    E[x] = a/(a+b)
    var[x] = a*b / ((a+b)**2 * (a+b+1))
    """
    def __init__(self, dim, a, b, E=None):
        Distribution.__init__(self, dim)

        # Initialise parameters
        a = s.ones(dim)*a
        b = s.ones(dim)*b
        self.params = { 'a':a, 'b':b }

        # Initialise expectations
        if E is None: 
            self.updateExpectations()
        else:
            self.expectations = { 'E':s.ones(dim)*E }

        # Check that dimensionalities match
        self.CheckDimensionalities()

    def updateExpectations(self):
        a, b = self.params['a'], self.params['b']
        E = s.divide(a,a+b)
        lnE = special.digamma(a) - special.digamma(a+b)
        lnEInv = special.digamma(b) - special.digamma(a+b) # expectation of ln(1-X)
        self.expectations = { 'E':E, 'lnE':lnE, 'lnEInv':lnEInv }

if __name__ == "__main__":
    a = Beta(dim=(10,20), a=1, b=1, E=3)
    MultivariateGaussian(dim=(1,10), mean=stats.norm.rvs(loc=0, scale=1, size=(10,)), cov=s.eye(10,10), E=None)
    exit()

