

from __future__ import division
from time import time
import cPickle as pkl
import scipy as s
import os
from sys import path
import scipy.special as special
import scipy.stats as stats
import numpy.linalg  as linalg


from scGFA.core.simulate import Simulate
from scGFA.core.BayesNet import BayesNet
from scGFA.core.utils import *
from scGFA.core.multiview_nodes import *
from scGFA.core.seeger_nodes import *
from scGFA.core.nodes import Constant_Node
from scGFA.core.nonsparse_updates import *

from scGFA.run.init_nodes import *

###################
## Generate data ##
###################

# Define dimensionalities
M = 3
N = 100
D = s.asarray([1000,1000,1000])
K = 6


# Simulate data
data = {}
tmp = Simulate(M=M, N=N, D=D, K=K)

data['Z'] = s.zeros((N,K))
data['Z'][:,0] = s.sin(s.arange(N)/(N/20))
data['Z'][:,1] = s.cos(s.arange(N)/(N/20))
data['Z'][:,2] = 2*(s.arange(N)/N-0.5)
data['Z'][:,3] = stats.norm.rvs(loc=0, scale=1, size=N)
data['Z'][:,4] = stats.norm.rvs(loc=0, scale=1, size=N)
data['Z'][:,5] = stats.norm.rvs(loc=0, scale=1, size=N)

data['alpha'] = s.zeros((M,K))
data['alpha'][0,:] = [1,1,1e6,1,1e6,1e6]
data['alpha'][1,:] = [1,1e6,1,1e6,1,1e6]
data['alpha'][2,:] = [1e6,1,1,1e6,1e6,1]

data['W'], _ = tmp.initW_ard(alpha=data['alpha'])
data['mu'] = [ s.zeros(D[m]) for m in xrange(M)]
# data['tau']= [ s.ones(D[m])*1000 for m in xrange(M) ]
data['tau']= [ stats.uniform.rvs(loc=1,scale=3,size=D[m]) for m in xrange(M) ]

missingness = 0.0
Y_gaussian = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
	likelihood="gaussian", missingness=missingness)
Y_poisson = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
	likelihood="poisson", missingness=missingness)
Y_bernoulli = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
	likelihood="bernoulli", missingness=missingness)
# Y_binomial = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
	# likelihood="binomial", missingness=0.0, min_trials=10, max_trials=50)

data = Y_gaussian

######################
## Initialise model ##
######################

# Define initial number of latent variables
K = 8

# Define model dimensionalities
dim = {}
dim["M"] = M
dim["N"] = N
dim["D"] = D
dim["K"] = K


##############################
## Define the model options ##
##############################

model_opts = {}
model_opts['likelihood'] = ['gaussian']* M 
model_opts['k'] = K


# Define priors
model_opts["priorZ"] = { 'mean':0., 'cov':s.eye(K) }
model_opts["priorAlpha"] = { 'a':[1e-14]*M, 'b':[1e-14]*M }
model_opts["priorW"] = { 'mean':[s.nan]*M, 'cov':[s.nan]*M }
model_opts["priorTau"] = { 'a':[1e-14]*M, 'b':[1e-14]*M }

# Define initialisation options
model_opts["initZ"] = { 'mean':"random", 'cov':s.eye(K), 'E':None, 'E2':None }
model_opts["initAlpha"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[100.]*M }
model_opts["initW"] = { 'mean':["random"]*M, 'cov':[s.eye(K)*x for x in model_opts["initAlpha"]['E']], 'E':[None]*M}
model_opts["initTau"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[100.]*M }

# Define covariates
model_opts['covariates'] = None

# Define schedule of updates
model_opts['schedule'] = ("W","Z","Alpha","Tau","Y")


#############################
## Define training options ##
#############################

train_opts = {}
train_opts['elbofreq'] = 1
train_opts['maxiter'] = 200
train_opts['tolerance'] = 1E-2
train_opts['forceiter'] = True
# train_opts['dropK'] = { "by_norm":0.05, "by_pvar":None, "by_cor":None } # IT IS NOT WORKING
train_opts['dropK'] = { "by_norm":None, "by_pvar":None, "by_cor":None }
train_opts['savefreq'] = s.nan
train_opts['savefolder'] = s.nan
train_opts['verbosity'] = 2
# train_opts['trials'] = 1
# cores = 1
# keep_best_run = False



####################
## Start training ##
####################

# runMultipleTrials(data_opts, model_opts, train_opts, cores, keep_best_run)
# runSingleTrial(data["Y"], model_opts, train_opts)

init = init_nonsparse(dim, data, model_opts["likelihood"])


init.initZ(pmean=model_opts["priorZ"]["mean"], pcov=model_opts["priorZ"]["cov"],
           qmean=model_opts["initZ"]["mean"], qcov=model_opts["initZ"]["cov"], qE=model_opts["initZ"]["E"],
           covariates=model_opts['covariates'])

init.initW(pmean=model_opts["priorW"]["mean"], pcov=model_opts["priorW"]["cov"],
           qmean=model_opts["initW"]["mean"], qcov=model_opts["initW"]["cov"], qE=model_opts["initW"]["E"])

init.initAlpha(pa=model_opts["priorAlpha"]['a'], pb=model_opts["priorAlpha"]['b'],
               qa=model_opts["initAlpha"]['a'], qb=model_opts["initAlpha"]['b'], qE=model_opts["initAlpha"]['E'])


init.initTau(pa=model_opts["priorTau"]['a'], pb=model_opts["priorTau"]['b'],
             qa=model_opts["initTau"]['a'], qb=model_opts["initTau"]['b'], qE=model_opts["initTau"]['E'])

init.initY()

# Define the markov blanket of each node
# print "Defining Markov Blankets...\n"
init.MarkovBlanket()

# Initialise Bayesian Network
# print "Initialising Bayesian network...\n"
net = BayesNet(dim=dim, trial=1, schedule=model_opts["schedule"], nodes=init.getNodes(), options=train_opts)

# Start training
net.iterate()

# print net.getNodes()["Z"].getParameters()["cov"]
# print net.getNodes()["Z"].getParameters()["cov"][0]
Z = net.getNodes()["Z"].getExpectations()["E"]
print s.absolute(corr(Z.T,Z.T))
# print s.diag(net.getNodes()["Z"].getParameters()["cov"][0])

# print "Number of factors:"
# print s.absolute(net.getNodes()["Z"].getExpectations()["E"]).mean(axis=0)
# print s.sum( s.absolute(net.getNodes()["Z"].getExpectations()["E"]).mean(axis=0) > 0.01)
