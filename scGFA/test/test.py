"""
Script to test missing values
"""

from __future__ import division
from time import time
import cPickle as pkl
import scipy as s
import os
import scipy.special as special
import scipy.stats as stats
import numpy.linalg  as linalg

# Import manually defined functions
from scGFA.core.simulate import Simulate
from scGFA.core.BayesNet import BayesNet
from scGFA.core.multiview_nodes import *
from scGFA.core.seeger_nodes import *
from scGFA.core.nodes import *
from scGFA.core.sparse_updates import *
from scGFA.core.utils import *
from scGFA.run.run_utils import *

###################
## Generate data ##
###################

# import numpy; numpy.random.seed(4)

# Define dimensionalities
M = 3
N = 50
D = s.asarray([500,500,500])
K = 6


## Simulate data  ##
data = {}
tmp = Simulate(M=M, N=N, D=D, K=K)

data['Z'] = s.zeros((N,K))
data['Z'][:,0] = s.sin(s.arange(N)/(N/20))
data['Z'][:,1] = s.cos(s.arange(N)/(N/20))
data['Z'][:,2] = 2*(s.arange(N)/N-0.5)
data['Z'][:,3] = stats.norm.rvs(loc=0, scale=1, size=N)
data['Z'][:,4] = stats.norm.rvs(loc=0, scale=1, size=N)
data['Z'][:,5] = stats.norm.rvs(loc=0, scale=1, size=N)

# data['Z'] = stats.norm.rvs(loc=0, scale=1, size=(N,K))


data['alpha'] = [ s.zeros(K,) for m in xrange(M) ]
data['alpha'][0] = [1,1,1e6,1,1e6,1e6]
data['alpha'][1] = [1,1e6,1,1e6,1,1e6]
data['alpha'][2] = [1e6,1,1,1e6,1e6,1]

theta = [ s.ones(K)*0.5 for m in xrange(M) ]
# theta = [ s.ones(K)*1.0 for m in xrange(M) ]
data['S'], data['W'], data['W_hat'], _ = tmp.initW_spikeslab(theta=theta, alpha=data['alpha'])

data['mu'] = [ s.zeros(D[m]) for m in xrange(M)]
data['tau']= [ stats.uniform.rvs(loc=1,scale=3,size=D[m]) for m in xrange(M) ]
# data['tau']= [ stats.uniform.rvs(loc=0.1,scale=3,size=D[m]) for m in xrange(M) ]

missingness = 0.
Y_gaussian = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
	likelihood="gaussian", missingness=missingness)
Y_poisson = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
	likelihood="poisson", missingness=missingness)
Y_bernoulli = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
	likelihood="bernoulli", missingness=missingness)
# Y_binomial = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
# 	likelihood="binomial", min_trials=10, max_trials=50, missingness=missingness)

data["Y"] = ( Y_gaussian[0], Y_gaussian[1], Y_gaussian[2] )
# data["Y"] = ( Y_bernoulli[0], Y_bernoulli[1], Y_bernoulli[2] )
# data["Y"] = ( Y_gaussian[0], Y_poisson[1], Y_bernoulli[2] )



#################################
## Initialise Bayesian Network ##
#################################

# Define initial number of latent variables
K = 20

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
model_opts['learnTheta'] = True
model_opts['k'] = K
model_opts['ardZ'] = False
model_opts['ardW'] = "extended"

####################################
## Define priors (P distribution) ##
####################################

# Latent Variables
if model_opts['ardZ']:
  model_opts["priorZ"] = { 'mean':0., 'var':s.nan }
  model_opts["priorAlphaZ"] = { 'a':1e-5, 'b':1e-5 }
else:
  model_opts["priorZ"] = { 'mean':0., 'var':1. }

# Weights
model_opts["priorSW"] = { 'Theta':[s.nan]*M, 'mean_S0':[s.nan]*M, 'var_S0':[s.nan]*M, 'mean_S1':[s.nan]*M, 'var_S1':[s.nan]*M } # Not required
if model_opts['ardW'] == "basic":
  model_opts["priorAlphaW"] = { 'a':[1e-5]*M, 'b':[1e-5]*M }
elif model_opts['ardW'] == "extended":
  model_opts["priorAlphaW"] = { 'a':[1e-5]*M, 'b':[1e-5]*M }
if model_opts['learnTheta']:
    model_opts["priorTheta"] = { 'a':[1.]*M, 'b':[1.]*M }

# Noise
model_opts["priorTau"] = { 'a':[1e-5]*M, 'b':[1e-5]*M }

#############################################
## Define initialisations (Q distribution) ##
#############################################

# Latent variables
model_opts["initZ"] = { 'mean':"orthogonal", 'var':1., 'E':None, 'E2':None }
if model_opts['ardZ']:
  model_opts["initAlphaZ"] = { 'a':s.nan, 'b':s.nan, 'E':1. }

# Weights
model_opts["initAlphaW"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[10.]*M }
model_opts["initSW"] = { 'Theta':[0.5]*M,
                          'mean_S0':[0.]*M, 'var_S0':model_opts["initAlphaW"]['E'],
                          # 'mean_S1':["random"]*M, 'var_S1':[1.]*M,
                          'mean_S1':[0.]*M, 'var_S1':[1.]*M,
                          'ES':[None]*M, 'EW_S0':[None]*M, 'EW_S1':[None]*M}
if model_opts['learnTheta']:
    model_opts["initTheta"] = { 'a':[1.]*M, 'b':[1.]*M, 'E':[None]*M }
else:
    model_opts["initTheta"] = { 'value':[args.initTheta]*M }

# Noise
model_opts["initTau"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[100.]*M }


# Define covariates
model_opts['covariates'] = None


#############################
## Define training options ##
#############################

train_opts = {}
# Define schedule of updates
# train_opts['schedule'] = ("Y","Tau","SW", "Z", "Clusters", "Theta", "Alpha")
train_opts['schedule'] = ("SW","Z","AlphaZ","AlphaW","Tau","Theta")
train_opts['elbofreq'] = 1
train_opts['maxiter'] = 3000
# train_opts['tolerance'] = 1E-2
train_opts['tolerance'] = 0.01
train_opts['forceiter'] = True
train_opts['drop'] = { "by_norm":None, "by_pvar":None, "by_cor":None, "by_r2":None }
train_opts['startdrop'] = 10
train_opts['freqdrop'] = 1
train_opts['savefreq'] = s.nan
train_opts['savefolder'] = s.nan
train_opts['verbosity'] = 2
train_opts['trials'] = 1
train_opts['cores'] = 1
keep_best_run = False



####################
## Start training ##
####################

data_opts = {}
data_opts["outfile"] = "/tmp/test.h5"
data_opts['view_names'] = ["g","p","b"]

# runMultipleTrials(data["Y"], data_opts, model_opts, train_opts, keep_best_run)
runSingleTrial(data["Y"], data_opts, model_opts, train_opts)
