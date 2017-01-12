

from __future__ import division
from time import time
import cPickle as pkl
import scipy as s
import os
from sys import path
import scipy.special as special
import scipy.stats as stats
import numpy.linalg  as linalg


path.insert(0,"../")
from simulate import Simulate
from BayesNet import BayesNet
from utils import *
from multiview_nodes import *
from seeger_nodes import Binomial_PseudoY_Node, Poisson_PseudoY_Node, Bernoulli_PseudoY_Node
from nodes import Constant_Node
from nonsparse_updates import Y_Node, Alpha_Node, W_Node, Tau_Node, Z_Node, Y_Node

###################
## Generate data ##
###################

# Define dimensionalities
M = 1
N = 100
D = s.asarray([100,])
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

data['alpha'] = s.zeros((M,K))
data['alpha'][0,:] = [1,1,1e6,1,1e6,1e6]
# data['alpha'][1,:] = [1,1e6,1,1e6,1,1e6]
# data['alpha'][2,:] = [1e6,1,1,1e6,1e6,1]

data['W'], _ = tmp.initW_ard(alpha=data['alpha'])
data['mu'] = [ s.zeros(D[m]) for m in xrange(M)]
# data['tau']= [ s.ones(D[m])*1000 for m in xrange(M) ]
data['tau']= [ stats.uniform.rvs(loc=1,scale=3,size=D[m]) for m in xrange(M) ]
Y_gaussian = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'], 
	likelihood="gaussian", missingness=0.00)
Y_poisson = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'], 
	likelihood="poisson", missingness=0.00)
Y_bernoulli = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'], 
	likelihood="bernoulli", missingness=0.0)
# Y_binomial = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'], 
	# likelihood="binomial", missingness=0.0, min_trials=10, max_trials=50)

# data["Y"] = ( Y_gaussian[0], Y_poisson[1], Y_bernoulli[2] )
# data["Y"] = Y_bernoulli
# data["Y"] = Y_poisson
# data["Y"] = Y_binomial
data["Y"] = Y_gaussian

# likelihood = ["gaussian","gaussian","gaussian"]
likelihood = ["gaussian"]
view_names = ["gaussian"]
# view_names = ["gaussian1","gaussian2","gaussian3"]

#################################
## Initialise Bayesian Network ##
#################################

# Define initial number of latent variables
K = 10

# Define model dimensionalities
dim = {}
dim["M"] = M
dim["N"] = N
dim["D"] = D
dim["K"] = K

###############
## Add nodes ##
###############

# Z without covariates (variational node)
Z_qmean = s.stats.norm.rvs(loc=0, scale=1, size=(N,K))
Z_qcov = s.repeat(s.eye(K)[None,:,:],N,0)
Z = Z_Node(dim=(N,K), pmean=s.nan, pcov=s.nan, qmean=Z_qmean, qcov=Z_qcov, qE=Z_qmean)

Z.updateExpectations()

# Z with covariates (variational node)
# Z_qmean = s.stats.norm.rvs(loc=0, scale=1, size=(N,K-1))
# Z_qmean = s.c_[ Z_qmean, s.asarray([True,False]*(N/2), dtype=s.float32) ]
# Z_qcov = s.repeat(s.eye(K)[None,:,:],N,0)
# Z = Z_Node(dim=(N,K), qmean=Z_qmean, qcov=Z_qcov)
# Z.updateExpectations()
# Z.setCovariates(idx=K-1)


# alpha (variational node)
alpha_list = [None]*M
alpha_pa = 1e-14
alpha_pb = 1e-14
# alpha_qa = s.ones(K)
alpha_qb = s.ones(K)
alpha_qE = s.ones(K)*100
for m in xrange(M):
	alpha_qa = alpha_pa + s.ones(K)*D[m]/2
	alpha_list[m] = Alpha_Node(dim=(K,), pa=alpha_pa, pb=alpha_pb, qa=alpha_qa[m], qb=alpha_qb, qE=alpha_qE)
alpha = Multiview_Variational_Node((K,)*M, *alpha_list)


# W (variational node)
W_list = [None]*M
for m in xrange(M):
	W_qmean = s.stats.norm.rvs(loc=0, scale=1, size=(D[m],K))
	W_qcov = s.repeat(a=s.eye(K)[None,:,:], repeats=D[m] ,axis=0)
	W_list[m] = W_Node(dim=(D[m],K), pmean=s.nan, pcov=s.nan, qmean=W_qmean, qcov=W_qcov, qE=W_qmean)
W = Multiview_Variational_Node(M, *W_list)


# tau/kappa (mixed node)
tau_list = [None]*M
for m in xrange(M):
	if likelihood[m] == "poisson":
		tmp = 0.25 + 0.17*s.amax(data["Y"][m],axis=0) 
		tau_list[m] = Observed_Local_Node(dim=(D[m],), value=tmp)
	elif likelihood[m] == "bernoulli":
		tmp = s.ones(D[m])*0.25 
		tau_list[m] = Observed_Local_Node(dim=(D[m],), value=tmp)
	elif likelihood[m] == "binomial":
		tmp = 0.25*s.amax(data["Y"]["tot"][m],axis=0)
		tau_list[m] = Observed_Local_Node(dim=(D[m],), value=tmp)
	elif likelihood[m] == "gaussian":
		tau_pa = 1e-14
		tau_pb = 1e-14
		tau_qa = tau_pa + s.ones(D[m])*N/2
		tau_qb = s.zeros(D[m])
		tau_qE = s.zeros(D[m]) + 100
		tau_list[m] = Tau_Node(dim=(D[m],), pa=tau_pa, pb=tau_pb, qa=tau_qa, qb=tau_qb, qE=tau_qE)
tau = Multiview_Mixed_Node(M,*tau_list)



# Y/Yhat (mixed node)
Y_list = [None]*M
for m in xrange(M):
	if likelihood[m] == "gaussian":
		Y_list[m] = Y_Node(dim=(N,D[m]), obs=data['Y'][m])
	elif likelihood[m] == "poisson":
		Y_list[m] = Poisson_PseudoY_Node(dim=(N,D[m]), obs=data['Y'][m], Zeta=None, E=None)
	elif likelihood[m] == "bernoulli":
		Y_list[m] = Bernoulli_PseudoY_Node(dim=(N,D[m]), obs=data['Y'][m], Zeta=None, E=None)
	elif likelihood[m] == "binomial":
		Y_list[m] = Binomial_PseudoY_Node(dim=(N,D[m]), tot=data['Y']["tot"][m], obs=data['Y']["obs"][m], Zeta=None, E=None)
Y = Multiview_Mixed_Node(M, *Y_list)


############################
## Define Markov Blankets ##
############################

Z.addMarkovBlanket(W=W, tau=tau, Y=Y)
for m in xrange(M):
	alpha.nodes[m].addMarkovBlanket(W=W.nodes[m])
	W.nodes[m].addMarkovBlanket(Z=Z, tau=tau.nodes[m], alpha=alpha.nodes[m], Y=Y.nodes[m])
	if likelihood[m] == "gaussian":
		Y.nodes[m].addMarkovBlanket(Z=Z, W=W.nodes[m], tau=tau.nodes[m])
		tau.nodes[m].addMarkovBlanket(W=W.nodes[m], Z=Z, Y=Y.nodes[m])
	else:
		Y.nodes[m].addMarkovBlanket(Z=Z, W=W.nodes[m], kappa=tau.nodes[m])


# Define update schedule
schedule = ["Y","W","Z","alpha","tau"]
nodes = { "W":W, "tau":tau, "Z":Z, "Y":Y, "alpha":alpha }

#############################
## Define training options ##
#############################

options = {}
options['maxiter'] = 500
options['tolerance'] = 1E-2
options['forceiter'] = False
options['elbofreq'] = 1
options['dropK'] = { 'by_norm':None, 'by_pvar':None, 'by_cor':None }
options['savefreq'] = options['maxiter']+1
options['savefolder'] = "/tmp/test"
options['verbosity'] = 2

####################
## Start training ##
####################

net = BayesNet(dim=dim, schedule=schedule, nodes=nodes, options=options)
net.iterate()

