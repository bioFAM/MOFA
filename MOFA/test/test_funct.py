"""
Script to test the model
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
from MOFA.core.simulate import Simulate
from MOFA.core.BayesNet import BayesNet
from MOFA.core.multiview_nodes import *
from MOFA.core.nongaussian_nodes import *
from MOFA.core.nodes import *
from MOFA.core.updates import *
from MOFA.core.utils import *
from MOFA.run.run_utils import *


###################
## Generate data ##
###################

# import numpy; numpy.random.seed(4)
def test(M=3, N=100, K=10, D=None):
    if D is None:
        D = s.asarray([1000]*M)
    elif isinstance(D, (int, long, float)):
        D = s.array([D]*M)
    else:
        assert len(D) == M, 'wrong shape for D'


    ## Simulate data  ##
    data = {}
    tmp = Simulate(M=M, N=N, D=D, K=K)

    # data['Z'] = s.zeros((N,K))
    # data['Z'][:,0] = s.sin(s.arange(N)/(N/20))
    # data['Z'][:,1] = s.cos(s.arange(N)/(N/20))
    # data['Z'][:,2] = 2*(s.arange(N)/N-0.5)
    # data['Z'][:,3] = stats.norm.rvs(loc=0, scale=1, size=N)
    # data['Z'][:,4] = stats.norm.rvs(loc=0, scale=1, size=N)
    # data['Z'][:,5] = stats.norm.rvs(loc=0, scale=1, size=N)
    data['Z'] = stats.norm.rvs(loc=0, scale=1, size=(N,K))


    data['alpha'] = [s.random.choice([1., 1e6], K) for m in xrange(M)]

    theta = [ s.ones(K)*1.0 for m in xrange(M) ]
    data['S'], data['W'], data['W_hat'], _ = tmp.initW_spikeslab(theta=theta, alpha=data['alpha'])

    data['mu'] = [ s.ones(D[m])*0. for m in xrange(M)]
    data['tau']= [ stats.uniform.rvs(loc=1,scale=3,size=D[m]) for m in xrange(M) ]
    # data['tau']= [ stats.uniform.rvs(loc=0.1,scale=3,size=D[m]) for m in xrange(M) ]

    missingness = 0.0
    missing_view = 0.0
    # Y_warp = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
    #   likelihood="warp", missingness=missingness, missing_view=missing_view)
    Y_gaussian = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
    	likelihood="gaussian", missingness=missingness, missing_view=missing_view)
    # Y_poisson = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
    	# likelihood="poisson", missingness=missingness, missing_view=missing_view)
    Y_bernoulli = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
    	likelihood="bernoulli", missingness=missingness, missing_view=missing_view)
    # Y_binomial = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
    # 	likelihood="binomial", min_trials=10, max_trials=50, missingness=missingness)

    likelihoods = np.random.choice(['gaussian', 'bernoulli'], M)
    data["Y"] = [None] * M
    for i in range(M):
        lik = likelihoods[i]
        if lik == 'gaussian':
            data["Y"][i] = Y_gaussian[i]
        elif lik == 'poisson':
            data["Y"][i] = Y_poisson[i]
        elif lik == 'bernoulli':
            data["Y"][i] = Y_bernoulli[i]
    data["Y"] = tuple(data["Y"])
    #

    # ( Y_gaussian[0], Y_gaussian[1], Y_gaussian[2] )
    # data["Y"] = ( Y_gaussian[0], Y_bernoulli[1], Y_bernoulli[2] )
    # data["Y"] = ( Y_bernoulli[0], Y_bernoulli[1], Y_bernoulli[2] )
    # data["Y"] = ( Y_gaussian[0], Y_poisson[1], Y_bernoulli[2] )
    # data["Y"] = ( Y_warp[0], )
    # data["Y"] = ( Y_warp[0], Y_warp[1], Y_warp[2] )
    # data["Y"] = ( Y_gaussian[0], )

    ##################
    ## Data options ##
    ##################


    data_opts = {}
    data_opts["outfile"] = "/tmp/test.h5"
    data_opts['view_names'] = ['foo'] * M
    data_opts['center'] = [True]*M
    data_opts['scale_views'] = [False]*M
    data_opts['covariates'] = None

    # Center the data
    if data_opts['center'][m]:
      data[m] = (data["Y"][m] - data["Y"][m].mean())

    # Scale the views
    # NOT WORKING YET
    # if data_opts['scale_views'][m]:
    #   print type(data["Y"][m])
    #   print s.std(data["Y"][m], axis=0).sum()
    #   print data["Y"][m] / s.std(data["Y"][m], axis=0).sum()
    #   exit()
    #   data["Y"][m] = data["Y"][m].multiply(1./s.std(data["Y"][m], axis=0).sum())
    #   exit()

    #################################
    ## Initialise Bayesian Network ##
    #################################

    # Define initial number of latent variables
    # take the real one

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
    # define likelihoods
    model_opts['likelihood'] = likelihoods

    # Define type of covariates
    model_opts['nonsparse_covariates'] = None
    if data_opts['covariates'] is not None:
      dim["K"] += data_opts['covariates'].shape[1]
      K = dim["K"]

    # Define whether to learn means
    model_opts["learnMean"] = False
    if model_opts["learnMean"]:
      # MODIFY THIS
      data_opts['covariates'] = s.ones((dim["N"],1))
      dim["K"] += data_opts['covariates'].shape[1]
      K = dim["K"]
      model_opts['nonsparse_covariates'] = [True]

    # Define initial number of factors
    model_opts['k'] = K

    # Define sparsities
    model_opts['ardZ'] = False
    model_opts['ardW'] = "mk"

    # Define for which factors to learn Theta
    model_opts['learnTheta'] = 0.*s.ones((M,K))

    # Define schedule of updates
    # model_opts['schedule'] = ["Y","SW","Z","AlphaW","Theta","Tau"]
    model_opts['schedule'] = ['Y', "SW","Z","AlphaW","Theta","Tau"]

    ####################################
    ## Define priors (P distribution) ##
    ####################################

    # Latent Variables
    if model_opts['ardZ']:
      model_opts["priorZ"] = { 'mean':s.zeros((N,K)), 'var':s.nan }
      model_opts["priorAlphaZ"] = { 'a':s.ones(K)*1e-5, 'b':s.ones(K)*1e-5 }
    else:
      model_opts["priorZ"] = { 'mean':s.zeros((N,K)), 'var':s.ones((K,)) }

    # Weights
    model_opts["priorSW"] = { 'Theta':[s.nan]*M, 'mean_S0':[s.nan]*M, 'var_S0':[s.nan]*M, 'mean_S1':[s.nan]*M, 'var_S1':[s.nan]*M } # Not required
    if model_opts['ardW'] == "m":
      model_opts["priorAlphaW"] = { 'a':[1e-3]*M, 'b':[1e-3]*M }
    elif model_opts['ardW'] == "k":
      model_opts["priorAlphaW"] = { 'a':s.ones(K)*1e-3, 'b':s.ones(K)*1e-3 }
    elif model_opts['ardW'] == "mk":
      model_opts["priorAlphaW"] = { 'a':[s.ones(K)*1e-3]*M, 'b':[s.ones(K)*1e-3]*M }

    # Theta
    model_opts["priorTheta"] = { 'a':[s.ones(K) for m in xrange(M)], 'b':[s.ones(K) for m in xrange(M)] }
    for m in xrange(M):
      for k in xrange(K):
        if model_opts['learnTheta'][m,k]==0:
          model_opts["priorTheta"]["a"][m][k] = s.nan
          model_opts["priorTheta"]["b"][m][k] = s.nan


    # Noise
    model_opts["priorTau"] = { 'a':[s.ones(D[m])*1e-3 for m in xrange(M)], 'b':[s.ones(D[m])*1e-3 for m in xrange(M)] }


    #############################################
    ## Define initialisations (Q distribution) ##
    #############################################

    # Latent variables
    # model_opts["initZ"] = { 'mean':"orthogonal", 'var':s.ones((N,K)), 'E':None, 'E2':None }
    model_opts["initZ"] = { 'mean':"random", 'var':s.ones((N,K)), 'E':None, 'E2':None }
    if model_opts['ardZ']:
      model_opts["initAlphaZ"] = { 'a':s.nan, 'b':s.nan, 'E':s.ones(K)*100 }

    # Weights
    if model_opts['ardW'] == "m":
      model_opts["initAlphaW"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[1.]*M }
    elif model_opts['ardW'] == "k":
      model_opts["initAlphaW"] = { 'a':s.nan*s.ones(K), 'b':s.nan*s.ones(K), 'E':s.ones(K) }
    elif model_opts['ardW'] == "mk":
      model_opts["initAlphaW"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[1.*s.ones(K) for m in xrange(M)] }
    model_opts["initSW"] = { 'Theta':[s.ones((D[m],K)) for m in xrange(M)],
                              'mean_S0':[s.zeros((D[m],K)) for m in xrange(M)],
                              'var_S0':[s.ones((D[m],K)) for m in xrange(M)],
                              'mean_S1':[s.zeros((D[m],K)) for m in xrange(M)], # (TO-DO) allow also random
                              # 'mean_S1':[s.ones((D[m],K)) for m in xrange(M)],
                              'var_S1':[s.ones((D[m],K)) for m in xrange(M)],
                              'ES':[None]*M, 'EW_S0':[None]*M, 'EW_S1':[None]*M }

    # Theta
    model_opts["initTheta"] = { 'a':[s.ones(K) for m in xrange(M)], 'b':[s.ones(K) for m in xrange(M)], 'E':[s.zeros(K)*s.nan for m in xrange(M)] }
    for m in xrange(M):
      for k in xrange(K):
        if model_opts['learnTheta'][m,k]==0.:
          model_opts["initTheta"]["a"][m][k] = s.nan
          model_opts["initTheta"]["b"][m][k] = s.nan
          # model_opts["initTheta"]["E"][m][k] = 0.5 # CHANGEE THISSSSsSSSSSS
          model_opts["initTheta"]["E"][m][k] = 1.0


    # Noise
    model_opts["initTau"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[s.ones(D[m])*100 for m in xrange(M)] }


    ######################################################
    ## Modify priors and initialisations for covariates ##
    ######################################################

    if data_opts['covariates'] is not None:
      idx = xrange(data_opts['covariates'].shape[1])

      ## Prior distributions (P) ##
      # Latent variables
      if model_opts['ardZ']:
        model_opts["priorZ"]["mean"][:,idx] = s.nan
        model_opts["priorAlphaZ"]["a"][idx] = s.nan
        model_opts["priorAlphaZ"]["b"][idx] = s.nan
      else:
        model_opts["priorZ"]["var"][:,idx] = s.nan

      ## Variational distributions (Q) ##
      # Latent variables
      # model_opts["initZ"]["mean"][:,idx] = model_opts["covariates"]
      model_opts["initZ"]["var"][:,idx] = 0.
      if model_opts['ardZ']:
            model_opts["initAlphaZ"]["E"][idx] = s.nan

      ###########################
      ## Non-sparse covariates ##
      ###########################

      if any(model_opts['nonsparse_covariates']):
        idx = s.where(model_opts['nonsparse_covariates'])
        model_opts['learnTheta'][:,idx] = 0.

        # WHERE IS THE LOOP OVER M?????

        ## Prior distributions (P) ##
        # Theta
        model_opts["priorTheta"]['a'][m][idx] = s.nan
        model_opts["priorTheta"]['b'][m][idx] = s.nan

        ## Variational distributions (Q) ##
        for m in range(M):
          # Weights
          model_opts["initSW"]["Theta"][m][:,idx] = 1.
          # Theta
          model_opts["initTheta"]["a"][m][idx] = s.nan
          model_opts["initTheta"]["b"][m][idx] = s.nan
          model_opts["initTheta"]["E"][m][idx] = 1.

    #############################
    ## Define training options ##
    #############################

    train_opts = {}
    train_opts['elbofreq'] = 1
    train_opts['maxiter'] = 10000
    # train_opts['tolerance'] = 1E-2
    train_opts['tolerance'] = 0.1
    train_opts['forceiter'] = False
    train_opts['drop'] = { "by_norm":None, "by_pvar":None, "by_cor":None, "by_r2":0. }
    train_opts['startdrop'] = 10
    train_opts['freqdrop'] = 1
    train_opts['savefreq'] = s.nan
    train_opts['savefolder'] = s.nan
    train_opts['verbosity'] = 2
    train_opts['trials'] = 1
    train_opts['cores'] = 1


    ####################
    ## Start training ##
    ####################

    keep_best_run = False
    # runMultipleTrials(data["Y"], data_opts, model_opts, train_opts, keep_best_run)

    time1 = time()
    runSingleTrial(data["Y"], data_opts, model_opts, train_opts)
    time2 = time()

    # return {'time': (time2-time1)*1000, 'likelihoods': likelihoods}
    return (time2-time1)
