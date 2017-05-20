
"""
TO-DO:
- Init theta in SW init?
- Some variables are init as s.nan, that depend on the schedule
- properly handle differnet trials...
- Save training options does notw ork because of dict in dropK
- Warning if negative lower bound
- Multiple trials really decreasing elbo
"""

# Import required modules
import argparse
import pandas as pd
import scipy as s

# Import manual functions
from run_utils import *


# Read arguments
p = argparse.ArgumentParser( description='Run script for MOFA' )
p.add_argument( '--inFiles',           type=str, nargs='+', required=True,                  help='Input data files (including extension)' )
p.add_argument( '--outFile',           type=str, required=True,                             help='Output data file (hdf5 format)' )
p.add_argument( '--likelihoods',       type=str, nargs='+', required=True,                  help='Likelihoods (bernoulli, gaussian, poisson, binomial)')
p.add_argument( '--views',             type=str, nargs='+', required=True,                  help='View names')
p.add_argument( '--covariatesFile',    type=str,                                            help='Input covariate data' )
p.add_argument( '--schedule',          type=str, nargs="+", required=True,                  help='Update schedule' )
p.add_argument( '--learnTheta',        action='store_true',                                 help='learn the sparsity parameter from the spike and slab?' )
p.add_argument( '--learnMean',         action='store_true',                                 help='learn the feature mean?' )
p.add_argument( '--initTheta',         type=float, default=0.5 ,                            help='hyperparameter theta in case that learnTheta is set to False')
p.add_argument( '--tolerance',         type=float, default=0.1 ,                            help='tolerance for convergence (deltaELBO)')
p.add_argument( '--startDrop',         type=int, default=5 ,                                help='First iteration to start dropping factors')
p.add_argument( '--freqDrop',          type=int, default=1 ,                                help='Frequency for dropping factors')
p.add_argument( '--maskAtRandom',      type=float,nargs="+", default=None,                  help='Fraction of data to mask per view')
p.add_argument( '--maskNSamples',      type=int,nargs="+", default=None,                    help='Number of patients to mask per view')
p.add_argument( '--nostop',            action='store_true',                                 help='Do not stop even when convergence criterion is met?' )
p.add_argument( '--ardZ',              action='store_true',                                 help='Automatic Relevance Determination on the latent variables?' )
p.add_argument( '--ardW',              type=str, default="mk" ,                             help=' "m" = Per view, "k" = Per factor, "mk" = Per view and factor' )
p.add_argument( '--dropNorm',          type=float, default=None ,                           help='Threshold to drop latent variables based on norm of the latent variable' )
p.add_argument( '--dropR2',            type=float, default=None ,                           help='Threshold to drop latent variables based on coefficient of determination' )
p.add_argument( '-n', '--ntrials',     type=int, default=1,                                 help='Number of trials' )
p.add_argument( '-c', '--ncores',      type=int, default=1,                                 help='Number of cores' )
p.add_argument( '-i', '--iter',        type=int, default=10,                                help='Number of iterations' )
p.add_argument( '-lb', '--elbofreq',   type=int, default=1,                                 help='Frequency of computation of ELBO' )
p.add_argument( '-k', '--factors',     type=int, default=10,                                help='Number of latent variables')
p.add_argument( '-v', '--verbose',     action='store_true',                                 help='More detailed log messages')
args = p.parse_args()
    
#############################
## Define the data options ##
#############################

data_opts = {}

data_opts['input_files'] = args.inFiles
data_opts['outfile'] = args.outFile
data_opts['view_names'] = args.views
if args.learnMean:
  data_opts['center'] = [False for l in args.likelihoods]
else:
  data_opts['center'] = [True if l=="gaussian" else False for l in args.likelihoods]
data_opts['rownames'] = 0
data_opts['colnames'] = 0
# data_opts['rownames'] = None
# data_opts['colnames'] = None
data_opts['delimiter'] = " "

M = len(data_opts['input_files'])

if args.maskAtRandom is not None:
  data_opts['maskAtRandom'] = args.maskAtRandom
else:
  data_opts['maskAtRandom'] = [0]*M 

if args.maskNSamples is not None:
  data_opts['maskNSamples'] = args.maskNSamples
else:
  data_opts['maskNSamples'] = [0]*M 

# Sanity checks
assert M == len(data_opts['view_names']), "Length of view names and input files does not match"
assert M == len(data_opts['maskAtRandom']), "Length of MaskAtRandom and input files does not match"
assert M == len(data_opts['maskNSamples']), "Length of MaskAtRandom and input files does not match"

###############
## Load data ##
###############

# Load observations
data = loadData(data_opts)
N = data[0].shape[0]
D = [data[m].shape[1] for m in xrange(M)]

# Load covariates
if args.covariatesFile is not None:
    data_opts['covariates'] = pd.read_csv(args.covariatesFile, delimiter=" ", header=0)
    args.factors += data_opts['covariates'].shape[1]
else:
    data_opts['covariates'] = None

# If we want to learn the mean, we add a covariate where the latent variable is all 1s
if args.learnMean:
  if data_opts['covariates'] is not None:
    data_opts['covariates'].insert(0, "mean", s.ones(N,))
  else:
    data_opts['covariates'] = s.ones((N,1))
  args.factors += 1

##############################
## Define the model options ##
##############################

model_opts = {}

# Define number of latent factors
K = model_opts['k'] = args.factors

# Define likelihoods
model_opts['likelihood'] = args.likelihoods
assert M==len(model_opts['likelihood']), "Please specify one likelihood for each view"

# Define whether to learn the feature-wise means
model_opts["learnMean"] = args.learnMean

# Define the sparsity levels
model_opts['ardZ'] = args.ardZ
model_opts['ardW'] = args.ardW
if args.learnTheta:
  model_opts['learnTheta'] = s.ones((M,K))
else:
  model_opts['learnTheta'] = s.zeros((M,K))

# Define schedule of updates
model_opts['schedule'] = args.schedule
if model_opts["ardZ"]: 
  assert "AlphaZ" in args.schedule, "AlphaZ should be in the update schedule"
else:
  assert "AlphaZ" not in args.schedule, "AlphaZ should not be in the update schedule"


####################################
## Define priors (P distribution) ##
####################################

# Latent Variables
model_opts["priorZ"] = { 'mean':s.zeros((N,K)) }
if model_opts['ardZ']:
  model_opts["priorZ"]['var'] = s.ones((N,K))*s.nan
  model_opts["priorAlphaZ"] = { 'a':s.ones(K)*1e-3, 'b':s.ones(K)*1e-3 }
else:
  model_opts["priorZ"]['var'] = s.ones((N,K))*1. 

# Weights
model_opts["priorSW"] = { 'Theta':[s.nan]*M, 'mean_S0':[s.nan]*M, 'var_S0':[s.nan]*M, 'mean_S1':[s.nan]*M, 'var_S1':[s.nan]*M } # Not required
if model_opts['ardW'] == "m":
  model_opts["priorAlphaW"] = { 'a':[1e-3]*M, 'b':[1e-3]*M }
elif model_opts['ardW'] == "mk":
  model_opts["priorAlphaW"] = { 'a':[s.ones(K)*1e-3]*M, 'b':[s.ones(K)*1e-3]*M }
elif model_opts['ardW'] == "k":
  model_opts["priorAlphaW"] = { 'a':s.ones(K)*1e-3, 'b':s.ones(K)*1e-3 }

# Theta
model_opts["priorTheta"] = { 'a':[s.ones(K) for m in xrange(M)], 'b':[s.ones(K) for m in xrange(M)] }
for m in xrange(M):
  for k in xrange(K):
    if model_opts['learnTheta'][m,k]==0:
      model_opts["priorTheta"]["a"][m][k] = s.nan
      model_opts["priorTheta"]["b"][m][k] = s.nan

# Noise
model_opts["priorTau"] = { 'a':[s.ones(D[m])*1e-3 for m in xrange(M)], 'b':[s.ones(D[m])*1e-3 for m in xrange(M)] }


##############################################
## Define initialisations of Q distribution ##
##############################################

# Latent variables 
model_opts["initZ"] = { 'mean':"random", 'var':s.ones((N,K)), 'E':None, 'E2':None }
if model_opts['ardZ']:
  model_opts["initAlphaZ"] = { 'a':s.nan, 'b':s.nan, 'E':s.ones(K)*1e2 }

# Noise
model_opts["initTau"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[s.ones(D[m])*100 for m in xrange(M)] }

# ARD of weights
if model_opts['ardW'] == "m":
  # model_opts["initAlphaW"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[10.]*M }
  model_opts["initAlphaW"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[ K*D[m]/(data[m].std(axis=0)**2 - 1./model_opts["initTau"]["E"][m]).sum() for m in xrange(M) ] } 
elif model_opts['ardW'] == "mk":
  model_opts["initAlphaW"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[s.ones(K)*100. for m in xrange(M)] }
elif model_opts['ardW'] == "k":
  model_opts["initAlphaW"] = { 'a':s.nan*s.ones(K), 'b':s.nan*s.ones(K), 'E':s.ones(K)*100. }

# Weights
model_opts["initSW"] = { 
  # 'Theta':[0.5*s.ones((D[m],K)) for m in xrange(M)],
  'Theta':[args.initTheta*s.ones((D[m],K)) for m in xrange(M)],
  'mean_S0':[s.zeros((D[m],K)) for m in xrange(M)],
  # 'var_S0':[s.ones((D[m],K))*model_opts["initAlphaW"]['E'][m] for m in xrange(M)], # IS THIS ACTUALLY USED?
  'var_S0':[s.nan*s.ones((D[m],K)) for m in xrange(M)], # IS THIS ACTUALLY USED?
  'mean_S1':[s.zeros((D[m],K)) for m in xrange(M)], # (TO-DO) allow also random
  # 'mean_S1':[stats.norm.rvs(loc=0, scale=1, size=(D[m],K)) for m in xrange(M)],
  'var_S1':[s.ones((D[m],K)) for m in xrange(M)],
  'ES':[None]*M, 'EW_S0':[None]*M, 'EW_S1':[None]*M # It will be calculated from the parameters
}

# Theta
# model_opts["initTheta"] = { 'a':[s.ones(K)*1e-5 for m in xrange(M)], 'b':[s.ones(K)*1e-5 for m in xrange(M)], 'E':[s.nan*s.zeros(K) for m in xrange(M)] } # NOT SURE HOW TO INITIALISE THIS, WHY IS EXPECTATION NOT INITIALISED?
model_opts["initTheta"] = { 'a':[s.ones(K) for m in xrange(M)], 'b':[s.ones(K) for m in xrange(M)], 'E':[s.nan*s.zeros(K) for m in xrange(M)] } 
for m in xrange(M):
  for k in xrange(K):
    if model_opts['learnTheta'][m,k]==0.:
      model_opts["initTheta"]["a"][m][k] = s.nan
      model_opts["initTheta"]["b"][m][k] = s.nan
      model_opts["initTheta"]["E"][m][k] = args.initTheta

##############################
## Learn feature-wise means ##
##############################

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


if model_opts["learnMean"]:
  model_opts['learnTheta'][:,0] = 0.

  for m in range(M): 
    # Weights
    if args.likelihoods[m]=="gaussian":
      model_opts["initSW"]["mean_S1"][m][:,0] = data[m].mean(axis=0)
      model_opts["initSW"]["var_S1"][m][:,0] = 1e-5
    # Theta
    model_opts["initSW"]["Theta"][m][:,0] = 1.
    model_opts["priorTheta"]['a'][m][0] = s.nan
    model_opts["priorTheta"]['b'][m][0] = s.nan
    model_opts["initTheta"]["a"][m][0] = s.nan
    model_opts["initTheta"]["b"][m][0] = s.nan
    model_opts["initTheta"]["E"][m][0] = 1.

#################################
## Define the training options ##
#################################

train_opts = {}

# Maximum number of iterations
train_opts['maxiter'] = args.iter

# Lower bound computation frequency
train_opts['elbofreq'] = args.elbofreq

# Save model (to finish)
train_opts['savefreq'] = s.nan
train_opts['savefolder'] = s.nan

# Verbosity
train_opts['verbosity'] = 2

# Criteria to drop latent variables while training
train_opts['drop'] = { "by_norm":args.dropNorm, "by_pvar":None, "by_cor":None, "by_r2":args.dropR2 }
train_opts['startdrop'] = args.startDrop
train_opts['freqdrop'] = args.freqDrop

# Tolerance level for convergence
train_opts['tolerance'] = args.tolerance

# Do no stop even when convergence criteria is met?
train_opts['forceiter'] = args.nostop

# Number of trials
train_opts['trials'] = args.ntrials

# Number of cores
train_opts['cores'] = args.ncores

#####################
## Train the model ##
#####################

# Keep the trial with the highest lower bound?
keep_best_run = False

# Go!
# runSingleTrial(data, data_opts, model_opts, train_opts, seed=None)
runMultipleTrials(data, data_opts, model_opts, train_opts, keep_best_run)

