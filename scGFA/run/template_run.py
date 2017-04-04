
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
p.add_argument( '--initTheta',         type=float, default=0.5 ,                            help='hyperparameter theta in case that learnTheta is set to False')
p.add_argument( '--tolerance',         type=float, default=0.1 ,                            help='tolerance for convergence (deltaELBO)')
p.add_argument( '--startDrop',         type=int, default=5 ,                                help='First iteration to start dropping factors')
p.add_argument( '--freqDrop',         type=int, default=1 ,                                 help='Frequency for dropping factors')
p.add_argument( '--maskAtRandom',  	   type=float,nargs="+", default=None,                  help='Fraction of data to mask per view')
p.add_argument( '--maskNSamples',      type=int,nargs="+", default=None,                    help='Number of patients to mask per view')
p.add_argument( '--nostop',            action='store_true',                                 help='Do not stop even when convergence criterion is met?' )
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
data_opts['center'] = [True if l=="gaussian" else False for l in args.likelihoods]
data_opts['rownames'] = 0
data_opts['colnames'] = 0
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

# pprint(data_opts)
# print "\n"

# Load data
data = loadData(data_opts)

##############################
## Define the model options ##
##############################


model_opts = {}

# Define likelihoods
model_opts['likelihood'] = args.likelihoods
assert M==len(model_opts['likelihood']), "Please specify one likelihood for each view"

# Define number of latent variables
model_opts['k'] = args.factors

# Define whether to learn the sparsity parameter of the spike and slab
model_opts['learnTheta'] = args.learnTheta

# Define priors
model_opts["priorZ"] = { 'mean':0., 'var':1. }
# model_opts["priorAlpha"] = { 'a':[1e-5]*M, 'b':[1e-5]*M }
model_opts["priorAlpha"] = { 'a':[1e-14]*M, 'b':[1e-14]*M }
model_opts["priorSW"] = { 'Theta':[s.nan]*M, 'mean_S0':[s.nan]*M, 'var_S0':[s.nan]*M, 'mean_S1':[s.nan]*M, 'var_S1':[s.nan]*M }
# model_opts["priorTau"] = { 'a':[1e-5]*M, 'b':[1e-5]*M }
model_opts["priorTau"] = { 'a':[1e-14]*M, 'b':[1e-14]*M }
if model_opts['learnTheta']:
    model_opts["priorTheta"] = { 'a':[1.]*M, 'b':[1.]*M }

# Define initialisations
model_opts["initZ"] = { 'mean':"random", 'var':1., 'E':None, 'E2':None }
# model_opts["initAlpha"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[100.]*M }
model_opts["initAlpha"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[1.]*M }
model_opts["initSW"] = { 'Theta':[0.5]*M,
                          'mean_S0':[0.]*M, 'var_S0':model_opts["initAlpha"]['E'],
                          # 'mean_S1':["random"]*M, 'var_S1':[1.]*M,
                          'mean_S1':[0.]*M, 'var_S1':[1.]*M,
                          'ES':[None]*M, 'EW_S0':[None]*M, 'EW_S1':[None]*M}
# model_opts["initTau"] = { 'a':[1.,1.,None], 'b':[1.,1.,None], 'E':[100.,100.,None] }
model_opts["initTau"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[100.]*M }

if model_opts['learnTheta']:
    model_opts["initTheta"] = { 'a':[1.]*M, 'b':[1.]*M, 'E':[None]*M }
else:
    model_opts["initTheta"] = { 'value':[args.initTheta]*M }


# Define schedule of updates
model_opts['schedule'] = args.schedule

# pprint(model_opts)
# print "\n"

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
train_opts['drop'] = { "by_norm":None, "by_pvar":None, "by_cor":None, "by_r2":None }
train_opts['startdrop'] = args.startDrop
train_opts['freqdrop'] = args.freqDrop

# Tolerance level for convergence
train_opts['tolerance'] = args.tolerance

# Do no stop even when convergence criteria is met?
train_opts['forceiter'] = args.nostop

# Covariates
if args.covariatesFile is not None:
    model_opts['covariates'] = pd.read_csv(args.covariatesFile, delimiter=" ", header=0, index_col=0)
else:
    model_opts['covariates'] = None

# Number of trials
train_opts['trials'] = args.ntrials

# Number of cores
train_opts['cores'] = args.ncores

# Keep the trial with the highest lower bound?
keep_best_run = False

# pprint(train_opts)
# print "\n"

# Go!
runMultipleTrials(data, data_opts, model_opts, train_opts, keep_best_run)
# runMultipleTrials(data_opts, model_opts, train_opts, keep_best_run)

