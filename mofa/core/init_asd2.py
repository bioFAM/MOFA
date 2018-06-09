import argparse
import pandas as pd
import scipy as s
import sys
from time import sleep

from .build_model import *

"""
TO-DO: 
- ADD DETAILED EXPLANATION OF ARGUMENTS
- ADd sanity checks
- Add print messages
- Within each method, check that the pipeline order is met

Pipeline
(1) Parse data options
(2) Parse train options or parse model options
(3) Parse data processing options
(4) Load the data or define priors or define variational
(5) Train

CHANGE PARSE TO SET
SET PRIORS AND VARIATIONAL DIST IN BUILD_MODEL, ONLY BIOFAM BRANCH
"""

class entry_point():
  def __init__(self):
    self.print_banner()
    self.dimensionalities = {}

  def print_banner(self):
    """ Method to print the MOFA banner """

    banner = """
    ###########################################################
    ###                 __  __  ___  _____ _                ### 
    ###                |  \/  |/ _ \|  ___/ \               ### 
    ###                | |\/| | | | | |_ / _ \              ### 
    ###                | |  | | |_| |  _/ ___ \             ### 
    ###                |_|  |_|\___/|_|/_/   \_\            ### 
    ###                                                     ###
    ########################################################### """

    print(banner)
    sleep(2)

  def set_data_options(self, 
    inFiles, outFile, views, 
    delimiter=" ", header_cols=False, header_rows=False
    ):
    """ Parse I/O data options """

    # TO-DO: sanity checks 
    # - Check that input file exists
    # - Check that output directory exists: warning if not
    # TO-DO: Verbosity, print messages 

    self.data_opts = {}

    # I/O
    if type(inFiles) is not list:
      inFiles = [inFiles]
    self.data_opts['input_files'] = inFiles
    self.data_opts['outfile'] = outFile
    self.data_opts['delimiter'] = delimiter

    # View names
    if type(views) is not list:
      views = [views]
    self.data_opts['view_names'] = views

    # Headers
    if header_rows is True:
      self.data_opts['rownames'] = 0
    else:
      self.data_opts['rownames'] = None

    if header_cols is True:
      self.data_opts['colnames'] = 0
    else:
      self.data_opts['colnames'] = None


    self.dimensionalities["M"] = len(self.data_opts['input_files'])

  def set_train_options(self, 
    iter=5000, elbofreq=1, startSparsity=100, tolerance=0.01, 
    startDrop=1, freqDrop=1, dropR2=0, nostop=False, verbose=False, seed=None
    ):
    """ Parse training options """

    # TO-DO: verbosity, print more messages

    self.train_opts = {}

    # Maximum number of iterations
    self.train_opts['maxiter'] = int(iter)

    # Lower bound computation frequency
    self.train_opts['elbofreq'] = int(elbofreq)

    # Verbosity
    self.train_opts['verbose'] = verbose

    # Criteria to drop latent variables while training
    self.train_opts['drop'] = { "by_r2":float(dropR2) }
    self.train_opts['startdrop'] = int(startDrop)
    self.train_opts['freqdrop'] = int(freqDrop)
    if (dropR2>0): print("\nDropping factors with minimum threshold of {0}% variance explained".format(dropR2))


    # Tolerance level for convergence
    self.train_opts['tolerance'] = float(tolerance)

    # Do no stop even when convergence criteria is met
    self.train_opts['forceiter'] = nostop

    # Iteration to activate spike and slab sparsity
    self.train_opts['startSparsity'] = int(startSparsity)

    # Seed
    if seed is None:
      seed = 0
    self.train_opts['seed'] = int(seed)

  def set_model_options(self, factors, likelihoods, schedule=None, sparsity=True, learnIntercept=False):
    """ Parse model options """

    # TO-DO: SANITY CHECKS 

    self.model_opts = {}

    # Define initial number of latent factors
    K = self.dimensionalities["K"] = self.model_opts['factors'] = int(factors)
    M = self.dimensionalities["M"]
    
    # Define likelihoods
    self.model_opts['likelihoods'] = likelihoods
    if type(self.model_opts['likelihoods']) is not list:
      self.model_opts['likelihoods'] = [self.model_opts['likelihoods']]

    assert len(self.model_opts['likelihoods'])==M, "Please specify one likelihood for each view"
    assert set(self.model_opts['likelihoods']).issubset(set(["gaussian","bernoulli","poisson"]))

    # Define whether to learn the feature-wise means
    self.model_opts["learnIntercept"] = learnIntercept
    if learnIntercept: 
      self.model_opts['factors'] += 1
      self.dimensionalities["K"] += 1
      K += 1

    # Define for which factors and views should we learn 'theta', the sparsity of the factor
    if type(sparsity) is bool:
      # self.model_opts['sparsity'] = True
      self.model_opts['sparsity'] = [s.ones(K) for m in range(M)]
    # elif type(sparsity) is list:
    #   self.model_opts['sparsity'] = True
    #   assert len(sparsity)==M, "--sparsity has to be a binary vector with length number of views"
    #   self.model_opts['sparsity'] = [ sparsity[m]*s.ones(K) for m in range(M) ]
    else:
       print("--sparsity has to be either 1 or 0 or a binary vector with length number of views")
       exit()

    # Define schedule of updates
    if schedule is not None:
      print("\nWarning... we recommend using the default training schedule\n")
      self.model_opts['schedule'] = schedule
    else:
      self.model_opts['schedule'] = ( "Y", "SW", "Z", "AlphaW", "Theta", "Tau" )

  def set_dataprocessing_options(self, 
    center_features=False, scale_features=False, scale_views=False, 
    maskAtRandom=None, maskNSamples=None, RemoveIncompleteSamples=False
    ):

    """ Parse data processing options """

    # TO-DO: more verbose messages


    # Sanity checks
    M = self.dimensionalities["M"]
    assert len(self.data_opts['view_names'])==M, "Length of view names and input files does not match"

    # Data processing: center features
    if center_features is True:
      self.data_opts['center_features'] = [ True if l=="gaussian" else False for l in self.model_opts["likelihoods"] ]
    else:
      if not self.model_opts["learnIntercept"]: print("\nWarning... you are not centering the data and not learning the mean...\n")
      self.data_opts['center_features'] = [ False for l in self.model_opts["likelihoods"] ]

    # Data processing: scale views
    if scale_views is True:
      self.data_opts['scale_views'] = [ True if l=="gaussian" else False for l in self.model_opts["likelihoods"] ]
    else:
      self.data_opts['scale_views'] = [ False for l in self.model_opts["likelihoods"] ]

    # Data processing: scale features
    if scale_features:
      assert data_opts['scale_views'] is False, "Scale either entire views or features, not both"
      self.data_opts['scale_features'] = [ True if l=="gaussian" else False for l in self.model_opts["likelihoods"] ]
    else:
      self.data_opts['scale_features'] = [ False for l in self.model_opts["likelihoods"] ]


    # Data processing: mask values
    if maskAtRandom is not None:
      self.data_opts['maskAtRandom'] = data_opts['maskAtRandom']
      assert len(self.data_opts['maskAtRandom'])==M, "Length of MaskAtRandom and input files does not match"
    else:
      self.data_opts['maskAtRandom'] = [0]*M

    if maskNSamples is not None:
      self.data_opts['maskNSamples'] = data_opts['maskNSamples']
      assert len(self.data_opts['maskNSamples'])==M, "Length of MaskAtRandom and input files does not match"
    else:
      self.data_opts['maskNSamples'] = [0]*M

    # Remove incomplete samples?
    self.data_opts['RemoveIncompleteSamples'] = RemoveIncompleteSamples

  def load_data(self):
    """ Load the data """

    # Load observations
    self.data = loadData(self.data_opts)

    # Remove samples with missing views
    if self.data_opts['RemoveIncompleteSamples']:
      self.data = removeIncompleteSamples(self.data)

    # Calculate dimensionalities
    M = self.dimensionalities["M"]
    N = self.dimensionalities["N"] = self.data[0].shape[0]
    D = self.dimensionalities["D"] = [self.data[m].shape[1] for m in range(M)]
    
    # Load covariates (NOT IMPLEMENTED)
    # if self.data_opts['covariatesFile'] is not None:

    #   self.data_opts['covariates'] = pd.read_csv(self.data_opts['covariatesFile'], delimiter=" ", header=None).as_matrix()
    #   print("Loaded covariates from " + self.data_opts['covariatesFile'] + "with shape " + str(self.data_opts['covariates'].shape) + "...")

    #   # Scale covariates
    #   self.data_opts['scale_covariates'] = self.data_opts['scale_covariates']
    #   if len(self.data_opts['scale_covariates']) == 1 and self.data_opts['covariates'].shape[1] > 1:
    #     self.data_opts['scale_covariates'] = self.data_opts['scalecovariates'][0] * s.ones(self.data_opts['covariates'].shape[1])
    #   elif type(self.data_opts['scale_covariates'])==list:
    #     assert len(self.data_opts['scale_covariates']) == self.data_opts['covariates'].shape[1], "'scale_covariates' has to be the same length as the number of covariates"
    #   self.data_opts['scale_covariates'] = [ bool(x) for x in self.data_opts['scale_covariates'] ]

    #   # Add the number of covariates to the total number of factors
    #   self.model_opts['factors'] += self.data_opts['covariates'].shape[1]

    #   # Parse covariates
    #   self.parse_covariates()

    # else:
    #   self.data_opts['scale_covariates'] = False
    #   self.data_opts['covariates'] = None

    self.data_opts['covariates'] = None
    self.data_opts['scale_covariates'] = False

  def define_priors(self):
    """ Define priors of the model"""

    N = self.dimensionalities["N"]
    K = self.dimensionalities["K"]
    M = self.dimensionalities["M"]
    D = self.dimensionalities["D"]
    
    # Latent Variables
    self.model_opts["priorZ"] = { 'mean':s.zeros((N,K)) }
    self.model_opts["priorZ"]['var'] = s.ones((K,))*1.

    # Weights
    self.model_opts["priorSW"] = { 'Theta':[s.nan]*M, 'mean_S0':[s.nan]*M, 'var_S0':[s.nan]*M, 'mean_S1':[s.nan]*M, 'var_S1':[s.nan]*M } # Not required
    self.model_opts["priorAlphaW"] = { 'a':[s.ones(K)*1e-14]*M, 'b':[s.ones(K)*1e-14]*M }

    # Theta
    self.model_opts["priorTheta"] = { 'a':[s.ones(K,) for m in range(M)], 'b':[s.ones(K,) for m in range(M)] }
    for m in range(M):
      for k in range(K):
        if self.model_opts['sparsity'][m][k]==0:
          self.model_opts["priorTheta"]["a"][m][k] = s.nan
          self.model_opts["priorTheta"]["b"][m][k] = s.nan

    # Tau
    self.model_opts["priorTau"] = { 'a':[s.ones(D[m])*1e-14 for m in range(M)], 'b':[s.ones(D[m])*1e-14 for m in range(M)] }

  def initialise_variational(self, initTheta=1.):
    """ Initialise variational distributions of the model"""

    N = self.dimensionalities["N"]
    K = self.dimensionalities["K"]
    M = self.dimensionalities["M"]
    D = self.dimensionalities["D"]

    # Latent variables
    self.model_opts["initZ"] = { 'mean':"random", 'var':s.ones((K,)), 'E':None, 'E2':None }

    # Tau
    self.model_opts["initTau"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[s.ones(D[m])*100 for m in range(M)] }

    # ARD of weights
    self.model_opts["initAlphaW"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[s.ones(K)*1. for m in range(M)] }

    # Theta
    self.model_opts["initTheta"] = { 'a':[s.ones(K,) for m in range(M)], 'b':[s.ones(K,) for m in range(M)], 'E':[s.nan*s.zeros((D[m],K)) for m in range(M)] }
    if type(initTheta) is float:
      self.model_opts['initTheta']['E'] = [s.ones((D[m],K))*initTheta for m in range(M)]
    # elif type(self.model_opts["initTheta"]) == list:
    #   assert len(self.model_opts["initTheta"]) == M, "--initTheta has to be a binary vector with length number of views"
    #   self.model_opts['initTheta']['E']= [ self.model_opts["initTheta"][m]*s.ones((D[m],K)) for m in range(M) ]
    else:
       print("Error: 'initTheta' must be a float")
       exit()

    for m in range(M):
      for k in range(K):
        if self.model_opts['sparsity'][m][k]==0.:
          self.model_opts["initTheta"]["a"][m][k] = s.nan
          self.model_opts["initTheta"]["b"][m][k] = s.nan

    # Weights
    self.model_opts["initSW"] = { 
      'Theta':[ self.model_opts['initTheta']['E'][m] for m in range(M)],
      'mean_S0':[s.zeros((D[m],K)) for m in range(M)],
      'var_S0':[s.nan*s.ones((D[m],K)) for m in range(M)],
      'mean_S1':[s.zeros((D[m],K)) for m in range(M)],
      # 'mean_S1':[stats.norm.rvs(loc=0, scale=1, size=(D[m],K)) for m in range(M)],
      'var_S1':[s.ones((D[m],K)) for m in range(M)],
      'ES':[None]*M, 'EW_S0':[None]*M, 'EW_S1':[None]*M # It will be calculated from the parameters
    }

  def parse_covariates(self):
    """ Parse covariates """

    K = self.dimensionalities["K"]
    M = self.dimensionalities["M"]

    # Define idx for covariate factors
    idx = range(self.data_opts['covariates'].shape[1])

    # Ignore mean and variance in the prior distribution of Z
    # self.model_opts["priorZ"]["mean"][idx] = s.nan
    self.model_opts["priorZ"]["var"][idx] = s.nan

    # Ignore variance in the variational distribution of Z
    # The mean has been initialised to the covariate values
    self.model_opts["initZ"]["var"][idx] = 0.

  def parse_intercept(self):
    """ Parse intercept factor """

    K = self.dimensionalities["K"]
    M = self.dimensionalities["M"]
    N = self.dimensionalities["N"]

    # If we want to learn the intercept, we add a constant covariate of 1s
    if self.model_opts['learnIntercept']:
      if self.data_opts['covariates'] is not None:
        self.data_opts['covariates'] = s.insert(self.data_opts['covariates'], obj=0, values=1, axis=1)
        self.data_opts['scale_covariates'].insert(0,False)
      else:
        self.data_opts['covariates'] = s.ones((N,1))
        self.data_opts['scale_covariates'] = [False]

      # Parse intercept
      # self.model_opts['factors'] += 1
      # self.dimensionalities["K"] += 1

      # Remove sparsity from the Intercept factor
      # TO-DO: CHECK THAT THE MODEL IS ALREADY NOT SPARSE
      # TO-DO: RECHECK THIS, ITS UGLY
      # stop if not self.model_opts["learnIntercept"] == TRUE

      for m in range(M):

        # Weights
        if self.model_opts["likelihoods"][m]=="gaussian":
          self.model_opts["initSW"]["mean_S1"][m][:,0] = self.data[m].mean(axis=0)
          self.model_opts["initSW"]["var_S1"][m][:,0] = 1e-5

        # Theta
        self.model_opts['sparsity'][m][0] = 0.
        self.model_opts["initSW"]["Theta"][m][:,0] = 1.
        self.model_opts["priorTheta"]['a'][m][0] = s.nan
        self.model_opts["priorTheta"]['b'][m][0] = s.nan
        self.model_opts["initTheta"]["a"][m][0] = s.nan
        self.model_opts["initTheta"]["b"][m][0] = s.nan
        self.model_opts["initTheta"]["E"][m][:,0] = 1.

  def train_model(self):
    """ Train the model """

    sys.stdout.flush()
    runMOFA(self.data, self.data_opts, self.model_opts, self.train_opts, self.train_opts['seed'])
    sys.stdout.flush()

