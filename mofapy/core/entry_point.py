import argparse
import pandas as pd
import scipy as s
import sys
from time import sleep

from .build_model import *
from .utils import *

"""
Entry point for the MOFA model

The order should be:
(1) Set the data (set_data)
(1) Set model options (set_model_options)
(2) Set data options (set_data_options)
(3) Parse the data (parse_data)
(4) Set train options (set_train_options)
(5) define priors and intiialisations (define_priors and define_init)
(6) parse_intercept
(7) Train the model (train_model)
(8) Save the model (save_model)
"""

class entry_point():
  def __init__(self):
    self.print_banner()
    self.dimensionalities = {}

    self.data_opts = {}
    self.model_opts = {}
    self.train_opts = {}

    # Covariates are depreciated
    self.data_opts['covariates'] = None
    self.data_opts['scale_covariates'] = False

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
    ########################################################### 
   \n 
    """

    print(banner)
    # sleep(2)

  def set_data(self, data):
    """ Method to set the data 

    PARAMETERS
    ----------
    data: several options:
    - a dictionary where each key is the view names and the object is a numpy array or a pandas data frame
    - a list where each element is a numpy array or a pandas data frame
    """

    # Sanity check
    if isinstance(data, dict):
      data = list(data.values())
    elif isinstance(data, list):
      pass
    else:
      print("Error: Data not recognised")
      sys.stdout.flush()
      exit()

    assert s.all([ isinstance(data[m],s.ndarray) or isinstance(data[m],pd.DataFrame) for m in range(len(data)) ]), "Error, input data is not a numpy.ndarray"

    # Verbose message
    for m in range(len(data)):
      print("Loaded view %d with %d samples and %d features..." % (m, data[m].shape[0], data[m].shape[1]))

    # Doing QC on the data
    self.data = qcData(data)

    # Save dimensionalities
    self.dimensionalities["M"] = len(self.data)
    self.dimensionalities["N"] = self.data[0].shape[0]
    self.dimensionalities["D"] = [self.data[m].shape[1] for m in range(len(self.data))]


  def parse_data(self):
    """ Method to parse the data """

    # Parsing the data: centering, scaling, etc.
    self.parsed_data = parseData(self.data, self.data_opts)

    # Remove samples with missing views
    if self.data_opts['RemoveIncompleteSamples']:
      self.data = removeIncompleteSamples(self.data)
      self.parsed_data = removeIncompleteSamples(self.parsed_data)
      self.dimensionalities["N"] = self.parsed_data[0].shape[0] # Update dimensionalities

  def set_train_options(self, iter=5000, elbofreq=1, startSparsity=100, tolerance=0.01, 
    startDrop=5, freqDrop=1, endDrop=9999, dropR2=0, nostop=False, verbose=False, seed=None
    ):
    """ Set training options """


    # Maximum number of iterations
    self.train_opts['maxiter'] = int(iter)

    # Lower bound computation frequency
    self.train_opts['elbofreq'] = int(elbofreq)

    # Verbosity
    self.train_opts['verbose'] = verbose

    # Criteria to drop latent variables while training
    self.train_opts['drop'] = { "by_r2":float(dropR2), "by_norm":1e-10 }
    self.train_opts['startdrop'] = int(startDrop)
    self.train_opts['freqdrop'] = int(freqDrop)
    self.train_opts['enddrop'] = int(endDrop)
    # print("\nDropping factors with minimum threshold of {0}% variance explained".format(dropR2*100))


    # Tolerance level for convergence
    self.train_opts['tolerance'] = float(tolerance)

    # Do no stop even when convergence criteria is met
    self.train_opts['forceiter'] = nostop

    # Iteration to activate spike and slab sparsity
    self.train_opts['startSparsity'] = int(startSparsity)
    if hasattr(self, 'model_opts'):
      if self.model_opts["sparsity"] is False:  self.train_opts['startSparsity'] = 999999999

    # Define schedule of updates
    self.train_opts['schedule'] = ( "Y", "SW", "Z", "Alpha", "Theta", "Tau" )

    # Seed
    if seed is None:
      seed = 0
    self.train_opts['seed'] = int(seed)

  def set_model_options(self, factors, likelihoods, sparsity=True, learnIntercept=False):
    """ Set model options 

        PARAMETERS
        ----------
        factors: int
          initial number of factors
        likelihoods: character or list of characters
          likelihood per view. Choose from 'gaussian', 'poisson' and 'bernoulli'
        sparsity: bool
          Fraction of values to mask at random
        learnIntercept: bool
          Number of samples to mask at random
    """

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

    # Define whether to use spike and slab sparsity or not
    if sparsity:
      self.model_opts['sparsity_bool'] = True
      self.model_opts['sparsity'] = [s.ones(K) for m in range(M)]
    else:
      self.model_opts['sparsity_bool'] = False
      print("\nWarning... sparsity is desactivated, we recommend using it\n")
      self.model_opts['sparsity'] = [s.zeros(K) for m in range(M)]
      if hasattr(self, 'train_opts'): self.train_opts['startSparsity'] = 999999999

  def set_data_options(self, view_names=None, center_features=True, scale_features=False, scale_views=False, 
    maskAtRandom=None, maskNSamples=None, RemoveIncompleteSamples=False
    ):

    """ Set data processing options

        PARAMETERS
        ----------
        center_features: bool
          center the features to zero mean?
        scale_features: bool
          scale the features to zero mean?
        scale_views: bool
          scale the features to unit variance?
        maskAtRandom: numeric
          Fraction of values to mask at random
        maskNSamples: bool
          Number of samples to mask at random
        RemoveIncompleteSamples: bool
          Remove samples that are not profiled for all omics
    """
    # Sanity checks
    assert "likelihoods" in self.model_opts, "Likelihoods not found in model options"

    # Sanity checks
    M = self.dimensionalities["M"]
    if view_names is None:
      self.data_opts['view_names'] = ["view_%d" % m for m in range(1,M+1)] 
    else:
      if isinstance(view_names,str): view_names = [ view_names ]
      assert len(view_names)==M, "Length of view names and number of views do not match"
      self.data_opts['view_names'] = view_names

    # Data processing: center features
    if center_features is True:
      self.data_opts['center_features'] = [ True if l=="gaussian" else False for l in self.model_opts["likelihoods"] ]
    else:
      if not self.model_opts["learnIntercept"]: print("\nWarning... you are not centering the data and not learning the intercept. The model is not going to work...\n")
      self.data_opts['center_features'] = [ False for l in self.model_opts["likelihoods"] ]

    # Data processing: scale views
    if scale_views is True:
      print("Warning: you are scaling the Gaussian views to unit variance. As long as the scale of the different views is not massively different, the model does not require this")
      self.data_opts['scale_views'] = [ True if l=="gaussian" else False for l in self.model_opts["likelihoods"] ]
    else:
      self.data_opts['scale_views'] = [ False for l in self.model_opts["likelihoods"] ]

    # Data processing: scale features
    if scale_features:
      print("Warning: you are scaling Gaussian the features to unit variance. This is only recommended if the difference in variances between features are driven by technical and not biological reasons")
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
    # self.model_opts["priorAlpha"] = { 'a':[s.ones(K)*1e-14]*M, 'b':[s.ones(K)*1e-14]*M }
    self.model_opts["priorAlpha"] = { 'a':[s.ones(K)*1e-5]*M, 'b':[s.ones(K)*1e-5]*M }

    # Theta
    self.model_opts["priorTheta"] = { 'a':[s.ones(K,) for m in range(M)], 'b':[s.ones(K,) for m in range(M)] }
    for m in range(M):
      for k in range(K):
        if self.model_opts['sparsity'][m][k]==0:
          self.model_opts["priorTheta"]["a"][m][k] = s.nan
          self.model_opts["priorTheta"]["b"][m][k] = s.nan

    # Tau
    # self.model_opts["priorTau"] = { 'a':[s.ones(D[m])*1e-14 for m in range(M)], 'b':[s.ones(D[m])*1e-14 for m in range(M)] }
    self.model_opts["priorTau"] = { 'a':[s.ones(D[m])*1e-5 for m in range(M)], 'b':[s.ones(D[m])*1e-5 for m in range(M)] }

  def define_init(self, initTheta=1.):
    """ Define Initialisations of the model

        PARAMETERS
        ----------
        initTheta flaot
          initialisation for theta. Default is 1. (no sparsity)
    """

    N = self.dimensionalities["N"]
    K = self.dimensionalities["K"]
    M = self.dimensionalities["M"]
    D = self.dimensionalities["D"]

    # Latent variables
    self.model_opts["initZ"] = { 'mean':"random", 'var':s.ones((K,)), 'E':None, 'E2':None }

    # Tau
    self.model_opts["initTau"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[s.ones(D[m])*100 for m in range(M)] }

    # ARD of weights
    self.model_opts["initAlpha"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[s.ones(K)*1. for m in range(M)] }

    # Theta
    self.model_opts["initTheta"] = { 'a':[s.ones(K,) for m in range(M)], 'b':[s.ones(K,) for m in range(M)], 'E':[s.nan*s.zeros(K,) for m in range(M)] }
    if type(initTheta) is float:
      self.model_opts['initTheta']['E'] = [s.ones(K,)*initTheta for m in range(M)]
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
      'Theta':[ s.repeat(self.model_opts['initTheta']['E'][m][None,:],self.dimensionalities["D"][m],0) for m in range(M)],
      'mean_S0':[s.zeros((D[m],K)) for m in range(M)],
      'var_S0':[s.nan*s.ones((D[m],K)) for m in range(M)],
      'mean_S1':[s.zeros((D[m],K)) for m in range(M)],
      # 'mean_S1':[stats.norm.rvs(loc=0, scale=1, size=(D[m],K)) for m in range(M)],
      'var_S1':[s.ones((D[m],K)) for m in range(M)],
      'ES':[None]*M, 'EW_S0':[None]*M, 'EW_S1':[None]*M # It will be calculated from the parameters
    }

  def parse_covariates(self):
    """ Parse covariates """

    print("Covariates are not functional")
    exit()

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

    # Sanity checks
    # TO-DO: CHECK THAT MODEL_OPTS AND DATA_OPTS ARE PROPERLY DEFINED

    K = self.dimensionalities["K"]
    M = self.dimensionalities["M"]
    N = self.dimensionalities["N"]

    # If we want to learn the intercept, we add a constant covariate of 1s
    if self.model_opts['learnIntercept']:
      if self.data_opts['covariates'] is not None:
        self.data_opts['covariates'] = s.insert(self.data_opts['covariates'], obj=0, values=1., axis=1)
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
        # if self.model_opts["likelihoods"][m]=="gaussian":
        self.model_opts["initSW"]["mean_S1"][m][:,0] = s.nanmean(self.parsed_data[m], axis=0)
        self.model_opts["initSW"]["var_S1"][m][:,0] = 1e-10

        # Theta
        self.model_opts['sparsity'][m][0] = 0.
        self.model_opts["initSW"]["Theta"][m][:,0] = 1.
        self.model_opts["priorTheta"]['a'][m][0] = s.nan
        self.model_opts["priorTheta"]['b'][m][0] = s.nan
        self.model_opts["initTheta"]["a"][m][0] = s.nan
        self.model_opts["initTheta"]["b"][m][0] = s.nan
        self.model_opts["initTheta"]["E"][m][0] = 1.

  def train_model(self):
    """ Train the model """

    # Sanity checks
    # TO-DO: CHECK THAT DATA_OPTS, TRAIN_OPTS AND MODEL_OPTS ARE PROPERLY DEFINED
    assert hasattr(self, 'data'), "Data has to be defined before training the model"

    print("\nLikelihoods are defined as:")
    for a,b in zip(self.data_opts['view_names'], self.model_opts['likelihoods']):
      print("\t%s: %s" % (a,b))

    sys.stdout.flush()
    self.model = runMOFA(self.parsed_data, self.data_opts, self.model_opts, self.train_opts, self.train_opts['seed'])
    sys.stdout.flush()

  def save_model(self, outfile, sample_names=None, feature_names=None):
    """ Save the model """

    # Sanity checks
    assert hasattr(self, 'data'), "Data has to be defined before training the model"
    assert hasattr(self, 'model'), "No trained model found"

    # Create output directory
    if not os.path.isdir(os.path.dirname(outfile)):
        print("Output directory does not exist, creating it...")
        os.makedirs(os.path.dirname(outfile))

    # Verbose messages
    print("Saving model in %s...\n" % outfile)
    if sample_names is None:
      print("Sample names not provided...")
    if feature_names is None:
      print("Feature names not provided...")

    # Save the model
    saveModel(self.model, 
      data=self.data,      # Uncentered data
      outfile=outfile, 
      view_names=self.data_opts['view_names'],
      sample_names=sample_names, 
      feature_names=feature_names, 
      train_opts=self.train_opts, 
      model_opts=self.model_opts
    )

    sys.stdout.flush()
