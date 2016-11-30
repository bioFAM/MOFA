"""
"""

# Import required modules
import argparse
import os
import pandas as pd
import scipy as s
from pprint import pprint
from sys import path
from joblib import Parallel, delayed

# Import manual functions
path.insert(0,"../")
from init_nodes import *
from BayesNet import BayesNet


# Function to load the data
def loadData(data_opts):
    Y = list()
    for m in xrange(len(data_opts['input_files'])):
        file = data_opts['input_files'][m]

        # Read file (with row and column names)
        tmp = pd.read_csv(file, delimiter=data_opts["delimiter"], header=data_opts["colnames"], index_col=data_opts["rownames"])

        # Center the data
        if data_opts['center'][m]: 
            tmp = (tmp - tmp.mean())

        Y.append(tmp)
    return Y

# Function to run a single trial of the model
def runSingleTrial(data, model_opts, train_opts, seed=None):

    # set the seed
    s.random.seed(seed)

    ######################
    ## Define the model ##
    ######################

    # Define dimensionalities
    M = len(data)
    N = data[0].shape[0]
    D = s.asarray([ data[m].shape[1] for m in xrange(M) ])
    K = model_opts["k"]

    dim = {'M':M, 'N':N, 'D':D, 'K':K }

    # Define and initialise the nodes

    init = init_scGFA(dim, data, model_opts["likelihood"])

    init.initSW(ptheta=model_opts["prior_SW"]["theta"], pmean=model_opts["prior_SW"]["mean"], pvar=model_opts["prior_SW"]["var"],
                qtheta=model_opts["init_SW"]["theta"],  qmean=model_opts["init_SW"]["mean"],  qvar=model_opts["init_SW"]["var"])

    init.initZ(type="random", pmean=model_opts["init_Z"]["mean"], pvar=model_opts["init_Z"]["var"])
    exit()
    init.initAlpha(pa=model_opts["prior_alpha"]['a'], pb=model_opts["prior_alpha"]['b'], 
                   qa=model_opts["init_alpha"]['a'], qb=model_opts["init_alpha"]['b'], qE=model_opts["init_alpha"]['E'])
    init.initTau(pa=model_opts["prior_tau"]['a'], pb=model_opts["prior_tau"]['b'], 
                 qa=model_opts["init_tau"]['a'], qb=model_opts["init_tau"]['b'], qE=model_opts["init_tau"]['E'])

    if model_opts['learn_theta']:
        init.initThetaLearn(pa=model_opts["prior_theta"]['a'], pb=model_opts["prior_theta"]['b'], 
                            qa=model_opts["init_theta"]['a'], qb=model_opts["init_theta"]['b'], qE=model_opts["init_theta"]['E'])

    else:
        init.initThetaConst(value=model_opts["init_theta"]['value'])

    init.initY()

    init.MarkovBlanket()


    ##################################
    ## Add the nodes to the network ##
    ##################################

    # Initialise Bayesian Network
    net = BayesNet(dim=dim)

    # Initialise sparse model
    net.addNodes(Theta=init.Theta, SW=init.SW, tau=init.Tau, Z=init.Z, Y=init.Y, alpha=init.Alpha)

    # Add training schedule
    net.setSchedule(model_opts["schedule"])

    # Add training options
    net.options = train_opts

    ####################
    ## Start training ##
    ####################

    net.iterate()

    return net

# Function to run multiple trials of the model
def runMultipleTrials(data_opts, model_opts, train_opts, cores):
    
    # Create the output folders
    if not os.path.exists(train_opts['outdir']):
        os.makedirs(train_opts['outdir'])
    if not os.path.exists(os.path.join(train_opts['outdir'],"data")):
        os.makedirs(os.path.join(train_opts['outdir'],"data"))
    if not os.path.exists(os.path.join(train_opts['outdir'],"model")):
            os.makedirs(os.path.join(train_opts['outdir'],"model"))
    if not os.path.exists(os.path.join(train_opts['outdir'],"stats")):
            os.makedirs(os.path.join(train_opts['outdir'],"stats"))
    if not os.path.exists(os.path.join(train_opts['outdir'],"opts")):
            os.makedirs(os.path.join(train_opts['outdir'],"opts"))

    ###################
    ## Load the data ##
    ###################

    data = loadData(data_opts)

    ###################
    ## Run the model ##
    ###################

    trials = Parallel(n_jobs=cores)(delayed(runSingleTrial)(data,model_opts,train_opts) for i in xrange(train_opts['trials']))

    #########################
    ## Analyse the results ##
    #########################


if __name__ == '__main__':

    # Define the data options
    data_opts = {}
    data_opts['input_files'] = \
    (
        "/Users/ricard/git/britta/processed_data/joined/expr.txt",
        "/Users/ricard/git/britta/processed_data/joined/met1.txt",
        "/Users/ricard/git/britta/processed_data/joined/mut.txt"
    )    
    # (
        # "/tmp/test0.txt",
        # "/tmp/test1.txt"
    # )
    data_opts['center'] = (True,True,False)
    data_opts['view_names'] = ("Expression","Methylation","Mutation")
    # data_opts['rownames'] = None
    data_opts['rownames'] = 0
    # data_opts['colnames'] = None
    data_opts['colnames'] = 0
    data_opts['delimiter'] = "\t"
    
    # Define the model options
    model_opts = {}
    model_opts['likelihood'] = ("gaussian","gaussian","bernoulli")
    model_opts['learn_theta'] = True
    model_opts['k'] = 10
    

    # Define priors
    M = len(data_opts['view_names'])
    model_opts["prior_Z"] = { 'mean':[0.]*M, 'var':[1.]*M }
    model_opts["prior_alpha"] = { 'a':[1e-14]*M, 'b':[1e-14]*M }
    model_opts["prior_SW"] = { 'theta':[.5]*M, 'mean':[0.]*M, 'var':[1.]*M }
    model_opts["prior_tau"] = { 'a':[1e-14]*M, 'var':[1e-14]*M }
    if model_opts['learn_theta']: 
        model_opts["prior_theta"] = { 'a':[1]*M, 'b':[1]*M }

    # Define initialisation options
    model_opts["init_Z"] = { 'mean':[0.]*M, 'var':[1.]*M }
    model_opts["init_alpha"] = { 'a':[1e-14]*M, 'b':[1e-14]*M, 'E':100. }
    model_opts["init_SW"] = { 'theta':[0.5]*M, 'mean':[0.]*M, 'var':[1.]*M }
    model_opts["init_tau"] = { 'a':[1e-14]*M, 'var':[1e-14]*M, 'E':100. }
    if model_opts['learn_theta']: 
        model_opts["init_theta"] = { 'a':[1.]*M, 'b':[1.]*M, 'E':[0.5]*M }
    else:
        model_opts["init_theta"] = { 'value':[0.5]*M }


    # Define schedule of updates
    model_opts['schedule'] = ["Y","SW","Z","alpha","tau","theta"]

    # Define the training options
    train_opts = {}
    train_opts['maxiter'] = 50
    train_opts['elbofreq'] = 1
    train_opts['outdir'] = "/tmp/out" 
    train_opts['savefreq'] = 10 
    train_opts['savefolder'] = "/tmp/tmp"
    train_opts['trials'] = 1
    train_opts['verbosity'] = 1
    train_opts['dropK'] = True
    train_opts['dropK_threshold'] = 0.01
    train_opts['forceiter'] = False
    train_opts['tolerance'] = 0.01

    # Define the number of cores
    cores = 1

    # Go!
    runMultipleTrials(data_opts, model_opts, train_opts, cores)
