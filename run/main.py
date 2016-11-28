"""
Script to run a single trial of scGFA 
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


# Function to run a single trial of the model
def runSingleTrial(data, model_opts, train_opts, seed=None):

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
    if model_opts["sparse"]:
        init = init_scGFA(dim,data,model_opts["likelihood"])
        init.initSW(S_ptheta=0.5) 
    else:
        init = init_GFA(dim,data,model_opts["likelihood"])
        init.initW() 
    init.initZ(type="random")
    init.initAlpha(pa=1e-5, pb=1e-5, qb=1., qE=100.)
    init.initTau(pa=1e-5, pb=1e-5, qb=1., qE=100.)
    init.initZeta()
    init.initY()
    init.MarkovBlanket()


    ##################################
    ## Add the nodes to the network ##
    ##################################

    # Initialise Bayesian Network
    net = BayesNet(dim=dim)

    # Initialise sparse model
    if model_opts["sparse"]:
        net.addNodes(Zeta=init.Zeta, SW=init.SW, tau=init.Tau, Z=init.Z, Y=init.Y, alpha=init.Alpha)
        # this si wrong, make general
        schedule = ["Zeta","Y","SW","Z","alpha","tau"]

    # Initialise non-sparse model
    # else:
        # net.addNodes(W=init.W, tau=init.Tau, Z=init.Z, Y=init.Y, alpha=init.Alpha)
        # schedule = ["W","Z","alpha","tau"]

    # Add training schedule
    net.setSchedule(schedule)

    # Add training options
    net.options = train_opts

    ####################
    ## Start training ##
    ####################

    net.iterate()

    return net

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
    model_opts['sparse'] = True
    model_opts['k'] = 10
    
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
