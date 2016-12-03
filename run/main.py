"""
"""

# Import required modules
import argparse
import os
from time import time
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
def runSingleTrial(data, model_opts, train_opts, seed=None, trial=1):

    # set the seed
    if seed is None:
        seed = int(round(time()*1000)%1e6)
    # s.random.seed(seed)
    print "Running trial number %d with seed %d\n" % (trial,seed)

    ######################
    ## Define the model ##
    ######################

    # Define dimensionalities
    M = len(data)
    N = data[0].shape[0]
    D = s.asarray([ data[m].shape[1] for m in xrange(M) ])
    K = model_opts["k"]

    dim = {'M':M, 'N':N, 'D':D, 'K':K }

    ## Define and initialise the nodes ##

    init = init_scGFA(dim, data, model_opts["likelihood"])

    init.initSW(ptheta=model_opts["priorSW"]["Theta"], pmean_S0=model_opts["priorSW"]["mean_S0"], pvar_S0=model_opts["priorSW"]["var_S0"], pmean_S1=model_opts["priorSW"]["mean_S1"], pvar_S1=model_opts["priorSW"]["var_S1"],
                qtheta=model_opts["initSW"]["Theta"], qmean_S0=model_opts["initSW"]["mean_S0"], qvar_S0=model_opts["initSW"]["var_S0"], qmean_S1=model_opts["initSW"]["mean_S1"], qvar_S1=model_opts["initSW"]["var_S1"],
                qEW_S0=model_opts["initSW"]["EW_S0"], qEW_S1=model_opts["initSW"]["EW_S1"], qES=model_opts["initSW"]["ES"])

    init.initZ(pmean=model_opts["priorZ"]["mean"], pvar=model_opts["priorZ"]["var"],
               qmean=model_opts["initZ"]["mean"], qvar=model_opts["initZ"]["var"], qE=model_opts["initZ"]["E"], qE2=model_opts["initZ"]["E2"])

    init.initAlpha(pa=model_opts["priorAlpha"]['a'], pb=model_opts["priorAlpha"]['b'], 
                   qa=model_opts["initAlpha"]['a'], qb=model_opts["initAlpha"]['b'], qE=model_opts["initAlpha"]['E'])

    init.initTau(pa=model_opts["priorTau"]['a'], pb=model_opts["priorTau"]['b'], 
                 qa=model_opts["initTau"]['a'], qb=model_opts["initTau"]['b'], qE=model_opts["initTau"]['E'])


    if model_opts['learnTheta']:
        init.initThetaLearn(pa=model_opts["priorTheta"]['a'], pb=model_opts["priorTheta"]['b'],
                            qa=model_opts["initTheta"]['a'],  qb=model_opts["initTheta"]['b'], qE=model_opts["initTheta"]['E'])
    else:
        init.initThetaConst(value=model_opts["initTheta"]['value'])

    init.initY()

    # Define the markov blanket of each node
    init.MarkovBlanket()

    # Initialise expectations of the required nodes
    init.initExpectations("Theta")

    ##################################
    ## Add the nodes to the network ##
    ##################################

    # Initialise Bayesian Network
    # net = BayesNet(dim=dim, trial=trial)
    net = BayesNet(dim=dim, trial=trial, schedule=model_opts["schedule"], nodes=init.getNodes(), options=train_opts)

    # Add nodes to the network
    # net.addNodes(Theta=init.Theta, SW=init.SW, Tau=init.Tau, Z=init.Z, Y=init.Y, Alpha=init.Alpha)
    
    # Add training schedule
    # net.setSchedule(model_opts["schedule"])

    # print net.schedule
    # print net.nodes

    # Add training options
    # net.options = train_opts

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

    trials = Parallel(n_jobs=cores, backend="threading")(
        delayed(runSingleTrial)(data,model_opts,train_opts,None,i) for i in xrange(1,train_opts['trials']+1))

    # trials = [0]*train_opts['trials']
    # for i in xrange(1,train_opts['trials']+1):
    #     trials[i-1] = runSingleTrial(data,model_opts,train_opts,None,i)
    # print "Finished"

    #########################
    ## Analyse the results ##
    #########################

    # Select the trial with the best lower bound
    print trials
    lb = map(lambda x: x.getTrainingStats()["elbo"][-1], trials)
    print lb
    print s.argmax(lb)

    # Save the results

if __name__ == '__main__':

    # Define the data options
    data_opts = {}
    data_opts['input_files'] = \
    (
        "/Users/ricard/git/britta/processed_data/joined/expr.txt",
    #     "/Users/ricard/git/britta/processed_data/joined/met1.txt",
    #     "/Users/ricard/git/britta/processed_data/joined/mut.txt"
    )    
    # data_opts['center'] = (True,True,False)
    data_opts['center'] = (True,)
    data_opts['view_names'] = ("Expression",)
    # data_opts['rownames'] = None
    data_opts['rownames'] = 0
    # data_opts['colnames'] = None
    data_opts['colnames'] = 0
    data_opts['delimiter'] = "\t"
    
    # Define the model options
    model_opts = {}
    # model_opts['likelihood'] = ("gaussian","gaussian","bernoulli")
    model_opts['likelihood'] = ("gaussian",)
    model_opts['learnTheta'] = True
    model_opts['k'] = 10
    

    # Define priors
    M = len(data_opts['view_names'])
    model_opts["priorZ"] = { 'mean':0., 'var':1. }
    model_opts["priorAlpha"] = { 'a':[1e-14]*M, 'b':[1e-14]*M }
    model_opts["priorSW"] = { 'Theta':[s.nan]*M, 'mean_S0':[s.nan]*M, 'var_S0':[s.nan]*M, 'mean_S1':[s.nan]*M, 'var_S1':[s.nan]*M }
    model_opts["priorTau"] = { 'a':[1e-14]*M, 'b':[1e-14]*M }
    if model_opts['learnTheta']: model_opts["priorTheta"] = { 'a':[1.]*M, 'b':[1.]*M }

    # Define initialisation options
    model_opts["initZ"] = { 'mean':0., 'var':1., 'E':"random", 'E2':1. }
    model_opts["initAlpha"] = { 'a':[1e-14]*M, 'b':[1e-14]*M, 'E':[100.]*M }
    model_opts["initSW"] = { 'Theta':[0.5]*M, 
                              'mean_S0':[0.]*M, 'var_S0':model_opts["initAlpha"]['E'], 
                              'mean_S1':[0.]*M, 'var_S1':[1.]*M,
                              'ES':[None]*M, 'EW_S0':[None]*M, 'EW_S1':[None]*M}
    model_opts["initTau"] = { 'a':[1.,1.,None], 'b':[1.,1.,None], 'E':[100.,100.,None] }

    if model_opts['learnTheta']: 
        model_opts["initTheta"] = { 'a':[1.]*M, 'b':[1.]*M, 'E':[0.5]*M }
    else:
        model_opts["initTheta"] = { 'value':[0.5]*M }


    # Define schedule of updates
    model_opts['schedule'] = ("Y","SW","Z","Alpha","Tau","Theta")

    # Define the training options
    train_opts = {}
    train_opts['maxiter'] = 3
    train_opts['elbofreq'] = 1
    train_opts['outdir'] = "/tmp/out" 
    train_opts['savefreq'] = 5
    train_opts['savefolder'] = "/tmp/tmp"
    train_opts['verbosity'] = 2
    train_opts['dropK'] = {}
    train_opts['dropK']['by_norm'] = 0.01
    train_opts['dropK']['by_pvar'] = None
    train_opts['dropK']['by_cor'] = 0.80
    train_opts['forceiter'] = False
    train_opts['tolerance'] = 0.01

    # Define the number of trials and cores
    train_opts['trials'] = 2
    cores = 2

    # Go!
    runMultipleTrials(data_opts, model_opts, train_opts, cores)
