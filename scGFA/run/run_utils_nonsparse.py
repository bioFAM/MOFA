
import scipy as s
from sys import path
from time import time
import pandas as pd
import numpy as np
from joblib import Parallel, delayed

from init_nodes_nonsparse import *
from scGFA.core.BayesNet import BayesNet

def loadData(data_opts, verbose=True):

    print "\n"
    print "#"*18
    print "## Loading data ##"
    print "#"*18
    print "\n"

    Y = list()
    for m in xrange(len(data_opts['input_files'])):
        file = data_opts['input_files'][m]

        # Read file (with row and column names)
        Y.append( pd.read_csv(file, delimiter=data_opts["delimiter"], header=data_opts["colnames"], index_col=data_opts["rownames"]) )
        print "Loaded %s with dim (%d,%d)..." % (file, Y[m].shape[0], Y[m].shape[1])

        # Center the data
        if data_opts['center'][m]:
            Y[m] = (Y[m] - Y[m].mean())
    return Y

# Function to run a single trial of the model
def runSingleTrial(data, data_opts, model_opts, train_opts, seed=None, trial=1, verbose=False):

    # set the seed
    if seed is None:
        seed = int(round(time()*1000)%1e6)
    s.random.seed(seed)

    print "\n"
    print "#"*45
    print "## Running trial number %d with seed %d ##" % (trial,seed)
    print "#"*45
    print "\n"

    ###########################
    ## Perform sanity checks ##
    ###########################

    if not os.path.isdir(os.path.dirname(data_opts["outfile"])):
        print "Error: Output directory does not exist"
        exit()

    # If it doesnt exist, create the output folder
    # outdir = os.path.dirname(data_opts['outfile'])
    # if not os.path.exists(outdir): os.makedirs(outdir)

    ####################
    ## Parse the data ##
    ####################

    # Mask
    if 'maskAtRandom' in data_opts or 'maskNSamples' in data_opts:
        if any(data_opts['maskAtRandom']) or any(data_opts['maskNSamples']):
            data = maskData(data, data_opts)

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

    # if verbose: print "Initialising nodes...\n"

    init = init_nonsparse(dim, data, model_opts["likelihood"], seed=seed)

    init.initZ(pmean=model_opts["priorZ"]["mean"], pcov=model_opts["priorZ"]["cov"],
               qmean=model_opts["initZ"]["mean"], qcov=model_opts["initZ"]["cov"], qE=model_opts["initZ"]["E"],
               covariates=model_opts['covariates'])

    init.initW(pmean=model_opts["priorW"]["mean"], pcov=model_opts["priorW"]["cov"],
               qmean=model_opts["initW"]["mean"], qcov=model_opts["initW"]["cov"], qE=model_opts["initW"]["E"])

    init.initAlpha(pa=model_opts["priorAlpha"]['a'], pb=model_opts["priorAlpha"]['b'],
               qa=model_opts["initAlpha"]['a'], qb=model_opts["initAlpha"]['b'], qE=model_opts["initAlpha"]['E'])

    init.initTau(pa=model_opts["priorTau"]['a'], pb=model_opts["priorTau"]['b'],
                 qa=model_opts["initTau"]['a'], qb=model_opts["initTau"]['b'], qE=model_opts["initTau"]['E'])

    init.initY()

    # Define the markov blanket of each node
    # print "Defining Markov Blankets...\n"
    init.MarkovBlanket()

    ##################################
    ## Add the nodes to the network ##
    ##################################

    # Initialise Bayesian Network
    # print "Initialising Bayesian network...\n"
    net = BayesNet(dim=dim, trial=trial, schedule=model_opts["schedule"], nodes=init.getNodes(), options=train_opts)

    ####################
    ## Start training ##
    ####################

    # print "Starting training...\n"

    net.iterate()

    #####################
    ## Finish training ##
    #####################

    return net

# Function to run multiple trials of the model
def runMultipleTrials(data, data_opts, model_opts, train_opts, keep_best_run, verbose=True):

    #########################
    ## Run parallel trials ##
    #########################

    trained_models = Parallel(n_jobs=train_opts['cores'])(
        delayed(runSingleTrial)(data,data_opts,model_opts,train_opts,None,i) for i in xrange(1,train_opts['trials']+1))

    print "\n"
    print "#"*43
    print "## Training finished, processing results ##"
    print "#"*43
    print "\n"

    #####################
    ## Process results ##
    #####################


    # Select the trial with the best lower bound or keep all models
    if train_opts['trials'] > 1:
        if keep_best_run:
            lb = map(lambda x: x.getTrainingStats()["elbo"][-1], trained_models)
            save_models = [ trials[s.argmax(lb)] ]
            outfiles = [ data_opts['outfile'] ]
        else:
            save_models = trained_models
            tmp = os.path.splitext(data_opts['outfile'])
            outfiles = [ tmp[0]+"_"+str(t)+tmp[1]for t in xrange(train_opts['trials']) ]
    else:
        save_models = trained_models
        outfiles = [ data_opts['outfile'] ]

    # Save the results
    sample_names = data[0].index.tolist()
    feature_names = [  data[m].columns.values.tolist() for m in xrange(len(data)) ]
    for t in xrange(len(save_models)):
        print "Saving model %d in %s...\n" % (t,outfiles[t])
        saveModel(save_models[t], outfile=outfiles[t], view_names=data_opts['view_names'],
            sample_names=sample_names, feature_names=feature_names, train_opts=train_opts, model_opts=model_opts)