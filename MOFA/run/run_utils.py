
import scipy as s
from sys import path
from time import time
import pandas as pd
import numpy as np
from joblib import Parallel, delayed

from init_nodes import *
from MOFA.run.init_nodes import *
from MOFA.core.BayesNet import BayesNet
from MOFA.core.utils import *

"""
To-do: initialuse MuZ
"""

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

    init = initModel(dim, data, model_opts["likelihood"], seed=seed)

    # Latent variables
    init.initZ(pmean=model_opts["priorZ"]["mean"], pvar=model_opts["priorZ"]["var"],
               qmean=model_opts["initZ"]["mean"], qvar=model_opts["initZ"]["var"], qE=model_opts["initZ"]["E"], qE2=model_opts["initZ"]["E2"],
               covariates=data_opts['covariates'])

    # Sparse weights
    init.initSW(ptheta=model_opts["priorSW"]["Theta"], pmean_S0=model_opts["priorSW"]["mean_S0"], pvar_S0=model_opts["priorSW"]["var_S0"], pmean_S1=model_opts["priorSW"]["mean_S1"], pvar_S1=model_opts["priorSW"]["var_S1"],
                qtheta=model_opts["initSW"]["Theta"], qmean_S0=model_opts["initSW"]["mean_S0"], qvar_S0=model_opts["initSW"]["var_S0"], qmean_S1=model_opts["initSW"]["mean_S1"], qvar_S1=model_opts["initSW"]["var_S1"],
                qEW_S0=model_opts["initSW"]["EW_S0"], qEW_S1=model_opts["initSW"]["EW_S1"], qES=model_opts["initSW"]["ES"])

    # ARD on latent variables
    if model_opts["ardZ"]:
        init.initAlphaZ(pa=model_opts["priorAlphaZ"]['a'], pb=model_opts["priorAlphaZ"]['b'],
                       qa=model_opts["initAlphaZ"]['a'], qb=model_opts["initAlphaZ"]['b'], qE=model_opts["initAlphaZ"]['E'])

    # ARD on weights
    if model_opts["ardW"] == "m":
        init.initAlphaW_m(pa=model_opts["priorAlphaW"]['a'], pb=model_opts["priorAlphaW"]['b'],
                              qa=model_opts["initAlphaW"]['a'], qb=model_opts["initAlphaW"]['b'], qE=model_opts["initAlphaW"]['E'])
    if model_opts["ardW"] == "k":
        init.initAlphaW_k(pa=model_opts["priorAlphaW"]['a'], pb=model_opts["priorAlphaW"]['b'],
                              qa=model_opts["initAlphaW"]['a'], qb=model_opts["initAlphaW"]['b'], qE=model_opts["initAlphaW"]['E'])
    elif model_opts["ardW"] == "mk":
        init.initAlphaW_mk(pa=model_opts["priorAlphaW"]['a'], pb=model_opts["priorAlphaW"]['b'],
                                 qa=model_opts["initAlphaW"]['a'], qb=model_opts["initAlphaW"]['b'], qE=model_opts["initAlphaW"]['E'])

    # Precision of the zero-mean normally-distributed noise
    init.initTau(pa=model_opts["priorTau"]['a'], pb=model_opts["priorTau"]['b'],
                 qa=model_opts["initTau"]['a'], qb=model_opts["initTau"]['b'], qE=model_opts["initTau"]['E'])

    # Sparsity on the weights
    if len(s.unique(model_opts['learnTheta'])) == 1:
        # All are infered
        if s.unique(model_opts['learnTheta'])==1.:
            init.initThetaLearn(pa=model_opts["priorTheta"]['a'], pb=model_opts["priorTheta"]['b'],
                qa=model_opts["initTheta"]['a'],  qb=model_opts["initTheta"]['b'], qE=model_opts["initTheta"]['E'])
        # None are infered
        elif s.unique(model_opts['learnTheta'])==0.:
            init.initThetaConst(value=model_opts["initTheta"]['E'])
    # Some are infered
    else:
        init.initThetaMixed(pa=model_opts["priorTheta"]['a'], pb=model_opts["priorTheta"]['b'],
            qa=model_opts["initTheta"]['a'],  qb=model_opts["initTheta"]['b'], qE=model_opts["initTheta"]['E'],
            learnTheta=model_opts['learnTheta'])

    init.initY()

    # Define the markov blanket of each node
    nodes = init.getNodes()
    if model_opts["ardZ"]:
        nodes["Z"].addMarkovBlanket(SW=nodes["SW"], Tau=nodes["Tau"], Y=nodes["Y"], Alpha=nodes["AlphaZ"])
        nodes["AlphaZ"].addMarkovBlanket(Z=nodes["Z"])
    else:
        nodes["Z"].addMarkovBlanket(SW=nodes["SW"], Tau=nodes["Tau"], Y=nodes["Y"])
    nodes["Theta"].addMarkovBlanket(SW=nodes["SW"])
    nodes["AlphaW"].addMarkovBlanket(SW=nodes["SW"])
    nodes["SW"].addMarkovBlanket(Z=nodes["Z"], Tau=nodes["Tau"], Alpha=nodes["AlphaW"], Y=nodes["Y"], Theta=nodes["Theta"])
    nodes["Y"].addMarkovBlanket(Z=nodes["Z"], SW=nodes["SW"], Tau=nodes["Tau"])
    nodes["Tau"].addMarkovBlanket(Z=nodes["Z"], SW=nodes["SW"], Y=nodes["Y"])


    # Initialise expectations of the required nodes
    # TO-DO we have to do something with this.....
    init.initExpectations("Theta")

    ##################################
    ## Add the nodes to the network ##
    ##################################

    # Initialise Bayesian Network
    net = BayesNet(dim=dim, trial=trial, schedule=model_opts["schedule"], nodes=init.getNodes(), options=train_opts)

    ####################
    ## Start training ##
    ####################

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

    # trained_models = Parallel(n_jobs=train_opts['cores'], backend="threading")(
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
