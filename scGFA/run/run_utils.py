
import scipy as s
from sys import path
from time import time
import pandas as pd
import numpy as np
from joblib import Parallel, delayed

from init_nodes import *
from scGFA.run.init_nodes import *
from scGFA.core.BayesNet import BayesNet

# def pprint(d, indent=0):
#     for key, value in d.iteritems():
#         print '\t' * indent + str(key)
#         if isinstance(value, dict):
#             pprint(value, indent+1)
#         else:
#             print '\t' * (indent+1) + str(value)

# Function to load the data

def maskData(data, data_opts):

    print "Masking data..."
    
    for m in xrange(len(data)):

        # Mask values at random
        D = data[m].shape[1]
        N = data[m].shape[0]
        p2Mask = data_opts['maskAtRandom'][m]
        if p2Mask != 0:
            idxMask = np.zeros(N*D)
            idxMask[:int(round(N*D*p2Mask))]  = 1
            np.random.shuffle(idxMask)
            idxMask=np.reshape(idxMask, [N, D])
            data[m] = data[m].mask(idxMask==1)

        # Mask samples in a complete view
        Nsamples2Mask = data_opts['maskNSamples'][m]
        if Nsamples2Mask != 0:
            idxMask = np.random.choice(N, size=Nsamples2Mask, replace = False)
            tmp = data[m].copy()
            tmp.ix[idxMask, :] = pd.np.nan
            data[m] = tmp

    return data

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

    init = init_sparse(dim, data, model_opts["likelihood"], seed=seed)

    init.initZ(pmean=model_opts["priorZ"]["mean"], pvar=model_opts["priorZ"]["var"],
               qmean=model_opts["initZ"]["mean"], qvar=model_opts["initZ"]["var"], qE=model_opts["initZ"]["E"], qE2=model_opts["initZ"]["E2"],
               covariates=model_opts['covariates'])


    init.initSW(ptheta=model_opts["priorSW"]["Theta"], pmean_S0=model_opts["priorSW"]["mean_S0"], pvar_S0=model_opts["priorSW"]["var_S0"], pmean_S1=model_opts["priorSW"]["mean_S1"], pvar_S1=model_opts["priorSW"]["var_S1"],
                qtheta=model_opts["initSW"]["Theta"], qmean_S0=model_opts["initSW"]["mean_S0"], qvar_S0=model_opts["initSW"]["var_S0"], qmean_S1=model_opts["initSW"]["mean_S1"], qvar_S1=model_opts["initSW"]["var_S1"],
                qEW_S0=model_opts["initSW"]["EW_S0"], qEW_S1=model_opts["initSW"]["EW_S1"], qES=model_opts["initSW"]["ES"])

    init.initAlpha(pa=model_opts["priorAlpha"]['a'], pb=model_opts["priorAlpha"]['b'],
                   qa=model_opts["initAlpha"]['a'], qb=model_opts["initAlpha"]['b'], qE=model_opts["initAlpha"]['E'])


    init.initTau(pa=model_opts["priorTau"]['a'], pb=model_opts["priorTau"]['b'],
                 qa=model_opts["initTau"]['a'], qb=model_opts["initTau"]['b'], qE=model_opts["initTau"]['E'])

    init.initClusters()

    if model_opts['learnTheta']:
        init.initThetaLearn(pa=model_opts["priorTheta"]['a'], pb=model_opts["priorTheta"]['b'],
                            qa=model_opts["initTheta"]['a'],  qb=model_opts["initTheta"]['b'], qE=model_opts["initTheta"]['E'])
    else:
        init.initThetaConst(value=model_opts["initTheta"]['value'])

    init.initY()

    # Define the markov blanket of each node
    # print "Defining Markov Blankets...\n"
    init.MarkovBlanket()

    # Initialise expectations of the required nodes
    # TO-DO we have to do something with this.....
    init.initExpectations("Theta")

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

    # If it doesnt exist, create the output folder
    outdir = os.path.dirname(data_opts['outfile'])
    if not os.path.exists(outdir): os.makedirs(outdir)

    ###################
    ## Load the data ##
    ###################

    # data = loadData(data_opts, verbose)

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