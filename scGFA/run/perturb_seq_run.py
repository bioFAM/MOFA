
# Import required modules
import argparse
import os
from time import time
import pandas as pd
import scipy as s
from sys import path
from joblib import Parallel, delayed
from socket import gethostname

# Import manual functions
from init_nodes import *
from scGFA.core.BayesNet import BayesNet


def pprint(d, indent=0):
    for key, value in d.iteritems():
        print '\t' * indent + str(key)
        if isinstance(value, dict):
            pprint(value, indent+1)
        else:
            print '\t' * (indent+1) + str(value)

# Function to load the data
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

        tmp = pd.read_csv(file, delimiter=data_opts["delimiter"])
        print "Loaded %s with dim (%d,%d)..." % (file, tmp.shape[0], tmp.shape[1])

        # transpose dataframe
        if data_opts['transpose'][m]:
            tmp = tmp.transpose()
        # Center the data
        if data_opts['center'][m]:
            tmp = (tmp - tmp.mean())

        Y.append(tmp)
    return Y

# Function to run a single trial of the model
def runSingleTrial(data, model_opts, train_opts, seed=None, trial=1, verbose=False):

    # set the seed
    if seed is None:
        seed = int(round(time()*1000)%1e6)
    # s.random.seed(seed)

    print "\n"
    print "#"*45
    print "## Running trial number %d with seed %d ##" % (trial,seed)
    print "#"*45
    print "\n"


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

    if verbose: print "Initialising nodes...\n"

    init = init_scGFA(dim, data, model_opts["likelihood"], seed=seed)

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

    # needs to init cluster nodes too even when there is no kown clusters
    init.initClusters()


    if model_opts['learnTheta']:
        init.initThetaLearn(pa=model_opts["priorTheta"]['a'], pb=model_opts["priorTheta"]['b'],
                            qa=model_opts["initTheta"]['a'],  qb=model_opts["initTheta"]['b'], qE=model_opts["initTheta"]['E'])
    else:
        init.initThetaConst(value=model_opts["initTheta"]['value'])

    init.initY()

    # Define the markov blanket of each node
    print "Defining Markov Blankets...\n"
    init.MarkovBlanket()

    # Initialise expectations of the required nodes
    # TO-DO we have to do something with this.....
    init.initExpectations("Theta")

    ##################################
    ## Add the nodes to the network ##
    ##################################

    # Initialise Bayesian Network
    print "Initialising Bayesian network...\n"
    net = BayesNet(dim=dim, trial=trial, schedule=model_opts["schedule"], nodes=init.getNodes(), options=train_opts)

    ####################
    ## Start training ##
    ####################

    print "Starting training...\n"

    net.iterate()

    #####################
    ## Finish training ##
    #####################

    return net

# Function to run multiple trials of the model
def runMultipleTrials(data_opts, model_opts, train_opts, cores, keep_best_run, verbose=True):

    # If it doesnt exist, create the output folder
    outdir = os.path.dirname(train_opts['outfile'])
    if not os.path.exists(outdir): os.makedirs(outdir)

    ###################
    ## Load the data ##
    ###################

    data = loadData(data_opts, verbose)

    #########################
    ## Run parallel trials ##
    ########################

    seed = None

    trained_models = []
    for i in range(train_opts['trials']):
        trained_models.append(runSingleTrial(data,model_opts,train_opts,seed,i+1,verbose))
    # trained_models = Parallel(n_jobs=cores, backend="threading")(
    #     delayed(runSingleTrial)(data,model_opts,train_opts,seed,i,verbose) for i in xrange(1,train_opts['trials']+1))

    print "\n"
    print "#"*43
    print "## Training finished, processing results ##"
    print "#"*43
    print "\n"

    #####################
    ## Process results ##
    #####################


    # Select the trial with the best lower bound or keep all models
    if keep_best_run:
        lb = map(lambda x: x.getTrainingStats()["elbo"][-1], trained_models)
        save_models = [ trials[s.argmax(lb)] ]
        outfiles = [ train_opts['outfile'] ]
    else:
        save_models = trained_models
        tmp = os.path.splitext(train_opts['outfile'])
        outfiles = [ tmp[0]+str(t)+tmp[1]for t in xrange(train_opts['trials']) ]

    # Save the results
    sample_names = data[0].index.tolist()
    feature_names = [  data[m].columns.values.tolist() for m in xrange(len(data)) ]
    for t in xrange(len(save_models)):
        print "Saving model %d in %s...\n" % (t,outfiles[t])
        saveModel(save_models[t], outfile=outfiles[t], view_names=data_opts['view_names'],
            sample_names=sample_names, feature_names=feature_names)



if __name__ == '__main__':

    #############################
    ## Define the data options ##
    #############################

    data_opts = {}

    # if 'Kvothe' in gethostname():
    #     base_folder = "/Users/ricard/git/gastrulation/expr/scGFA/expr/filt_data"
    # elif 'yoda' in gethostname():
    #     base_folder = ""
    # else:
    #     print "Computer not recognised"
    #     exit()
    base_folder = '/homes/arnol/multi_view_FA/perturb_seq/data/concatenated_bmdc/'

    data_opts['view_names'] = ( "gene_expression", "guide")

    data_opts['input_files'] = [base_folder+'/dc_both_filt_fix_tp10k.txt',
                                base_folder+'/all_guides.txt']
    M = len(data_opts['input_files'])
    data_opts['center'] = [True]*M
    # data_opts['rownames'] = 0
    # data_opts['colnames'] = 0
    data_opts['delimiter'] = "\t"
    data_opts['transpose'] = [True]*M

    # pprint(data_opts)
    # print "\n"

    ##############################
    ## Define the model options ##
    ##############################

    model_opts = {}
    model_opts['likelihood'] = ["gaussian"]*M
    model_opts['learnTheta'] = False
    model_opts['k'] = 10


    # Define priors
    model_opts["priorZ"] = { 'mean':0., 'var':1. }
    model_opts["priorAlpha"] = { 'a':[1e-5]*M, 'b':[1e-5]*M }
    model_opts["priorSW"] = { 'Theta':[s.nan]*M, 'mean_S0':[s.nan]*M, 'var_S0':[s.nan]*M, 'mean_S1':[s.nan]*M, 'var_S1':[s.nan]*M }
    model_opts["priorTau"] = { 'a':[1e-5]*M, 'b':[1e-5]*M }
    if model_opts['learnTheta']:
        model_opts["priorTheta"] = { 'a':[1.]*M, 'b':[1.]*M }

    # Define initialisation options
    model_opts["initZ"] = { 'mean':"random", 'var':1., 'E':None, 'E2':None }
    model_opts["initAlpha"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[100.]*M }
    model_opts["initSW"] = { 'Theta':[0.5]*M,
                              'mean_S0':[0.]*M, 'var_S0':model_opts["initAlpha"]['E'],
                              'mean_S1':["random"]*M, 'var_S1':[1.]*M,
                              'ES':[None]*M, 'EW_S0':[None]*M, 'EW_S1':[None]*M}
    # model_opts["initTau"] = { 'a':[1.,1.,None], 'b':[1.,1.,None], 'E':[100.,100.,None] }
    model_opts["initTau"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[100.]*M }

    if model_opts['learnTheta']:
        model_opts["initTheta"] = { 'a':[1.]*M, 'b':[1.]*M, 'E':[None]*M }
    else:
        model_opts["initTheta"] = { 'value':[0.5]*M }


    # Define schedule of updates
    model_opts['schedule'] = ("SW","Z","Alpha","Tau","Theta")

    # pprint(model_opts)
    # print "\n"

    #################################
    ## Define the training options ##
    #################################

    train_opts = {}
    train_opts['maxiter'] = 300
    train_opts['elbofreq'] = 1
    # if 'Kvothe' in gethostname():
    #     train_opts['outfile'] = "/Users/ricard/git/gastrulation/expr/scGFA/expr/out/singleview.hdf5"
    # elif 'yoda' in gethostname():
    #     train_opts['outfile'] = ""
    train_opts['outfile'] = '/Users/damienarnol1/Documents/local/pro/PhD/perturb_seq/output/res.h5'
    train_opts['savefreq'] = s.nan
    train_opts['savefolder'] = s.nan
    train_opts['verbosity'] = 2
    train_opts['dropK'] = { "by_norm":0.01, "by_pvar":None, "by_cor":0.80 }
    train_opts['forceiter'] = True
    train_opts['tolerance'] = 0.01

    # model_opts['covariates'] = pd.read_csv("%s/covariates.txt" % base_folder, delimiter="\t", header=0, index_col=0)
    model_opts['covariates'] = None

    # Define the number of trials and cores
    train_opts['trials'] = 1
    cores = 1
    keep_best_run = False

    # pprint(data_opts)
    # print "\n"

    # print "Cores: %d " % cores
    # print "\n"

    # Go!
    runMultipleTrials(data_opts, model_opts, train_opts, cores, keep_best_run)
