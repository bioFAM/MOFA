
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


# Import manual functions
path.insert(0,"../")
from init_nodes import *
from BayesNet import BayesNet

"""
alpha with datavar
to-do: check order updates
initialise pca
save output properly
center the data
somehow store metadata in the model
train with pandas dataframe
"""

def get_args():
    """ Function to parse the arguments """
    parser = argparse.ArgumentParser()

    # I/O options
    parser.add_argument("-m","--met_inputfiles", type=str, nargs="?", help="Input methylation files")
    parser.add_argument("-e","--expr_inputfiles", type=str, nargs="+", help="Input expression files")
    parser.add_argument("-o","--outdir", type=str, help="Output folder to save trained model")
    parser.add_argument("-sparse","--sparse", type=bool, help="Use element-wise sparsity (spike and slab)?", default=False)

    # Model options
    parser.add_argument("-s","--view_names", type=str, nargs="+", help="Names of the views")
    parser.add_argument("-l","--likelihood", type=str, nargs="+", help="Likelihood for each view")
    parser.add_argument("-k", type=int, help="Initial number of latent variables", default=100)
    # parser.add_argument("-n","--n_trials", type=int, help="Number of trials", default=1)
    parser.add_argument("-i","--iter", type=int, help="Number of iterations", default=100)
    parser.add_argument("-v","--verbosity", type=int, help="Verbosity level (0,1,2)", default=1)
    parser.add_argument("-dropK","--dropK", type=bool, help="Drop inactive latent variables?", default=True)
    parser.add_argument("-dropK_threshold","--dropK_threshold", type=float, help="DropK threshold (specify)", default=0.01)
    parser.add_argument("-f","--forceiter", type=bool, help="Keep iterating after convergence criterion is met?", default=True)
    # parser.add_argument("-c","--cores", type=int, help="Number of cores", default=1)
    parser.add_argument("-savefreq","--savefreq", type=int, help="Frequency of saving a temporary copy of the model", default=100000)
    parser.add_argument("-savefolder","--savefolder", type=str, help="Folder to save the temporary copies of the model", default="/tmp/scGFA")
    parser.add_argument("-elbofreq","--elbofreq", type=int, help="Frequency of lower bound calculation", default=1)
    
    # Read arguments
    args = vars(parser.parse_args())

    return args


def main(options):
    # - options: dictionary with the options (given in the arguments)

    # Create the output folders
    # print "Creating output folders...\n"
    if not os.path.exists(model_options['outdir']):
        os.makedirs(model_options['outdir'])
    if not os.path.exists(os.path.join(model_options['outdir'],"data")):
        os.makedirs(os.path.join(model_options['outdir'],"data"))
    if not os.path.exists(os.path.join(model_options['outdir'],"model")):
            os.makedirs(os.path.join(model_options['outdir'],"model"))

    ###################
    ## Load the data ##
    ###################

    # print "Loading the data...\n"
    # options['expr_inputfiles'] = ['/Users/ricard/data/scMT/expr/processed/filt/tmp/e_matrix.txt']
    # options['outdir'] = '/Users/ricard/git/scGFA/scMT/tmp'
    # options['savefolder'] = '/tmp/scGFA'
    # options['K'] = 10
    # options["likelihood"] = ["gaussian"]

    # Load expression data
    e = list()
    # e_meta = list()
    for file in options["expr_inputfiles"]:
        tmp = pd.read_csv(file, sep=' ', header=0, index_col=0)
        # e.append(tmp.as_matrix())
        e.append(tmp)
        # asd = {'genes':tmp.index , 'samples':tmp.columns}
        # print asd
        # exit()
        # e_meta.append()

    # Collect everything into a single list
    # Y = e+m
    data = e

    #####################
    ## Filter the data ##
    #####################

    ######################
    ## Define the model ##
    ######################

    # Define dimensionalities
    M = len(data)
    N = data[0].shape[0]
    D = s.asarray([ data[m].shape[1] for m in xrange(M) ])
    K = options["k"]

    dim = {'M':M, 'N':N, 'D':D, 'K':K }

    # Define and initialise the nodes
    if options["sparse"]:
        init = init_scGFA(dim,data,options["likelihood"])
        init.initSW(S_ptheta=0.5) 
    else:
        init = init_GFA(dim,data,options["likelihood"])
        init.initW() 
    init.initZ(type="random")
    init.initAlpha(pa=1e-14, pb=1e-14, qb=1., qE=1.)
    init.initTau(pa=1e-14, pb=1e-14, qb=1., qE=1.)
    init.initY()
    init.MarkovBlanket()

    ##################################
    ## Add the nodes to the network ##
    ##################################

    # Initialise Bayesian Network
    net = BayesNet(dim=dim)

    if options["sparse"]:
        net.addNodes(SW=init.SW, tau=init.Tau, Z=init.Z, Y=init.Y, alpha=init.Alpha)
        schedule = ["SW","Z","alpha","tau"]
    else:
        net.addNodes(W=init.W, tau=init.Tau, Z=init.Z, Y=init.Y, alpha=init.Alpha)
        schedule = ["W","Z","alpha","tau"]
    net.setSchedule(schedule)

    #############################
    ## Define training options ##
    #############################

    opt = {}
    opt['maxiter'] = options["iter"]
    opt['tolerance'] = 1E-2
    opt['forceiter'] = options["forceiter"]
    opt['elbofreq'] = options["elbofreq"]
    opt['dropK'] = options["dropK"]
    opt['dropK_threshold'] = options["dropK_threshold"]
    opt['savefreq'] = options["savefreq"]
    opt['savefolder'] = options["savefolder"]
    opt['verbosity'] = options["verbosity"]
    net.options = opt


    ####################
    ## Start training ##
    ####################

    net.iterate()

    ##################
    ## Save results ##
    ##################

    saveModel(net, outdir=os.path.join(model_options['outdir'],"model"), compress=True)

    pass

if __name__ == '__main__':

    # Get model options specified by the arguments
    model_options = get_args()
    print "Model options:"
    pprint(model_options)
    print "\n\n"

    # Go!
    main(model_options)
