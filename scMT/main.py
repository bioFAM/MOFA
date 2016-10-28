
"""
Script to run scGFA
"""

# Import required modules
import argparse
import os

# Import manual functions

"""
Format of input files
alpha with datavar
to-do: check order updates
initialise pca
save output properly
input data
"""

def get_args():
    """ Function to parse the arguments """
    parser = argparse.ArgumentParser()

    # I/O options
    parser.add_argument("-m","-met_inputfiles", type=str, nargs="?", help="Input methylation files")
    parser.add_argument("-e","-expr_inputfiles", type=str, nargs="?", help="Input expression files")
    parser.add_argument("-o","-outdir", type=str, help="Output folder to save trained model")

    # Model options
    parser.add_argument("-view_names", type=str, nargs="?", help="Names of the views")
    parser.add_argument("-K", type=int, help="Initial number of latent variables", default=100)
    parser.add_argument("-n","--n_trials", type=int, help="Number of trials", default=1)
    parser.add_argument("-i","--iter", type=int, help="Number of iterations", default=100)
    parser.add_argument("-v","--verbosity", type=int, help="Verbosity level (0,1,2)", default=1)

    # 
    parser.add_argument("-c","--cores", type=int, help="Number of cores", default=1)
    
    # Read arguments
    args = vars(parser.parse_args())

    return args


def main(options):
    # - options: dictionary with the options (given in the arguments)

    # Create the output folders
    if not os.path.exists(model_options['outdir']):
        os.makedirs(model_options['outdir'])
    if not os.path.exists(os.path.join(model_options['outdir'],"data")):
        os.makedirs(os.path.join(model_options['outdir'],"data"))
    if not os.path.exists(os.path.join(model_options['outdir'],"model")):
            os.makedirs(os.path.join(model_options['outdir'],"model"))

    ###################
    ## Load the data ##
    ###################

    print "Loading the data..."

    # Load expression data
    e = list()
    for file in options["expr_inputfiles"]:
        e.append( np.load(file(file,"rb")) )

    # Collect everything into a single list
    Y = e+m

    #####################
    ## Filter the data ##
    #####################

    ######################
    ## Define the model ##
    ######################

    # Define dimensionalities
    dim = {}
    M = 
    N = 
    D = 
    K = 

    # Define the nodes



    ############################
    ## Define Markov Blankets ##
    ############################

    Z.addMarkovBlanket(W=W, tau=tau, Y=Y)
    for m in xrange(dim["M"]):
        alpha.nodes[m].addMarkovBlanket(W=W.nodes[m])
        W.nodes[m].addMarkovBlanket(Z=Z, tau=tau.nodes[m], alpha=alpha.nodes[m], Y=Y.nodes[m])
        if m in M_gaussian:
            Y.nodes[m].addMarkovBlanket(Z=Z, W=W.nodes[m], tau=tau.nodes[m])
            tau.nodes[m].addMarkovBlanket(W=W.nodes[m], Z=Z, Y=Y.nodes[m])
        else:
            Zeta.nodes[m].addMarkovBlanket(Z=Z, W=W.nodes[m])
            Y.nodes[m].addMarkovBlanket(Z=Z, W=W.nodes[m], kappa=tau.nodes[m], zeta=Zeta.nodes[m])



if __name__ == '__main__':

    # Get model options specified by the arguments
    model_options = get_args()
    print "Model options:"
    pprint(model_options)
    print "\n\n"

    # Go!
    main(model_options)
