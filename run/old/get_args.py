def get_args():
    """ Function to parse the arguments """
    parser = argparse.ArgumentParser()

    # I/O options
    parser.add_argument("-Y","--inputfiles", type=str, nargs="?", help="Input view files")
    parser.add_argument("-o","--outdir", type=str, help="Output folder to save trained model")
    parser.add_argument("-savefreq","--savefreq", type=int, help="Frequency of saving a temporary copy of the model", default=100000)
    parser.add_argument("-savefolder","--savefolder", type=str, help="Folder to save the temporary copies of the model", default="/tmp/scGFA")

    # Data options
    parser.add_argument("-center","--center", type=bool, help="Center the data (column-wise)", default=True)
    parser.add_argument("-s","--view_names", type=str, nargs="+", help="Names of the views")
    parser.add_argument("-l","--likelihood", type=str, nargs="+", help="Likelihood for each view")
    
    # Training options
    parser.add_argument("-n","--trials", type=int, help="Number of trials", default=1)
    parser.add_argument("-i","--maxiter", type=int, help="Number of iterations", default=100)
    parser.add_argument("-v","--verbosity", type=int, help="Verbosity level (0,1,2)", default=1)
    parser.add_argument("-dropK","--dropK", type=bool, help="Drop inactive latent variables?", default=True)
    parser.add_argument("-dropK_threshold","--dropK_threshold", type=float, help="DropK threshold (specify)", default=0.01)
    parser.add_argument("-f","--forceiter", type=bool, help="Keep iterating after convergence criterion is met?", default=True)
    parser.add_argument("-c","--cores", type=int, help="Number of cores", default=1)
    parser.add_argument("-elbofreq","--elbofreq", type=int, help="Frequency of lower bound calculation", default=1)
    
    # Model options
    parser.add_argument("-sparse","--sparse", type=bool, help="Use element-wise sparsity (spike and slab)?", default=False)
    parser.add_argument("-k", type=int, help="Initial number of latent variables", default=100)
    
    # Read arguments
    args = vars(parser.parse_args())

    return args
