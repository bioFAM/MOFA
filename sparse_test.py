"""
Script to test the spike and slab updates with gaussian and non-gaussian likelihoods
"""

# from __future__ import division
from time import time
import cPickle as pkl
import scipy as s
import os
import scipy.special as special
import scipy.stats as stats
import numpy.linalg  as linalg

# Import manually defined functions
from simulate import Simulate
from BayesNet import BayesNet
from multiview_nodes import *
from seeger_nodes import Binomial_PseudoY_Node, Poisson_PseudoY_Node, Bernoulli_PseudoY_Node, Zeta_Node
from local_nodes import Local_Node, Observed_Local_Node
from sparse_updates import Y_Node, Alpha_Node, SW_Node, Tau_Node, Z_Node
from utils import *

import numpy as np

from joblib import Parallel, delayed
import pdb

N_tests = 200

sparsity_truth = s.zeros([3 * N_tests, 8])
sparsity_inferred_opt = s.zeros([3 * N_tests, 12])
sparsity_inferred_no_opt = s.zeros([3 * N_tests, 12])

alpha_opt = s.zeros([3 * N_tests, 12])
alpha_no_opt = s.zeros([3 * N_tests, 12])

# choose whether to put the theta prior factor/view-wise or just view wise
factor_wise = True
ARD_bool = True

"""
Define function to run one single test. Results are saved in global variables
"""
def run_test(test_ix):
    print test_ix
    ###################
    ## Generate data ##
    ###################

    # Define dimensionalities
    M = 3
    N = 100
    D = s.asarray([100, 100, 100])
    K = 6


    ## Simulate data  ##
    data = {}
    tmp = None
    tmp = Simulate(M=M, N=N, D=D, K=K)

    data['Z'] = s.zeros((N, K))
    data['Z'][:, 0] = s.sin(s.arange(N) / (N / 20))
    data['Z'][:, 1] = s.cos(s.arange(N) / (N / 20))
    data['Z'][:, 2] = 2 * (s.arange(N) / N - 0.5)
    data['Z'][:, 3] = stats.norm.rvs(loc=0, scale=1, size=N)
    data['Z'][:, 4] = stats.norm.rvs(loc=0, scale=1, size=N)
    data['Z'][:, 5] = stats.norm.rvs(loc=0, scale=1, size=N)

    # Add a known covariate
    # data['Z'] = s.c_[ data['Z'], s.asarray([True,False]*(N/2), dtype=s.float32) ]
    # covariate = 6

    # data['alpha'] = s.zeros((M,K))
    # data['alpha'][0,:] = [1,1,1e6,1,1e6,1e6]
    # data['alpha'][1,:] = [1,1e6,1,1e6,1,1e6]
    # data['alpha'][2,:] = [1e6,1,1,1e6,1e6,1]

    data['alpha'] = [s.zeros(K, ) for m in xrange(M)]
    data['alpha'][0] = [1, 1, 1e6, 1, 1e6, 1e6]
    data['alpha'][1] = [1, 1e6, 1, 1e6, 1, 1e6]
    data['alpha'][2] = [1e6, 1, 1, 1e6, 1e6, 1]

    # theta = [ s.ones(K)*0.5 for m in xrange(M) ]
    # theta = [ s.ones(K)*0.5 for m in xrange(M) ]
    # if not factor_wise, initialisation should be done as it is, and then
    # rtansformed inside the code to be of the right dimension
    # TODO to change, initialisation and that should be it
    if factor_wise:
        theta_priors = s.random.uniform(0, 1, M*K)
        theta_priors = s.reshape(theta_priors, [M, K])
        theta = [theta_priors[m, :] for m in xrange(M)]
        sparsity_truth[3 * test_ix:(3 * test_ix + 3), 0] = test_ix
        sparsity_truth[3 * test_ix:(3 * test_ix + 3), 1] = [0,1,2]
        sparsity_truth[3 * test_ix:(3 * test_ix + 3), 2:] = theta_priors
        # TODO store the truth in sparsity_truth
    else:
        theta_priors = s.random.uniform(0, 1, 3)
        sparsity_truth[3 * test_ix:(3 * test_ix + 3), 0] = test_ix
        sparsity_truth[3 * test_ix:(3 * test_ix + 3), 1] = [0,1,2]
        sparsity_truth[3 * test_ix:(3 * test_ix + 3), 2:] = s.reshape(theta_priors, [3,1])
        theta = [s.ones(K) * theta_priors[m] for m in xrange(M)]

    data['S'], data['W'], data['W_hat'], _ = tmp.initW_spikeslab(theta=theta, alpha=data['alpha'])

    data['mu'] = [s.zeros(D[m]) for m in xrange(M)]
    data['tau'] = [stats.uniform.rvs(loc=1, scale=5, size=D[m]) for m in xrange(M)]
    Y_gaussian = None
    Y_gaussian = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
                                  likelihood="gaussian", missingness=0.05)
    # Y_poisson = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'], likelihood="poisson")
    # Y_bernoulli = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'], likelihood="bernoulli")
    # Y_binomial = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'], likelihood="binomial",
    # min_trials=10, max_trials=50)

    # data["Y"] = ( Y_gaussian[0], Y_poisson[1], Y_bernoulli[2] )
    # data["Y"] = Y_bernoulli
    # data["Y"] = Y_poisson
    # data["Y"] = Y_binomial
    data["Y"] = Y_gaussian

    # M_bernoulli = [2]
    # M_poisson = [1]
    # M_gaussian = [0]
    # M_binomial = []

    M_bernoulli = []
    M_poisson = []
    M_gaussian = [0, 1, 2]
    M_binomial = []

    # M_bernoulli = []
    # M_poisson = []
    # M_gaussian = []
    # M_binomial = [0,1,2]

    #################################
    ## Initialise Bayesian Network ##
    #################################
    for optimise_theta_bool in [True, False]:
        net = BayesNet(nodes={}, schedule=())

        # Define initial number of latent variables
        K = 10

        # Define model dimensionalities
        net.dim["M"] = M
        net.dim["N"] = N
        net.dim["D"] = D
        net.dim["K"] = K

        ###############
        ## Add nodes ##
        ###############

        # Z without covariates (variational node)
        Z_pmean = 0.
        Z_pvar = 1.
        Z_qmean = s.stats.norm.rvs(loc=0, scale=1, size=(N, K))
        Z_qvar = s.ones((N, K))
        Z = Z_Node(dim=(N, K), pmean=Z_pmean, pvar=Z_pvar, qmean=Z_qmean, qvar=Z_qvar)

        # Z with covariates (variational node)
        # Z_pmean = 0.
        # Z_pvar = 1.
        # Z_qmean = s.stats.norm.rvs(loc=0, scale=1, size=(N,K-1))
        # Z_qmean = s.c_[ Z_qmean, s.asarray([True,False]*(N/2), dtype=s.float32) ]
        # Z_qvar = s.ones((N,K))
        # Z = Z_Node(dim=(N,K), pmean=Z_pmean, pvar=Z_pvar, qmean=Z_qmean, qvar=Z_qvar)
        # Z.setCovariates(idx=K-1)

        # alpha (variational node)
        alpha_list = [None] * M
        alpha_pa = 1e-14
        alpha_pb = 1e-14
        alpha_qb = s.ones(K)
        alpha_qE = s.ones(K) * 100
        for m in xrange(M):
            alpha_qa = alpha_pa + s.ones(K) * D[m] / 2
            alpha_list[m] = Alpha_Node(dim=(K,), pa=alpha_pa, pb=alpha_pb, qa=alpha_qa[m], qb=alpha_qb, qE=alpha_qE)
            alpha_list[m].updateExpectations()
        alpha = Multiview_Variational_Node((K,) * M, *alpha_list)


        # W (variational node)
        SW_list = [None] * M
        # TODO make sure dimension of ptheta fine everywhere inside
        # TODO intitialisaton here didnt change, only inside -> changed to right dimensions
        S_ptheta = 0.5
        for m in xrange(M):
            S_qtheta = s.ones((D[m], K)) * S_ptheta
            W_qmean = s.stats.norm.rvs(loc=0, scale=1, size=(D[m], K))
            W_qvar = s.ones((D[m], K))
            SW_list[m] = SW_Node(dim=(D[m], K), ptheta=S_ptheta, qtheta=S_qtheta,
                                 qmean=W_qmean, qvar=W_qvar,
                                 optimise_theta_bool = optimise_theta_bool,
                                 pi_opt_per_factor = factor_wise)
        SW = Multiview_Variational_Node(M, *SW_list)

        # tau/kappa (mixed node)
        tau_list = [None] * M
        for m in xrange(M):
            if m in M_poisson:
                tmp = 0.25 + 0.17 * s.amax(data["Y"][m], axis=0)
                tau_list[m] = Observed_Local_Node(dim=(D[m],), value=tmp)
            elif m in M_bernoulli:
                tmp = s.ones(D[m]) * 0.25
                tau_list[m] = Observed_Local_Node(dim=(D[m],), value=tmp)
            elif m in M_binomial:
                tmp = 0.25 * s.amax(data["Y"]["tot"][m], axis=0)
                tau_list[m] = Observed_Local_Node(dim=(D[m],), value=tmp)
            elif m in M_gaussian:
                tau_pa = 1e-14
                tau_pb = 1e-14
                tau_qa = tau_pa + s.ones(D[m]) * N / 2
                tau_qb = s.zeros(D[m])
                tau_qE = s.zeros(D[m]) + 100
                tau_list[m] = Tau_Node(dim=(D[m],), pa=tau_pa, pb=tau_pb, qa=tau_qa, qb=tau_qb, qE=tau_qE)
        tau = Multiview_Mixed_Node(M, *tau_list)


        # zeta (local node)
        # not initialised since it is the first update
        Zeta_list = [None] * M
        for m in xrange(M):
            if m not in M_gaussian:
                Zeta_list[m] = Zeta_Node(dim=(N, D[m]), initial_value=None)
            else:
                Zeta_list[m] = None
        Zeta = Multiview_Local_Node(M, *Zeta_list)


        # Y/Yhat (mixed node)
        Y_list = [None] * M
        for m in xrange(M):
            if m in M_gaussian:
                Y_list[m] = Y_Node(dim=(N, D[m]), obs=data['Y'][m])
            elif m in M_poisson:
                Y_list[m] = Poisson_PseudoY_Node(dim=(N, D[m]), obs=data['Y'][m], E=None)
            elif m in M_bernoulli:
                Y_list[m] = Bernoulli_PseudoY_Node(dim=(N, D[m]), obs=data['Y'][m], E=None)
            elif m in M_binomial:
                Y_list[m] = Binomial_PseudoY_Node(dim=(N, D[m]), tot=data['Y']["tot"][m], obs=data['Y']["obs"][m], E=None)
        Y = Multiview_Mixed_Node(M, *Y_list)


        ############################
        ## Define Markov Blankets ##
        ############################

        Z.addMarkovBlanket(SW=SW, tau=tau, Y=Y)
        for m in xrange(M):
            alpha.nodes[m].addMarkovBlanket(SW=SW.nodes[m])
            SW.nodes[m].addMarkovBlanket(Z=Z, tau=tau.nodes[m], alpha=alpha.nodes[m], Y=Y.nodes[m])
            if m in M_gaussian:
                Y.nodes[m].addMarkovBlanket(Z=Z, SW=SW.nodes[m], tau=tau.nodes[m])
                tau.nodes[m].addMarkovBlanket(SW=SW.nodes[m], Z=Z, Y=Y.nodes[m])
            else:
                Zeta.nodes[m].addMarkovBlanket(Z=Z, W=SW.nodes[m])
                Y.nodes[m].addMarkovBlanket(Z=Z, W=SW.nodes[m], kappa=tau.nodes[m], zeta=Zeta.nodes[m])

        ##################################
        ## Update required expectations ##
        ##################################

        SW.updateExpectations()
        Z.Q.updateExpectations()

        ##################################
        ## Add the nodes to the network ##
        ##################################

        net.addNodes(Zeta=Zeta, SW=SW, tau=tau, Z=Z, Y=Y, alpha=alpha)

        # Define update schedule
        if ARD_bool:
            schedule = ["Zeta", "Y", "SW", "Z", "alpha", "tau"]
        else:
            schedule = ["Zeta", "Y", "SW", "Z", "tau"]
        net.setSchedule(schedule)

        #############################
        ## Define training options ##
        #############################

        options = {}
        options['maxiter'] = 100
        options['tolerance'] = 1E-2
        options['forceiter'] = True
        # options['elbofreq'] = options['maxiter']+1
        options['elbofreq'] = 1
        options['dropK'] = True
        options['dropK_threshold'] = 0.01
        options['savefreq'] = options['maxiter'] + 1
        options['savefolder'] = "/tmp/test"
        options['verbosity'] = 0
        net.options = options


        ####################
        ## Start training ##
        ####################

        net.iterate()

        ##################
        ## Save results ##
        ##################


        outdir = "/tmp/test"

        if not os.path.exists(outdir):
            os.makedirs(outdir)
        if not os.path.exists(os.path.join(outdir, "data")):
            os.makedirs(os.path.join(outdir, "data"))
        if not os.path.exists(os.path.join(outdir, "model")):
            os.makedirs(os.path.join(outdir, "model"))
        if not os.path.exists(os.path.join(outdir, "stats")):
            os.makedirs(os.path.join(outdir, "stats"))
        if not os.path.exists(os.path.join(outdir, "opts")):
            os.makedirs(os.path.join(outdir, "opts"))


        if optimise_theta_bool:
            # save the Pi parameters
            sparsity_inferred_opt[3 * test_ix:(3 * test_ix + 3), 0] = test_ix
            sparsity_inferred_opt[3 * test_ix:(3 * test_ix + 3), 1] = [0, 1, 2]
            sparsity_inferred_opt_tmp = np.vstack([sw.P_theta for sw in SW_list])
            sparsity_inferred_opt[3 * test_ix:(3 * test_ix + 3), 2:(sparsity_inferred_opt_tmp.shape[1]+2)] = sparsity_inferred_opt_tmp
            # sparsity_inferred_opt[3 * test_ix:(3 * test_ix + 3), 2] = [sw.P_theta for sw in SW_list]

            # save the alpha parameters
            alpha_opt[3 * test_ix:(3 * test_ix + 3), 0] = test_ix
            alpha_opt[3 * test_ix:(3 * test_ix + 3), 1] = [0, 1, 2]
            alpha_tmp = np.vstack([alpha.getExpectations()['E'] for alpha in alpha_list])
            alpha_opt[3 * test_ix:(3 * test_ix + 3), 2:(alpha_tmp.shape[1]+2)] = alpha_tmp

        else:
            # save the Pi parameters
            sparsity_inferred_no_opt[3 * test_ix:(3 * test_ix + 3), 0] = test_ix
            sparsity_inferred_no_opt[3 * test_ix:(3 * test_ix + 3), 1] = [0, 1, 2]
            sparsity_inferred_no_opt_tmp = np.vstack([sw.P_theta for sw in SW_list])
            sparsity_inferred_no_opt[3 * test_ix:(3 * test_ix + 3), 2:(sparsity_inferred_no_opt_tmp.shape[1]+2)] = sparsity_inferred_no_opt_tmp
            # sparsity_inferred_no_opt[3 * test_ix:(3 * test_ix + 3), 2] = [sw.optimise_theta() for sw in SW_list]

            # save the alpha parameters
            alpha_no_opt[3 * test_ix:(3 * test_ix + 3), 0] = test_ix
            alpha_no_opt[3 * test_ix:(3 * test_ix + 3), 1] = [0, 1, 2]
            alpha_tmp = np.vstack([alpha.getExpectations()['E'] for alpha in alpha_list])
            alpha_no_opt[3 * test_ix:(3 * test_ix + 3), 2:(alpha_tmp.shape[1]+2)] = alpha_tmp



        # Save the model parameters and expectations
        # print "\nSaving model..."
        # file_name = os.path.join(outdir, "model")
        # saveModel(net, outdir=file_name, only_first_moments=True)
        #
        # # Save training statistics
        # print "\nSaving training stats..."
        # file_name = os.path.join(outdir, "stats")
        # opts = saveTrainingStats(model=net, outdir=file_name)
        #
        # # Save training options
        # print "\nSaving training opts..."
        # file_name = os.path.join(outdir, "opts")
        # opts = saveTrainingOpts(model=net,
        #                         outdir=file_name)  # saving real values of sparsity and infered values
        #
        # del net, Zeta, SW, tau, Z, Y, alpha


"""
Running all tests in parallel
"""
# Parallel(n_jobs=8)(delayed(run_test)(ix) for ix in range(0, N_tests))

"""
Running tests sequentiallly
"""
for ix in range(0,N_tests):
    run_test(ix)

# test_dir = '/Users/damienarnol1/Documents/local/pro/PhD/scGFA_all/tests/model_100/'
# test_dir = '/Users/damienarnol1/Documents/local/pro/PhD/scGFA_all/tests/test_100_per_factor_What/'
test_dir = '/tmp/'

if factor_wise:
    factors = ' '.join(['factor' + str(ix) for ix in range(10)])
    file_header = 'run view ' + factors
else:
    file_header = 'run factor res'

# if factor_wise:
#     factors = ' '.join(['factor' + str(ix) for ix in range(6)])
#     file_header = 'run view ' + factors
# else:
#     file_header = 'run view res'

s.savetxt(test_dir + 'sparsity_truth.txt',
          sparsity_truth,
          header = file_header, comments='')

# if factor_wise:
#     factors = ' '.join(['inferred_pi' + str(ix) for ix in range(10)])
#     file_header = 'run view ' + factors
# else:
#     file_header = 'run factor true_pi'
s.savetxt(test_dir + 'sparsity_inf_opt.txt',
          sparsity_inferred_opt,
          header = file_header, comments='')

# if factor_wise:
#     factors = ' '.join(['estimated_pi' + str(ix) for ix in range(10)])
#     file_header = 'run view ' + factors
# else:
#     file_header = 'run factor true_pi'
s.savetxt(test_dir + 'sparsity_inf_no_opt.txt',
          sparsity_inferred_no_opt,
          header = file_header, comments='')

# writing alpha parameters into a file
factors = ' '.join(['factor' + str(ix) for ix in range(10)])
file_header = 'run view ' + factors
s.savetxt(test_dir + 'alpha_no_opt.txt',
          alpha_no_opt,
          header = file_header, comments='')

factors = ' '.join(['factor' + str(ix) for ix in range(10)])
file_header = 'run view ' + factors
s.savetxt(test_dir + 'alpha_opt.txt',
          alpha_opt,
          header = file_header, comments='')
