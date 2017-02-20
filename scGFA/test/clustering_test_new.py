from __future__ import division
import scipy as s
from sys import path
import scipy.stats as stats
import pandas as pd


# Import manually defined functions
from scGFA.core.simulate import Simulate
from scGFA.core.BayesNet import BayesNet
from scGFA.core.multiview_nodes import *
from scGFA.core.nodes import Constant_Node
from scGFA.core.seeger_nodes import Binomial_PseudoY_Node, Poisson_PseudoY_Node, Bernoulli_PseudoY_Node
from scGFA.core.sparse_updates import Y_Node, Alpha_Node, SW_Node, Tau_Node, Z_Node, Theta_Node, Theta_Constant_Node, Cluster_Node_Gaussian
from scGFA.core.utils import *

from scGFA.core.mixed_nodes import Mixed_Theta_Nodes
import pandas as pd


def run_test(use_annotations, seed=None, swap_factor=False):
    if seed is not None:
        np.random.seed(seed)

    """
    define dimensions
    """
    M = 2
    N = 300
    D = s.asarray([150, 200])
    K = 2

    """
    define clusters
    """
    generative_clusters = np.random.choice([0,1], N)

    """
    define data
    """
    data = {}
    tmp = Simulate(M=M, N=N, D=D, K=K)

    """
    define Zs
    """
    # cluster parameters
    means = [4., -4.]
    sigmas = [1.3, 1.3]
    z_score = (means[0] - means[1])/s.sqrt(sigmas[0]**2. + sigmas[1]**2.)

    data['Z'] = s.zeros((N,K))

    cluster_values = s.array([[4,0], [-4,0]])
    tmp_Z_1 = stats.norm.rvs(loc=means[0], scale=sigmas[0], size=N)
    tmp_Z_2 = stats.norm.rvs(loc=means[1], scale=sigmas[1], size=N)
    data['Z'][:,0] = tmp_Z_1 * generative_clusters + tmp_Z_2 *(1-generative_clusters)
    #
    tmp_Z_3 = stats.norm.rvs(loc=0, scale=0.3, size=N)
    tmp_Z_4 = stats.norm.rvs(loc=0, scale=0.3, size=N)
    data['Z'][:,1] = tmp_Z_3 * generative_clusters + tmp_Z_4 *(1-generative_clusters)



    """
    define alphas
    """
    data['alpha'] = [ s.zeros(K,) for m in xrange(M) ]
    data['alpha'][0] = [1,1]
    data['alpha'][1] = [10,10]


    """
    define annotations
    """
    # factor annotated
    annotation_strength = 1e-3
    annotation_theta1 = s.random.choice([annotation_strength, 1.-annotation_strength], p=[0.8, 0.2], size = D[0])
    annotation_theta2 = s.random.choice([annotation_strength, 1.-annotation_strength], p=[0.8, 0.2], size = D[0])
    theta_view_1 = s.concatenate((annotation_theta1[:,None], annotation_theta2[:, None]), axis=1)

    annotation_theta1_2 = s.random.choice([annotation_strength, 1.-annotation_strength], p=[0.8, 0.2], size = D[1])
    annotation_theta2_2 = s.random.choice([annotation_strength, 1.-annotation_strength], p=[0.8, 0.2], size = D[1])
    theta_view_2 = s.concatenate((annotation_theta1_2[:,None], annotation_theta2_2[:, None]), axis=1)

    theta = [theta_view_1, theta_view_2]

    """
    generate data
    """
    data['S'], data['W'], data['W_hat'], _ = tmp.initW_spikeslab(theta=theta, alpha=data['alpha'], annotation=True)
    data['mu'] = [ s.zeros(D[m]) for m in xrange(M)]

    # data['tau']= [ stats.uniform.rvs(loc=1,scale=3,size=D[m]) for m in xrange(M) ]
    data['tau']= [ s.ones(D[m])*2 for m in xrange(M) ]

    Y_gaussian = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
    	likelihood="gaussian", missingness=0.00)

    data["Y"] = Y_gaussian

    likelihood = ["gaussian", "gaussian"]

    view_names = ["foo", "bar"]


    #################################
    ## Initialise Bayesian Network ##
    #################################

    # Define initial number of latent variables
    K = 2

    dim = {}
    dim["M"] = M
    dim["N"] = N
    dim["D"] = D
    dim["K"] = K

    ###############
    ## Add nodes ##
    ###############

    # Z without covariates (variational node)
    Z_pmean = 0.
    Z_pvar = 1.
    # Z_qmean = s.stats.norm.rvs(loc=0, scale=1, size=(N,K))
    Z_qmean = s.zeros((N,K))
    Z_qvar = s.ones((N,K))
    Z = Z_Node(dim=(N,K), pmean=Z_pmean, pvar=Z_pvar, qmean=Z_qmean, qvar=Z_qvar)


    cluster_q_var =1.
    cluster_q_mean =0.
    cluster_p_var =1.
    cluster_p_mean =0.
    # clusters known
    clusters = generative_clusters
    Cluster_Node = Cluster_Node_Gaussian(cluster_p_mean, cluster_p_var, cluster_q_mean,
    					  cluster_q_var, clusters, K)

    # alpha (variational node)
    alpha_list = [None]*M
    alpha_pa = 1e-14
    alpha_pb = 1e-14
    alpha_qb = s.ones(K)
    alpha_qa = s.ones(K)
    alpha_qE = s.ones(K)*100
    for m in xrange(M):
    	alpha_list[m] = Alpha_Node(dim=(K,), pa=alpha_pa, pb=alpha_pb, qa=alpha_qa, qb=alpha_qb, qE=alpha_qE)
    	alpha_qa = alpha_pa + s.ones(K)*D[m]/2
    	# alpha_list[m] = Alpha_Node(dim=(K,), pa=alpha_pa, pb=alpha_pb, qa=alpha_qa[m], qb=alpha_qb, qE=alpha_qE)
    Alpha = Multiview_Variational_Node((K,)*M, *alpha_list)


    # SW (variational node)
    SW_list = [None]*M
    P_theta = s.nan
    P_mean_S1 = s.nan
    P_mean_S0 = s.nan
    P_var_S1 = s.nan
    P_var_S0 = s.nan
    for m in xrange(M):
    	Q_theta = s.ones((D[m],K))*0.5
    	Q_mean_S1 = s.stats.norm.rvs(loc=0, scale=1, size=(D[m],K))
    	Q_mean_S0 = 0.
    	Q_var_S0 = alpha_list[m].getExpectation()
    	Q_var_S1 = s.ones((D[m],K))
    	# SW_list[m] = SW_Node(dim=(D[m],K), pmean=W_pmean, pvar=W_pvar, ptheta=S_ptheta, qmean=W_qmean, qvar=W_qvar, qtheta=S_qtheta)
    	SW_list[m] = SW_Node(dim=(D[m],K),
            pmean_S0=P_mean_S0, pmean_S1=P_mean_S0, pvar_S0=P_var_S0, pvar_S1=P_var_S1, ptheta=P_theta,
            qmean_S0=Q_mean_S0, qmean_S1=Q_mean_S1, qvar_S0=Q_var_S0, qvar_S1=Q_var_S1, qtheta=Q_theta, qEW_S0=None, qEW_S1=None, qES=None)
    SW = Multiview_Variational_Node(M, *SW_list)

    # tau/kappa (mixed node)
    tau_list = [None]*M
    for m in xrange(M):
    	if likelihood[m] == "poisson":
    		tmp = 0.25 + 0.17*s.amax(data["Y"][m],axis=0)
    		tau_list[m] = Constant_Node(dim=(D[m],), value=tmp)
    	elif likelihood[m] == "bernoulli":
    		tmp = s.ones(D[m])*0.25
    		tau_list[m] = Constant_Node(dim=(D[m],), value=tmp)
    	elif likelihood[m] == "binomial":
    		tmp = 0.25*s.amax(data["Y"]["tot"][m],axis=0)
    		tau_list[m] = Constant_Node(dim=(D[m],), value=tmp)
    	elif likelihood[m] == "gaussian":
    		tau_pa = 1e-3
    		tau_pb = 1e-3
    		tau_qa = tau_pa + s.ones(D[m])*N/2
    		tau_qb = s.nan
    		tau_qE = s.zeros(D[m]) + 100
    		tau_list[m] = Tau_Node(dim=(D[m],), pa=tau_pa, pb=tau_pb, qa=tau_qa, qb=tau_qb, qE=tau_qE)
    Tau = Multiview_Mixed_Node(M,*tau_list)

    # Y/Yhat (mixed node)
    Y_list = [None]*M
    for m in xrange(M):
    	if likelihood[m] == "gaussian":
    		Y_list[m] = Y_Node(dim=(N,D[m]), value=data['Y'][m])
    	elif likelihood[m] == "poisson":
    		Y_list[m] = Poisson_PseudoY_Node(dim=(N,D[m]), obs=data['Y'][m], Zeta=None, E=None)
    	elif likelihood[m] == "bernoulli":
    		Y_list[m] = Bernoulli_PseudoY_Node(dim=(N,D[m]), obs=data['Y'][m], Zeta=None, E=None)
    	elif likelihood[m] == "binomial":
    		Y_list[m] = Binomial_PseudoY_Node(dim=(N,D[m]), tot=data['Y']["tot"][m], obs=data['Y']["obs"][m], Zeta=None, E=None)
    Y = Multiview_Mixed_Node(M, *Y_list)

    # Theta
    theta_pa = 1.
    theta_pb = 1.
    theta_qb = 1.
    theta_qa = 1.
    theta_qE = None

    Theta_list = [None for i in range(M)]
    if use_annotations:
        annotated_node = Theta_Constant_Node((D[0],K), theta[0], N)
        Theta_list[0] = annotated_node
    else:
        Theta_list[0] = Theta_Node(dim=(K,), pa=theta_pa, pb=theta_pb, qa=theta_qa, qb=theta_qb, qE=theta_qE)

    Theta_list[1] = Theta_Node(dim=(K,), pa=theta_pa, pb=theta_pb, qa=theta_qa, qb=theta_qb, qE=theta_qE)

    if use_annotations:
        Theta = Multiview_Mixed_Node(M, *Theta_list)
    else:
        Theta = Multiview_Variational_Node(M, *Theta_list)


    ############################
    ## Define Markov Blankets ##
    ############################

    Z.addMarkovBlanket(SW=SW, Tau=Tau, Y=Y, Cluster=Cluster_Node)
    Cluster_Node.addMarkovBlanket(Z=Z)
    for m in xrange(M):
    	Theta.nodes[m].addMarkovBlanket(SW=SW.nodes[m])
    	Alpha.nodes[m].addMarkovBlanket(SW=SW.nodes[m])
    	SW.nodes[m].addMarkovBlanket(Z=Z, Tau=Tau.nodes[m], Alpha=Alpha.nodes[m], Y=Y.nodes[m], Theta=Theta.nodes[m])
    	if likelihood[m] == "gaussian":
    		Y.nodes[m].addMarkovBlanket(Z=Z, SW=SW.nodes[m], Tau=Tau.nodes[m])
    		Tau.nodes[m].addMarkovBlanket(SW=SW.nodes[m], Z=Z, Y=Y.nodes[m])
    	else:
    		Y.nodes[m].addMarkovBlanket(Z=Z, W=SW.nodes[m], kappa=Tau.nodes[m])

    # Update required expectations
    Z.updateExpectations()
    Theta.updateExpectations()

    schedule = ["Z", "SW","Alpha","Tau","Theta","Cluster"]
    nodes = { "Theta":Theta, "SW":SW, "Tau":Tau, "Z":Z, "Y":Y, "Alpha":Alpha,
     		  "Cluster": Cluster_Node}

    #############################
    ## Define training options ##
    #############################

    options = {}
    options['maxiter'] = 500
    options['tolerance'] = 1E-2
    options['forceiter'] = False
    # options['elbofreq'] = options['maxiter']+1
    options['elbofreq'] = 1
    options['dropK'] = { 'by_norm':None, 'by_pvar':None, 'by_cor':None }
    options['savefreq'] = options['maxiter']+1
    options['savefolder'] = "/tmp/test"
    options['verbosity'] = -1


    ####################
    ## Start training ##
    ####################

    net = BayesNet(dim=dim, schedule=schedule, nodes=nodes, options=options)
    net.iterate()

    """
    Build results DataFrame
    """
    # inferred  clusters
    tmp = Cluster_Node.getParameters()
    rescaling_factor = _rescaling_factor(data['W'], SW, K)
    # tmp = test_results['Cluster'].Q.getExpectations()['E'] / rescaling_factor

    # rescale values


    return Cluster_Node, z_score, rescaling_factor

def _rescaling_factor(data_W, res_W, n_k):
    # rescale learnt weight matrix
    n_views = len(res_W.getExpectations())
    rescaling_factor = np.zeros([n_views, n_k])
    for view in range(n_views):
        W = res_W.getExpectations()[view]['ESW']
        W_norm  = (W**2.).sum(axis=0)
        true_W_norm= (data_W[view]**2.).sum(axis=0)
        rescaling_factor[view, :] = np.sqrt(true_W_norm/W_norm)

    # rescaling_factor = np.sqrt(true_W_norm/W_norm)

    return rescaling_factor
