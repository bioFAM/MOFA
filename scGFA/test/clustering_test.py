import sparse_test
import numpy as np
import scipy.stats as ss

def run_test():
    factor_ranks = np.zeros([100, 3])
    np.random.seed(1)
    for i in range(10):
        test_results = sparse_test.run_test()
        # approximative method to align the learnt latent variables with truth
        # rescale learnt weight matrix
        W = test_results['W'].getExpectations()[0]['ESW']
        W_norm  = (W**2.).sum(axis=0)
        true_W_norm= (test_results['data_W'][0]**2.).sum(axis=0)

        rescaling_factor = np.sqrt(true_W_norm/W_norm)
        # W = W * rescaling_factor

        # permute latent variables so that the latent variable-wise norm computed from


        # rescale cluster mean terms
        # print 'rescaling ', rescaling_factor
        tmp = test_results['Cluster'].Q.getExpectations()['E'] / rescaling_factor
        print tmp
        # factor_ranks[i, :] = ss.rankdata((0.5*abs(tmp)).sum(axis=0))
        # factor_ranks[i,:] = tmp

    # print factor_ranks
    # print factor_ranks.sum(axis=0)

    # import pdb; pdb.set_trace()

def rescale_clustering_values(test_results):
    # rescale learnt weight matrix
    W = test_results['W'].getExpectations()[0]['ESW']
    W_norm  = (W**2.).sum(axis=0)
    true_W_norm= (test_results['data_W'][0]**2.).sum(axis=0)

    rescaling_factor = np.sqrt(true_W_norm/W_norm)

    # rescale cluster mean terms
    # print 'rescaling ', rescaling_factor
    tmp = test_results['Cluster'].Q.getExpectations()['E'] / rescaling_factor
    return tmp
