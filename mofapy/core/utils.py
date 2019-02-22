from __future__ import division
from time import sleep
from copy import deepcopy

import numpy as np
import numpy.ma as ma
import os
import h5py
import sys

"""
Module to define util functions

TO-DO: 
- Create a proper class for saving the models
- Create a proper class for loading and parsing the data
- Move the util math functions into a math.py
"""


def removeIncompleteSamples(data):
    """ Method to remove samples with missing views 

    PARAMETERS
    ----------
    data: list  
    """
    print("\nRemoving incomplete samples...")

    M = len(data)
    N = data[0].shape[0]
    samples_to_remove = []
    for n in range(N):
        for m in range(M):
            if np.all(np.isnan(data[m][n,:])):
                samples_to_remove.append(n)
                break

    if len(samples_to_remove) > 0:
        print("A total of " + str(len(samples_to_remove)) + " sample(s) have at least a missing view and will be removed")

    data_filt = [None]*M
    samples_to_keep = np.setdiff1d(range(N),samples_to_remove)
    for m in range(M):
        data_filt[m] = data[m][samples_to_keep,:]

    return data_filt

def maskData(data, data_opts):
    """ Method to mask values of the data, 
    It is mainly to test missing values and to evaluate imputation
    
    PARAMETERS
    ----------
    data_opts: dic
    
    """
    print("Not functional")
    sys.stdout.flush()
    exit()

    print("Masking data with the following options:")
    print("at random:")
    print(data_opts['maskAtRandom'])
    print("full cases:")
    print(data_opts['maskNSamples'])


    for m in range(len(data)):

        # Mask values at random
        D = data[m].shape[1]
        N = data[m].shape[0]
        p2Mask = data_opts['maskAtRandom'][m]
        if p2Mask != 0:
            idxMask = np.zeros(N*D)
            idxMask[:int(round(N*D*p2Mask))] = 1
            np.random.shuffle(idxMask)
            idxMask = np.reshape(idxMask, [N, D])
            data[m] = data[m].mask(idxMask==1)

        # Mask samples in a complete view
        Nsamples2Mask = data_opts['maskNSamples'][m]
        if Nsamples2Mask != 0:
            idxMask = np.random.choice(N, size=Nsamples2Mask, replace = False)
            # idxMask = np.arange(Nsamples2Mask)
            # print idxMask
            tmp = data[m].copy()
            tmp.ix[idxMask,:] = pd.np.nan
            data[m] = tmp

    return data

def parseData(data, data_opts):
    """ Method to do parse the data
    
    TO-DO: CHECK IF ANY SAMPLE HAS MISSING VALUES IN ALL VIEWS 

    PARAMETERS
    ----------
    data: list of numpy arrays or pandas dataframes
    """
    M = len(data)
    parsed_data = deepcopy(data)
    for m in range(M):
        # Convert to float32
        # parsed_data[m] = parsed_data[m].astype(pd.np.float32)
        parsed_data[m] = parsed_data[m].astype(np.float32)

        # For some reason, reticulate stores missing values in integer matrices as -2147483648
        parsed_data[m][parsed_data[m] == -2147483648] = np.nan
        # Center the features
        if data_opts['center_features'][m]:
            print("Centering features for view " + data_opts["view_names"][m] + "...")
            parsed_data[m] = parsed_data[m] - np.nanmean(parsed_data[m],axis=0)

        # Scale the views to unit variance
        if data_opts['scale_views'][m]:
            print("Scaling view " + data_opts["view_names"][m] + " to unit variance...")
            parsed_data[m] = parsed_data[m] / np.nanstd(parsed_data[m])

        # Scale the features to unit variance
        if data_opts['scale_features'][m]:
            print("Scaling features for view " + data_opts["view_names"][m] + " to unit variance...")
            parsed_data[m] = parsed_data[m] / np.nanstd(parsed_data[m], axis=0, )

    print("\nAfter parsing the data:")
    for m in range(M): print("view %s has %d samples and %d features..." % (data_opts["view_names"][m], parsed_data[m].shape[0], parsed_data[m].shape[1]))

    return parsed_data

def qcData(data):
    """ Method to do quality control on the data
    
    TO-DO: CHECK IF ANY SAMPLE HAS MISSING VALUES IN ALL VIEWS 

    PARAMETERS
    ----------
    data: list of numpy arrays or pandas dataframes
    """

    M = len(data)

    # Check that the dimensions match
    if len(set([data[m].shape[0] for m in range(M)])) != 1:
        if all([data[m].shape[1] for m in range(M)]):
            print("\nWarning: columns seem to be the shared axis, transposing the data...")
            for m in range(M): data[m] = data[m].T
        else:
            print("\nError: Dimensionalities do not match, aborting. Make sure that either columns or rows are shared!")
            sys.stdout.flush()
            exit()


    # Sanity checks on the data
    print ("\n" +"#"*46)
    print("## Doing sanity checks and parsing the data ##")
    print ("#"*46 + "\n")
    for m in range(M):

        # Detect features with complete missing values
        nas = np.isnan(data[m]).mean(axis=0)
        if np.any(nas==1.):
            print("Error: %d features(s) on view %d have missing values in all samples, please remove them before running the model." % ( (nas==1.).sum(), m) )
            sys.stdout.flush()
            exit()

        # Detect features with no variance
        var = data[m].std(axis=0) 
        if np.any(var==0.):
            print("Error: %d features(s) on view %d have zero variance. Please, remove lowly variable features for each omic separately, they can cause numerical issues." % ( (var==0.).sum(),m) )
            sys.stdout.flush()
            exit()

    return data

def dotd(A, B, out=None):
    """Diagonal of :math:`\mathrm A\mathrm B^\intercal`.
    If ``A`` is :math:`n\times p` and ``B`` is :math:`p\times n`, it is done in :math:`O(pn)`.
    Args:
        A (array_like): Left matrix.
        B (array_like): Right matrix.
        out (:class:`numpy.ndarray`, optional): copy result to.
    Returns:
        :class:`numpy.ndarray`: Resulting diagonal.
    """
    A = ma.asarray(A, float)
    B = ma.asarray(B, float)
    if A.ndim == 1 and B.ndim == 1:
        if out is None:
            return ma.dot(A, B)
        return ma.dot(A, B, out)

    if out is None:
        out = ma.empty((A.shape[0], ), float)

    out[:] = ma.sum(A * B.T, axis=1)
    return out

def nans(shape, dtype=float):
    """ Method to create an array filled with missing values """
    a = np.empty(shape, dtype)
    a.fill(np.nan)
    return a

def corr(A,B):
    """ Method to efficiently compute correlation coefficients between two matrices 
    
    PARMETERS
    ---------
    A: np array
    B: np array
    """

    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:,None]
    B_mB = B - B.mean(1)[:,None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1);
    ssB = (B_mB**2).sum(1);

    # Finally get corr coeff
    return np.dot(A_mA,B_mB.T)/np.sqrt(np.dot(ssA[:,None],ssB[None]))

# NOT HERE
def logdet(X):
    return np.log(np.linalg.det(X))
    # UC = np.linalg.cholesky(X)
    # return 2*sum(np.log(np.diag(UC)))


# NOT HERE
def ddot(d, mtx, left=True):
    """Multiply a full matrix by a diagonal matrix.
    This function should always be faster than dot.

    Input:
      d -- 1D (N,) array (contains the diagonal elements)
      mtx -- 2D (N,N) array
      left: is the diagonal matrix on the left or on the right of the product?

    Output:
      ddot(d, mts, left=True) == dot(diag(d), mtx)
      ddot(d, mts, left=False) == dot(mtx, diag(d))
    """
    if left:
        return (d*mtx.T).T
    else:
        return d*mtx


def saveParameters(model, hdf5, view_names=None):
    """ Method to save the parameters of the model in an hdf5 file
    
    PARAMETERS
    ----------
    model: a BayesNet instance
    hdf5: 
    view_names
    """
    
    # Get nodes from the model
    nodes = model.getNodes()

    # Create groups
    param_grp = hdf5.create_group("parameters")

    # Iterate over nodes
    for node in nodes:

        # Collect node parameters
        parameters = nodes[node].getParameters()

        # Create node subgroup
        node_subgrp = param_grp.create_group(node)

        # Multi-view nodes
        if type(parameters) == list:
            # Loop through the views
            for m in range(len(parameters)):
                if view_names is not None:
                    tmp = view_names[m]
                else:
                    tmp = "%d" % m
                # Create subsubgroup for the view
                view_subgrp = node_subgrp.create_group(tmp)
                # Loop through the parameters of the view
                if parameters[m] is not None:
                    # Variational nodes
                    if type(parameters[m]) == dict:
                        for param_name in parameters[m].keys():
                            if parameters[m][param_name] is not None:
                                view_subgrp.create_dataset(param_name, data=parameters[m][param_name].T)
                    # Non-variational nodes (no distributions)
                    elif type(parameters[m]) == np.ndarray:
                           view_subgrp.create_dataset("value", data=parameters[m].T)

        # Single-view nodes
        else:
            for param_name in parameters.keys():
                node_subgrp.create_dataset("%s" % (param_name), data=parameters[param_name].T)
    pass

def saveExpectations(model, hdf5, view_names=None):
    """ Method to save the expectations of the model in an hdf5 file
    
    PARAMETERS
    ----------
    model: a BayesNet instance
    hdf5: 
    view_names:
    """
    # Get nodes from the model
    nodes = model.getNodes()

    exp_grp = hdf5.create_group("expectations")

    # Iterate over nodes
    for node in nodes:

        # Collect node expectations
        expectations = nodes[node].getExpectations()


        # Multi-view nodes
        if type(expectations) == list:

            # Create subgroup for the node
            node_subgrp = exp_grp.create_group(node)

            # Iterate over views
            for m in range(len(expectations)):
                if view_names is not None:
                    view = view_names[m]
                else:
                    view = "%d" % m

                # Collect expectations
                exp = expectations[m]["E"]
                if exp  is not None:
                    if type(exp) == ma.core.MaskedArray:
                        tmp = ma.filled(exp, fill_value=np.nan)
                        node_subgrp.create_dataset(view, data=tmp.T)
                    else:
                        node_subgrp.create_dataset(view, data=exp.T)

        # Single-view nodes
        else:
            exp_grp.create_dataset(node, data=expectations["E"].T)

def saveTrainingStats(model, hdf5):
    """ Method to save the training statistics in an hdf5 file
    
    PARAMETERS
    ----------
    model: a BayesNet instance
    hdf5: 
    """
    stats = model.getTrainingStats()
    stats_grp = hdf5.create_group("training_stats")
    stats_grp.create_dataset("activeK", data=stats["activeK"])
    stats_grp.create_dataset("elbo", data=stats["elbo"])
    stats_grp.create_dataset("elbo_terms", data=stats["elbo_terms"].T)
    stats_grp['elbo_terms'].attrs['colnames'] = [a.encode('utf8') for a in stats["elbo_terms"].columns.values]

def saveTrainingOpts(opts, hdf5):
    """ Method to save the training options in an hdf5 file
    
    PARAMETERS
    ----------
    opts:
    hdf5: 
    """
    # Remove dictionaries from the options
    for k,v in opts.copy().items():
        if type(v)==dict:
            for k1,v1 in v.items():
                opts[str(k)+"_"+str(k1)] = v1
            opts.pop(k)

    if 'schedule' in opts.keys():
        del opts['schedule']

    # Create HDF5 data set
    hdf5.create_dataset("training_opts", data=np.array(list(opts.values()), dtype=np.float))
    hdf5['training_opts'].attrs['names'] = np.asarray(list(opts.keys())).astype('S')

def saveModelOpts(opts, hdf5):
    """ Method to save the model options in an hdf5 file
    
    PARAMETERS
    ----------
    opts: model options
    hdf5: h5py.File instance
    """
    opts_interest = ["learnIntercept","likelihoods","sparsity_bool","factors"]
    opts = dict((k, opts[k]) for k in opts_interest)
    opts["sparsity"] = opts.pop("sparsity_bool")
    opts["numFactors"] = opts.pop("factors")
    grp = hdf5.create_group('model_opts')
    for k,v in opts.items():
        grp.create_dataset(k, data=np.asarray(v).astype('S'))
    grp[k].attrs['names'] = np.asarray(list(opts.keys())).astype('S')

def saveTrainingData(model, hdf5, data, view_names=None, sample_names=None, feature_names=None, likelihoods=None):
    """ Method to save the training data in an hdf5 file
    
    PARAMETERS
    ----------
    model: a BayesNet instance
    hdf5: h5py.File instance
    view_names: list with view names as characters. If None, no view names will be saved
    sample_names: list with sample names as characters. If None, no sample names will be saved
    feature_names: list with feature names as characters. If None, no feature names will be saved
    """
    # data = model.getTrainingData()

    # Make sure that missing values are masked
    # for i in range(len(data)):
    #     data[i].data[data[i].mask] = np.nan

    data_grp = hdf5.create_group("data")
    samples_grp = hdf5.create_group("samples")
    features_grp = hdf5.create_group("features")
    intercept_grp = hdf5.create_group("intercept")

    if sample_names is not None:
        samples_grp.create_dataset("samples", data=np.array(sample_names, dtype='S50'))

    # if likelihoods is not None:
    #     data_grp.attrs['likelihood'] = np.array(likelihoods, dtype='S50')

    for m in range(len(data)):
        view = view_names[m] if view_names is not None else str(m)
        data_grp.create_dataset(view, data=data[m].T)
        intercept_grp.create_dataset(view, data=np.nanmean(data[m],axis=0))
        # if likelihoods[m] is "gaussian":
        # else:
            # intercept_grp.create_dataset(view, data=model.getNodes()["Y"].nodes[m].means)
        if feature_names is not None:
            features_grp.create_dataset(view, data=np.array(feature_names[m], dtype='S50'))
        

def saveModel(model, data, outfile, train_opts, model_opts, view_names=None, sample_names=None, feature_names=None):
    """ Method to save the model in an hdf5 file
    
    PARAMETERS
    ----------
    TO-FILL....
    """

    # QC checks
    assert model.trained == True, "Model is not trained yet"
    if view_names is not None:
        assert len(np.unique(view_names)) == len(view_names), 'View names must be unique'
    if sample_names is not None:
        assert len(np.unique(sample_names)) == len(sample_names), 'Sample names must be unique'
    if feature_names is not None:
        for x in feature_names: assert len(np.unique(x)) == len(x), 'Feature names must be unique'

    # Create output directory
    if not os.path.isdir(os.path.dirname(outfile)):
        print("Output directory does not exist, creating it...")
        os.makedirs(os.path.dirname(outfile))

    # For some reason h5py orders the datasets alphabetically, so we have to sort the likelihoods accordingly
    idx = sorted(range(len(view_names)), key=lambda k: view_names[k])
    tmp = [model_opts["likelihoods"][idx[m]] for m in range(len(model_opts["likelihoods"]))]
    model_opts["likelihoods"] = tmp

    # Open HDF5 handler
    hdf5 = h5py.File(outfile,'w')

    # Save expectations
    saveExpectations(model,hdf5,view_names)

    # Save parameters
    # saveParameters(model,hdf5,view_names)

    # Save training statistics
    saveTrainingStats(model,hdf5)

    # Save training options
    saveTrainingOpts(train_opts,hdf5)

    # Save model options
    saveModelOpts(model_opts,hdf5)

    # Save training data
    saveTrainingData(model, hdf5, data, view_names, sample_names, feature_names, model_opts["likelihoods"])

    # Close HDF5 file
    hdf5.close()
