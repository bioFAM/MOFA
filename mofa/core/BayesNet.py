
"""
This module is used to define the class containing the entire Bayesian Network,
and the corresponding attributes/methods to train the model, set algorithmic options, calculate lower bound, etc.

To-do:
- default values for schedule and options
"""

from __future__ import division
from time import time
import os
import scipy as s
import pandas as pd
import sys

from .variational_nodes import Variational_Node
from .utils import corr, nans



class BayesNet(object):
    def __init__(self, dim, nodes, options):
        """ Initialisation of a Bayesian network

        PARAMETERS
        ----------
        dim: dict
            keyworded dimensionalities, ex. {'N'=10, 'M'=3, ...}
        nodes: dict
            dictionary with all nodes where the keys are the name of the node and the values are instances the 'Node' class
        options: dict
            training options, such as schedule of updates, maximum number of iterations, training options, etc.
        """

        self.dim = dim
        self.nodes = nodes
        self.options = options

        # Training flag
        self.trained = False

    def removeInactiveFactors(self, by_norm=0, by_r2=0):
        """Method to remove inactive factors

        PARAMETERS
        ----------
        by_norm: float
            threshold to shut down factors based on the norm of the latent variable
            CURRENTLY NOT IMPLEMENTED
        by_r2: float
            threshold to shut down factors based on the coefficient of determination
        """
        drop_dic = {}

        # Shut down based on norm of latent variable vectors
        if by_norm > 0:
            Z = self.nodes["Z"].getExpectation()
            drop_dic["by_norm"] = s.where((Z**2).mean(axis=0) < by_norm)[0]
            if len(drop_dic["by_norm"]) > 0:
                print("\n...A Factor has a norm smaller than {0}, dropping it and recomputing ELBO...\n".format(by_norm))
                drop_dic["by_norm"] = [ s.random.choice(drop_dic["by_norm"]) ]

        # Shut down based on coefficient of determination with respect to the residual variance
        #   Advantages: it takes into account both weights and latent variables, is based on how well the model fits the data
        #   Disadvantages: slow, doesnt work with non-gaussian data
        if by_r2>0:
            Z = self.nodes['Z'].getExpectation()
            Y = self.nodes["Y"].getExpectation()
            W = self.nodes["SW"].getExpectation()
            all_r2 = s.zeros([self.dim['M'], self.dim['K']])
            initial_k = 0
            for m in range(self.dim['M']):

                # Fetch the mask for missing vlaues
                mask = self.nodes["Y"].getNodes()[m].getMask()

                # Fetch the data
                Ym = Y[m].data.copy()

                # If there is an intercept term, regress it out
                if s.all(Z[:,0]==1.):
                    Ym -= W[m][:,0]
                    all_r2[:,0] = 1.
                    initial_k = 1

                # Subtract mean and compute sum of squares (denominator)
                # Ym[mask] = 0.
                # Ym -= s.mean(Ym, axis=0)
                # Ym[mask] = 0.


                # Compute R2
                SS = (Ym**2.).sum()
                for k in range(initial_k,self.dim['K']):
                    Ypred_mk = s.outer(Z[:,k], W[m][:,k])
                    Ypred_mk[mask] = 0.
                    Res = ((Ym - Ypred_mk)**2.).sum()
                    all_r2[m,k] = 1. - Res/SS

            # Select factor to remove. If multiple, then just pick one at random.
            drop_dic["by_r2"] = s.where( (all_r2>by_r2).sum(axis=0) == 0)[0] 
            # drop_dic["by_r2"] = s.where( ((all_r2)>by_r2+1e-16).sum(axis=0) == 0)[0] 
            if len(drop_dic["by_r2"]) > 0: 
                print("\n...A Factor explains less than {0}% of variance, dropping it and recomputing ELBO...\n".format(by_r2*100))
                drop_dic["by_r2"] = [ s.random.choice(drop_dic["by_r2"]) ] # drop one factor at a time

        # Drop the factors
        drop = s.unique(s.concatenate(list(drop_dic.values())))
        if len(drop) > 0:
            for node in self.nodes.keys():
                self.nodes[node].removeFactors(drop)
            self.dim['K'] -= len(drop)

        # TO-DO: THIS ALSO COUNTS COVARIATES
        if self.dim['K']==0:
            print("Shut down all components, no structure found in the data.")
            sys.stdout.flush()
            exit()

        pass

    def iterate(self):
        """Method to start iterating and updating the variables using the VB algorithm"""

        # Define some variables to monitor training
        nodes = list(self.getVariationalNodes().keys())
        elbo = pd.DataFrame(data = nans((self.options['maxiter'], len(nodes)+1 )), columns = nodes+["total"] )
        activeK = nans((self.options['maxiter']))
        
        # Precompute updates
        for n in self.nodes:
            self.nodes[n].precompute()

        # Start training
        for i in range(self.options['maxiter']):
            t = time();

            # Remove inactive latent variables
            if (i>=self.options["startdrop"]) and (i%self.options['freqdrop'])==0 and i<=self.options["enddrop"]:
                self.removeInactiveFactors(**self.options['drop'])
                activeK[i] = self.dim["K"]

            # Update node by node, with E and M step merged
            for node in self.options["schedule"]:
                if node=="Theta" and i<self.options['startSparsity']:
                    continue
                self.nodes[node].update()

            # if i==self.options['startSparsity'] and "Theta" in self.options["schedule"]:
            #     print("\n...Activating sparsity, recomputing ELBO...\n")

            # Calculate Evidence Lower Bound
            if (i+1) % self.options['elbofreq'] == 0 and activeK[i]==activeK[i-1] and i!=self.options['startSparsity']:
                elbo.iloc[i] = self.calculateELBO()

                # Print first iteration
                if i==0:
                    print("Iteration 1: time=%.2f ELBO=%.2f, Factors=%d" % (time()-t,elbo.iloc[i]["total"], (~self.nodes["Z"].covariates).sum() ))
                    if self.options['verbose']:
                        print("".join([ "%s=%.2f  " % (k,v) for k,v in elbo.iloc[i].drop("total").iteritems() ]) + "\n")

                else:
                    # Check convergence using the ELBO
                    delta_elbo = elbo.iloc[i]["total"]-elbo.iloc[i-self.options['elbofreq']]["total"]

                    # Print ELBO monitoring
                    if not s.isnan(delta_elbo):
                        print("Iteration %d: time=%.2f ELBO=%.2f, deltaELBO=%.4f, Factors=%d" % (i+1, time()-t, elbo.iloc[i]["total"], delta_elbo, (~self.nodes["Z"].covariates).sum() ))
                    if self.options['verbose']:
                        print("".join([ "%s=%.2f  " % (k,v) for k,v in elbo.iloc[i].drop("total").iteritems() ]) + "\n")
                    if delta_elbo<0 and i!=self.options['startSparsity'] and self.options['verbose']: print("Warning, lower bound is decreasing..."); print('\a')

                    # Assess convergence
                    if (0 <= delta_elbo < self.options['tolerance']) and (not self.options['forceiter']):
                        activeK = activeK[:(i+1)]
                        elbo = elbo[:(i+1)]
                        print ("Converged!\n")
                        break

            # Do not calculate lower bound
            else:
                print("Iteration %d: time=%.2f, Factors=%d" % (i+1,time()-t,(~self.nodes["Z"].covariates).sum()))

            # Flush (we need this to print when running on the cluster)
            sys.stdout.flush()

        # Finish by collecting the training statistics
        self.train_stats = { 'activeK':activeK, 'elbo':elbo["total"].values, 'elbo_terms':elbo.drop("total",1) }
        self.trained = True

    def getParameters(self, *nodes):
        """Method to collect all parameters of a given set of nodes (all by default)

        PARAMETERS
        ----------
        nodes: iterable
            name of the nodes
        """

        if len(nodes) == 0: nodes = self.nodes.keys()
        params = {}
        for node in nodes:
            tmp = self.nodes[node].getParameters()
            if tmp != None: params[node] = tmp
        return params

    def getExpectations(self, only_first_moments=False, *nodes):
        """Method to collect all expectations of a given set of nodes (all by default)
        
        PARAMETERS
        ----------
        only_first_moments: bool
            get only first moments?
        nodes: list
            name of the nodes
        """

        if len(nodes) == 0: nodes = self.nodes.keys()
        expectations = {}
        for node in nodes:
            if only_first_moments:
                tmp = self.nodes[node].getExpectation()
            else:
                tmp = self.nodes[node].getExpectations()
            expectations[node] = tmp
        return expectations

    def getNodes(self):
        """ Method to return all nodes """
        return self.nodes

    def getVariationalNodes(self):
        """ Method to return all variational nodes """
        return { k:v for k,v in self.nodes.items() if isinstance(v,Variational_Node) }

    def getTrainingStats(self):
        """ Method to return training statistics """
        return self.train_stats

    def getTrainingOpts(self):
        """ Method to return training options """
        return self.options

    def getTrainingData(self):
        """ Method to return training data """
        return self.nodes["Y"].getValues()

    def calculateELBO(self, *nodes):
        """Method to calculate the Evidence Lower Bound of the model

        PARAMETERS
        ----------
        nodes: iterable
            list/tuple with the name of the nodes
        """
        if len(nodes) == 0: nodes = self.getVariationalNodes().keys()
        elbo = pd.Series(s.zeros(len(nodes)+1), index=list(nodes)+["total"])
        for node in nodes:
            elbo[node] = float(self.nodes[node].calculateELBO())
            elbo["total"] += elbo[node]
        return elbo
