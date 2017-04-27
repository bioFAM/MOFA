
"""
Module to initalise the nodes
"""

import scipy as s
import scipy.stats as stats
from sys import path
import sklearn.decomposition

from scGFA.core.nodes import *
from scGFA.core.multiview_nodes import *
from scGFA.core.seeger_nodes import *
from scGFA.core.sparse_updates import *


# General class to initialise a MOFA model
class initModel(object):
    def __init__(self, dim, data, lik, seed=None):
        # Inputs:
        #  dim (dic): keyworded dimensionalities
        #    N for the number of samples, M for the number of views
        #    K for the number of latent variables, D for the number of features (per view, so it is a list)
        #  data (list of ndarrays of length M): observed data
        #  lik (list): likelihood type for each view
        self.data = data
        self.lik = lik
        self.N = dim["N"]
        self.K = dim["K"]
        self.M = dim["M"]
        self.D = dim["D"]

        self.nodes = {}

        s.random.seed(seed)

    def getNodes(self):
        return { k:v for (k,v) in self.nodes.iteritems()}


# Class to iniailise the (sparse) MOFAmodel
class init_sparse(initModel):
    def __init__(self, *args, **kwargs):
        super(init_sparse, self).__init__(*args, **kwargs)

    # def __init__(self, dim, data, lik, seed=None):
        # initModel.__init__(self, dim, data, lik, seed)

    def initZ(self, pmean, pvar, qmean, qvar, qE=None, qE2=None, covariates=None):
        # Method to initialise the latent variables
        # covariates (nd array): matrix of covariates with dimensions (nsamples,ncovariates)

        # Initialise first moment
        if qmean is not None:
            if isinstance(qmean,str):
                if qmean == "random": # Random initialisation of latent variables
                    qmean = stats.norm.rvs(loc=0, scale=1, size=(self.N,self.K))

                elif qmean == "orthogonal": # Latent variables are initialised randomly but ensuring orthogonality
                    pca = sklearn.decomposition.PCA(n_components=self.K, copy=True, whiten=True)
                    pca.fit(stats.norm.rvs(loc=0, scale=1, size=(self.N,9999)).T)
                    qmean = pca.components_.T

                elif qmean == "pca": # Latent variables are initialised from PCA in the concatenated matrix
                    pca = sklearn.decomposition.PCA(n_components=self.K, copy=True, whiten=True)
                    pca.fit(s.concatenate(self.data,axis=0).T)
                    qmean = pca.components_.T

            elif isinstance(qmean,s.ndarray):
                assert qmean.shape == (self.N,self.K)

            elif isinstance(qmean,(int,float)):
                qmean = s.ones((self.N,self.K)) * qmean

            else:
                print "Wrong initialisation for Z"
                exit()

        # Add covariates
        if covariates is not None:
            idx_covariates = s.array(range(covariates.shape[1]))
            qmean[:,idx_covariates] = covariates
        else:
            idx_covariates = None

        # Initialise the node
        # self.Z = Constant_Node(dim=(self.N,self.K), value=qmean)
        self.Z = Z_Node(dim=(self.N,self.K),
                        pmean=s.ones((self.N,self.K))*pmean,
                        pvar=s.ones((self.N,self.K))*pvar,
                        qmean=s.ones((self.N,self.K))*qmean,
                        qvar=s.ones((self.N,self.K))*qvar,
                        qE=qE, qE2=qE2,
                        idx_covariates=idx_covariates)
        self.nodes["Z"] = self.Z

    def initSW(self, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0, qEW_S1, qES):
        # Method to initialise the spike-slab variable (product of bernoulli and gaussian variables)
        SW_list = [None]*self.M
        for m in xrange(self.M):

           # Initialise first moment
            if isinstance(qmean_S1[m],str):
                if qmean_S1[m] == "random":
                    qmean_S1[m] = stats.norm.rvs(loc=0, scale=1, size=(self.D[m],self.K))
                else:
                    print "%s initialisation not implemented for SW" % qmean_S1[m]
                    exit()
            elif isinstance(qmean_S1[m],s.ndarray):
                assert qmean_S1[m].shape == (self.D[m],self.K), "Wrong dimensionality"
            elif isinstance(qmean_S1[m],(int,float)):
                qmean_S1[m] = s.ones((self.D[m],self.K)) * qmean_S1[m]
            else:
                print "Wrong initialisation for SW"
                exit()
            SW_list[m] = SW_Node(
                dim=(self.D[m],self.K),

                ptheta=ptheta[m],
                pmean_S0=pmean_S0[m],
                pvar_S0=pvar_S0[m],
                pmean_S1=pmean_S1[m],
                pvar_S1=pvar_S1[m],

                qtheta=qtheta[m],
                qmean_S0=qmean_S0[m],
                qvar_S0=qvar_S0[m],
                qmean_S1=qmean_S1[m],
                qvar_S1=qvar_S1[m],

                qES=qES[m],
                qEW_S0=qEW_S0[m],
                qEW_S1=qEW_S1[m],
                )

        self.SW = Multiview_Variational_Node(self.M, *SW_list)
        self.nodes["SW"] = self.SW

    def initAlphaW_extended(self, pa, pb, qa, qb, qE):
        # Method to initialise the precision of the group-wise ARD prior
        # Inputs:
        #  pa (float): 'a' parameter of the prior distribution
        #  pb (float): 'b' parameter of the prior distribution
        #  qb (float): initialisation of the 'b' parameter of the variational distribution
        #  qE (float): initial expectation of the variational distribution
        alpha_list = [None]*self.M
        for m in xrange(self.M):
            alpha_list[m] = AlphaW_Node(dim=(self.K,), pa=pa[m], pb=pb[m], qa=qa[m], qb=qb[m], qE=qE[m])
            # alpha_list[m] = Constant_Node(dim=(self.K,), value=qE[m])
            # alpha_list[m].factors_axis = 0
        self.AlphaW = Multiview_Variational_Node(self.M, *alpha_list)
        # self.AlphaW = Multiview_Constant_Node(self.M, *alpha_list)
        self.nodes["AlphaW"] = self.AlphaW

    def initAlphaW_basic(self, pa, pb, qa, qb, qE):
        # Method to initialise the precision of the group-wise ARD prior
        # Inputs:
        #  pa (float): 'a' parameter of the prior distribution
        #  pb (float): 'b' parameter of the prior distribution
        #  qb (float): initialisation of the 'b' parameter of the variational distribution
        #  qE (float): initial expectation of the variational distribution
        alpha_list = [None]*self.M
        for m in xrange(self.M):
            alpha_list[m] = AlphaW_Node(dim=(1,), pa=pa[m], pb=pb[m], qa=qa[m], qb=qb[m], qE=qE[m])
            # alpha_list[m] = Constant_Node(dim=(1,), value=qE[m])
        self.AlphaW = Multiview_Variational_Node(self.M, *alpha_list)
        # self.AlphaW = Multiview_Constant_Node(self.M, *alpha_list)
        self.nodes["AlphaW"] = self.AlphaW

    def initAlphaZ(self, pa, pb, qa, qb, qE, covariates=None):
        # Method to initialise the precision of the latent variable ARD prior
        # Inputs:
        #  pa (float): 'a' parameter of the prior distribution
        #  pb (float): 'b' parameter of the prior distribution
        #  qb (float): initialisation of the 'b' parameter of the variational distribution
        #  qE (float): initial expectation of the variational distribution

        self.AlphaZ = AlphaZ_Node(dim=(self.K,), pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        # self.AlphaZ = Constant_Node(dim=(self.K,), value=qE)
        # self.AlphaZ.factors_axis = 0
        self.nodes["AlphaZ"] = self.AlphaZ

    def initTau(self, pa, pb, qa, qb, qE):
        # Method to initialise the precision of the noise
        # Inputs:
        #  pa (float): 'a' parameter of the prior distribution
        #  pb (float): 'b' parameter of the prior distribution
        #  qb (float): initialisation of the 'b' parameter of the variational distribution
        #  qE (float): initial expectation of the variational distribution
        tau_list = [None]*self.M
        for m in xrange(self.M):
            if self.lik[m] == "poisson":
                tmp = 0.25 + 0.17*s.amax(self.data[m],axis=0)
                tau_list[m] = Constant_Node(dim=(self.D[m],), value=tmp)
            elif self.lik[m] == "bernoulli":
                # tmp = s.ones(self.D[m])*0.25
                # tau_list[m] = Constant_Node(dim=(self.D[m],), value=tmp)
                tau_list[m] = Tau_Jaakkola(dim=(self.D[m],), value=1)
            elif self.lik[m] == "binomial":
                tmp = 0.25*s.amax(self.data["tot"][m],axis=0)
                tau_list[m] = Constant_Node(dim=(self.D[m],), value=tmp)
            elif self.lik[m] == "gaussian":
                tau_list[m] = Tau_Node(dim=(self.D[m],), pa=pa[m], pb=pb[m], qa=qa[m], qb=qb[m], qE=qE[m])
        self.Tau = Multiview_Mixed_Node(self.M,*tau_list)
        self.nodes["Tau"] = self.Tau

    def initY(self):
        # Method to initialise the observed data
        Y_list = [None]*self.M
        for m in xrange(self.M):
            if self.lik[m]=="gaussian":
                Y_list[m] = Y_Node(dim=(self.N,self.D[m]), value=self.data[m])
            elif self.lik[m]=="poisson":
                Y_list[m] = Poisson_PseudoY(dim=(self.N,self.D[m]), obs=self.data[m], E=None)
            elif self.lik[m]=="bernoulli":
                # Y_list[m] = Bernoulli_PseudoY(dim=(self.N,self.D[m]), obs=self.data[m], E=None)
                tmp = stats.norm.rvs(loc=0, scale=1, size=(self.N,self.D[m]))
                Y_list[m] =  Bernoulli_PseudoY_Jaakkola(dim=(self.N,self.D[m]), obs=self.data[m], E=tmp)
            elif self.lik[m]=="binomial":
                Y_list[m] = Binomial_PseudoY(dim=(self.N,self.D[m]), tot=data["tot"][m], obs=data["obs"][m], E=None)
        self.Y = Multiview_Mixed_Node(self.M, *Y_list)
        self.nodes["Y"] = self.Y

    def initTheta(self, pa, pb, qa, qb, qE, learnTheta):
        # Method to initialie a general theta node
        # Inputs:
        #  pa (float): 'a' parameter of the prior distribution
        #  pb (float): 'b' parameter of the prior distribution
        #  qb (float): initialisation of the 'b' parameter of the variational distribution
        #  qE (float): initial expectation of the variational distribution
        #  
        Theta_list = [None] * self.M
        for m in xrange(self.M):
            Klearn = learnTheta[m,:]==1.
            Kconst = learnTheta[m,:]==0.
            Theta_list[m] = Mixed_Theta_Nodes(
                LearnTheta=Theta_Node(dim=(s.sum(Klearn),), pa=pa[m][Klearn], pb=pb[m][Klearn], qa=qa[m][Klearn], qb=qb[m][Klearn], qE=qE[m][Klearn]),
                ConstTheta=Theta_Constant_Node(dim=((self.D[m],s.sum(Kconst),)), value=s.repeat(qE[m][Kconst][None,:], self.D[m], 0), N_cells=1.),
                idx=learnTheta[m,:]
            )
        # WHAT IS THIS???
        # self.Theta = Multiview_Mixed_Theta_Nodes(self.M, *Theta_list)
        # self.Theta = Multiview_Variational_Node(self.M, *Theta_list)
        self.Theta = Multiview_Mixed_Node(self.M, *Theta_list)
        self.nodes["Theta"] = self.Theta

    def initThetaLearn(self, pa, pb, qa, qb, qE):
        # Method to initialise the theta node
        Theta_list = [None] * self.M
        for m in xrange(self.M):
            Theta_list[m] = Theta_Node(dim=(self.K,), pa=pa[m], pb=pb[m], qa=qa[m], qb=qb[m], qE=qE[m])
        self.Theta = Multiview_Variational_Node(self.M, *Theta_list)
        self.nodes["Theta"] = self.Theta

    def initThetaConst(self, value):
        # Method to initialise the theta node
        Theta_list = [None] * self.M
        for m in xrange(self.M):
            Theta_list[m] = Theta_Constant_Node(dim=(self.K,), value=value[m], N_cells=1.)
        self.Theta = Multiview_Constant_Node(self.M, *Theta_list)
        self.nodes["Theta"] = self.Theta

    def initClusters(self, clusters=None, pmean=0, pvar=1, qmean=0, qvar=1, qE=None):
        if clusters is None:
            clusters = s.zeros(self.N, int)
        self.Clusters = Cluster_Node_Gaussian(pmean, pvar, qmean, qvar, clusters, self.K)
        # self.Clusters = Constant_Node(pmean, pvar, qmean, qvar, clusters, self.K)
        self.nodes['Clusters'] = self.Clusters

    def initExpectations(self, *nodes):
        # Method to initialise the expectations of some nodes
        for node in nodes:
            self.nodes[node].updateExpectations()

    def MarkovBlanket(self):
        # Method to define the markov blanket

        if 'AlphaZ' in self.nodes:
            if 'Cluster' in self.nodes:
                self.Z.addMarkovBlanket(SW=self.SW, Tau=self.Tau, Y=self.Y, Alpha=self.AlphaZ, Cluster=self.Clusters)
                self.AlphaZ.addMarkovBlanket(Z=self.Z, Cluster=self.Clusters)
                self.Clusters.addMarkovBlanket(Z=self.Z, Alpha=self.AlphaZ)
            else:
                self.Z.addMarkovBlanket(SW=self.SW, Tau=self.Tau, Y=self.Y, Alpha=self.AlphaZ)
                self.AlphaZ.addMarkovBlanket(Z=self.Z)
        else:
            if 'Cluster' in self.nodes:
                self.Clusters.addMarkovBlanket(Z=self.Z)
                self.Z.addMarkovBlanket(SW=self.SW, Tau=self.Tau, Y=self.Y, Cluster=self.Clusters)
            else:
                self.Z.addMarkovBlanket(SW=self.SW, Tau=self.Tau, Y=self.Y)

        for m in xrange(self.M):
            self.Theta.nodes[m].addMarkovBlanket(SW=self.SW.nodes[m])
            self.AlphaW.nodes[m].addMarkovBlanket(SW=self.SW.nodes[m])
            self.SW.nodes[m].addMarkovBlanket(Z=self.Z, Tau=self.Tau.nodes[m], Alpha=self.AlphaW.nodes[m], Y=self.Y.nodes[m], Theta=self.Theta.nodes[m])
            if self.lik[m]=="gaussian":
                self.Y.nodes[m].addMarkovBlanket(Z=self.Z, SW=self.SW.nodes[m], Tau=self.Tau.nodes[m])
                self.Tau.nodes[m].addMarkovBlanket(SW=self.SW.nodes[m], Z=self.Z, Y=self.Y.nodes[m])
            elif self.lik[m]=="bernoulli":
                self.Tau.nodes[m].addMarkovBlanket(Z=self.Z, SW=self.SW.nodes[m], Y=self.Y.nodes[m]) # Jaakkola
                self.Y.nodes[m].addMarkovBlanket(Z=self.Z, SW=self.SW.nodes[m])
            else:
                self.Y.nodes[m].addMarkovBlanket(Z=self.Z, SW=self.SW.nodes[m], kappa=self.Tau.nodes[m])
