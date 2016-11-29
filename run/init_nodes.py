
"""
Module to initalise the nodes
"""

import scipy as s
import scipy.stats as stats
from sys import path
import sklearn.decomposition

path.insert(0,"../")
from multiview_nodes import *
from local_nodes import *
from seeger_nodes import *

from sparse_updates import *

# General class to initialise a Group Factor Analysis model
class initModel(object):
    def __init__(self, dim, data, lik):
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


# Class to iniailise the (sparse) Group Factor Analysis model
class init_scGFA(initModel):
    def __init__(self, dim, data, lik):
        initModel.__init__(self, dim, data, lik)

    def initZ(self, pmean, pvar, type="random"):
        # Method to initialise the latent variables
        # Inputs:
        #  type (str): random, orthogonal, pca

        # Initialise the mean of the Q distribution
        if type == "random":
            qmean = stats.norm.rvs(loc=0, scale=1, size=(self.N,self.K))
        elif type == "orthogonal":
            print "Not implemented"
            exit()
        elif type == "pca":
            pca = sklearn.decomposition.PCA(n_components=self.K, copy=True, whiten=True)
            tmp = s.concatenate(self.data,axis=0).T
            pca.fit(tmp)
            qmean = pca.components_.T

        # Initialise the variance of the Q distribution
        qvar = s.ones((self.N,self.K))

        # Initialise the node
        self.Z = Z_Node(dim=(self.N,self.K), pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qmean)
        # self.Z.updateExpectations()


        def __init__(self, dim, qmean, qvar, ptheta, qtheta):

    def initSW(self, ptheta, pmean, pvar, qtheta, qmean, qvar):
        # Method to initialise the spike-slab variable (product of bernoulli and gaussian variables)
        # TO-DO: RIGHT NOW PMEAN AND PVAR ARE NOT USED
        SW_list = [None]*self.M
        for m in xrange(self.M):
            qtheta = s.ones((self.D[m],self.K))*qtheta
            qmean = stats.norm.rvs(loc=qmean, scale=qvar, size=(self.D[m],self.K))
            qvar = s.ones((self.D[m],self.K))*qvar
            SW_list[m] = SW_Node(dim=(self.D[m],self.K), ptheta=ptheta, qtheta=qtheta, qmean=qmean, qvar=qvar)
        self.SW = Multiview_Variational_Node(self.M, *SW_list)

    def initAlpha(self, pa, pb, qa, qb, qE):
        # Method to initialise the precision of the group-wise ARD prior
        # Inputs:
        #  pa (float): 'a' parameter of the prior distribution
        #  pb (float): 'b' parameter of the prior distribution
        #  qb (float): initialisation of the 'b' parameter of the variational distribution
        #  qE (float): initial expectation of the variational distribution
        alpha_list = [None]*self.M
        for m in xrange(self.M):
            alpha_list[m] = Alpha_Node(dim=(self.K,), pa=pa, pb=pb, qa=qa, qb=s.ones(self.K)*qb, qE=s.ones(self.K)*qE)
        self.Alpha = Multiview_Variational_Node((self.K,)*self.M, *alpha_list)


    def initTau(self, pa, pb, qb, qE):
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
                tau_list[m] = Observed_Local_Node(dim=(self.D[m],), value=tmp)
            elif self.lik[m] == "bernoulli":
                tmp = s.ones(self.D[m])*0.25 
                tau_list[m] = Observed_Local_Node(dim=(self.D[m],), value=tmp)
            elif self.lik[m] == "binomial":
                tmp = 0.25*s.amax(self.data["tot"][m],axis=0)
                tau_list[m] = Observed_Local_Node(dim=(self.D[m],), value=tmp)
            elif self.lik[m] == "gaussian":
                qa = pa + self.N/2
                tau_list[m] = Tau_Node(dim=(self.D[m],), pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.Tau = Multiview_Mixed_Node(self.M,*tau_list)


    def initY(self):
        # Method to initialise the observed data
        Y_list = [None]*self.M
        for m in xrange(self.M):
            if self.lik[m]=="gaussian":
                Y_list[m] = Y_Node(dim=(self.N,self.D[m]), obs=self.data[m])
            elif self.lik[m]=="poisson":
                Y_list[m] = Poisson_PseudoY_Node(dim=(self.N,self.D[m]), obs=self.data[m], E=None)
            elif self.lik[m]=="bernoulli":
                Y_list[m] = Bernoulli_PseudoY_Node(dim=(self.N,self.D[m]), obs=self.data[m], E=None)
            elif self.lik[m]=="binomial":
                Y_list[m] = Binomial_PseudoY_Node(dim=(self.N,self.D[m]), tot=data["tot"][m], obs=data["obs"][m], E=None)
        self.Y = Multiview_Mixed_Node(self.M, *Y_list)


    def initThetaLearn(self, ptheta, pa, pb, qa., qb):
        # Method to initialise the theta node
        # TO-DO: ADD ANNOTATIONS
        Theta_list = [None] * M
        learn_theta = True
        for m in xrange(M):
            if learn_theta:
                Theta_list[m] = Theta_Node_No_Annotation((self.K,))
            else:
        Theta = Multiview_Mixed_Node(M, *Theta_list)

                Theta_list[m] = Constant_Node((D[m],K),0.5)



    def MarkovBlanket(self):
        # Method to define the markov blanket
        self.Z.addMarkovBlanket(SW=self.SW, tau=self.Tau, Y=self.Y)
        for m in xrange(self.M):
            self.Alpha.nodes[m].addMarkovBlanket(SW=self.SW.nodes[m])
            self.SW.nodes[m].addMarkovBlanket(Z=self.Z, tau=self.Tau.nodes[m], alpha=self.Alpha.nodes[m], Y=self.Y.nodes[m])
            if self.lik[m] is "gaussian":
                self.Y.nodes[m].addMarkovBlanket(Z=self.Z, SW=self.SW.nodes[m], tau=self.Tau.nodes[m])
                self.Tau.nodes[m].addMarkovBlanket(SW=self.SW.nodes[m], Z=self.Z, Y=self.Y.nodes[m])
            else:
                self.Zeta.nodes[m].addMarkovBlanket(Z=self.Z, W=self.SW.nodes[m])
                self.Y.nodes[m].addMarkovBlanket(Z=self.Z, W=self.SW.nodes[m], kappa=self.Tau.nodes[m], zeta=self.Zeta.nodes[m])
        # Update expectations of SW (we need to do it here because it requires the markov blanket to be defined )
        self.SW.updateExpectations()
