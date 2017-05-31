from __future__ import division
import numpy.ma as ma
import numpy as np

from time import time

# Import manually defined functions
from variational_nodes import *
from utils import *
from nodes import Constant_Node
from mixed_nodes import Mixed_Theta_Nodes

import scipy.special as special

"""
Module to define the nodes and the corresponding updates of the model

Current nodes:
    Y_Node: observed data
    SW_Node: spike and slab weights
    Tau_Node: precision of the zero-mean normally distributed noise
    AlphaW_Node_basic: ARD precision on the weights, per view
    AlphaW_Node_extended: ARD precision on the weights, per factor (and view)
    AlphaZ_Node: ARD precision on the latent variables, per factor
    Z_Node: latent variables
    Theta_Node: learning sparsity parameter of the spike and slab
    Theta_Constant_Node: fixed sparsity parameter of the spike and slab
    MuZ_Node: learning mean for latent variables, allowing cluster-specific

All nodes belong to the 'Variational_Node' class. They share the following
    methods:
    - precompute: precompute some terms to speed up the calculations
    - calculateELBO: calculate evidence lower bound using current estimates of expectations/params
    - getParameters: return current parameters
    - getExpectations: return current expectations
    - updateParameters: update parameters using current estimates of expectations
    - updateExpectations: update expectations using current estimates of parameters

    attributes:
    - markov_blanket: dictionary that defines the set of nodes that are in the markov blanket of the current node
    - Q: an instance of Distribution() which contains the specification of the variational distribution
    - P: an instance of Distribution() which contains the specification of the prior distribution
    - dim: dimensionality of the node
"""

class Y_Node(Constant_Variational_Node):
    def __init__(self, dim, value):
        Constant_Variational_Node.__init__(self, dim, value)

        # Create a boolean mask of the data to hidden missing values
        if type(self.value) != ma.MaskedArray:
            self.mask()

        # Precompute some terms
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        self.N = self.dim[0] - ma.getmask(self.value).sum(axis=0)
        self.D = self.dim[1]
        self.likconst = -0.5*s.sum(self.N)*s.log(2.*s.pi)

    def mask(self):
        # Mask the observations if they have missing values
        self.value = ma.masked_invalid(self.value)

    def calculateELBO(self):
        # Calculate evidence lower bound
        # We use the trick that the update of Tau already contains the Gaussian likelihod.
        # However, it is important that the lower bound is calculated after the update of Tau is performed
        tauQ_param = self.markov_blanket["Tau"].getParameters("Q")
        tauP_param = self.markov_blanket["Tau"].getParameters("P")
        tau_exp = self.markov_blanket["Tau"].getExpectations()
        lik = self.likconst + 0.5*s.sum(self.N*(tau_exp["lnE"])) - s.dot(tau_exp["E"],tauQ_param["b"]-tauP_param["b"])
        return lik

class Tau_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        super(Tau_Node,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        self.D = self.dim[0]
        self.lbconst = s.sum(self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']))

    def updateParameters(self):

        # Collect expectations from other nodes
        Y = self.markov_blanket["Y"].getExpectation()
        tmp = self.markov_blanket["SW"].getExpectations()
        SW,SWW = tmp["E"], tmp["ESWW"]
        Ztmp = self.markov_blanket["Z"].getExpectations()
        Z,ZZ = Ztmp["E"],Ztmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # Calculate terms for the update
        term1 = s.square(Y).sum(axis=0).data
        term2 = 2*(Y*s.dot(Z,SW.T)).sum(axis=0).data
        term3 = ma.array(ZZ.dot(SWW.T), mask=ma.getmask(Y)).sum(axis=0)
        SWZ = ma.array(SW.dot(Z.T), mask=ma.getmask(Y).T)
        term4 = dotd(SWZ, SWZ.T) - ma.array(s.dot(s.square(Z),s.square(SW).T),mask=ma.getmask(Y)).sum(axis=0)
        tmp = term1 - term2 + term3 + term4

        # Perform updates of the Q distribution
        Qa = Pa + (Y.shape[0] - ma.getmask(Y).sum(axis=0))/2.
        Qb = Pb + tmp/2.

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations from current node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.expectations['E'], self.Q.expectations['lnE']

        # Do the calculations
        lb_p = self.lbconst + s.sum((Pa-1.)*QlnE) - s.sum(Pb*QE)
        lb_q = s.sum(Qa*s.log(Qb)) + s.sum((Qa-1.)*QlnE) - s.sum(Qb*QE) - s.sum(special.gammaln(Qa))

        return lb_p - lb_q

class AlphaW_Node_m(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        super(AlphaW_Node_m,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

    def updateParameters(self):

        # Collect expectations from other nodes
        tmp = self.markov_blanket["SW"].getExpectations()
        ES,EWW = tmp["ES"],tmp["EWW"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # Perform updates
        Qa = Pa + 0.5*ES.shape[0]*ES.shape[1]
        Qb = Pb + 0.5*EWW.sum()

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.getExpectations()['E'], self.Q.getExpectations()['lnE']

        # Do the calculations
        lb_p = (Pa*s.log(Pb)).sum() - special.gammaln(Pa).sum() + ((Pa-1.)*QlnE).sum() - (Pb*QE).sum()
        lb_q = (Qa*s.log(Qb)).sum() - special.gammaln(Qa).sum() + ((Qa-1.)*QlnE).sum() - (Qb*QE).sum()

        return lb_p - lb_q

class AlphaW_Node_mk(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        # Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super(AlphaW_Node_mk,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        # self.lbconst = self.K * ( self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )
        # self.lbconst = s.sum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) )
        self.factors_axis = 0

    def updateParameters(self):

        # Collect expectations from other nodes
        tmp = self.markov_blanket["SW"].getExpectations()
        ES,EWW = tmp["ES"],tmp["EWW"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # Perform updates
        Qa = Pa + 0.5*ES.shape[0]
        Qb = Pb + 0.5*EWW.sum(axis=0)

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.getExpectations()['E'], self.Q.getExpectations()['lnE']

        # Do the calculations
        lb_p = (Pa*s.log(Pb)).sum() - special.gammaln(Pa).sum() + ((Pa-1.)*QlnE).sum() - (Pb*QE).sum()
        lb_q = (Qa*s.log(Qb)).sum() - special.gammaln(Qa).sum() + ((Qa-1.)*QlnE).sum() - (Qb*QE).sum()

        return lb_p - lb_q

class AlphaW_Node_k(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        # Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super(AlphaW_Node_k,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        # self.lbconst = self.K * ( self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )
        # self.lbconst = s.sum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) )
        self.factors_axis = 0

    def updateParameters(self):

        # Collect expectations from other nodes
        tmp = self.markov_blanket["SW"].getExpectations()

        # Concatenate
        M = len(tmp)
        ES = s.concatenate([tmp[m]["ES"] for m in xrange(M)],axis=0)
        EWW = s.concatenate([tmp[m]["EWW"] for m in xrange(M)],axis=0)

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # Perform updates
        Qa = Pa + 0.5*ES.shape[0]
        Qb = Pb + 0.5*EWW.sum(axis=0)

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.getExpectations()['E'], self.Q.getExpectations()['lnE']

        # Do the calculations
        lb_p = (Pa*s.log(Pb)).sum() - special.gammaln(Pa).sum() + ((Pa-1.)*QlnE).sum() - (Pb*QE).sum()
        lb_q = (Qa*s.log(Qb)).sum() - special.gammaln(Qa).sum() + ((Qa-1.)*QlnE).sum() - (Qb*QE).sum()

        return lb_p - lb_q

class SW_Node(BernoulliGaussian_Unobserved_Variational_Node):
    # TOO MANY ARGUMENTS, SHOULD WE USE **KWARGS AND *KARGS ONLY?
    # def __init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0=None, qEW_S1=None, qES=None):
    def __init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0=None, qEW_S1=None, qES=None):
        super(SW_Node,self).__init__(dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0, qEW_S1, qES)
        self.precompute()

    def precompute(self):
        self.D = self.dim[0]
        self.factors_axis = 1

    def updateParameters(self):
        # Collect expectations from other nodes
        Ztmp = self.markov_blanket["Z"].getExpectations()
        Z,ZZ = Ztmp["E"],Ztmp["E2"]
        tau = self.markov_blanket["Tau"].getExpectation()
        Y = self.markov_blanket["Y"].getExpectation()
        alpha = self.markov_blanket["Alpha"].getExpectation().copy()
        thetatmp = self.markov_blanket['Theta'].getExpectations()
        theta_lnE, theta_lnEInv  = thetatmp['lnE'], thetatmp['lnEInv']

        # Collect parameters and expectations from P and Q distributions of this node
        SW = self.Q.getExpectations()["E"]
        Q = self.Q.getParameters()
        Qmean_S1, Qvar_S1, Qtheta = Q['mean_S1'], Q['var_S1'], Q['theta']

        # Check dimensions of Theta and and expand if necessary
        if theta_lnE.shape != Qmean_S1.shape:
            theta_lnE = s.repeat(theta_lnE[None,:],Qmean_S1.shape[0],0)
        if theta_lnEInv.shape != Qmean_S1.shape:
            theta_lnEInv = s.repeat(theta_lnEInv[None,:],Qmean_S1.shape[0],0)
        # Check dimensions of Tau and and expand if necessary
        if tau.shape != Y.shape:
            tau = s.repeat(tau[None,:], Y.shape[0], axis=0)
        tau = ma.masked_where(ma.getmask(Y), tau)
        # Check dimensions of Alpha and and expand if necessary
        if alpha.shape[0] == 1:
            alpha = s.repeat(alpha[:], self.dim[1], axis=0)

        # Update each latent variable in turn
        for k in xrange(self.dim[1]):
            # Calculate intermediate steps
            term1 = (theta_lnE-theta_lnEInv)[:,k]
            term2 = 0.5*s.log(alpha[k])
            term3 = 0.5*s.log(ma.dot(ZZ[:,k],tau) + alpha[k])
            term4_tmp1 = ma.dot((tau*Y).T,Z[:,k]).data
            term4_tmp2 = ( tau * s.dot((Z[:,k]*Z[:,s.arange(self.dim[1])!=k].T).T, SW[:,s.arange(self.dim[1])!=k].T) ).sum(axis=0)
            term4_tmp3 = ma.dot(ZZ[:,k].T,tau) + alpha[k]
            term4 = 0.5*s.divide((term4_tmp1-term4_tmp2)**2,term4_tmp3)
            # Update S
            # NOTE there could be some precision issues in S --> loads of 1s in result
            Qtheta[:,k] = 1./(1.+s.exp(-(term1+term2-term3+term4)))
            # Update W
            Qvar_S1[:,k] = s.divide(1,term4_tmp3)
            Qmean_S1[:,k] = Qvar_S1[:,k]*(term4_tmp1-term4_tmp2)

            # Update Expectations for the next iteration
            SW[:,k] = Qtheta[:,k] * Qmean_S1[:,k]

        # Save updated parameters of the Q distribution
        # self.Q.setParameters(mean_S0=s.zeros((self.D,self.dim[1])), var_S0=s.repeat(1/alpha[None,:],self.D,0), mean_S1=Qmean_S1, var_S1=Qvar_S1, theta=Qtheta )
        self.Q.setParameters(mean_S0=s.zeros((self.D,self.dim[1])), var_S0=s.zeros((1,)), mean_S1=Qmean_S1, var_S1=Qvar_S1, theta=Qtheta )

    def calculateELBO(self):

        # Collect parameters and expectations
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
        S,WW = Qexp["ES"], Qexp["EWW"]
        Qvar = Qpar['var_S1']
        theta = self.markov_blanket['Theta'].getExpectations()

        # Get ARD sparsity or prior variance
        if "Alpha" in self.markov_blanket:
            alpha = self.markov_blanket['Alpha'].getExpectations().copy()
            if alpha["E"].shape[0] == 1:
                alpha["E"] = s.repeat(alpha["E"][:], self.dim[1], axis=0)
                alpha["lnE"] = s.repeat(alpha["lnE"][:], self.dim[1], axis=0)
        else:
            print "Not implemented"
            exit()

        # Calculate ELBO for W
        lb_pw = (self.D*alpha["lnE"].sum() - s.sum(alpha["E"]*WW))/2.
        lb_qw = -0.5*self.dim[1]*self.D - 0.5*(S*s.log(Qvar) + (1.-S)*s.log(1./alpha["E"])).sum() # IS THE FIRST CONSTANT TERM CORRECT???
        lb_w = lb_pw - lb_qw

        # Calculate ELBO for S
        lb_ps = S*theta['lnE'] + (1.-S)*theta['lnEInv']
        lb_qs = S*s.log(S) + (1.-S)*s.log(1.-S)
        lb_ps[s.isnan(lb_ps)] = 0.
        lb_qs[s.isnan(lb_qs)] = 0.
        lb_s = s.sum(lb_ps) - s.sum(lb_qs)

        return lb_w + lb_s

class Theta_Node(Beta_Unobserved_Variational_Node):
    """
    This class comtain a Theta node associate to factors for which
    we dont have annotations.

    The inference is done per view and factor, so the dimension of the node is the
    number of non-annotated factors

    the updateParameters function needs to know what factors are non-annotated in
    order to choose from the S matrix
    """

    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        # Beta_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super(Theta_Node,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        self.factors_axis = 0
        self.Ppar = self.P.getParameters()

    def updateParameters(self, factors_selection=None):
        # factors_selection (np array or list): indices of factors that are non-annotated

        # Collect expectations from other nodes
        S = self.markov_blanket['SW'].getExpectations()["ES"]

        # Precompute terms
        # TODO check that it is ok with dimensions of S !! (because S is a pointer, so might be messed up)

        if factors_selection is not None:
            tmp1 = S[:,factors_selection].sum(axis=0)
        else:
            tmp1 = S.sum(axis=0)

        # Perform updates
        Qa = self.Ppar['a'] + tmp1
        Qb = self.Ppar['b'] + S.shape[0]-tmp1

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):

        # Collect parameters and expectations
        Qpar,Qexp = self.getParameters(), self.getExpectations()
        Pa, Pb, Qa, Qb = self.Ppar['a'], self.Ppar['b'], Qpar['a'], Qpar['b']
        QE, QlnE, QlnEInv = Qexp['E'], Qexp['lnE'], Qexp['lnEInv']

        # minus cross entropy of Q and P
        lb_p = ma.masked_invalid( (Pa-1.)*QlnE + (Pb-1.)*QlnEInv - special.betaln(Pa,Pb) ).sum()

        # minus entropy of Q
        lb_q = ma.masked_invalid( (Qa-1.)*QlnE + (Qb-1.)*QlnEInv - special.betaln(Qa,Qb) ).sum()

        return lb_p - lb_q

class Theta_Constant_Node(Constant_Variational_Node):
    """
    Dimensions of Theta_Constant_Node should be (D[m], K)
    """
    def __init__(self, dim, value, N_cells=1):
        super(Theta_Constant_Node, self).__init__(dim, value)
        self.N_cells = N_cells
        self.precompute()

    def precompute(self):
        self.E = self.value
        # TODO this is wrong with missing values -> need to correct N_cells to account for the cells in which a given gene is missing 
        self.lnE = self.N_cells * s.log(self.value)
        self.lnEInv = self.N_cells * s.log(1.-self.value)

    def getExpectations(self):
        return { 'E':self.E, 'lnE':self.lnE, 'lnEInv':self.lnEInv }

    def removeFactors(self, idx, axis=0):
        # Ideally we want this node to use the removeFactors defined in Node()
        # but the problem is that we also need to update the "expectations", so i need
        # to call precompute()
        self.value = s.delete(self.value, idx, axis)
        self.precompute()
        self.updateDim(axis=axis, new_dim=self.dim[axis]-len(idx))

class Z_Node(UnivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None, qE2=None, idx_covariates=None):
        super(Z_Node,self).__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)
        self.precompute()

        # Define indices for covariates
        if idx_covariates is not None:
            self.covariates[idx_covariates] = True

    def precompute(self):
        # Precompute terms to speed up computation
        self.N = self.dim[0]
        self.covariates = np.zeros(self.dim[1], dtype=bool)
        self.factors_axis = 1

    def getLvIndex(self):
        # Method to return the index of the latent variables (without covariates)
        latent_variables = np.array(range(self.dim[1]))
        if any(self.covariates):
            # latent_variables = np.delete(latent_variables, latent_variables[self.covariates])
            latent_variables = latent_variables[~self.covariates]
        return latent_variables

    def updateParameters(self):

        # Collect expectations from the markov blanket
        Y = self.markov_blanket["Y"].getExpectation()
        SWtmp = self.markov_blanket["SW"].getExpectations()
        tau = self.markov_blanket["Tau"].getExpectation()
        latent_variables = self.getLvIndex() # excluding covariates from the list of latent variables

        if "Mu" in self.markov_blanket:
            Mu = self.markov_blanket['Mu'].getExpectation()
        else:
            Mu = self.P.getParameters()["mean"]

        if "Alpha" in self.markov_blanket:
            Alpha = self.markov_blanket['Alpha'].getExpectation()
            Alpha = s.repeat(Alpha[None,:], self.N, axis=0)
        else:
            Alpha = 1./self.P.getParameters()["var"]

        # Check dimensionality of Tau and expand if necessary (for Jaakola's bound only)
        for m in xrange(len(Y)):
            if tau[m].shape != Y[m].shape:
                tau[m] = s.repeat(tau[m].copy()[None,:], self.N, axis=0)
                tau[m] = ma.masked_where(ma.getmask(Y[m]), tau[m])

        # Collect parameters from the P and Q distributions of this node
        Q = self.Q.getParameters().copy()
        Qmean, Qvar = Q['mean'], Q['var']

        Qvar_copy = Qvar.copy()
        Qmean_copy = Qmean.copy()

        ### START LOOP ###
        M = len(Y)
        for k in latent_variables:
            tmp = s.zeros((self.N,))
            for m in xrange(M):
                if k in latent_variables:
                    Qvar[:,k] += ma.dot(tau[m],SWtmp[m]["ESWW"][:,k])
                tmp += ma.dot(tau[m]*(Y[m] - s.dot( Qmean[:,s.arange(self.dim[1])!=k] , SWtmp[m]["E"][:,s.arange(self.dim[1])!=k].T )), SWtmp[m]["E"][:,k])
            if k in latent_variables:
                Qvar[:,k] = 1./(Alpha[:,k]+Qvar[:,k])
            Qmean[:,k] = Qvar[:,k] * (  Alpha[:,k]*Mu[:,k] + tmp )
        ### END LOOP ###

        # ### START CONCATENATE ###
        # # Qvar = Qvar_copy
        # # Qmean = Qmean_copy
        # # Concatenate multi-view nodes to avoid looping over M (maybe its not a good idea)
        # M = len(Y)
        # Y = ma.concatenate([Y[m] for m in xrange(M)],axis=1)
        # SW = s.concatenate([SWtmp[m]["E"]for m in xrange(M)],axis=0)
        # SWW = s.concatenate([SWtmp[m]["ESWW"] for m in xrange(M)],axis=0)
        # tau = ma.concatenate([tau[m] for m in xrange(M)],axis=1)


        # # Update variance
        # # Qvar_copy = Qvar.copy()
        # CHECK THAT HTIS IS CORRECT WITH ALPHA BEING A K-VECTOR
        # Qvar = 1./(Alpha + ma.dot(tau,SWW))

        # # restoring values of the variance for the covariates
        # if any(self.covariates):
        #     Qvar[self.covariates] = Qvar_copy[self.covariates]

        # # Update mean
        # for k in latent_variables:
        #     Qmean[:,k] = Qvar[:,k] * (  Alpha*Mu[:,k] +
        #                                 ma.dot(tau*(Y - s.dot( Qmean[:,s.arange(self.dim[1])!=k] , SW[:,s.arange(self.dim[1])!=k].T )), SW[:,k])  )
        # ### END CONCATENATE ###

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean=Qmean, var=Qvar)

    def calculateELBO(self):
        # Collect parameters and expectations of current node
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
        Qmean, Qvar = Qpar['mean'], Qpar['var']
        QE, QE2 = Qexp['E'],Qexp['E2']

        if "Mu" in self.markov_blanket:
            PE, PE2 = self.markov_blanket['Mu'].getExpectations()['E'], self.markov_blanket['Mu'].getExpectations()['E2']
        else:
            PE, PE2 = self.P.getParameters()["mean"], s.zeros((self.N,self.dim[1]))

        if "Alpha" in self.markov_blanket:
            Alpha = self.markov_blanket['Alpha'].getExpectations().copy() # Notice that this Alpha is the ARD prior on Z, not on W.
            Alpha["E"] = s.repeat(Alpha["E"][None,:], self.N, axis=0)
            Alpha["lnE"] = s.repeat(Alpha["lnE"][None,:], self.N, axis=0)
        else:
            Alpha = { 'E':1./self.P.getParameters()["var"], 'lnE':s.log(1./self.P.getParameters()["var"]) }

        # This ELBO term contains only cross entropy between Q and P,and entropy of Q. So the covariates should not intervene at all
        latent_variables = self.getLvIndex()
        Alpha["E"], Alpha["lnE"] = Alpha["E"][:,latent_variables], Alpha["lnE"][:,latent_variables]
        Qmean, Qvar = Qmean[:, latent_variables], Qvar[:, latent_variables]
        PE, PE2 = PE[:, latent_variables], PE2[:, latent_variables]
        QE, QE2 = QE[:, latent_variables], QE2[:, latent_variables]

        # compute term from the exponential in the Gaussian
        tmp1 = 0.5*QE2 - PE*QE + 0.5*PE2
        tmp1 = -(tmp1 * Alpha['E']).sum()

        # compute term from the precision factor in front of the Gaussian
        tmp2 = 0.5*Alpha["lnE"].sum()

        lb_p = tmp1 + tmp2
        # lb_q = -(s.log(Qvar).sum() + self.N*self.dim[1])/2. # I THINK THIS IS WRONG BECAUSE SELF.DIM[1] ICNLUDES COVARIATES
        lb_q = -(s.log(Qvar).sum() + self.N*len(latent_variables))/2.

        return lb_p-lb_q

class AlphaZ_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        # Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super(AlphaZ_Node,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        self.factors_axis = 0

    def updateParameters(self):

        # Collect expectations from other nodes
        Ztmp = self.markov_blanket["Z"].getExpectations()
        Z, ZZ = Ztmp["E"], Ztmp["E2"]

        if 'Mu' in self.markov_blanket:
            Mutmp = self.markov_blanket['Mu'].getExpectations()
            MuE, MuE2 = Mutmp['E'], Mutmp['E2']
        else:
            MuE = self.markov_blanket['Z'].P.getParameters()["mean"]
            MuE2 = s.zeros(Z.shape)

        # Collect parameters from the P distributions of this node
        P = self.P.getParameters()
        Pa, Pb = P['a'], P['b']

        # Perform updates
        Qa = Pa + 0.5*Z.shape[0]
        Qb = Pb + 0.5*ZZ.sum(axis=0)
        Qb -= 2*(Z*MuE).sum(axis=0)
        Qb += MuE2.sum(axis=0)

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.expectations['E'], self.Q.expectations['lnE']

        # This ELBO term contains only cross entropy between Q and P,and entropy of Q. So the covariates should not intervene at all
        latent_variables = self.markov_blanket["Z"].getLvIndex()
        Pa, Pb, Qa, Qb = Pa[latent_variables], Pb[latent_variables], Qa[latent_variables], Qb[latent_variables]
        QE, QlnE = QE[latent_variables], QlnE[latent_variables]

        # Do the calculations
        # lb_p = s.nansum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) ) + s.nansum((Pa-1)*QlnE) - s.nansum(Pb*QE)
        lb_p = (Pa*s.log(Pb)).sum() - (special.gammaln(Pa)).sum() + ((Pa-1)*QlnE).sum() - (Pb*QE).sum()
        lb_q = (Qa*s.log(Qb)).sum() + ((Qa-1)*QlnE).sum() - (Qb*QE).sum() - special.gammaln(Qa).sum()

        return lb_p - lb_q

class MuZ_Node(UnivariateGaussian_Unobserved_Variational_Node):
    """ """
    def __init__(self, pmean, pvar, qmean, qvar, clusters, n_Z, cluster_dic=None, qE=None, qE2=None):
        # compute dim from numbers of clusters (n_clusters * Z)
        self.clusters = clusters
        self.N = len(self.clusters)
        self.n_clusters = len(np.unique(clusters))
        dim = (self.n_clusters, n_Z)
        self.factors_axis = 1
        super(Cluster_Node_Gaussian, self).__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)

    def getExpectations(self):
        # reshape the values to N_samples * N_factors and return
        QExp = self.Q.getExpectations()
        expanded_expectation = QExp['E'][self.clusters, :]
        expanded_E2 = QExp['E2'][self.clusters, :]
        # do we need to expand the variance as well -> not used I think
        return {'E': expanded_expectation , 'E2': expanded_E2}

    def updateParameters(self):
        Ppar = self.P.getParameters()
        Z = self.markov_blanket['Z'].Q.getExpectation()

        if "Alpha" in self.markov_blanket:
            Alpha = self.markov_blanket['Alpha'].getExpectation().copy() # Notice that this Alpha is the ARD prior on Z, not on W.
            Alpha = s.repeat(Alpha[None,:], self.N, axis=0)
        else:
            Alpha = 1./self.markov_blanket['Z'].P.getParameters()["var"]

        Qmean, Qvar = self.Q.getParameters()['mean'], self.Q.getParameters()['var']
        ZTauMean = Z * Alpha

        # TODO merge two loops when sure it's clean
        # update of the variance
        for c in range(self.n_clusters):
            mask = (self.clusters == c)
            tmp = (Alpha[mask, :]).sum(axis=0)
            Qvar[c,:] = tmp
        Qvar += 1./Ppar['var']
        Qvar = 1./Qvar

        # update of the mean
        for c in range(self.n_clusters):
            mask = (self.clusters == c)
            tmp = (ZTauMean[mask, :]).sum(axis=0)
            Qmean[c,:] = tmp
        Qmean = Qmean + Ppar['mean']/Ppar['var']
        Qmean *= Qvar

        self.Q.setParameters(mean=Qmean, var=Qvar)

    def calculateELBO(self):
        PParam = self.P.getParameters()
        PVar, Pmean = PParam['var'], PParam['mean']

        QExp = self.Q.getExpectations()
        QE2, QE = QExp['E2'], QExp['E']

        Qvar = self.Q.getParameters()['var']

        # Cluster terms corersponding to covariates should not intervene
        # filtering the covariates out
        latent_variables = self.markov_blanket['Z'].getLvIndex()
        PVar, Pmean = PVar[:, latent_variables], Pmean[:, latent_variables]
        QE2, QE = QE2[:, latent_variables], QE[:, latent_variables]
        Qvar = Qvar[:, latent_variables]

        # minus cross entropy
        tmp = -(0.5 * s.log(PVar)).sum()
        tmp2 = - ((0.5/PVar) * (QE2 - 2.*QE*Pmean + Pmean**2.)).sum()

        # entropy of Q
        tmp3 = 0.5 * (s.log(Qvar)).sum()
        tmp3 += 0.5 * self.dim[0] * len(latent_variables)

        return tmp + tmp2 + tmp3
