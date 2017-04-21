from __future__ import division
import numpy.linalg  as linalg
# from numpy_sugar.linalg import dotd
from numpy_sugar.ma import dotd
import numpy.ma as ma
import numpy as np

# Import manually defined functions
from variational_nodes import *
from utils import *
from nodes import Constant_Node
from mixed_nodes import Mixed_Theta_Nodes

import scipy.special as special

"""
###################################################
## Updates for the Sparse Group Factor Analysis ##
###################################################

Extensions with respect to the Gaussian Group Factor Analysis:
- Element-wise spike and slab

(Derivation of equations can be found in file XX)

Current nodes:
    Y_Node: observed data
    SW_Node: spike and slab weights
    Tau_Node: precision of the noise
    Alpha_Node: ARD precision
    Z_Node: latent variables

Each node is a Variational_Node() class with the following main variables:
    Important methods:
    - precompute: precompute some terms to speed up the calculations
    - calculateELBO: calculate evidence lower bound using current estimates of expectations/params
    - getParameters: return current parameters
    - getExpectations: return current expectations
    - updateParameters: update parameters using current estimates of expectations
    - updateExpectations: update expectations using current estimates of parameters

    Important attributes:
    - markov_blanket: dictionary that defines the set of nodes that are in the markov blanket of the current node
    - Q: an instance of Distribution() which contains the specification of the variational distribution
    - P: an instance of Distribution() which contains the specification of the prior distribution
    - dim: dimensionality of the node

"""

class Y_Node(Constant_Variational_Node):
    """
    Observed data node
    """
    def __init__(self, dim, value):
        Constant_Variational_Node.__init__(self, dim, value)

        # Create a boolean mask of the data to hidden missing values
        if type(self.value) != ma.MaskedArray:
            self.mask()

        # Precompute some terms
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        # self.N = self.dim[0]
        self.N = self.dim[0] - ma.getmask(self.value).sum(axis=0)
        self.D = self.dim[1]
        # self.likconst = -0.5*self.N*self.D*s.log(2*s.pi)
        self.likconst = -0.5*s.sum(self.N)*s.log(2.*s.pi)

    def mask(self):
        # Mask the observations if they have missing values
        self.value = ma.masked_invalid(self.value)

    def calculateELBO(self):
        # Calculate evidence lower bound using the trick that the update of Tau already contains the Gaussian likelihod 
        tauQ_param = self.markov_blanket["Tau"].getParameters("Q")
        tauP_param = self.markov_blanket["Tau"].getParameters("P")
        tau_exp = self.markov_blanket["Tau"].getExpectations()
        lik = self.likconst + 0.5*s.sum(self.N*(tau_exp["lnE"])) - s.dot(tau_exp["E"],tauQ_param["b"]-tauP_param["b"])
        return lik

# class Z_Node(UnivariateGaussian_Unobserved_Variational_Node):
#     def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None, qE2=None, idx_covariates=None):
#         super(Z_Node,self).__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)
#         self.precompute()

#         # Define indices for covariates
#         if idx_covariates is not None:
#             self.covariates[idx_covariates] = True

#     def precompute(self):
#         # Precompute terms to speed up computation
#         self.N = self.dim[0]
#         self.covariates = np.zeros(self.dim[1], dtype=bool)
#         self.factors_axis = 1

#     def getLvIndex(self):
#         # Method to return the index of the latent variables (without covariates)
#         latent_variables = np.array(range(self.dim[1]))
#         if any(self.covariates):
#             latent_variables = np.delete(latent_variables, latent_variables[self.covariates])
#         return latent_variables

#     def updateParameters(self):

#         # Collect expectations from the markov blanket
#         Y = self.markov_blanket["Y"].getExpectation()
#         SWtmp = self.markov_blanket["SW"].getExpectations()
#         tau = self.markov_blanket["Tau"].getExpectation()
#         Mu = self.markov_blanket['Cluster'].getExpectation()

#         # Check dimensionality of Tau and expand if necessary (for Jaakola's bound only)
#         for m in xrange(len(Y)):
#             if tau[m].shape != Y[m].shape:
#                 tau[m] = s.repeat(tau[m][None,:], self.N, axis=0)

#         # Collect parameters from the P and Q distributions of this node
#         P,Q = self.P.getParameters(), self.Q.getParameters()
#         Pvar, Qmean, Qvar = P['var'], Q['mean'], Q['var']

#         # Concatenate multi-view nodes to avoid looping over M (maybe its not a good idea)
#         M = len(Y)
#         Y = ma.concatenate([Y[m] for m in xrange(M)],axis=1)
#         SW = s.concatenate([SWtmp[m]["E"]for m in xrange(M)],axis=0)
#         SWW = s.concatenate([SWtmp[m]["ESWW"] for m in xrange(M)],axis=0)
#         tau = s.concatenate([tau[m] for m in xrange(M)],axis=1)

#         # Update variance
#         Qvar_copy = Qvar.copy()
#         Qvar = 1./(1./Pvar + s.dot(tau,SWW))
#         # print "start test:"
#         # print 1./(1./Pvar + s.dot(tau,SWW))
#         # print "end test:"

#         # restoring values of the variance for the covariates
#         if any(self.covariates):
#             Qvar[:, self.covariates] = Qvar_copy[:, self.covariates]

#         # Update mean
#         latent_variables = self.getLvIndex() # excluding covariates from the list of latent variables
#         for k in latent_variables:
#             Qmean[:,k] = Qvar[:,k] * (  1./Pvar[:,k]*Mu[:,k] + ma.dot(
#                 tau*(Y - s.dot( Qmean[:,s.arange(self.dim[1])!=k] , SW[:,s.arange(self.dim[1])!=k].T )),
#                 SW[:,k])  )

#         # Save updated parameters of the Q distribution
#         self.Q.setParameters(mean=Qmean, var=Qvar)

#     def calculateELBO(self):
#         # Collect parameters and expectations of current node
#         Ppar,Qpar,Qexp = self.P.getParameters(), self.Q.getParameters(), self.Q.getExpectations()
#         Pvar, Qmean, Qvar = Ppar['var'], Qpar['mean'], Qpar['var']
#         PE, PE2 = self.markov_blanket['Cluster'].getExpectations()['E'], self.markov_blanket['Cluster'].getExpectations()['E2']
#         QE, QE2 = Qexp['E'],Qexp['E2']

#         # This ELBO term contains only cross entropy between Q and P,and entropy of Q. So the covariates should not intervene at all
#         latent_variables = self.getLvIndex()
#         Pvar, Qmean, Qvar = Pvar[:, latent_variables], Qmean[:, latent_variables], Qvar[:, latent_variables]
#         PE, PE2 = PE[:, latent_variables], PE2[:, latent_variables]
#         QE, QE2 = QE[:, latent_variables],QE2[:, latent_variables]

#         # compute term from the exponential in the Gaussian
#         tmp1 = 0.5*QE2 - PE*QE + 0.5*PE2
#         tmp1 = -(tmp1/Pvar).sum()

#         # compute term from the precision factor in front of the Gaussian
#         tmp2 = - (s.log(Pvar)/2.).sum()

#         lb_p = tmp1 + tmp2
#         lb_q = - (s.log(Qvar).sum() + self.N*self.dim[1])/2.

#         return lb_p-lb_q

class Tau_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        # Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super(Tau_Node,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        self.D = self.dim[0]
        # self.lbconst = s.sum(self.D*(self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a'])))
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
        term1 = (Y**2.).sum(axis=0).data
        term2 = (2.*Y*s.dot(Z,SW.T)).sum(axis=0).data
        term3 = ma.array(ZZ.dot(SWW.T), mask=ma.getmask(Y)).sum(axis=0)
        SWZ = ma.array(SW.dot(Z.T), mask=ma.getmask(Y).T)
        # term4 = s.diag(ma.dot( SWZ,SWZ.T )) - ma.array(s.dot(Z**2,(SW**2).T),mask=ma.getmask(Y)).sum(axis=0) # VERY SLOW
        term4 = dotd(SWZ, SWZ.T) - ma.array(s.dot(Z**2,(SW**2).T),mask=ma.getmask(Y)).sum(axis=0)
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

class AlphaW_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        # Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super(AlphaW_Node,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        # self.lbconst = self.K * ( self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )
        # self.lbconst = s.sum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) )
        self.factors_axis = 0

    def updateParameters(self):

        # Collect expectations from other nodes
        tmp = self.markov_blanket["SW"].getExpectations()
        ES,EWW,ESWW = tmp["ES"],tmp["EWW"],tmp["ESWW"]

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
        QE, QlnE = self.Q.expectations['E'], self.Q.expectations['lnE']

        # Do the calculations
        lb_p = s.sum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) ) + s.sum((Pa-1)*QlnE) - s.sum(Pb*QE)
        lb_q = s.sum(Qa*s.log(Qb)) + s.sum((Qa-1)*QlnE) - s.sum(Qb*QE) - s.sum(special.gammaln(Qa))

        return lb_p - lb_q

class SW_Node(BernoulliGaussian_Unobserved_Variational_Node):
    # TOO MANY ARGUMENTS, SHOULD WE USE **KWARGS AND *KARGS ONLY?
    # def __init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0=None, qEW_S1=None, qES=None):
    def __init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0=None, qEW_S1=None, qES=None):
        super(SW_Node,self).__init__(dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0, qEW_S1, qES)
        self.precompute()

    def precompute(self):
        self.D = self.dim[0]
        # self.K = self.dim[1]
        self.factors_axis = 1

    def updateParameters(self):

        # Collect expectations from other nodes
        Ztmp = self.markov_blanket["Z"].getExpectations()
        Z,ZZ = Ztmp["E"],Ztmp["E2"]
        tau = self.markov_blanket["Tau"].getExpectation()
        Y = self.markov_blanket["Y"].getExpectation()
        alpha = self.markov_blanket["Alpha"].getExpectation()
        thetatmp = self.markov_blanket['Theta'].getExpectations() # TODO make general in mixed nodw
        theta_lnE, theta_lnEInv  = thetatmp['lnE'], thetatmp['lnEInv']

        # Collect parameters and expectations from P and Q distributions of this node
        SW = self.Q.getExpectations()["E"]
        Q = self.Q.getParameters()
        Qmean_S1, Qvar_S1, Qtheta = Q['mean_S1'], Q['var_S1'], Q['theta']

        # Check dimensions of Theta and Tau and expand if necessary
        # I THINK THIS SHOULD NOT BE HERE, BUT WE COULDNT COME UP WITH A SOLUTION..
        if theta_lnE.shape != Qmean_S1.shape:
            theta_lnE = s.repeat(theta_lnE[None,:],Qmean_S1.shape[0],0)
        if theta_lnEInv.shape != Qmean_S1.shape:
            theta_lnEInv = s.repeat(theta_lnEInv[None,:],Qmean_S1.shape[0],0)
        if tau.shape != Y.shape:
            tau = s.repeat(tau[None,:], Y.shape[0], axis=0)

        # Update each latent variable in turn
        for k in xrange(self.dim[1]):

            # Calculate intermediate steps
            term1 = (theta_lnE-theta_lnEInv)[:,k]
            term2 = 0.5*s.log(alpha[k])
            term3 = 0.5*s.log(ma.dot(ZZ[:,k],tau) + alpha[k])
            term4_tmp1 = ma.dot((tau*Y).T,Z[:,k]).data
            tmp = tau * s.dot((Z[:,k]*Z[:,s.arange(self.dim[1])!=k].T).T, SW[:,s.arange(self.dim[1])!=k].T)
            term4_tmp2 = ma.array(tmp, mask=ma.getmask(Y)).sum(axis=0)
            term4_tmp3 = ma.dot(ZZ[:,k].T,tau) + alpha[k]
            term4 = 0.5*s.divide((term4_tmp1-term4_tmp2)**2,term4_tmp3)

            # Update S
            # NOTE there could be some precision issues in S --> loads of 1s in result
            Qtheta[:,k] = 1/(1+s.exp(-(term1+term2-term3+term4)))

            # Update W
            Qvar_S1[:,k] = s.divide(1,term4_tmp3)
            Qmean_S1[:,k] = Qvar_S1[:,k]*(term4_tmp1-term4_tmp2)

            # Update Expectations for the next iteration
            SW[:,k] = Qtheta[:,k] * Qmean_S1[:,k]

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean_S0=s.zeros((self.D,self.dim[1])), var_S0=s.repeat(1/alpha[None,:],self.D,0),
                             mean_S1=Qmean_S1, var_S1=Qvar_S1, theta=Qtheta )

    def calculateELBO(self):

        # Collect parameters and expectations
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
        S,WW = Qexp["ES"], Qexp["EWW"]
        Qvar = Qpar['var_S1']
        alpha = self.markov_blanket["Alpha"].getExpectations()
        theta = self.markov_blanket['Theta'].getExpectations()


        # Calculate ELBO for W
        lb_pw = (self.D*alpha["lnE"].sum() - s.sum(alpha["E"]*WW))/2
        # lb_qw = -0.5*self.dim[1]*self.D - 0.5*s.log(S*Qvar + ((1-S)/alpha["E"])).sum()
        lb_qw = -0.5*self.dim[1]*self.D - 0.5*(S*s.log(Qvar) + (1-S)*s.log(1/alpha["E"])).sum()
        lb_w = lb_pw - lb_qw

        # Calculate ELBO for S
        lb_ps_tmp = S*theta['lnE'] + (1.-S)*theta['lnEInv']
        lb_qs_tmp = S*s.log(S) + (1.-S)*s.log(1.-S)

        # Ignore NAs
        lb_ps_tmp[s.isnan(lb_ps_tmp)] = 0.
        lb_qs_tmp[s.isnan(lb_qs_tmp)] = 0.

        lb_ps = s.sum(lb_ps_tmp)
        lb_qs = s.sum(lb_qs_tmp)
        lb_s = lb_ps - lb_qs

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
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
        Pa, Pb, Qa, Qb = self.Ppar['a'], self.Ppar['b'], Qpar['a'], Qpar['b']
        QE, QlnE, QlnEInv = Qexp['E'], Qexp['lnE'], Qexp['lnEInv']

        # minus cross entropy of Q and P
        tmp1 = (Pa-1.)*QlnE + (Pb-1.)*QlnEInv - special.betaln(Pa,Pb)
        lb_p = tmp1.sum()

        # minus entropy of Q
        tmp2 = (Qa-1.)*QlnE + (Qb-1.)*QlnEInv - special.betaln(Qa,Qb)
        lb_q = tmp2.sum()

        return lb_p - lb_q

class Theta_Constant_Node(Constant_Variational_Node):
    """
    Dimensions of Theta_Constant_Node should be (D[m], K)
    """
    def __init__(self, dim, value, N_cells):
        super(Theta_Constant_Node, self).__init__(dim, value)
        self.N_cells = N_cells
        self.precompute()

    def precompute(self):
        self.E = self.value
        self.lnE = self.N_cells * s.log(self.value)
        self.lnEInv = self.N_cells * s.log(1.-self.value)

    def getExpectations(self):
        return { 'E':self.E, 'lnE':self.lnE, 'lnEInv':self.lnEInv }

    def removeFactors(self, idx, axis=1):
        # Ideally we want this node to use the removeFactors defined in Node()
        # but the problem is that we also need to update the "expectations", so i need
        # to call precompute()
        self.value = s.delete(self.value, idx, axis)
        self.precompute()
        self.updateDim(axis=axis, new_dim=self.dim[axis]-len(idx))



####################
## ARD prior on Z ##
####################

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
            latent_variables = np.delete(latent_variables, latent_variables[self.covariates])
        return latent_variables

    def updateParameters(self):

        # Collect expectations from the markov blanket
        Y = self.markov_blanket["Y"].getExpectation()
        SWtmp = self.markov_blanket["SW"].getExpectations()
        tau = self.markov_blanket["Tau"].getExpectation()
        Alpha = self.markov_blanket['Alpha'].getExpectation() # Notice that this Alpha is the ARD prior on Z, not on W.

        # Check dimensionality of Tau and expand if necessary (for Jaakola's bound only)
        for m in xrange(len(Y)):
            if tau[m].shape != Y[m].shape:
                tau[m] = s.repeat(tau[m][None,:], self.N, axis=0)

        # Collect parameters from the P and Q distributions of this node
        # P,Q = self.P.getParameters(), self.Q.getParameters()
        Q = self.Q.getParameters()
        # Pvar, Qmean, Qvar = P['var'], Q['mean'], Q['var']
        Qmean, Qvar = Q['mean'], Q['var']

        # Concatenate multi-view nodes to avoid looping over M (maybe its not a good idea)
        M = len(Y)
        Y = ma.concatenate([Y[m] for m in xrange(M)],axis=1)
        SW = s.concatenate([SWtmp[m]["E"]for m in xrange(M)],axis=0)
        SWW = s.concatenate([SWtmp[m]["ESWW"] for m in xrange(M)],axis=0)
        tau = s.concatenate([tau[m] for m in xrange(M)],axis=1)

        # Update variance
        Qvar_copy = Qvar.copy()
        Qvar = 1./(Alpha + s.dot(tau,SWW))

        # restoring values of the variance for the covariates
        if any(self.covariates):
            Qvar[:, self.covariates] = Qvar_copy[:, self.covariates]

        # Update mean
        latent_variables = self.getLvIndex() # excluding covariates from the list of latent variables
        for k in latent_variables:
            Qmean[:,k] = Qvar[:,k] * ( ma.dot(tau*(Y - s.dot( Qmean[:,s.arange(self.dim[1])!=k] , SW[:,s.arange(self.dim[1])!=k].T )), SW[:,k])  )

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean=Qmean, var=Qvar)

    def calculateELBO(self):
        # Collect parameters and expectations of current node
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
        Qmean, Qvar = Qpar['mean'], Qpar['var']
        Alpha = self.markov_blanket['Alpha'].getExpectations()
        QE, QE2 = Qexp['E'],Qexp['E2']

        # This ELBO term contains only cross entropy between Q and P,and entropy of Q. So the covariates should not intervene at all
        latent_variables = self.getLvIndex()
        Alpha["E"], Alpha["lnE"] = Alpha["E"][latent_variables], Alpha["lnE"][latent_variables]
        Qmean, Qvar = Qmean[:, latent_variables], Qvar[:, latent_variables]
        QE, QE2 = QE[:, latent_variables], QE2[:, latent_variables]

        lb_p = (self.N*Alpha["lnE"].sum() - s.sum(Alpha["E"]*QE2))/2.
        lb_q = - (s.log(Qvar).sum() + self.N*self.dim[1])/2.

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

        # Collect parameters from the P distributions of this node
        P = self.P.getParameters()
        Pa, Pb = P['a'], P['b']

        # Perform updates
        Qa = Pa + 0.5*Z.shape[0]
        Qb = Pb + 0.5*ZZ.sum(axis=0)

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.expectations['E'], self.Q.expectations['lnE']

        # Do the calculations
        lb_p = s.sum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) ) + s.sum((Pa-1)*QlnE) - s.sum(Pb*QE)
        lb_q = s.sum(Qa*s.log(Qb)) + s.sum((Qa-1)*QlnE) - s.sum(Qb*QE) - s.sum(special.gammaln(Qa))

        return lb_p - lb_q
