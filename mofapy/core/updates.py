from __future__ import division
import numpy.ma as ma
import numpy as np
import warnings
from time import time
import scipy.special as special

# Import manually defined functions
from .variational_nodes import *
from .utils import *
from .nodes import Constant_Node
from .mixed_nodes import Mixed_Theta_Nodes


warnings.filterwarnings('ignore')

"""
Module to define the nodes and the corresponding updates of the model
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

        # Create a boolean mask of the data to hide missing values
        if type(self.value) != ma.MaskedArray:
            self.mask()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        self.N = self.dim[0] - ma.getmask(self.value).sum(axis=0)
        self.likconst = -0.5*s.sum(self.N)*s.log(2.*s.pi)
        self.means = self.value.mean(axis=0).data

    def mask(self):
        # Mask the observations if they have missing values
        self.value = ma.masked_invalid(self.value)

    def getMask(self):
        return ma.getmask(self.value)

    def calculateELBO(self):
        # Calculate evidence lower bound
        # We use the trick that the update of Tau already contains the Gaussian likelihod.
        # However, it is important that the lower bound is calculated after the update of Tau is performed
        tauQ_param = self.markov_blanket["Tau"].getParameters("Q")
        tauP_param = self.markov_blanket["Tau"].getParameters("P")
        tau_exp = self.markov_blanket["Tau"].getExpectations(expand=False)

        # Important: this assumes that the Tau update has been done beforehand
        lik = self.likconst + 0.5*s.sum(self.N*(tau_exp["lnE"])) - s.dot(tau_exp["E"],tauQ_param["b"]-tauP_param["b"])
        return lik

class Tau_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        super(Tau_Node,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

    def precompute(self):

        # constant update of Qa
        Y = self.markov_blanket["Y"].getExpectation()
        # mask = ma.getmask(Y)
        mask = self.markov_blanket["Y"].getMask()
        # Y = Y.data

        Qa = self.P.getParameters()['a'] + (Y.shape[0] - mask.sum(axis=0))/2.
        self.Q.params['a'] = Qa

        # constant lower bound term
        self.lbconst = s.sum(self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']))

    def updateParameters(self):

        # Collect expectations from other nodes
        Y = self.markov_blanket["Y"].getExpectation()
        # mask = ma.getmask(Y)
        mask = self.markov_blanket["Y"].getMask()
        
        Wtmp = self.markov_blanket["SW"].getExpectations()
        Ztmp = self.markov_blanket["Z"].getExpectations()
        SW, SWW = Wtmp["E"], Wtmp["ESWW"]
        Z, ZZ = Ztmp["E"], Ztmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # Mask matrices
        Y = Y.data
        Y[mask] = 0.

        # Calculate temporary terms for the update
        ZW = Z.dot(SW.T)
        ZW[mask] = 0.

        # Calculate terms for the update
        term1 = s.square(Y).sum(axis=0)

        term2 = ZZ.dot(SWW.T)
        term2[mask] = 0
        term2 = term2.sum(axis=0)

        term3 = np.dot(np.square(Z),np.square(SW).T)
        term3[mask] = 0.
        term3 = -term3.sum(axis=0)
        term3 += np.square(ZW).sum(axis=0)

        ZW *= Y  # Warning for developers: ZW becomes ZWY
        term4 = 2.*(ZW.sum(axis=0))

        tmp = term1 + term2 + term3 - term4

        # Perform updates of the Q distribution
        Qb = Pb + tmp/2.

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=self.Q.params['a'], b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations from current node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.expectations['E'], self.Q.expectations['lnE']

        # Do the calculations
        lb_p = self.lbconst + s.sum((Pa-1.)*QlnE) - s.sum(Pb*QE)
        lb_q = s.sum(Qa*s.log(Qb)) + s.sum((Qa-1.)*QlnE) - s.sum(Qb*QE) - s.sum(special.gammaln(Qa))

        return lb_p - lb_q

    def getExpectations(self, expand=True):
        QExp = self.Q.getExpectations()
        if expand:
            N = self.markov_blanket['Z'].dim[0]
            expanded_E = s.repeat(QExp['E'][None, :], N, axis=0)
            expanded_lnE = s.repeat(QExp['lnE'][None, :], N, axis=0)
            return {'E': expanded_E, 'lnE': expanded_lnE}
        else:
            return QExp

    def getExpectation(self, expand=True):
        QExp = self.getExpectations(expand)
        return QExp['E']

class Alpha_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        # Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super(Alpha_Node,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

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

    def getExpectations(self, expand=False):
        QExp = self.Q.getExpectations()
        QExp['E'] = QExp['E']
        QExp['lnE'] = QExp['lnE']
        if expand:
            D = self.markov_blanket['SW'].dim[0]
            expanded_E = s.repeat(QExp['E'][None, :], D, axis=0)
            expanded_lnE = s.repeat(QExp['lnE'][None, :], D, axis=0)
            return {'E': expanded_E, 'lnE': expanded_lnE}
        else:
            return QExp

    def getExpectation(self, expand=False):
        QExp = self.getExpectations(expand)
        return QExp['E']

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
    def __init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0=None, qEW_S1=None, qES=None):
        super(SW_Node,self).__init__(dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0, qEW_S1, qES)

    def precompute(self):
        self.factors_axis = 1

    def updateParameters(self):
        # Collect expectations from other nodes
        Ztmp = self.markov_blanket["Z"].getExpectations()
        Z,ZZ = Ztmp["E"],Ztmp["E2"]
        tau = self.markov_blanket["Tau"].getExpectation()
        Y = self.markov_blanket["Y"].getExpectation()
        alpha = self.markov_blanket["Alpha"].getExpectation(expand=False)
        thetatmp = self.markov_blanket["Theta"].getExpectations()
        theta_lnE, theta_lnEInv  = thetatmp['lnE'], thetatmp['lnEInv']
        mask = self.markov_blanket["Y"].getMask()
        # mask = ma.getmask(Y)

        # Collect parameters and expectations from P and Q distributions of this node
        SW = self.Q.getExpectations()["E"]
        Q = self.Q.getParameters()
        Qmean_S1, Qvar_S1, Qtheta = Q['mean_S1'], Q['var_S1'], Q['theta']

        # Mask matrices
        # Ymean = Y.mean(axis=0)
        Y = Y.data
        Y[mask] = 0.
        tau[mask] = 0.

        # precompute terms used for all factors
        tauYT = (tau*Y).T

        for k in range(self.dim[1]):
            # Calculate intermediate steps
            term1 = (theta_lnE-theta_lnEInv)[k]
            term2 = 0.5*s.log(alpha[k])
            term3 = 0.5*s.log(s.dot(ZZ[:,k], tau) + alpha[k])

            term4_tmp1 = s.dot(tauYT,Z[:,k])

            term4_tmp2_1 = SW[:,s.arange(self.dim[1])!=k].T
            term4_tmp2_2 = (Z[:,k]*Z[:,s.arange(self.dim[1])!=k].T).T
            term4_tmp2 = s.dot(term4_tmp2_2, term4_tmp2_1)
            term4_tmp2 *= tau
            term4_tmp2 = term4_tmp2.sum(axis=0)

            term4_tmp3 = s.dot(ZZ[:,k].T,tau) + alpha[k]

            term4 = 0.5*s.divide(s.square(term4_tmp1-term4_tmp2),term4_tmp3)

            # Update S
            # NOTE there could be some precision issues in S --> loads of 1s in result
            Qtheta[:,k] = 1./(1.+s.exp(-(term1+term2-term3+term4)))

            # Update W
            Qvar_S1[:,k] = 1./term4_tmp3
            Qmean_S1[:,k] = Qvar_S1[:,k]*(term4_tmp1-term4_tmp2)

            # Update Expectations for the next iteration
            SW[:,k] = Qtheta[:,k] * Qmean_S1[:,k]

        # Save updated parameters of the Q distribution
        # self.Q.setParameters(mean_S0=0., var_S0=1./alpha, mean_S1=Qmean_S1, var_S1=Qvar_S1, theta=Qtheta )
        self.Q.setParameters(mean_S0=s.zeros((self.dim[0],self.dim[1])), var_S0=s.repeat(1./alpha[None,:],self.dim[0],0), mean_S1=Qmean_S1, var_S1=Qvar_S1, theta=Qtheta )

    def calculateELBO(self):

        # Collect parameters and expectations
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
        S,WW = Qexp["ES"], Qexp["EWW"]
        Qvar = Qpar['var_S1']
        theta = self.markov_blanket['Theta'].getExpectations()
        alpha = self.markov_blanket["Alpha"].getExpectations(expand=False)

        # Calculate ELBO for W
        lb_pw = 0.5*(self.dim[0]*alpha["lnE"].sum() - s.sum(alpha["E"]*WW))
        lb_qw = -0.5*self.dim[1]*self.dim[0] - 0.5*(S*s.log(Qvar) + (1.-S)*s.log(1./alpha["E"])).sum()
        lb_w = lb_pw - lb_qw

        # Calculate ELBO for S
        # TO-DO: CHECK THAT THE BROADCASTING IS CORRECT FOR THETA
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

    def precompute(self):
        self.factors_axis = 0
        self.Ppar = self.P.getParameters()

    def updateParameters(self, factors_selection=None):
        # factors_selection (np array or list): indices of factors that are non-annotated

        # Collect expectations from other nodes
        S = self.markov_blanket['SW'].getExpectations()["ES"]

        # Precompute terms
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
        # TO-DO: IS THIS CORRECT? DO WE NEED TO SUM OVER D?
        
        # Collect parameters and expectations
        Qpar,Qexp = self.getParameters(), self.getExpectations()
        Pa, Pb, Qa, Qb = self.Ppar['a'], self.Ppar['b'], Qpar['a'], Qpar['b']
        QE, QlnE, QlnEInv = Qexp['E'], Qexp['lnE'], Qexp['lnEInv']

        # minus cross entropy of Q and P
        lb_p = (Pa-1.)*QlnE + (Pb-1.)*QlnEInv - special.betaln(Pa,Pb)
        lb_p[np.isnan(lb_p)] = 0

        # minus entropy of Q
        lb_q = (Qa-1.)*QlnE + (Qb-1.)*QlnEInv - special.betaln(Qa,Qb)
        lb_q[np.isnan(lb_q)] = 0

        return lb_p.sum() - lb_q.sum()

class Theta_Constant_Node(Constant_Variational_Node):
    """
    Dimensions of Theta_Constant_Node should be (D[m], K)
    """
    def __init__(self, dim, value):
        super(Theta_Constant_Node, self).__init__(dim, value)

    def precompute(self):
        self.E = self.value
        self.lnE = s.log(self.value)
        self.lnEInv = s.log(1.-self.value)

    def getExpectations(self):
        return { 'E':self.E, 'lnE':self.lnE, 'lnEInv':self.lnEInv }

    def removeFactors(self, idx, axis=1):
        # Ideally we want this node to use the removeFactors defined in Node()
        # but the problem is that we also need to update the "expectations", so i need
        # to call precompute()
        self.value = s.delete(self.value, idx, axis)
        self.precompute()
        self.updateDim(axis=axis, new_dim=self.dim[axis]-len(idx))

class Z_Node(UnivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None, qE2=None, idx_covariates=None):
        super(Z_Node,self).__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)

        self.covariates = np.zeros(self.dim[1], dtype=bool)

        # Define indices for covariates
        if idx_covariates is not None:
            self.covariates[idx_covariates] = True

    def precompute(self):
        # Precompute terms to speed up computation
        self.N = self.dim[0]
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
        # mask = [ ma.getmask(Y[m]) for m in range(len(Y)) ]
        mask = self.markov_blanket["Y"].getMask()

        # Collect parameters from the prior or expectations from the markov blanket
        Alpha = 1./self.P.getParameters()["var"]

        # Mask data
        for m in range(len(Y)):
            # Mask tau
            tau[m][mask[m]] = 0.
            # Mask Y
            Y[m] = Y[m].data
            Y[m][mask[m]] = 0.

        # Collect parameters from the P and Q distributions of this node
        Q = self.Q.getParameters().copy()
        Qmean, Qvar = Q['mean'], Q['var']

        M = len(Y)
        for k in latent_variables:
            foo = s.zeros((self.N,))
            bar = s.zeros((self.N,))
            for m in range(M):
                foo += np.dot(tau[m], SWtmp[m]["ESWW"][:,k])

                bar_tmp1 = SWtmp[m]["E"][:,k]

                bar_tmp2 = - s.dot(Qmean[:, s.arange(self.dim[1]) != k], SWtmp[m]["E"][:, s.arange(self.dim[1]) != k].T)
                bar_tmp2 += Y[m]
                bar_tmp2 *= tau[m]
                bar += np.dot(bar_tmp2, bar_tmp1)

            Qvar[:,k] = 1./(Alpha[:,k]+foo)
            Qmean[:,k] = Qvar[:,k] * bar

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean=Qmean, var=Qvar)

    def calculateELBO(self):
        # Collect parameters and expectations of current node
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
        Qmean, Qvar = Qpar['mean'], Qpar['var']
        QE, QE2 = Qexp['E'],Qexp['E2']
        Alpha = { 'E':1./self.P.getParameters()["var"], 'lnE':s.log(1./self.P.getParameters()["var"]) }

        # This ELBO term contains only cross entropy between Q and P,and entropy of Q. So the covariates should not intervene at all
        latent_variables = self.getLvIndex()
        Alpha["E"], Alpha["lnE"] = Alpha["E"][:,latent_variables], Alpha["lnE"][:,latent_variables]
        Qmean, Qvar = Qmean[:, latent_variables], Qvar[:, latent_variables]
        QE, QE2 = QE[:, latent_variables], QE2[:, latent_variables]

        # compute term from the exponential in the Gaussian
        QE2 *= 0.5
        tmp1 = -(QE2 * Alpha['E']).sum()

        # compute term from the precision factor in front of the Gaussian
        tmp2 = 0.5*Alpha["lnE"].sum()

        lb_p = tmp1 + tmp2
        lb_q = -0.5*(s.log(Qvar).sum() + self.N*len(latent_variables))

        return lb_p-lb_q
