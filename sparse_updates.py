from __future__ import division
import numpy.linalg  as linalg
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
        self.likconst = -0.5*s.sum(self.N)*s.log(2*s.pi)

    def mask(self):
        # Mask the observations if they have missing values
        self.value = ma.masked_invalid(self.value)

    def calculateELBO(self):
        tauQ_param = self.markov_blanket["Tau"].getParameters("Q")
        tauP_param = self.markov_blanket["Tau"].getParameters("P")
        tau_exp = self.markov_blanket["Tau"].getExpectations()
        lik = self.likconst + 0.5*s.sum(self.N*(tau_exp["lnE"])) - s.dot(tau_exp["E"],tauQ_param["b"]-tauP_param["b"])
        return lik


class Z_Node(UnivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None, qE2=None, idx_covariates=None):
        # UnivariateGaussian_Unobserved_Variational_Node.__init__(self, dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE)
        super(Z_Node,self).__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)
        self.precompute()

        if idx_covariates is not None:
            self.covariates[idx_covariates] = True

    def precompute(self):
        self.N = self.dim[0]
        self.covariates = np.zeros(self.dim[1], dtype=bool)
        self.factors_axis = 1

    def updateParameters(self):
        # TODO check what is needed from the new node (exp or param) and how
        # Collect expectations from other nodes
        # TO DO: MAKE THIS FASTER
        Y = self.markov_blanket["Y"].getExpectation()
        SWtmp = self.markov_blanket["SW"].getExpectations()
        tau = self.markov_blanket["Tau"].getExpectation()

        ClustPrior = self.markov_blanket['Cluster']
        Mu = ClustPrior.getExpectations()['E']

        M = len(Y)
        Y = ma.concatenate([Y[m] for m in xrange(M)],axis=1)
        SW = s.concatenate([SWtmp[m]["ESW"]for m in xrange(M)],axis=0)
        SWW = s.concatenate([SWtmp[m]["ESWW"] for m in xrange(M)],axis=0)
        tau = s.concatenate([tau[m] for m in xrange(M)],axis=0)

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Q = self.Q.getParameters()
        Pvar, Qmean = P['var'], Q['mean']
        Qmean =  Q['mean']

        # Variance
        # POSSIBLE MISTAKE: THE PLUS ONE HERE? OR IS THIS RELATED TO PVAR?
        tmp = (tau*SWW.T).sum(axis=1)
        tmp = s.repeat(tmp[None,:],self.N,0)
        tmp += 1./Pvar  # adding the prior precision to the updated precision
        Qvar = 1./tmp

        # Mean
        if any(self.covariates):
            covariates = Qmean[:,self.covariates]

        for k in xrange(self.dim[1]):
            tmp1 = SW[:,k]*tau
            tmp2 = Y - s.dot( Qmean[:,s.arange(self.dim[1])!=k] , SW[:,s.arange(self.dim[1])!=k].T )
            tmp3 = ma.dot(tmp2,tmp1)
            tmp3 += 1./Pvar[:,k] * Mu[:,k]
            Qmean[:,k] = Qvar[:,k] * tmp3

        # Do not update the latent variables associated with known covariates
        if any(self.covariates):
            Qmean[:,self.covariates] = covariates

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean=Qmean, var=Qvar)

    def calculateELBO(self):
        # TODO see what's left in the ELBO here and what should be moved to the new node
        # Collect parameters and expectations
        Ppar,Qpar,Qexp = self.P.getParameters(), self.Q.getParameters(), self.Q.getExpectations()
        Pvar, Qmean, Qvar = Ppar['var'], Qpar['mean'], Qpar['var']
        PE, PE2 = self.markov_blanket['Cluster'].getExpectations()['E'], self.markov_blanket['Cluster'].getExpectations()['E2']

        QE,QE2 = Qexp['E'],Qexp['E2']

        # compute term from the exponential in the Gaussian
        tmp1 = 0.5*QE2 - PE*QE + 0.5*PE2
        tmp1 = -(tmp1/Pvar).sum()

        # compute term from the precision factor in front of the Gaussian (TODO should be computed only once)
        tmp2 = - (s.log(Pvar)/2.).sum()

        lb_p = tmp1 + tmp2
        lb_q = - (s.log(Qvar).sum() + self.N*self.dim[1])/2.

        return lb_p-lb_q

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
        SW,SWW = tmp["ESW"], tmp["ESWW"]
        Ztmp = self.markov_blanket["Z"].getExpectations()
        Z,ZZ = Ztmp["E"],Ztmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # Calculate terms for the update
        term1 = (Y**2).sum(axis=0).data
        term2 = 2*(Y*s.dot(Z,SW.T)).sum(axis=0).data
        term3 = (ZZ.dot(SWW.T)).sum(axis=0)
        term4 = s.diag(s.dot( SW.dot(Z.T), Z.dot(SW.T) )) - s.dot(Z**2,(SW**2).T).sum(axis=0)
        tmp = term1 - term2 + term3 + term4

        # Perform updates of the Q distribution
        Qa = Pa + (Y.shape[0] - ma.getmask(Y).sum(axis=0))/2
        Qb = Pb + tmp/2

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.expectations['E'], self.Q.expectations['lnE']

        # Do the calculations
        lb_p = self.lbconst + s.sum((Pa-1)*QlnE) - s.sum(Pb*QE)
        lb_q = s.sum(Qa*s.log(Qb)) + s.sum((Qa-1)*QlnE) - s.sum(Qb*QE) - s.sum(special.gammaln(Qa))
        return lb_p - lb_q

class Alpha_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        # Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super(Alpha_Node,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        # self.K = self.dim[0]
        # self.lbconst = self.K * ( self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )
        self.lbconst = s.sum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) )
        self.factors_axis = 0

    def updateParameters(self):

        # Collect expectations from other nodes
        tmp = self.markov_blanket["SW"].getExpectations()
        ES,EWW,ESWW = tmp["ES"],tmp["EWW"],tmp["ESWW"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # ARD prior on ???
        # Qa = Pa + 0.5*EWW.shape[0]
        # Qb = Pb + 0.5*EWW.sum(axis=0)

        # ARD prior on ????
        # Qa = Pa + 0.5*ES.sum(axis=0)
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
        lb_p = self.lbconst + s.sum((Pa-1)*QlnE) - s.sum(Pb*QE)
        lb_q = s.sum(Qa*s.log(Qb)) + s.sum((Qa-1)*QlnE) - s.sum(Qb*QE) - s.sum(special.gammaln(Qa))

        return lb_p - lb_q

class SW_Node(BernoulliGaussian_Unobserved_Variational_Node):
    # SHOULD WE USE **KWARGS AND *KARGS ONLY?
    # def __init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0=None, qEW_S1=None, qES=None):
    def __init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0=None, qEW_S1=None, qES=None):
        super(SW_Node,self).__init__(dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0, qEW_S1, qES)
        # BernoulliGaussian_Unobserved_Variational_Node.__init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0, qEW_S1, qES)
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
        thetatmp = self.markov_blanket['Theta'].getExpectations() # TODO make general in mixed node
        theta_lnE, theta_lnEInv  = thetatmp['lnE'], thetatmp['lnEInv']
        # theta_E = thetatmp['E']

        # Collect parameters and expectations from P and Q distributions of this node
        SW = self.Q.getExpectations()["ESW"]
        Q = self.Q.getParameters()
        Qmean_S1, Qvar_S1, Qtheta = Q['mean_S1'], Q['var_S1'], Q['theta']

        # check dimensions of theta and expand if necessary
        # I THINK THIS SHOULD NOT BE HERE...
        if theta_lnE.shape != Qmean_S1.shape:
            theta_lnE = s.repeat(theta_lnE[None,:],Qmean_S1.shape[0],0)
        if theta_lnEInv.shape != Qmean_S1.shape:
            theta_lnEInv = s.repeat(theta_lnEInv[None,:],Qmean_S1.shape[0],0)

        all_term1 = theta_lnE - theta_lnEInv
        # all_term1 = s.log(theta_E/(1.-theta_E))

        for k in xrange(self.dim[1]):

            term1 = all_term1[:,k]
            term2 = 0.5*s.log(s.divide(alpha[k],tau))
            term3 = 0.5*s.log(s.sum(ZZ[:,k]) + s.divide(alpha[k],tau))
            term41 = ma.dot(Y.T,Z[:,k]).data
            term42 = s.dot( SW[:,s.arange(self.dim[1])!=k] , (Z[:,k]*Z[:,s.arange(self.dim[1])!=k].T).sum(axis=1) )
            term43 = s.sum(ZZ[:,k]) + s.divide(alpha[k],tau)
            term4 = 0.5*tau * s.divide((term41-term42)**2,term43)

            # Update S
            Qtheta[:,k] = 1/(1+s.exp(-(term1+term2-term3+term4)))

            # Update W
            Qmean_S1[:,k] = s.divide(term41-term42,term43)
            Qvar_S1[:,k] = s.divide(1,tau*term43)

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
        lb_qw = -0.5*self.dim[1]*self.D - 0.5*s.log(S*Qvar + ((1-S)/alpha["E"])).sum()
        lb_w = lb_pw - lb_qw

        # Calculate ELBO for S
        # Slower = 0.00001
        # Supper = 0.99999
        # S[S<Slower] = Slower
        # S[S>Supper] = Supper

        lb_ps = s.sum( S*theta['lnE'] + (1-S)*theta['lnEInv'])

        lb_qs_tmp = S*s.log(S) + (1-S)*s.log(1-S)
        lb_qs_tmp[s.isnan(lb_qs_tmp)] = 0

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

        # Collect expectations from other nodes
        S = self.markov_blanket['SW'].getExpectations()["ES"]


        # Precompute terms
        if factors_selection is not None:
            S = S[:,factors_selection]
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

        # D = self.markov_blanket['SW'].getDimensions()[0]

        # minus cross entropy of Q and P
        tmp1 = (Pa-1.)*QlnE + (Pb-1.)*QlnEInv - special.betaln(Pa,Pb)
        lb_p = tmp1.sum()

        # minus entropy of Q
        tmp2 = (Qa-1.)*QlnE + (Qb-1.)*QlnEInv - special.betaln(Qa,Qb)
        lb_q = tmp2.sum()

        return lb_p - lb_q

class Theta_Constant_Node(Constant_Variational_Node):
    def __init__(self, dim, value):
        super(Theta_Constant_Node, self).__init__(dim, value)
        self.precompute()

    def precompute(self):
        self.E = self.value
        self.lnE = s.log(self.value)
        self.lnEInv = s.log(1-self.value)

    def getExpectations(self):
        return { 'E':self.E, 'lnE':self.lnE, 'lnEInv':self.lnEInv }

    def removeFactors(self, idx, axis=0):
        # Ideally we want this node to use the removeFactors defined in Node()
        # but the problem is that we also need to update the "expectations", so i need
        # to call precompute()
        self.value = s.delete(self.value, idx, axis)
        self.precompute()
        self.updateDim(axis=axis, new_dim=self.dim[axis]-len(idx))

# we need a list of list of index to define the clusters (dictionary ? could also
# be a vector of length N_samples with cluster index)
# there should be a defalut when there is no clusters
# TODO do updatdes, but before that implement test to check 'pipes' are ok
class Cluster_Node_Gaussian(UnivariateGaussian_Unobserved_Variational_Node):
    # TODO need to implement droping a latent variable

    """ """
    def __init__(self, pmean, pvar, qmean, qvar, clusters, n_Z, cluster_dic=None, qE=None, qE2=None):
        # compute dim from numbers of clusters (n_clusters * Z)
        self.clusters = clusters
        self.n_clusters = len(np.unique(clusters))
        dim = (self.n_clusters, n_Z)
        super(Cluster_Node_Gaussian, self).__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)


    def getExpectations(self):
        # reshape the values to N_samples * N_factors and return
        QExp = self.Q.getExpectations()
        expanded_expectation = QExp['E'][self.clusters, :]
        expanded_E2 = QExp['E2'][self.clusters, :]
        # do we need to expand the variance as well ?
        return {'E': expanded_expectation , 'E2': expanded_E2}

    def updateParameters(self):

        Ppar = self.P.getParameters()
        ZQPar = self.markov_blanket['Z'].Q.getParameters()
        Qmean, Qvar = self.Q.getParameters()['mean'], self.Q.getParameters()['var']

        ZTau = 1./ZQPar['var']
        ZTauMean = ZQPar['mean']/ZQPar['var']

        # update of the variance
        for c in range(self.n_clusters):
            mask = (self.clusters == c)
            tmp = (ZTau[mask, :]).sum(axis=0)
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
        import pdb; pdb.set_trace()
        PParam = self.P.getParameters()
        PVar, Pmean = PParam['var'], PParam['mean']

        QExp = self.Q.getExpectations()
        QE2, QE = QExp['E2'], QExp['E']

        Qvar = self.Q.getParameters()['var']

        # minus cross entropy
        tmp = -(0.5 * s.log(PVar)).sum()
        tmp2 = - ((0.5/PVar) * (QE2 - 2.*QE*Pmean + Pmean**2.)).sum()

        # entropy of Q
        tmp3 = 0.5 * (s.log(Qvar)).sum()
        tmp3 += 0.5 * self.dim[0] * self.dim[1]

        return tmp + tmp2 + tmp3


# TODO do we need this ?
class Cluster_Node_Constant(Constant_Variational_Node):
    """ """
    def __init__(self,  dim, pmean, pvar, qmean, qvar, qE=None, qE2=None):
        super(Cluster_Node, self).__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)


    def updateParameters(self):

        pass

    def calculateELBO(self):
        return 0
