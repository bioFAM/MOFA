from variational_nodes import Variational_Node
from nodes import Constant_Node
import scipy as s

import pdb

# Comment: Only differences between Mixed_Nodes and Multiview_Mixed_Node nodes are that :
#   1- nodes in mixed nodes are NOT all updated
#   2- some nodes dont have an expectation but a value instead
#   3- Mixed_Nodes return a matrix as an expectation, as the concatenation of
#      all nodes is done inside and not outside to provide the same interface as
#      a normal node

class Mixed_Theta_Nodes(Variational_Node, Constant_Node):
    """
    Class for a mixture of LearningTheta and constant ConstantTheta nodes.
    For a total of K factors, some (Klearn) will learn Theta whereas for the others (Kconst) it will be constant
        K = Klearn + Kconst
    """
    def __init__(self, LearnTheta, ConstTheta, idx):
        # Inputs:
        # - LearnTheta: Theta_Node with dimensions (Klearn,)
        # - ConstTheta: Theta_Constant_Node with dimensions (D,Kconst) or (Kconst,1) - NOT IMPLEMENTED YET -
        # - idx: list or numpy array indicating which factors are LearnTheta(idx=1. or idx=True) and which are ConstTheta(idx=0. or idx=False)
        self.constTheta = ConstTheta
        self.learnTheta = LearnTheta

        self.K = ConstTheta.dim[1] + LearnTheta.dim[0]
        self.D = ConstTheta.dim[0]

        self.idx = idx
        
    def addMarkovBlanket(self, **kargs):
        # SHOULD WE ALSO ADD MARKOV BLANKET FOR CONSTHTETA???
        self.learnTheta.addMarkovBlanket(**kargs)

    def getExpectations(self):

        # Get expectations from ConstTheta nodes (D,Kconst)
        Econst = self.constTheta.getExpectations().copy()

        # Get expectations from LearnTheta nodes and expand to (D,Kconst)
        Elearn = self.learnTheta.getExpectations().copy()
        Elearn["E"] = s.repeat(Elearn["E"][None,:], self.D, 0)
        Elearn["lnE"] = s.repeat(Elearn["lnE"][None,:], self.D, 0)
        Elearn["lnEInv"] = s.repeat(Elearn["lnEInv"][None,:], self.D, 0)

        # Concatenate expectations to (D,K)
        E = s.concatenate((Econst["E"], Elearn["E"]), axis=1)
        lnE = s.concatenate((Econst["lnE"], Elearn["lnE"]), axis=1)
        lnEInv = s.concatenate((Econst["lnEInv"], Elearn["lnEInv"]), axis=1)

        # Permute to the right order given by self.idx
        idx = s.concatenate((s.nonzero(1-self.idx)[0],s.where(self.idx)[0]), axis=0)
        E, lnE, lnEinv = E[:,idx], lnE[:,idx], lnEInv[:,idx]
        return dict({'E': E, 'lnE': lnE, 'lnEInv':lnEInv})

    def getExpectation(self):
        return self.getExpectations()['E']

    def updateExpectations(self):
        self.learnTheta.updateExpectations()

    def updateParameters(self):
        # the argument contains the indices of the non_annotated factors
        self.learnTheta.updateParameters(s.nonzero(self.idx)[0])

    def calculateELBO(self):
        return self.learnTheta.calculateELBO()

    def removeFactors(self, *idx):
        for i in idx:
            if self.idx[idx] == 1:
                self.learnTheta.removeFactors(s.where(i == s.nonzero(self.idx)[0])[0])
            else:
                self.constTheta.removeFactors(s.where(i == s.nonzero(1-self.idx)[0])[0])
            self.idx = self.idx[s.arange(self.K)!=i]
            self.K -= 1

