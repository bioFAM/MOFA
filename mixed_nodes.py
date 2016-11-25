from variational_nodes import Variational_Node
from constant_nodes import Constant_Node
import scipy as s

import pdb

# Comment: Only differences between Mixed_Nodes and Multiview_Mixed_Node nodes are that :
#   1- nodes in mixed nodes are NOT all updated
#   2- some nodes dont have an expectation but a value instead
#   3- Mixed_Nodes return a matrix as an expectation, as the concatenation of
#      all nodes is done inside and not outside to provide the same interface as
#      a normal node

# when NOT-constant, might not be same dimensions -> to expand somewhere
class Mixed_Theta_Nodes(Variational_Node, Constant_Node):
    """
    Class for a mixture of annotated Theta and non annotated Theta

    Annotated_theta is a Constant_Node of dimension k_1
    Non_Annotated_Theta is a Theta_Node_No_Annotation of diemnsion K - k_1

    When updating the non_annotated part of a Mixed_Theta_Nodes, we need to pass
    as an argument to the Theta_Node_No_Annotation the indices of the factors which
    are not annotated
    """
    def __init__(self, Annotated_Theta, Non_Annotated_Theta):
        self.annotated_theta = Annotated_Theta  # Constant_Node
        self.non_annotated_theta = Non_Annotated_Theta  # Beta_Node

        self.K = Annotated_Theta.dim[1] + Non_Annotated_Theta.Q.dim[0]

        self.annotated_factors_ix = range(0, Annotated_Theta.dim[1])
        self.non_annotated_factors_ix = range(Annotated_Theta.dim[1], self.K)

    def addMarkovBlanket(self, **kargs):
        self.non_annotated_theta.addMarkovBlanket(**kargs)

    def getExpectations(self):
        # get expectations or values for each node and return concatenated array
        values_annotated = self.annotated_theta.getValue()
        exp = self.non_annotated_theta.getExpectations()['E']

        # deal with different dimensions
        pdb.set_trace()  # TODO check dimensions here
        if exp.shape != values_annotated.shape:
            exp = s.repeat(exp[None, :], values_annotated.shape[0], 0)

        tmp = s.concatenate((values_annotated, exp), axis=1)
        return dict({'E': tmp})

    def updateExpectations(self):
        self.non_annotated_theta.updateExpectations()

    def updateParameters(self):
        # the argument contains the indices of the non_annotated factors
        self.non_annotated_theta.updateParameters(self.non_annotated_factors_ix)

    # TODO implement
    def calculateELBO(self):
        return self.non_annotated_theta.calculateELBO()

    def removeFactors(self, *ix):
        keep = s.setdiff1d(s.arange(self.K),idx)

        keep_annotated = s.intersect1d(keep, self.annotated_factors_ix)
        keep_non_annotated = s.intersect1d(keep, self.non_annotated_factors_ix)

        self.annotated_theta.removeFactors(keep_annotated)
        self.non_annotated_theta.removeFactors(keep_non_annotated)

        self.annotated_factors_ix = xrange(len(keep_annotated))
        self.non_annotated_factors_ix = xrange(len(keep_non_annotated)) + len(self.annotated_factors_ix)

        self.K = len(self.annotated_factors_ix) + len(self.non_annotated_factors_ix)
