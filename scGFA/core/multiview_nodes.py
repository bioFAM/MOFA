import scipy as s

from nodes import Node
from variational_nodes import Variational_Node

"""
Module to define multi-view nodes in a Bayesian network

All multiview nodes (either variational or not variational) have two main attributes:
- M: total number of views
- activeM: in some occasions a particular node is active in only a subset of the M views. 
    This variable keeps track of which views contain the node.
- nodes: the set of (single-view) nodes
"""

class Multiview_Node(Node):
    """
    General class for multiview nodes in a Bayesian network
    """
    def __init__(self, M, *nodes):
        # - M: number of views
        # - nodes: list of M instances (or children) of the 'Node' class or None if the node is not defined in a particular view
        self.M = M

        # Some nodes are active in only a subset of views, this variables keeps track of these views.
        self.activeM = [ m for m,node in enumerate(nodes) if node is not None]

        # Initialise nodes
        self.nodes = nodes

    def removeFactors(self,idx):
        # Method to remove factors from the node
        #   - idx (ndarray): indices of the factors to be removed
        for m in self.activeM: self.nodes[m].removeFactors(idx)

    def getNodes(self):
        # Method to get the nodes
        return self.nodes
        
    def getExpectation(self):
        # Method to get the first moments (expectation)
        return [ self.nodes[m].getExpectation() for m in self.activeM ]

    def getExpectations(self):
        # Method to get all relevant moments
        return [ self.nodes[m].getExpectations() for m in self.activeM ]

    def getParameters(self):
        # Method to get  the parameters
        return [ self.nodes[m].getParameters() for m in self.activeM ]

    def updateDim(self, axis, new_dim, m=None):
        # Method to update the dimensionality of the node
        #  axis (int)
        #  new_dim (int)
        #  m (iterable): views to update
        assert s.all(m in self.activeM), "Trying to update the dimensionality of a node that doesnt exist in a view"
        M = self.activeM if m is None else m
        for m in M: self.nodes[m].updateDim(axis,new_dim)

class Multiview_Variational_Node(Multiview_Node, Variational_Node):
    """
    General class for multiview variational nodes.
    """
    def __init__(self, M, *nodes):
        # nodes: list of M 'Node' instances
        Multiview_Node.__init__(self, M, *nodes)
        for node in nodes: assert isinstance(node, Variational_Node)

    def update(self):
        # Method to update both parameters and expectations
        for m in self.activeM:
            self.nodes[m].updateParameters()
            self.nodes[m].updateExpectations()
    def updateExpectations(self):
        # Method to update expectations using current estimates of the parameters
        for m in self.activeM: self.nodes[m].updateExpectations()
    def updateParameters(self):
        # Method to update parameters using current estimates of the expectations
        for m in self.activeM: self.nodes[m].updateParameters()
    def calculateELBO(self):
        # Method to calculate variational evidence lower bound
        lb = [ self.nodes[m].calculateELBO() for m in self.activeM ]
        return sum(lb)

class Multiview_Constant_Node(Multiview_Node):
    """
    General class for multiview local nodes
    """
    def __init__(self, M, *nodes):
        # nodes: list of M 'Node' instances
        Multiview_Node.__init__(self, M, *nodes)

    def getValues(self):
        # Method to retun the constant values
        return [ self.nodes[m].getValue() for m in self.activeM ]

class Multiview_Mixed_Node(Multiview_Constant_Node, Multiview_Variational_Node):
    """
    General Class for multiview nodes that contain both variational and constant nodes
    """
    def __init__(self, M, *nodes):
        # M: number of views
        # nodes: list of M 'Node' instances
        Multiview_Node.__init__(self, M, *nodes)

    def update(self):
        # Method to update values of the nodes
        for m in self.activeM: self.nodes[m].update()

    def calculateELBO(self):
        # Method to calculate variational evidence lower bound
        # The lower bound of a multiview node is the sum of the lower bound of its 
        # corresponding single view variational nodes
        lb = 0
        for m in self.activeM:
            if isinstance(self.nodes[m],Variational_Node):
                lb += self.nodes[m].calculateELBO()
        return lb

    # def getValues(self):
        # AT SOME POINT WE HAVE TO REMOVE THIS FUNCTION AND REPLACE IT BY VALUE
        # return [ self.nodes[m].getValue() for m in self.activeM ]
