import scipy as s

from nodes import Node
from local_nodes import Local_Node
from variational_nodes import Variational_Node
from constant_nodes import Constant_Node

"""
Module to define multi-view nodes in a Bayesian network

All multiview nodes (either variational or not variational) have two main attributes:
- M: total number of views
- idx: in some occasions a particular node is active in only a subset of the M views. This variable keeps track of which views contain the node.
"""

class Multiview_Node(Node):
    """
    General class for multiview nodes in a Bayesian network
    """
    def __init__(self, M, *nodes):
        # Input:
        # - M: number of views
        # - nodes: list of M instances (or children) of the 'Node' class or None if the node is not defined in a particular view
        self.M = M

        # Some nodes are active in only a subset of views, this variables keeps track of these views.
        self.idx = [ idx for idx,node in enumerate(nodes) if node is not None]

        # Initialise nodes
        self.nodes = nodes

    def removeFactors(self,*idx):
        # Remove factors/latent variables
        for m in self.idx: self.nodes[m].removeFactors(idx)

    def getExpectation(self):
        # Get the current estimate of the first moment
        return [ self.nodes[m].getExpectation() for m in self.idx ]

    def getExpectations(self):
        # Get current estimate of all relevant moments
        return [ self.nodes[m].getExpectations() for m in self.idx ]

    def getParameters(self):
        # Get current estimate of parameters
        return [ self.nodes[m].getParameters() for m in self.idx ]

# class Multiview_Local_Node(Multiview_Node,Local_Node):
class Multiview_Local_Node(Multiview_Node):
    """
    General class for multiview local nodes
    """
    def __init__(self, M, *nodes):
        # nodes: list of M 'Node' instances
        Multiview_Node.__init__(self, M, *nodes)

    def update(self):
        # Update parameters
        for m in self.idx: self.nodes[m].update()
        pass
    def getValue(self):
        # Get the current estimates of the values
        return [ self.nodes[m].getValue() for m in self.idx ]
    # def getExpectation(self):
        # Get the current estimates of the values
        # return self.getValue()

class Multiview_Variational_Node(Multiview_Node, Variational_Node):
    """
    General class for multiview variational nodes.
    """
    def __init__(self, M, *nodes):
        # nodes: list of M 'Node' instances
        Multiview_Node.__init__(self, M, *nodes)
        for node in nodes: assert isinstance(node, Variational_Node)

    def update(self):
        # Update both parameters and expectations
        for m in self.idx:
            self.nodes[m].updateParameters()
            self.nodes[m].updateExpectations()
        pass
    def updateExpectations(self):
        # Update expectations using current estimates of the parameters
        for m in self.idx: self.nodes[m].updateExpectations()
        pass
    def updateParameters(self):
        # Update parameters using current estimates of the expectations
        for m in self.idx: self.nodes[m].updateParameters()
        pass
    def calculateELBO(self):
        # Calculate variational evidence lower bound
        lb = [ self.nodes[m].calculateELBO() for m in self.idx ]
        return sum(lb)

class Multiview_Mixed_Node(Multiview_Local_Node, Multiview_Variational_Node):
    """
    General Class for multiview nodes that contain both variational and local nodes
    """
    def __init__(self, M, *nodes):
        # nodes: list of M 'Node' instances
        Multiview_Node.__init__(self, M, *nodes)

    def update(self):
        # Update values of the node
        for m in self.idx: self.nodes[m].update()
        pass

    def calculateELBO(self):
        # Calculate variational evidence lower bound
        lb = 0
        for m in self.idx:
            if isinstance(self.nodes[m],Variational_Node):
                lb += self.nodes[m].calculateELBO()
        return lb

    def getObservations(self):
        return [ self.nodes[m].getObservations() for m in self.idx ]
