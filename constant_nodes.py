from nodes import Node
import scipy as s

class Constant_Node(Node):
    """

    """
    def __init__(self, dim, value):
        self.dim = dim
        if isinstance(value,(int,float)):
            self.value = value * s.ones(dim)
        else:
            assert value.shape == dim, "dimensionality mismatch"
            self.value = value

    def getValue(self):
        return self.value

    def removeFactors(self, *idx):
        if len(self.dim) == 1:
            keep = s.setdiff1d(s.arange(self.dim[0]),idx)
            self.value = self.value[keep]
            self.dim = (len(self.value),)
        else:
            keep = s.setdiff1d(s.arange(self.dim[1]),idx)
            self.value = self.value[:, keep]
            self.dim = self.value.shape
