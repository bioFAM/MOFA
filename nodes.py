

class Node(object):
    """ 
    General class for a node in a Bayesian network
    """
    def __init__(self, dim):
    	self.dim = dim
    	pass
    def addMarkovBlanket(self, **kwargs):
    	# Function to define the Markov blanket of the node 
        self.markov_blanket = kwargs

# class Observed_Node(Node):
#     """ 
#     General class for an observed node in a Bayesian network
#     """
#     def __init__(self, dim, obs):
#     	Node.__init__(self,dim)
#     	self.obs = obs

# class Multiview_Node(Node):
# 	""" 
# 	General class for multiview nodes in a Bayesian network
# 	"""
# 	def __init__(self, M, *nodes):
# 		# dim: list of M tuples with the dimensionality of the node in each view, ex. [ (10,5), (20,5), ...]
# 		# nodes: list of M 'Node' instances
# 		# Node.__init__(self, dim)

# 		# Initialise dimensionalities
# 		# self.M = len(self.dim)
# 		self.M = M

# 		# Initialise nodes
# 		for node in nodes: assert isinstance(node,Node), "Nodes have to be instances of the general Node class"
# 		self.nodes = nodes