
# Class to iniailise the traditional GFA model
# Main differences: no element-wise sparsity and latent variables do not fully factorise
# PROBABLY THIS DOESNT WORK ANYMORE
class init_nonsparse(initModel):
    def __init__(self, *args, **kwargs):
        super(init_nonsparse, self).__init__(*args, **kwargs)

    def initZ(self, pmean, pcov, qmean, qcov, qE=None, covariates=None):
        # Method to initialise the latent variables
        # covariates (nd array): matrix of covariates with dimensions (nsamples,ncovariates)

        # Initialise first moment
        if qmean is not None:
            if isinstance(qmean,str):
                if qmean == "random":
                    qmean = stats.norm.rvs(loc=0, scale=1, size=(self.N,self.K))
                elif qmean == "orthogonal":
                    print "Not implemented"
                    exit()
                elif qmean == "pca":
                    pca = sklearn.decomposition.PCA(n_components=self.K, copy=True, whiten=True)
                    tmp = s.concatenate(self.data,axis=0).T
                    pca.fit(tmp)
                    qmean = pca.components_.T

            elif isinstance(qmean,s.ndarray):
                assert qmean.shape == (self.N,self.K)

            elif isinstance(qmean,(int,float)):
                qmean = s.ones((self.N,self.K)) * qmean

            else:
                print "Wrong initialisation for Z"
                exit()

        # Add covariates
        if covariates is not None:
            # idx_covariates = s.arange(self.K)[-covariates.shape[1]:]
            # putting the covariates first instead
            # import pdb; pdb.set_trace()
            idx_covariates = s.array(range(covariates.shape[1]))

            qmean[:,idx_covariates] = covariates
            qvar = s.ones((self.N,self.K))*qvar
            qvar[:, idx_covariates] = 0.
            # qmean = s.c_[ qmean, covariates ]
            # idx_covariates = s.arange(covariates.shape[1]) + self.K
        else:
            idx_covariates = None

        # Initialise the node
        self.Z = Z_Node(dim=(self.N,self.K),
                        pmean=pmean,
                        pcov=pcov,
                        qmean=qmean,
                        qcov=qcov,
                        qE=qE,
                        idx_covariates=idx_covariates)
        self.nodes["Z"] = self.Z

    def initW(self, pmean, pcov, qmean, qcov, qE=None):
        # Method to initialise weights
        W_list = [None]*self.M
        for m in xrange(self.M):

           # Initialise first moment
            if isinstance(qmean[m],str):
                if qmean[m] == "random":
                    qmean[m] = stats.norm.rvs(loc=0, scale=1, size=(self.D[m],self.K))
                else:
                    print "%s initialisation not implemented for W" % qmean[m]
                    exit()
            elif isinstance(qmean[m],s.ndarray):
                assert qmean[m].shape == (self.D[m],self.K), "Wrong dimensionality"
            elif isinstance(qmean[m],(int,float)):
                qmean[m] = s.ones((self.N,self.K)) * qmean[m]
            else:
                print "Wrong initialisation for W"
                exit()

            # Initalise the node
            W_list[m] = W_Node(
                dim=(self.D[m],self.K),
                pmean=pmean[m], pcov=pcov[m],
                qmean=qmean[m], qcov=qcov[m], qE=qE[m]
                )

        self.W = Multiview_Variational_Node(self.M, *W_list)
        self.nodes["W"] = self.W

    def initAlpha(self, pa, pb, qa, qb, qE):
        # Method to initialise the precision of the group-wise ARD prior
        # Inputs:
        #  pa (float): 'a' parameter of the prior distribution
        #  pb (float): 'b' parameter of the prior distribution
        #  qb (float): initialisation of the 'b' parameter of the variational distribution
        #  qE (float): initial expectation of the variational distribution
        alpha_list = [None]*self.M
        for m in xrange(self.M):
            alpha_list[m] = Alpha_Node(dim=(self.K,), pa=pa[m], pb=pb[m], qa=qa[m], qb=qb[m], qE=qE[m])
        self.Alpha = Multiview_Variational_Node((self.K,)*self.M, *alpha_list)
        self.nodes["Alpha"] = self.Alpha

    def initTau(self, pa, pb, qa, qb, qE):
        # Method to initialise the precision of the noise
        # Inputs:
        #  pa (float): 'a' parameter of the prior distribution
        #  pb (float): 'b' parameter of the prior distribution
        #  qb (float): initialisation of the 'b' parameter of the variational distribution
        #  qE (float): initial expectation of the variational distribution
        tau_list = [None]*self.M
        for m in xrange(self.M):
            if self.lik[m] == "poisson":
                tmp = 0.25 + 0.17*s.amax(self.data[m],axis=0)
                tau_list[m] = Constant_Node(dim=(self.D[m],), value=tmp)
            elif self.lik[m] == "bernoulli":
                tmp = s.ones(self.D[m])*0.25
                tau_list[m] = Constant_Node(dim=(self.D[m],), value=tmp)
            elif self.lik[m] == "binomial":
                tmp = 0.25*s.amax(self.data["tot"][m],axis=0)
                tau_list[m] = Constant_Node(dim=(self.D[m],), value=tmp)
            elif self.lik[m] == "gaussian":
                tau_list[m] = Tau_Node(dim=(self.D[m],), pa=pa[m], pb=pb[m], qa=qa[m], qb=qb[m], qE=qE[m])
        self.Tau = Multiview_Mixed_Node(self.M,*tau_list)
        self.nodes["Tau"] = self.Tau

    def initY(self):
        # Method to initialise the observed data
        Y_list = [None]*self.M
        for m in xrange(self.M):
            if self.lik[m]=="gaussian":
                Y_list[m] = Y_Node(dim=(self.N,self.D[m]), value=self.data[m])
            elif self.lik[m]=="poisson":
                Y_list[m] = Poisson_PseudoY_Node(dim=(self.N,self.D[m]), obs=self.data[m], E=None)
            elif self.lik[m]=="bernoulli":
                Y_list[m] = Bernoulli_PseudoY_Node(dim=(self.N,self.D[m]), obs=self.data[m], E=None)
            elif self.lik[m]=="binomial":
                Y_list[m] = Binomial_PseudoY_Node(dim=(self.N,self.D[m]), tot=data["tot"][m], obs=data["obs"][m], E=None)
        self.Y = Multiview_Mixed_Node(self.M, *Y_list)
        self.nodes["Y"] = self.Y


    def initExpectations(self, *nodes):
        # Method to initialise the expectations of some nodes
        for node in nodes:
            self.nodes[node].updateExpectations()

    def MarkovBlanket(self):
        # Method to define the markov blanket
        self.Z.addMarkovBlanket(W=self.W, Tau=self.Tau, Y=self.Y, )
        for m in xrange(self.M):
            self.Alpha.nodes[m].addMarkovBlanket(W=self.W.nodes[m])
            self.W.nodes[m].addMarkovBlanket(Z=self.Z, Tau=self.Tau.nodes[m], Alpha=self.Alpha.nodes[m], Y=self.Y.nodes[m])
            if self.lik[m]=="gaussian":
                self.Y.nodes[m].addMarkovBlanket(Z=self.Z, W=self.W.nodes[m], Tau=self.Tau.nodes[m])
                self.Tau.nodes[m].addMarkovBlanket(W=self.W.nodes[m], Z=self.Z, Y=self.Y.nodes[m])
            else:
                self.Y.nodes[m].addMarkovBlanket(Z=self.Z, W=self.W.nodes[m], kappa=self.Tau.nodes[m])
