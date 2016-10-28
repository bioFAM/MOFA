    def initZ(dim, type="random", Y=None):
        # Function to initialise the latent variables
        # - dim (dic): dictionary with the dimensionalities (N,K,D,M)
        #  - type (str): random, orthogonal, pca
        #  - Y (list of np arrays with dim (N,Dm)): observed data, only for pca

        if type == "random":
            Z_qmean = s.stats.norm.rvs(loc=0, scale=1, size=(dim["N"],dim["K"]))
        elif type == "orthogonal":
            print "Not implemented"
            exit()
        elif type == "pca":
            print "Not implemented"
            exit()

        Z_qcov = s.repeat(s.eye(dim["K"])[None,:,:],dim["N"],0)
        Z = Z_Node(dim=(dim["N"],dim["K"]), qmean=Z_qmean, qcov=Z_qcov)
        # Z.updateExpectations()
        return Z

    def initAlpha(dim, pa=1e-14, pb=1e-14, qb=1., qE=1.)
        # Function to initialise the latent variables
        # - dim (dic): dictionary with the dimensionalities (N,K,D,M)
        alpha_list = [None]*dim["M"]
        qb = s.ones(dim["K"])*qb
        qE = s.ones(dim["K"])*qE
        for m in xrange(dim["M"]):
            qa = pa + s.ones(dim["K"])*dim["D"][m]/2
            alpha_list[m] = Alpha_Node(dim=(dim["K"],), pa=pa, pb=pb, qa=qa, qb=s.ones(dim["K"])*qb, qE=s.ones(dim["K"])*qE)
        alpha = Multiview_Variational_Node((dim["K"],)*dim["M"], *alpha_list)
        return alpha

    def initW(dim, type):
        # Function to initialise the weights
        #  - type (str): random, pca
        W_list = [None]*dim["M"]
        for m in xrange(dim["M"]):
            if type == "random":
                W_qmean = s.stats.norm.rvs(loc=0, scale=1, size=(dim["D"][m],dim["K"]))
            elif type == "pca":
                print "Not implemented"
                exit()
            W_qcov = s.repeat(a=s.eye(dim["K"])[None,:,:], repeats=dim["D"][m] ,axis=0)
            W_list[m] = W_Node(dim=(dim["D"][m],dim["K"]), qmean=W_qmean, qcov=W_qcov, qE=W_qmean)
        W = Multiview_Variational_Node(dim["M"], *W_list)
        return W


    def initTau(dim, lik, pa=1e-14, pb=1e-14, qb=0., qE=100.):
        # Function to initialise the precision of the noise
        # - dim (dic): dictionary with the dimensionalities (N,K,D,M)
        # lik (list): likelihood of each view (poisson,gaussian,bernoulli,binomial)
        tau_list = [None]*dim["M"]
        for m in xrange(dim["M"]):
            if lik[m] == "poisson"
                tmp = 0.25 + 0.17*s.amax(data["Y"][m],axis=0) 
                tau_list[m] = Observed_Local_Node(dim=(dim["D"][m],), value=tmp)
            elif lik[m] == "bernoulli"
                tmp = s.ones(dim["D"][m])*0.25 
                tau_list[m] = Observed_Local_Node(dim=(dim["D"][m],), value=tmp)
            elif lik[m] == "binomial"
                tmp = 0.25*s.amax(data["Y"]["tot"][m],axis=0)
                tau_list[m] = Observed_Local_Node(dim=(dim["D"][m],), value=tmp)
            elif lik[m] == "gaussian":
                qa = pa + s.ones(dim["D"][m])*dim["N"]/2
                qb = s.zeros(dim["D"][m]) + qb
                qE = s.zeros(dim["D"][m]) + qE
                tau_list[m] = Tau_Node(dim=(dim["D"][m],), pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        tau = Multiview_Mixed_Node(dim["M"],*tau_list)

        return tau

    def initZeta(dim,lik):
        # Function to initialise the local variable zeta of the seeger approach
        # not initialised since it is the first update
        Zeta_list = [None]*dim["M"]
        for m in xrange(dim["M"]):
            if lik[m] is not "gaussian":
                Zeta_list[m] = Zeta_Node(dim=(dim["N"],dim["D"][m]), initial_value=None) 
            else:
                Zeta_list[m] = None
        Zeta = Multiview_Local_Node(dim["M"], *Zeta_list)
        return Zeta

    def initY(dim,lik):
        # Function to initialise the observed data

        # Y/Yhat (mixed node)
        Y_list = [None]*dim["M"]
        for m in xrange(dim["M"]):
            if m in M_gaussian:
                Y_list[m] = Y_Node(dim=(dim["N"],dim["D"][m]), obs=data['Y'][m])
            elif m in M_poisson:
                Y_list[m] = Poisson_PseudoY_Node(dim=(dim["N"],dim["D"][m]), obs=data['Y'][m], E=None)
            elif m in M_bernoulli:
                Y_list[m] = Bernoulli_PseudoY_Node(dim=(dim["N"],dim["D"][m]), obs=data['Y'][m], E=None)
            elif m in M_binomial:
                Y_list[m] = Binomial_PseudoY_Node(dim=(dim["N"],dim["D"][m]), tot=data['Y']["tot"][m], obs=data['Y']["obs"][m], E=None)
        Y = Multiview_Mixed_Node(dim["M"], *Y_list)
        return Y