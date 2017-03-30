        ########
        # print "in tau"
        Y = self.markov_blanket["Y"].getExpectation().data
        SW = self.markov_blanket["SW"].getExpectations()
        Z = self.markov_blanket["Z"].getExpectations()

        N = Y.shape[0]
        D = Y.shape[1]
        K = Z["E"].shape[1]
        lik = s.zeros(D)
        for d in xrange(D):
            foo = 0
            for n in xrange(N):
                tmp = 0
                if not s.isnan(Y[n,d]):
                    for k in xrange(K):
                        tmp += SW["ESWW"][d,k]*Z["E2"][n,k]
                        for j in xrange(k+1,K):
                            tmp += 2*SW["ESW"][d,k]*Z["E"][n,k]*SW["ESW"][d,j]*Z["E"][n,j]
                    foo += 0.5*(Y[n,d]**2 - 2*Y[n,d]*s.dot(SW["ESW"][d,:],Z["E"][n,:]) + tmp)
            lik[d] = foo

        #######



                # Unvectorised
        # tau = self.markov_blanket["Tau"].getExpectations()
        # SW = self.markov_blanket["SW"].getExpectations()
        # Z = self.markov_blanket["Z"].getExpectations()
        # Y = self.value.data
        # N = self.dim[0]
        # D = self.dim[1]
        # K=SW["ES"].shape[1]

        # lik = 0
        # for n in xrange(N):
        #     for d in xrange(D):
        #         tmp = 0
        #         if not s.isnan(Y[n,d]):
        #             for k in xrange(K):
        #                 tmp += SW["ES"][d,k]*SW["EWW"][d,k]*Z["E2"][n,k]
        #                 for j in xrange(k+1,K):
        #                     tmp += 2*SW["ESW"][d,k]*Z["E"][n,k] * SW["ESW"][d,j]*Z["E"][n,j]
        #             lik += -0.5*s.log(2.*s.pi) + 0.5*tau["lnE"][d] - 0.5*tau["E"][d]*(Y[n,d]**2 - 2*Y[n,d]*s.dot(SW["ESW"][d,:],Z["E"][n,:]) + tmp)
        # print lik
        # exit()


        ## NON-VECTORISED ##


        tmp = self.markov_blanket["Z"].getExpectations()
        Z,ZZ = tmp["E"],tmp["E2"]
        tau = self.markov_blanket["Tau"].getExpectation()
        Y = self.markov_blanket["Y"].getExpectation()
        alpha = self.markov_blanket["Alpha"].getExpectation()

        SW = self.Q.getExpectations()["ESW"].copy()
        Q = self.Q.getParameters()
        Qmean_S1, Qvar_S1, Qtheta = Q['mean_S1'].copy(), Q['var_S1'].copy(), Q['theta'].copy()

        thetatmp = self.markov_blanket['Theta'].getExpectations() # TODO make general in mixed nodw
        theta_lnE, theta_lnEInv  = thetatmp['lnE'], thetatmp['lnEInv']  
        if theta_lnE.shape != Qmean_S1.shape:
            theta_lnE = s.repeat(theta_lnE[None,:],Qmean_S1.shape[0],0)
        if theta_lnEInv.shape != Qmean_S1.shape:
            theta_lnEInv = s.repeat(theta_lnEInv[None,:],Qmean_S1.shape[0],0)

        all_term1 = theta_lnE - theta_lnEInv

        ## Non-vectorised ##
        # oldmean = self.Q.mean.copy()
        # oldtheta = self.Q.theta.copy()
        # oldSW = SW.copy()
        for d in xrange(self.dim[0]):
            for k in xrange(self.dim[1]):
                term1 = all_term1[d,k]
                term2 = 0.5*s.log(alpha[k]/tau[d])
                term3 = 0.5*s.log(s.sum(ZZ[:,k]) + alpha[k]/tau[d])
                term4 = 0
                tmp1 = 0
                for n in xrange(Y.shape[0]):
                    if not s.isnan(Y.data[n,d]):
                        tmp1 += Y[n,d]*Z[n,k]
                tmp2 = 0
                for j in xrange(self.dim[1]):
                    if j!=k:
                        for n in xrange(Y.shape[0]):
                            if not s.isnan(Y.data[n,d]): 
                                tmp2 += SW[d,j]*Z[n,k]*Z[n,j]
                foo = tmp1 - tmp2
                bar = s.sum(ZZ[:,k]) + alpha[k]/tau[d]
                term4 = 0.5*tau[d]*foo**2 / bar

                # Update S
                Qtheta[d,k] = 1/(1+s.exp(-(term1+term2-term3+term4)))

                # Update W
                Qmean_S1[d,k] = foo/bar
                Qvar_S1[d,k] = 1/(tau[d]*bar)

                # Update expectations
                SW[d,k] = Qtheta[d,k] * Qmean_S1[d,k]