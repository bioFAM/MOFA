        
###################################################
## Vectorised updates for SW without general tau ##
###################################################

# for k in xrange(self.dim[1]):

#     # Calculate intermediate steps
#     term1 = (theta_lnE-theta_lnEInv)[:,k]
#     term2 = 0.5*s.log(s.divide(alpha[k],tau))
#     term3 = 0.5*s.log(s.sum(ZZ[:,k]) + s.divide(alpha[k],tau))
#     term41 = ma.dot(Y.T,Z[:,k]).data
#     tmp = s.dot((Z[:,k]*Z[:,s.arange(self.dim[1])!=k].T).T, SW[:,s.arange(self.dim[1])!=k].T)
#     term42 = ma.array(tmp, mask=ma.getmask(Y)).sum(axis=0)
#     # term42 = s.dot( SW[:,s.arange(self.dim[1])!=k] , (Z[:,k]*Z[:,s.arange(self.dim[1])!=k].T).sum(axis=1) ) # THIS IS WRONG WITH NAs
#     term43 = s.sum(ZZ[:,k]) + s.divide(alpha[k],tau)
#     term4 = 0.5*tau * s.divide((term41-term42)**2,term43)

#     # Update S
#     # NOTE there could be some precision issues in S --> loads of 1s in result
#     Qtheta[:,k] = 1/(1+s.exp(-(term1+term2-term3+term4)))

#     # Update W
#     Qmean_S1[:,k] = s.divide(term41-term42,term43)
#     Qvar_S1[:,k] = s.divide(1,tau*term43)

#     # Update Expectations for the next iteration
#     SW[:,k] = Qtheta[:,k] * Qmean_S1[:,k]