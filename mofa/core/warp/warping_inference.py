# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 14:41:10 2015

@author: Rukawa
"""


import numpy
import scipy.optimize

from .warping_functions import Warping_functions


class Warping_inference(object):
    
    """
    Super class for the warping functions
    z = f(x)
    """

    def __init__(self,func_type='tanh',I=3):
        self.options = dict()
        self.options['maxiter'] = 50
        self.options['maxfun'] = 50
        self.options['disp'] = False # Whether to print covergence message.
        self.method = 'L-BFGS-B'
#        self.method = 'TNC'               
        
        self.entity = Warping_functions(func_type,I)
        
    def f(self,x,i_not=-1):
        """
        Transform x using all warping terms except the i^th term
        If i is -1, all warping terms are used. 
        """
        res = 0.0
        for l in range(self.entity.I):
            if l!=i_not: 
                res += self.entity.f(x,l)
        return x + res


    def f_prime(self,x,i_not=-1):
        """
        Gradient of f w.r.t. to x of all warping terms except the i^th term
        """
        res = 0.0
        for l in range(self.entity.I):
            if l!=i_not: 
                res += self.entity.f_prime(x,l)
        return 1. + res


    def update_parameters(self,X,F,tau,mat_mask): 
        """
        Update parameters
        X ~ N(X|F,tau) where F = UV' + ...

        NOTE: PARTIAL_F_VAL AND PRIME_VAL NOT SURE IF ITS IN THE LOOP OR NOT
        """
        for i in range(self.entity.I):
            for warp_var in self.entity.param.keys():
                partial_f_val = self.f(X,i_not=i)
                partial_f_prime_val = self.f_prime(X,i_not=i)
                self._optimise_parameter(warp_var,i,X,F,tau,mat_mask,partial_f_val,partial_f_prime_val) 
                
    def _optimise_parameter(self,warp_var,i,X,F,tau,mat_mask,partial_f_val,partial_f_prime_val):
        a_i_old = self.entity.param['a'][i].x
        b_i_old = self.entity.param['b'][i].x
        c_i_old = self.entity.param['c'][i].x
        
        res = scipy.optimize.minimize(self._f_ELBO, jac=self._f_ELBO_d_param, x0=numpy.random.random(), method=self.method,args=(warp_var,i,X,F,tau,mat_mask,partial_f_val,partial_f_prime_val),options=self.options)

        if res.success and (res.nit > 0) and (res.jac < 1E-4):
            pass
            # print warp_var + ' = ' + str(self.entity.param[warp_var][i].x)
            # print self._f_ELBO_d_param(self.entity.param['a'][i].x, "a", i, X, F, tau, mat_mask,partial_f_val,partial_f_prime_val)
        else:
            if warp_var == 'a':
                self.entity.param['a'][i].x = a_i_old                
            elif warp_var == 'b':   
                self.entity.param['b'][i].x = b_i_old            
            elif warp_var == 'c':   
                self.entity.param['c'][i].x = c_i_old
            # print warp_var


    def _f_ELBO(self,param,param_opt,i,X,F,tau,mat_mask,partial_f_val,partial_f_prime_val):

        if param_opt == 'a':
            # self.entity.param['a'][i].x = param
            self.entity.param['a'][i].x = param[0]
        elif param_opt == 'b':
            # self.entity.param['b'][i].x = param
            self.entity.param['b'][i].x = param[0]
        elif param_opt == 'c':
            # self.entity.param['c'][i].x = param
            self.entity.param['c'][i].x = param[0]

        Z = partial_f_val + self.entity.f(X,i)
        ans = numpy.log(partial_f_prime_val + self.entity.f_prime(X,i)) 
        ans -= 0.5*tau * ((Z - F)**2)
        func = numpy.sum(ans)
        # assert func == numpy.sum( mat_mask * ans)
        return -func

    def _f_ELBO_d_param(self,param,param_opt,i,X,F,tau,mat_mask,partial_f_val,partial_f_prime_val):

        if param_opt == 'a':
            # self.entity.param['a'][i].x = param
            self.entity.param['a'][i].x = param[0]
            param_class = self.entity.param['a'][i]
        elif param_opt == 'b':
            # self.entity.param['b'][i].x = param
            self.entity.param['b'][i].x = param[0]
            param_class = self.entity.param['b'][i]
        elif param_opt == 'c':
            # self.entity.param['c'][i].x = param
            self.entity.param['c'][i].x = param[0]
            param_class = self.entity.param['c'][i]

        inv_f_prime_val = 1.0/(partial_f_prime_val + self.entity.f_prime(X,i))        
        Z = partial_f_val + self.entity.f(X,i)
        tau_Z_F = tau * (Z-F)            
        
        grad = inv_f_prime_val * self.entity.f_prime(X,i,partial=param_opt) 
        grad -= tau_Z_F * self.entity.f(X,i,partial=param_opt)       
        grad *= param_class.f_prime()

        res = numpy.array([(-1.0) * numpy.mean(grad)])
        # assert res == numpy.array([(-1.0) * numpy.mean(mat_mask*grad)])

        return res


    def pdf(self,d,g,mu_z,var_z,space):
#        sd_z = numpy.sqrt(1.00/tau_z)
        sd_z = numpy.sqrt(var_z)
        z = numpy.linspace(mu_z-3*sd_z,mu_z+3*sd_z,1000)
        pdf = scipy.stats.norm.pdf(z,loc=mu_z,scale=sd_z)
        if space == 'observation':
            x = self.f_inv(z)
            pdf *= self.f_prime(x) 
            return x,pdf
        else:
            return z,pdf


    def f_inv(self,z,x=None,rate = 0.1, max_iters=50):
        """
        Inverse function transformation
        """
        if x is None:
            x = numpy.zeros_like(z).astype(float) # How to choose starting points !
        update = numpy.inf
        t = 0
        while t == 0 or (numpy.abs(update).sum() > (1e-5) and t < max_iters):
            update = rate*(self.f(x) - z) / self.f_prime(x)
            x -= update
            t += 1 
        if t == max_iters:
            print("WARNING: maximum number of iterations reached in f_inv ")
        return x

    def f_inv_gauss_hermite(self,mu_z,tau_z,deg=10):
        """
        Mean of the inverse function transformation under a Gaussian density
        """
        x,w = numpy.polynomial.hermite.hermgauss(deg)
        sqrt_2_over_tau_z = numpy.sqrt(2.0/tau_z)
        mu_x = numpy.zeros_like(mu_z).astype(float)
        for t in range(deg):
            new_arg = x[t]*sqrt_2_over_tau_z + mu_z 
            mu_x += w[t]*self.f_inv(new_arg)
        return (1.0/numpy.sqrt(numpy.pi)) * mu_x


    # def plot(self, xmin, xmax):
    #     """
    #     Plot the warping function
    #     """
    #     y = numpy.arange(xmin, xmax, 0.01)
    #     f_y = self.f(y)
    #     plt.figure()
    #     plt.plot(y, f_y)
    #     plt.xlabel('x')
    #     plt.ylabel('f(x)')
        # plt.title('Warping function')
        
    