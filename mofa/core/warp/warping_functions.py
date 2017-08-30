# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 14:59:50 2016

@author: Rukawa
"""
import numpy

from . import factor_functions

class Warping_functions(object):
    ''' Type of the warping functions e.g. tanh. '''
    def __init__(self,func_type,I):
        self.I = I
        self.factor_ft_type = dict()
        self.factor_ft_type['a'] = factor_functions.Log_exp
        self.factor_ft_type['b'] = factor_functions.Log_exp
        self.factor_ft_type['c'] = factor_functions.Identity
        self.param = dict()
        self.param['a'] = [[]] * self.I
        self.param['b'] = [[]] * self.I
        self.param['c'] = [[]] * self.I
        for i in range(self.I):
            self.param['a'][i] = self.factor_ft_type['a'](numpy.random.random())
            self.param['b'][i] = self.factor_ft_type['b'](numpy.random.random())
            self.param['c'][i] = self.factor_ft_type['c'](numpy.random.random())
        if func_type == 'tanh':
            self.func_type = _Tanh
            self.func_type_str = 'tanh'
        elif func_type == 'logistic':
            self.func_type = _Sigmoid
            self.func_type_str = 'logistic'

    def f(self,x,i,partial = None):
        linear = _Linear.f(x,self.param['b'][i].f(),self.param['c'][i].f())
        if partial is None:
            return self.param['a'][i].f() * self.func_type.f(linear)
        elif partial == 'data':
            return self.param['a'][i].f() * self.param['b'][i].f() * self.func_type.f_prime(linear)
        elif partial == 'a':
            return self.func_type.f(linear) 
        elif partial == 'b':
            return self.param['a'][i].f() * self.func_type.f_prime(linear) * _Linear.f(x,self.param['b'][i].f(),self.param['c'][i].f(),partial = 'b')
        elif partial == 'c':
            return self.param['a'][i].f() * self.func_type.f_prime(linear) * _Linear.f(x,self.param['b'][i].f(),self.param['c'][i].f(),partial = 'c')
    
    
    def f_prime(self,x,i,partial = None):
        linear = _Linear.f(x,self.param['b'][i].f(),self.param['c'][i].f())
        if partial is None:
            return self.f(x,i,partial = 'data')
        elif partial == 'a':
            return self.param['b'][i].f() * self.func_type.f_prime(linear)
        elif partial == 'b':
            return self.param['a'][i].f() * self.param['b'][i].f() * self.func_type.f_double_prime(linear) * _Linear.f(x,self.param['b'][i].f(),self.param['c'][i].f(),partial = 'b') + self.param['a'][i].f() * self.func_type.f_prime(linear)
        elif partial == 'c':
            return self.param['a'][i].f() * self.param['b'][i].f() * self.func_type.f_double_prime(linear) * _Linear.f(x,self.param['b'][i].f(),self.param['c'][i].f(),partial = 'c')

    
class _Linear(object):
    @staticmethod
    def f(x,b,c,partial=None):
        if partial is None:
            return b * (x + c)
#            return x
        elif partial == 'data':         
            return b
        elif partial == 'b':         
            return x + c
        elif partial == 'c':         
            return b

class _Tanh(object):
    @staticmethod
    def f(x):
        return numpy.tanh(x)
    @staticmethod
    def f_prime(x):
        return 1. - _Tanh.f(x)**2        
    @staticmethod
    def f_double_prime(x):
        return - 2 * _Tanh.f(x) * _Tanh.f_prime(x)

class _Sigmoid(object):
    @staticmethod
    def f(x):
        return 1./(1. + numpy.exp(-1.*x))
    @staticmethod
    def f_prime(x):
        return _Sigmoid.f(x) - _Sigmoid.f(x)**2        
    @staticmethod
    def f_double_prime(x):
        return _Sigmoid.f_prime(x) - 2*_Sigmoid.f(x)*_Sigmoid.f_prime(x)
        
        

