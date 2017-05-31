# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 11:41:30 2015

@author: Rukawa
"""

import numpy

class Identity(object):
    def __init__(self,x):
        self.x = x
        
    def f(self):
        return self.x
    def f_prime(self):
        # return 1.0
        return numpy.ones((1,))
    def __str__(self):
        return ''

        
class Log_exp(object):
    def __init__(self,x):
        self.x = x
    def f(self):
        return numpy.log(1. + numpy.exp(self.x))
    def f_prime(self):
        return 1. - numpy.exp(-self.x)

class Exp(object):
    def __init__(self,x):
        self.x = x
    def f(self):
        return numpy.exp(self.x)
    def f_prime(self):
        return numpy.exp(self.x)
