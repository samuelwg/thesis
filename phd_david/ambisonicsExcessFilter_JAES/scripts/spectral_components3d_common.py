#!/usr/bin/env python
from __future__ import division
from sympy import exp, I, integrate, special, symbols, simplify, factorial, cos, sin, binomial
import sympy
import operator
import sys
from parameters import *
import maz

sympy.var("w t", real=True) # angular velocity and time
y,z=symbols("y z".split(), positive=True)
x=symbols("x")
to,=symbols(["to"], positive=True) # to=R/c
wo,=symbols(["wo"], positive=True) # y= -I*w*to = -I*wo; wo = I*y
print wo



N=sympy.Integer

def R(n) : return reduce(operator.mul, xrange(n,1,-2),N(1))
def sum(v) : return reduce(operator.add, v, N(0))
def prod(v) : return reduce(operator.mul, v, N(1))
def f(n) : return reduce(operator.mul, xrange(n,0,-1), N(1))

w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12 = [ w**i for i in xrange(1,12+1) ]
t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12 = [ t**i for i in xrange(1,12+1) ]
# rn = n!!
r1, r3, r5, r7, r9, r11, r13, r15, r17, r19, r21 = [R(i) for i in xrange(1,22,2) ]
r0, r2, r4, r6, r8, r10, r12, r14, r16, r18, r20 = [R(i) for i in xrange(0,21,2) ]
# bn = 2^n
b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18 = [ sympy.Integer(2)**(i) for i in xrange(18+1) ]
# fn = factorial(n)
f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12 = [ f(i) for i in xrange(12+1) ]
# cn = cos(n*t)
c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12 = [ cos(N(i)*t) for i in xrange(12+1) ]
# sn = sin(n*t)
s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12 = [ sin(N(i)*t) for i in xrange(12+1) ]
# zn = w^2-n^2
z0, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12 = [ (w**2-N(i)**2) for i in xrange(12+1) ]

BN = binomial
Pl = special.polynomials.legendre



