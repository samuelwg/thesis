#!/usr/bin/env python
from __future__ import division
from sympy import exp, I, integrate, special, symbols, simplify, factorial, cos, sin, binomial
import sympy
import operator
import sys
from parameters import *
import maz

t0=R/c
sympy.var("w t x", real=True)
N=sympy.Integer

def sum(v) : return reduce(operator.add, v, N(0))
def prod(v) : return reduce(operator.mul, v, N(1))
def f(n) : return reduce(operator.mul, xrange(n,0,-1), N(1))
def ff(n) : return reduce(operator.mul, xrange(n,0,-2), N(1))

"""
"""


maz.check("Polynomial", "initial (%i)", [
	sum([
		ff(2*n-2*m-1)/ff(2*m-1)*x**(2*m-n-1)
		for m in xrange(1,n//2+1) ])
for n in xrange(10) ])
maz.check("Polynomial", "m'=m+1 (%i)", [
	sum([
		ff(2*n-2*m-3)/ff(2*m+1)*x**(2*m-n+1)
	for m in xrange(n//2) ])
for n in xrange(10) ])
maz.check("Polynomial", "m'=2*m (%i)", [
	sum([
		ff(2*n-m-3)/ff(m+1)*x**(m-n+1)
	for m in xrange(0,2*(n//2),2) ])
for n in xrange(10) ])
maz.check("Polynomial", "removing the // (%i)", [
	sum([
		ff(2*n-m-3)/ff(m+1)*x**(m-n+1)
	for m in xrange(0,n-1,2) ])
for n in xrange(10) ])
maz.check("Polynomial", "m'=m-1 // (%i)", [
	sum([
		ff(2*n-2-m)/ff(m)*x**(m-n)
	for m in xrange(1,n,2) ])
for n in xrange(10) ])

maz.pprint("Polynomial")



