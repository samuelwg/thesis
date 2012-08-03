#!/usr/bin/env python
from __future__ import division
from spectral_components3d_common import *

# this file is to develop in a general way prod(1/(x+a+i), (i,0,n-1)) = sum(A_i/(x+a+i, (i,0,n-1)

def integerRootProductory(n,r,wo) :
	return (
		sum([
			(r-j)
			* N(-1)**(j)
			/ f(j)
			/ f(n-j)
			/(wo+j-r)
		for j in xrange(n+1) ])
		)
	# the original expression
	return (
		wo / prod([
			(wo+(-r+j))
			for j in xrange(n+1)
		])
	)


if __name__ == "__main__" :
	cases = reduce(operator.add, [[(n,r) for r in xrange(0,n+1) ] for n in xrange(0,6) ], [])
	print cases
	maz.cases("PartialFraction", cases)
	maz.check("PartialFraction", "initial, case %s", [
		wo / prod([
			(wo+(-r+j))
			for j in xrange(n+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "divide by wo, case %s", [
		prod([
			1/(wo+(-r+j))
			for j in xrange(n+1) if j!=r
		])
	for n,r in cases ])
	# TODO: Order zero does not match!!!
	maz.check("PartialFraction", "partial fraction, case %s", [
		sum([
			prod([
				1/N(r-j+(-r+i))
				for i in xrange(n+1) if i!=r and i!=j
			])
			/(wo+(-r+j))
			for j in xrange(n+1) if j!=r
		])
	for n,r in cases ])
	maz.check("PartialFraction", "removing r in the inner fractions, case %s", [
		sum([
			prod([
				1/N(i-j)
				for i in xrange(n+1) if i!=r and i!=j
			])
			/(wo+(-r+j))
			for j in xrange(n+1) if j!=r
		])
	for n,r in cases ])
	maz.check("PartialFraction", "split outer sum by r=j, case %s", [
		sum([ # j<r
			prod([
				1/N(i-j)
				for i in xrange(n+1) if i!=r and i!=j
			])
			/(wo+(-r+j))
			for j in xrange(r)
		])
		+
		sum([ # j>r
			prod([
				1/N(i-j)
				for i in xrange(n+1) if i!=r and i!=j
			])
			/(wo+(-r+j))
			for j in xrange(r+1,n+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "split inner prods by i=j, case %s", [
		sum([ # j<r
			prod([ # i<j
				1/N(i-j)
				for i in xrange(j) if i!=r
			])
			*
			prod([ # i>j
				1/N(i-j)
				for i in xrange(j+1,n+1) if i!=r
			])
			/(wo+(-r+j))
			for j in xrange(r)
		])
		+
		sum([ # j>r
			prod([ # i<j
				1/N(i-j)
				for i in xrange(j) if i!=r
			])
			*
			prod([ # i>j
				1/N(i-j)
				for i in xrange(j+1,n+1) if i!=r
			])
			/(wo+(-r+j))
			for j in xrange(r+1,n+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "restrictions on r have no sense if r not in i range, case %s", [
		sum([ # j<r
			prod([ # i<j
				1/N(i-j)
				for i in xrange(j)
			])
			*
			prod([ # i>j
				1/N(i-j)
				for i in xrange(j+1,n+1) if i!=r
			])
			/(wo+(-r+j))
			for j in xrange(r)
		])
		+
		sum([ # j>r
			prod([ # i<j
				1/N(i-j)
				for i in xrange(j) if i!=r
			])
			*
			prod([ # i>j
				1/N(i-j)
				for i in xrange(j+1,n+1)
			])
			/(wo+(-r+j))
			for j in xrange(r+1,n+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "splitting on i==r when in range, case %s", [
		sum([ # j<r
			prod([ # i<j
				1/N(i-j)
				for i in xrange(j)
			])
			*
			prod([ # i>j i<r
				1/N(i-j)
				for i in xrange(j+1,r)
			])
			*
			prod([ # i>j i>r
				1/N(i-j)
				for i in xrange(r+1,n+1)
			])
			/(wo+(-r+j))
			for j in xrange(r)
		])
		+
		sum([ # j>r
			prod([ # i<j
				1/N(i-j)
				for i in xrange(r)
			])
			*
			prod([ # i<j
				1/N(i-j)
				for i in xrange(r+1,j)
			])
			*
			prod([ # i>j
				1/N(i-j)
				for i in xrange(j+1,n+1)
			])
			/(wo+(-r+j))
			for j in xrange(r+1,n+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "compact code, case %s", [
		sum([ # j<r
			1 /(wo+(-r+j))
			* prod([ 1/N(i-j) for i in xrange(j) ]) # i<j
			* prod([ 1/N(i-j) for i in xrange(j+1,r) ]) # i>j i<r
			* prod([ 1/N(i-j) for i in xrange(r+1,n+1) ]) # i>j i>r
			for j in xrange(r)
		])
		+
		sum([ # j>r
			1/(wo+(-r+j))
			* prod([ 1/N(i-j) for i in xrange(r) ]) # i<j i<r
			* prod([ 1/N(i-j) for i in xrange(r+1,j) ]) # i<j i>r
			* prod([ 1/N(i-j) for i in xrange(j+1,n+1) ]) # i>j
			for j in xrange(r+1,n+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "i'=i-j, case %s", [
		sum([ # j<r
			1 /(wo+(-r+j))
			* prod([ 1/N(i) for i in xrange(-j,0) ]) # i<j
			* prod([ 1/N(i) for i in xrange(1,r-j) ]) # i>j i<r
			* prod([ 1/N(i) for i in xrange(r-j+1,n-j+1) ]) # i>j i>r
			for j in xrange(r)
		])
		+
		sum([ # j>r
			1/(wo+(-r+j))
			* prod([ 1/N(i) for i in xrange(-j,r-j) ]) # i<j i<r
			* prod([ 1/N(i) for i in xrange(r-j+1,0) ]) # i<j i>r
			* prod([ 1/N(i) for i in xrange(1,n-j+1) ]) # i>j
			for j in xrange(r+1,n+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "making all productories positive, case %s", [
		sum([ # j<r
			1 /(wo+(-r+j))
			* N(-1)**(j)
			* prod([ 1/N(i) for i in xrange(1,j+1) ]) # i<j
			* prod([ 1/N(i) for i in xrange(1,r-j) ]) # i>j i<r
			* prod([ 1/N(i) for i in xrange(r-j+1,n-j+1) ]) # i>j i>r
			for j in xrange(r)
		])
		+
		sum([ # j>r
			1/(wo+(-r+j))
			* N(-1)**N(j-r-1+r)
			* prod([ 1/N(i) for i in xrange(j-r+1, j+1) ]) # i<j i<r
			* prod([ 1/N(i) for i in xrange(1,-r+j) ]) # i<j i>r
			* prod([ 1/N(i) for i in xrange(1,n-j+1) ]) # i>j
			for j in xrange(r+1,n+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "simplify sign, case %s", [
		sum([ # j<r
			1 /(wo+(-r+j))
			* N(-1)**(j)
			* prod([ 1/N(i) for i in xrange(1,j+1) ]) # i<j
			* prod([ 1/N(i) for i in xrange(1,r-j) ]) # i>j i<r
			* prod([ 1/N(i) for i in xrange(r-j+1,n-j+1) ]) # i>j i>r
			for j in xrange(r)
		])
		-
		sum([ # j>r
			1/(wo+(-r+j))
			* N(-1)**N(j)
			* prod([ 1/N(i) for i in xrange(j-r+1, j+1) ]) # i<j i<r
			* prod([ 1/N(i) for i in xrange(1,-r+j) ]) # i<j i>r
			* prod([ 1/N(i) for i in xrange(1,n-j+1) ]) # i>j
			for j in xrange(r+1,n+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "all 1 based productories, case %s", [
		sum([ # j<r
			1 /(wo+(-r+j))
			* N(-1)**(j)
			* prod([ 1/N(i) for i in xrange(1,j+1) ]) # i<j
			* prod([ 1/N(i) for i in xrange(1,r-j) ]) # i>j i<r
			* prod([ 1/N(i) for i in xrange(1,n-j+1) ]) # i>j i>r
			/ prod([ 1/N(i) for i in xrange(1,r-j+1) ]) # i>j i>r
			for j in xrange(r)
		])
		-
		sum([ # j>r
			1/(wo+(-r+j))
			* N(-1)**N(j)
			* prod([ 1/N(i) for i in xrange(1,j+1) ]) # i<j i<r
			* prod([ 1/N(i) for i in xrange(1,-r+j) ]) # i<j i>r
			* prod([ 1/N(i) for i in xrange(1,n-j+1) ]) # i>j
			/ prod([ 1/N(i) for i in xrange(1,j-r+1) ]) # i<j i<r
			for j in xrange(r+1,n+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "turn into factorials, case %s", [
		sum([ # j<r
			1 /(wo+(-r+j))
			* N(-1)**(j)
			/ f(j)
			/ f(r-j-1)
			/ f(n-j)
			* f(r-j)
			for j in xrange(r)
		])
		-
		sum([ # j>r
			1/(wo+(-r+j))
			* N(-1)**N(j)
			/ f(j)
			/ f(j-r-1)
			/ f(n-j)
			* f(j-r)
			for j in xrange(r+1,n+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "j'=-r+j, case %s", [
		sum([ # j<r
			1 /(wo+(-r+(j+r)))
			* N(-1)**((j+r))
			/ f((j+r))
			/ f(r-(j+r)-1)
			/ f(n-(j+r))
			* f(r-(j+r))
			for j in xrange(-r,0)
		])
		-
		sum([ # j>r
			1/(wo+(-r+(j+r)))
			* N(-1)**N((j+r))
			/ f((j+r))
			/ f((j+r)-r-1)
			/ f(n-(j+r))
			* f((j+r)-r)
			for j in xrange(1,n-r+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "simplify, case %s", [
		sum([
				1 /(wo+j)
				* N(-1)**(j+r)
				/ f(j+r)
				/ f(-j-1)
				/ f(n-j-r)
				* f(-j)
			if j<0
			else
				1/(wo+j)
				* N(-1)**N(j+r+1)
				/ f(j+r)
				/ f(j-1)
				/ f(n-j-r)
				* f(j)
			for j in xrange(-r,n-r+1) if j
		])
	for n,r in cases ])
	maz.check("PartialFraction", "extracting common factors, case %s", [
		sum([
				1 /(wo+j)
				* N(-1)**(j+r)
				/ f(j+r)
				/ f(n-j-r)
				* (
					1
					* f(-j)
					/ f(-j-1)
				if j<0
				else
					-1
					* f(j)
					/ f(j-1)
				)
			for j in xrange(-r,n-r+1) if j
		])
	for n,r in cases ])
	maz.check("PartialFraction", "remaining part is always -j, case %s", [
		sum([
				1 /(wo+j)
				* N(-1)**(j+r)
				/ f(j+r)
				/ f(n-j-r)
				* -j
			for j in xrange(-r,n-r+1) if j
		])
	for n,r in cases ])
	maz.check("PartialFraction", "multiplying by j allows removing condition, case %s", [
		sum([
				1 /(wo+j)
				* N(-1)**(j+r)
				/ f(j+r)
				/ f(n-j-r)
				* -j
			for j in xrange(-r,n-r+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "multiplying by j allows removing condition, case %s", [
		sum([
				-j
				* N(-1)**(j+r)
				/ f(j+r)
				/ f(n-j-r)
				/(wo+j)
			for j in xrange(-r,n-r+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "j'=j+r, case %s", [
		sum([
				(r-j)
				* N(-1)**(j)
				/ f(j)
				/ f(n-j)
				/(wo+j-r)
			for j in xrange(n+1)
		])
	for n,r in cases ])
	maz.check("PartialFraction", "initial, case %s", [
		wo / prod([
			(wo+(-r+j))
			for j in xrange(n+1)
		])
	for n,r in cases ])






