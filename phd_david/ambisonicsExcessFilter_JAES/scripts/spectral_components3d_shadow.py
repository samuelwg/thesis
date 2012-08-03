#!/usr/bin/env python
from __future__ import division
from spectral_components3d_common import *
from spectral_components3d_Pli import Pli

"""
Context (taken from spectral_components3d_full):
The problem finding the spectral domain expression of Hl,0,+:
sqrt(2*l+1)/2/to * ( # TODO: maybe some 1/2pi is needed for fourier but this makes H0(0) be 0dB
	integrate(Pl(t/to) * exp(-I*w*t), (t,-to,0)) +
	integrate(cos(t/to) * Pl(sin(t/to)) * exp(-I*w*t), (t,0,to*pi/2)) 
	)
	Substitution:
		t = to * z
		dt = to * dz
		z0=-1 z1=0
		z0=0, z1= pi/2
sqrt(2*l+1)/2 * (
	+ integrate(Pl(z)*exp(-1j*w*to*z), (z,-1, 0))             # direct part
	+ integrate(cos(z)*Pl(sin(z))*exp(-1j*w*to*z), (z,0,pi/2) # shadowed part
	)
	Substitution:
		y = -1j*w*to
sqrt(2*l+1)/2 * (
	+ integrate(Pl(z)*exp(y*z), (z,-1, 0))             # direct part
	+ integrate(cos(z)*Pl(sin(z))*exp(y*z), (z,0,pi/2) # shadowed part
	)
"""


"""
################################################################
# Shadowed part integral
################################################################
Table[Integrate[Cos[t] LegendreP[l, Sin[t]] Exp[-I w t], t] /Exp(-I t w) , {l, 0, 10}]
integrate(cos(z)*Pl(sin(z))*exp(y*z), (z,0,pi/2))
	Substitute:
	x= sin(z), dx = cos(z) dz
	x0 = 0, x1 = 1
= integrate(Pl(x)*exp(y*asin(x)), (x,0,1))
	g(x) = exp(y*asin(x))
= integrate( Pl(x) g(x), x )
	Ii=integrate ith on x
	di=ith derivative (d/dx)**i
= I1( Pl(x) g(x))
	# By parts: int(udv) +int(vdu) = uv
	dv = g(x) dx
	v = I1(g(x))
	u = Pl(x)
	du = d1(Pl(x)) dx
= I1(g(x))*Pl(x) - I1( I1(g(x)) d1(Pl(x)) )
	# By parts:
	v = d1(Pl(x))
	dv = d2(Pl(x)) dx
	du = I1(g(x)) dx
	u = I2(g(x))
= I1(g(x))*Pl(x) - I2(g(x))d1(Pl(x)) + I1( I2(g(x),x) d2(Pl(x)))
= ...
	# Pl is a polynomial just has l non-zero derivatives
= sum ( [ (-1)**i  * I(i+1)(g(x)) * di(Pl(x)) for i in xrange(l+1) ] )
"""

if __name__ == "__main__" :
	"""
	# Go for Ii(g(x)) the ith integral of g on x

	# But x is not convenient variable now we need to undo the variable subst on z
		What does it means to integrate on x expresed on z
			sin(z)=x
			dx=cos(z)dz
		integrate(f(x), x) 
		= integrate(f(sin(z))*cos(z),z)
		We need to multiply by cos(z) before each integral on z
		Recall that: cos(z) = (exp(+I*z)+exp(-I*z))/2
	# g(x) = exp(y*asin(x))
	# I1(g(x))
		= integrate(exp(y*asin(x)),x)
			z=asin(x)
		= integrate(exp(y*z)*cos(z),z)
		= integrate(exp(y*z)*(exp(I*z)+exp(-I*z))/2,z)
	"""

	""" # checked but commented out for speed

	cosz = (exp(I*z)+exp(-I*z))/N(2)
	maz.check("PreIntegral_1_g_x", "Multiply g by cos(z) (%i)",
	[
		exp(y*z)*cosz
	])
	maz.check("PreIntegral_1_g_x", "cos expressed in exponential (%i)",
	[
		exp(y*z)*(exp(I*z)+exp(-I*z))/2
	])
	maz.check("Integral_1_g_x","initial taking the cos multiplied formula (%i)",
	[
		integrate(maz.recall("PreIntegral_1_g_x")[0],z)
	])
	maz.check("Integral_1_g_x","splitting in two integrals 1+%i",
	[
		+ integrate(exp(y*z)*exp(+I*z),z)/2
		+ integrate(exp(y*z)*exp(-I*z),z)/2 
	])
	maz.check("Integral_1_g_x","integrated (%i)",
	[
		exp(z*(y+I))/2/(y+I) + exp(z*(y-I))/2/(y-I)
	])
	maz.check("Integral_1_g_x","1/2 factor (%i)",
	[
		1/2 * (exp(z*(y+I))/(y+I) + exp(z*(y-I))/(y-I))
	])
	maz.check("Integral_1_g_x","rewritten a n*I (%i)",
	[
		1/2 * (
			+ exp(z*(y+1*I))/(y+1*I)
			+ exp(z*(y-1*I))/(y-1*I)
		)
	])
	maz.check("Integral_1_g_x","exp(z*y) factor extract (%i)",
	[
		1/2 * exp(z*y) * (
			+ exp(+1*I*z)/(y+1*I)
			+ exp(-1*I*z)/(y-1*I)
		)
	])
	maz.check("Integral_1_g_x","fraction sum (%i)",
	[
		1/2 * exp(z*y) * (
			+ (y+1*I) *exp(-1*I*z)
			+ (y-1*I) *exp(+1*I*z)
		)/(y*y+1)
	])
	maz.check("Integral_1_g_x","moved (y^2+1) denominator (%i)",
	[
		1/2/(y*y+1) * exp(z*y) * (
			+ (y+1*I) *exp(-1*I*z)
			+ (y-1*I) *exp(+1*I*z)
		)
	])
	I1s9 = maz.recall("Integral_1_g_x")[0]

	# I2

	maz.check("Pre-Integral_2_g_x", "I1*cos(z) (%i)",
	[
		I1s9 * (exp(+I*z)+exp(-I*z)) /2
	])
	maz.check("Pre-Integral_2_g_x", "Substitute I1 (%i)",
	[
		1/2 * (
			+ exp(z*(y+1*I))/(y+1*I)
			+ exp(z*(y-1*I))/(y-1*I)
		) * (exp(+I*z)+exp(-I*z)) /2
	])
	maz.check("Pre-Integral_2_g_x", "Distribute cos(z) (%i)",
	[
		1/4 * (
			+ exp(z*(y+1*I))/(y+1*I) * (exp(+I*z)+exp(-I*z))
			+ exp(z*(y-1*I))/(y-1*I) * (exp(+I*z)+exp(-I*z))
		)
	])
	maz.check("Pre-Integral_2_g_x", "Separate terms (%i)",
	[
		1/4 * (
			+ exp(z*(y+1*I))/(y+1*I) * (exp(+I*z))
			+ exp(z*(y-1*I))/(y-1*I) * (exp(+I*z))
			+ exp(z*(y+1*I))/(y+1*I) * (exp(-I*z))
			+ exp(z*(y-1*I))/(y-1*I) * (exp(-I*z))
		)
	])
	maz.check("Pre-Integral_2_g_x", "Separate terms (%i)",
	[
		1/4 * (
			+ exp(z*(y+2*I))/(y+1*I)
			+ exp(z*(y-0*I))/(y-1*I)
			+ exp(z*(y+0*I))/(y+1*I)
			+ exp(z*(y-2*I))/(y-1*I)
		)
	])
	maz.check("Pre-Integral_2_g_x", "extre (%i)",
	[
		1/4 /(y*y+1) * (
			+ (y-1*I)*exp(z*(y+2*I))
			+ (y+1*I)*exp(z*(y-0*I))
			+ (y-1*I)*exp(z*(y+0*I))
			+ (y+1*I)*exp(z*(y-2*I))
		)
	])
	maz.check("Pre-Integral_2_g_x", "extre (%i)",
	[
		1/4 /(y*y+1) * (
			+ (y-1*I)*exp(z*y)*exp(+2*I*z)
			+ (y+1*I)*exp(z*y)*exp(-0*I*z)
			+ (y-1*I)*exp(z*y)*exp(+0*I*z)
			+ (y+1*I)*exp(z*y)*exp(-2*I*z)
		)
	])

	maz.check("Integral_2_g_x", "initial (%i)",
	[
		integrate(maz.recall("Pre-Integral_2_g_x")[0],z)
	])
	maz.check("Integral_2_g_x", "separate term integrals (%i)",
	[
		1/4 /(y*y+1) * (
			+ (y-1*I)*integrate(exp(z*y)*exp(+2*I*z),z)
			+ (y+1*I)*integrate(exp(z*y)*exp(-0*I*z),z)
			+ (y-1*I)*integrate(exp(z*y)*exp(+0*I*z),z)
			+ (y+1*I)*integrate(exp(z*y)*exp(-2*I*z),z)
		)
	])
	maz.check("Integral_2_g_x", "integrate (%i)",
	[
		1/4 /(y*y+1) * (
			+ (y-1*I)*exp(z*y)*exp(+2*I*z)/(y+2*I)
			+ (y+1*I)*exp(z*y)*exp(-0*I*z)/(y-0*I)
			+ (y-1*I)*exp(z*y)*exp(+0*I*z)/(y+0*I)
			+ (y+1*I)*exp(z*y)*exp(-2*I*z)/(y-2*I)
		)
	])
	maz.check("Integral_2_g_x", "extract exp(zy-2zj) common factors (%i)",
	[
		1/4 /(y*y+1) * exp(z*y) * exp(-2*I*z) * (
			+ (y-1*I)*exp(+4*I*z)/(y+2*I)
			+ (y+1*I)*exp(+2*I*z)/(y-0*I)
			+ (y-1*I)*exp(+2*I*z)/(y+0*I)
			+ (y+1*I)*exp(+0*I*z)/(y-2*I)
		)
	])
	maz.check("Integral_2_g_x", "reorder terms (%i)",
	[
		1/4 /(y*y+1) * exp(z*y) * exp(-2*I*z) * (
			+ (y-1*I)*exp(+4*I*z)/(y+2*I)
			+ (y+1*I)*exp(+2*I*z)/(y-0*I)
			+ (y-1*I)*exp(+2*I*z)/(y+0*I)
			+ (y+1*I)*exp(+0*I*z)/(y-2*I)
		)
	])
	maz.check("Integral_2_g_x", "extract 1/(y^2+0) factor (%i)",
	[
		1/4 /(y-0*I)/(y*y+1) * exp(z*y) * exp(-2*I*z) * (
			+ (y-1*I)*(y+0*I)*exp(+4*I*z)/(y+2*I)
			+ (y+1*I)*exp(+2*I*z)
			+ (y-1*I)*exp(+2*I*z)
			+ (y+1*I)*(y+0*I)*exp(+0*I*z)/(y-2*I)
		)
	])
	maz.check("Integral_2_g_x", "extract 1/(y^2+4) factor (%i)",
	[
		1/4 /(y+0)/(y*y+1)/(y*y+4) * exp(z*y) * exp(-2*I*z) * (
			+ (y-1*I)*(y+0*I)*(y-2*I)*exp(+4*I*z)
			+ (y+1*I)*(y+2*I)*(y-2*I)*exp(+2*I*z)
			+ (y-1*I)*(y+2*I)*(y-2*I)*exp(+2*I*z)
			+ (y+1*I)*(y+2*I)*(y+0*I)*exp(+0*I*z)
		)
	])
	maz.check("Integral_2_g_x", "join exps with same sign (%i)",
	[
		1/4 /y/(y*y+1)/(y*y+4) * exp(z*y) * exp(-2*I*z) * (
			+ 1*(y-0*I)*(y-1*I)*(y-2*I)*exp(+4*I*z)
			+ 2*(y+2*I)*(y-0*I)*(y-2*I)*exp(+2*I*z)
			+ 1*(y-0*I)*(y+1*I)*(y+2*I)*exp(+0*I*z)
		)
	])
	maz.check("Integral_2_g_x", "extract recover y based factors (%i)",
	[
		1/4 * exp(z*y) * exp(-2*I*z) * (
			+ 1 * exp(0*I*z) / (y-1*I) / (y-2*I)
			+ 2 * exp(2*I*z) / (y-1*I) / (y+1*I)
			+ 1 * exp(4*I*z) / (y+1*I) / (y+2*I)
			)
	])

	I2s8	= maz.recall("Integral_2_g_x")[0]

	# I3 

	I3p0	= I2s8 * (exp(+I*z)+exp(-I*z)) /2
	I3p1	= 1/8 * exp(z*y) * exp(-2*I*z) * (
			+ 1 * exp(0*I*z) / (y-1*I) / (y-2*I)
			+ 2 * exp(2*I*z) / (y-1*I) / (y+1*I)
			+ 1 * exp(4*I*z) / (y+1*I) / (y+2*I)
		) * (exp(+I*z)+exp(-I*z))
	assertEquivalent(I3p0, I3p1,"Pre Integral 3, step 1")


	I3p2	= 1/8 * exp(z*y) * exp(-2*I*z) * (
			+ 1 * exp(0*I*z) / (y-1*I) / (y-2*I) * (exp(+I*z)+exp(-I*z))
			+ 2 * exp(2*I*z) / (y-1*I) / (y+1*I) * (exp(+I*z)+exp(-I*z))
			+ 1 * exp(4*I*z) / (y+1*I) / (y+2*I) * (exp(+I*z)+exp(-I*z))
		)
	assertEquivalent(I3p1, I3p2,"Pre Integral 3, step 2")
	I3p3	= 1/8 * exp(z*y) * exp(-3*I*z) * (
			+ 1 * exp(2*I*z) / (y-1*I) / (y-2*I)
			+ 2 * exp(4*I*z) / (y-1*I) / (y+1*I)
			+ 1 * exp(6*I*z) / (y+1*I) / (y+2*I)
			+ 1 * exp(0*I*z) / (y-1*I) / (y-2*I)
			+ 2 * exp(2*I*z) / (y-1*I) / (y+1*I)
			+ 1 * exp(4*I*z) / (y+1*I) / (y+2*I)
		)
	assertEquivalent(I3p1, I3p2,"Pre Integral 3, step 3")
	I3p4	= 1/8 * exp(z*y) * exp(-3*I*z) * (
			+ 1 * exp(6*I*z) / (y+1*I) / (y+2*I)
			+ 3 * exp(4*I*z) / (y-1*I) / (y+2*I)
			+ 3 * exp(2*I*z) / (y+1*I) / (y-2*I)
			+ 1 * exp(0*I*z) / (y-1*I) / (y-2*I)
		)
	assertEquivalent(I3p3, I3p4,"Pre Integral 3, step 4")
	I3p5	= 1/8 * (
			+ 1 * exp(z*y) * exp(-3*I*z) * exp(6*I*z) / (y+1*I) / (y+2*I)
			+ 3 * exp(z*y) * exp(-3*I*z) * exp(4*I*z) / (y-1*I) / (y+2*I)
			+ 3 * exp(z*y) * exp(-3*I*z) * exp(2*I*z) / (y+1*I) / (y-2*I)
			+ 1 * exp(z*y) * exp(-3*I*z) * exp(0*I*z) / (y-1*I) / (y-2*I)
		)
	assertEquivalent(I3p4, I3p5,"Pre Integral 3, step 5")

	I3s0	= integrate(I3p5, z)

	I3s1	= 1/8 * (
			+ 1 * integrate(exp(z*y) * exp(-3*I*z) * exp(6*I*z) ,z) / (y+1*I) / (y+2*I)
			+ 3 * integrate(exp(z*y) * exp(-3*I*z) * exp(4*I*z) ,z) / (y-1*I) / (y+2*I)
			+ 3 * integrate(exp(z*y) * exp(-3*I*z) * exp(2*I*z) ,z) / (y+1*I) / (y-2*I)
			+ 1 * integrate(exp(z*y) * exp(-3*I*z) * exp(0*I*z) ,z) / (y-1*I) / (y-2*I)
		)
	assertEquivalent(I3s0, I3s1,"Integral 3, step 1")
	I3s2	= 1/8 * (
			+ 1 * exp(z*y) * exp(-3*I*z) * exp(6*I*z) / (y+1*I) / (y+2*I) / (y+3*I)
			+ 3 * exp(z*y) * exp(-3*I*z) * exp(4*I*z) / (y-1*I) / (y+2*I) / (y+1*I)
			+ 3 * exp(z*y) * exp(-3*I*z) * exp(2*I*z) / (y+1*I) / (y-2*I) / (y-1*I)
			+ 1 * exp(z*y) * exp(-3*I*z) * exp(0*I*z) / (y-1*I) / (y-2*I) / (y-3*I)
		)
	assertEquivalent(I3s1, I3s2,"Integral 3, step 2")
	I3s3	= 1/8 * exp(z*y) * exp(-3*I*z) * (
			+ 1 * exp(6*I*z) / (y+1*I) / (y+2*I) / (y+3*I)
			+ 3 * exp(4*I*z) / (y-1*I) / (y+2*I) / (y+1*I)
			+ 3 * exp(2*I*z) / (y+1*I) / (y-2*I) / (y-1*I)
			+ 1 * exp(0*I*z) / (y-1*I) / (y-2*I) / (y-3*I)
		)
	assertEquivalent(I3s2, I3s3,"Integral 3, step 3")


	# I4

	I4p0	= I3s3 * (exp(+I*z)+exp(-I*z)) /2
	I4p1	= 1/8 * exp(z*y) * exp(-3*I*z) * (
			+ 1 * exp(6*I*z) / (y+1*I) / (y+2*I) / (y+3*I)
			+ 3 * exp(4*I*z) / (y-1*I) / (y+2*I) / (y+1*I)
			+ 3 * exp(2*I*z) / (y+1*I) / (y-2*I) / (y-1*I)
			+ 1 * exp(0*I*z) / (y-1*I) / (y-2*I) / (y-3*I)
		) * (exp(+I*z)+exp(-I*z)) /2
	assertEquivalent(I4p0, I4p1,"Pre Integral 4, step 1")
	I4p2	= 1/16 * exp(z*y) * exp(-4*I*z) * (
			+ 1 * exp(8*I*z) / (y+1*I) / (y+2*I) / (y+3*I)
			+ 3 * exp(6*I*z) / (y-1*I) / (y+2*I) / (y+1*I)
			+ 3 * exp(4*I*z) / (y+1*I) / (y-2*I) / (y-1*I)
			+ 1 * exp(2*I*z) / (y-1*I) / (y-2*I) / (y-3*I)
			+ 1 * exp(6*I*z) / (y+1*I) / (y+2*I) / (y+3*I)
			+ 3 * exp(4*I*z) / (y-1*I) / (y+2*I) / (y+1*I)
			+ 3 * exp(2*I*z) / (y+1*I) / (y-2*I) / (y-1*I)
			+ 1 * exp(0*I*z) / (y-1*I) / (y-2*I) / (y-3*I)
		)
	assertEquivalent(I4p1, I4p2,"Pre Integral 4, step 2")
	I4p3	= 1/16 * exp(z*y) * exp(-4*I*z) * (
			+ 1 * exp(8*I*z) / (y+1*I) / (y+2*I) / (y+3*I)
			+ 3 * (y+3*I) * exp(6*I*z) / (y-1*I) / (y+2*I) / (y+1*I) / (y+3*I)
			+ 1 * (y-1*I) * exp(6*I*z) / (y+1*I) / (y+2*I) / (y+3*I) / (y-1*I)
			+ 3 * (y+2*I) * exp(4*I*z) / (y+1*I) / (y-2*I) / (y-1*I) / (y+2*I)
			+ 3 * (y-2*I) * exp(4*I*z) / (y-1*I) / (y+2*I) / (y+1*I) / (y-2*I)
			+ 1 * (y+1*I) * exp(2*I*z) / (y-1*I) / (y-2*I) / (y-3*I) / (y+1*I)
			+ 3 * (y-3*I) * exp(2*I*z) / (y+1*I) / (y-2*I) / (y-1*I) / (y-3*I)
			+ 1 * exp(0*I*z) / (y-1*I) / (y-2*I) / (y-3*I)
		)
	assertEquivalent(I4p2, I4p3,"Pre Integral 4, step 3")

	I4p4	= 1/16 * exp(z*y) * exp(-4*I*z) * (
			+ 1 * exp(8*I*z) / (y+1*I) / (y+2*I) / (y+3*I)
			+ 4 * exp(6*I*z) / (y-1*I) / (y+1*I) / (y+3*I)
			+ 6 * exp(4*I*z) / (y+1*I) / (y-2*I) / (y-1*I) / (y+2*I) * y
			+ 4 * exp(2*I*z) / (y-1*I) / (y-3*I) / (y+1*I)
			+ 1 * exp(0*I*z) / (y-1*I) / (y-2*I) / (y-3*I)
		)
	assertEquivalent(I4p3, I4p4,"Pre Integral 4, step 4")

	I4p5	= 1/16 * (
			+ 1 * exp(z*y) * exp(-4*I*z) * exp(8*I*z) / (y+1*I) / (y+2*I) / (y+3*I)
			+ 4 * exp(z*y) * exp(-4*I*z) * exp(6*I*z) / (y-1*I) / (y+1*I) / (y+3*I)
			+ 6 * exp(z*y) * exp(-4*I*z) * exp(4*I*z) / (y+1*I) / (y-2*I) / (y-1*I) / (y+2*I) * y
			+ 4 * exp(z*y) * exp(-4*I*z) * exp(2*I*z) / (y-1*I) / (y-3*I) / (y+1*I)
			+ 1 * exp(z*y) * exp(-4*I*z) * exp(0*I*z) / (y-1*I) / (y-2*I) / (y-3*I)
		)
	assertEquivalent(I4p4, I4p5,"Pre Integral 4, step 5")


	I4s0	= integrate(I4p5, z)
	I4s1	= 1/16 * (
			+ 1 * integrate(exp(z*y) * exp(-4*I*z) * exp(8*I*z) ,z) / (y+1*I) / (y+2*I) / (y+3*I)
			+ 4 * integrate(exp(z*y) * exp(-4*I*z) * exp(6*I*z) ,z) / (y-1*I) / (y+1*I) / (y+3*I)
			+ 6 * integrate(exp(z*y) * exp(-4*I*z) * exp(4*I*z) ,z) / (y+1*I) / (y-2*I) / (y-1*I) / (y+2*I) * y
			+ 4 * integrate(exp(z*y) * exp(-4*I*z) * exp(2*I*z) ,z) / (y-1*I) / (y-3*I) / (y+1*I)
			+ 1 * integrate(exp(z*y) * exp(-4*I*z) * exp(0*I*z) ,z) / (y-1*I) / (y-2*I) / (y-3*I)
		)
	assertEquivalent(I4s0, I4s1,"Integral 4, step 1")
	I4s2	= 1/16 * (
			+ 1 * exp(z*y) * exp(-4*I*z) * exp(8*I*z) / (y+1*I) / (y+2*I) / (y+3*I) / (y+4*I)
			+ 4 * exp(z*y) * exp(-4*I*z) * exp(6*I*z) / (y-1*I) / (y+1*I) / (y+3*I) / (y+2*I)
			+ 6 * exp(z*y) * exp(-4*I*z) * exp(4*I*z) / (y+1*I) / (y-2*I) / (y-1*I) / (y+2*I)
			+ 4 * exp(z*y) * exp(-4*I*z) * exp(2*I*z) / (y-1*I) / (y-3*I) / (y+1*I) / (y-2*I)
			+ 1 * exp(z*y) * exp(-4*I*z) * exp(0*I*z) / (y-1*I) / (y-2*I) / (y-3*I) / (y-4*I)
		)
	assertEquivalent(I4s1, I4s2,"Integral 4, step 2")
	I4s3	=  1/16 * exp(z*y) * exp(-4*I*z) * (
			+ 1 * exp(8*I*z) / (y+1*I) / (y+2*I) / (y+3*I) / (y+4*I)
			+ 4 * exp(6*I*z) / (y-1*I) / (y+1*I) / (y+3*I) / (y+2*I)
			+ 6 * exp(4*I*z) / (y+1*I) / (y-2*I) / (y-1*I) / (y+2*I)
			+ 4 * exp(2*I*z) / (y-1*I) / (y-3*I) / (y+1*I) / (y-2*I)
			+ 1 * exp(0*I*z) / (y-1*I) / (y-2*I) / (y-3*I) / (y-4*I)
		)
	assertEquivalent(I4s2, I4s3, "Integral 4, step 3")
	I4s4	=  1/16 * exp(z*y) * (
			+ 1 * (
				+ exp(-4*I*z) * (y+1*I) * (y+2*I) * (y+3*I) * (y+4*I)
				+ exp(+4*I*z) * (y-1*I) * (y-2*I) * (y-3*I) * (y-4*I)
				) / (y*y+1) / (y*y+4) / (y*y+9) / (y*y+16)
			+ 4 * (
				+ exp(-2*I*z) * (y+3*I) * (y+2*I)
				+ exp(+2*I*z) * (y-3*I) * (y-2*I)
				) / (y*y+1) / (y*y+4) / (y*y+9) 
			+ 6 * exp(+0*I*z)
				/ (y*y+1) / (y*y+4)
		)
	assertEquivalent(I4s3, I4s4, "Integral 4, step 4")
	"""

	# from the above we deduce the following general formular for In(g(y,z))
	def In(z,y,order) :
		return 1/2**order * exp(z*(y+order*I)) * (
			sum([
				BN(order,k) * y * exp(-2*k*I*z)
				/ prod( [ 
					(y-j*I) 
					for j in xrange(-order+k, +k+1)
				])
				for k in xrange(0,order+1)
			]))

	""" # checked but commented out for speed

	assertEquivalent(maz.recall("Integral_1_g_x")[0], In(z,y,1), "Integral 1, composed")
	assertEquivalent(maz.recall("Integral_2_g_x")[0], In(z,y,2), "Integral 2, composed")
	assertEquivalent(I2s8, In(z,y,2), "Integral 2, composed")
	assertEquivalent(I3s2, In(z,y,3), "Integral 3, composed")
	assertEquivalent(I4s3, In(z,y,4), "Integral 4, composed")
	"""

	maz.check("IntroducingWo", "y=-I*wo where wo = w to (%i)",
	[
		In(z,-I*wo,order)
		for order in xrange(5)
	])

	maz.check("IntroducingWo", "Inlining function (%i)",
	[
		1/2**order * exp(z*(-I*wo+order*I)) * (
			sum([
				BN(order,k) * -I*wo * exp(-2*k*I*z)
				/ prod( [ 
					(-I*wo-j*I) 
					for j in xrange(-order+k, +k+1)
				])
				for k in xrange(0,order+1)
			]))
		for order in xrange(5)
	])
	maz.check("IntroducingWo", "Common factor -I*wo (%i)",
	[
		1/2**order * -I*wo * exp(z*(-wo*I+order*I))
		* sum([
			BN(order,k) * exp(-2*k*I*z)
			/ prod([
				(-wo*I-j*I)
				for j in xrange(-order+k, +k+1)
			])
			for k in xrange(0,order+1)
			])
		for order in xrange(5)
	], skip=[4])
	maz.check("IntroducingWo", "common factor I in outside exp (%i)",
	[
		1/2**order * -I*wo * exp(-z*I*(wo-order)) * (
		sum([
			BN(order,k) * exp(-2*k*I*z)
			/ prod([
				(-wo*I-j*I)
				for j in xrange(-order+k, +k+1)
			])
			for k in xrange(0,order+1)
		]))
		for order in xrange(5)
	])
	maz.check("IntroducingWo", "common factor I inside the productory (%i)",
	[
		1/2**order * -I*wo * exp(-z*I*(wo-order)) * (
		sum([
			BN(order,k) * exp(-2*k*I*z)
			/ prod([
				-I*(wo+j)
				for j in xrange(-order+k, +k+1)
			])
			for k in xrange(0,order+1)
		]))
		for order in xrange(5)
	])
	maz.check("IntroducingWo", "Removing the I form the productory (%i)",
	[
		1/2**order * wo * exp(-z*I*(wo-order)) * (
		sum([
			BN(order,k) * exp(-2*k*I*z) * I**(order)
			/ prod([
				(wo+j)
				for j in xrange(-order+k, +k+1)
			])
			for k in xrange(0,order+1)
		]))
		for order in xrange(5)
	])
	maz.check("IntroducingWo", "extracting I^order from the sumatory (%i)",
	[
		wo * exp(-I*z*wo) * (I/N(2))**order * exp(I*z*order) * (
		sum([
			BN(order,k) * exp(-2*k*I*z)
			/ prod( [ 
				(wo+j)
				for j in xrange(-order+k, +k+1)
			])
			for k in xrange(0,order+1)
		]))
		for order in xrange(5)
	])
	maz.check("IntroducingWo", "Introducing exp(Izn) into the sumatory (%i)",
	[
		wo * exp(-I*z*wo) * (I/2)**order * (
		sum([
			BN(order,k) * exp((order-2*k)*I*z)
			/ prod([ 
				(wo+j)
				for j in xrange(-order+k, +k+1)
			])
			for k in xrange(0,order+1)
		]))
		for order in xrange(5)
	])
	maz.check("IntroducingWo", "j'=j-order (%i)",
	[
		wo * exp(-I*z*wo) * (I/2)**order * (
		sum([
			BN(order,k) * exp((order-2*k)*I*z)
			/ prod([ 
				(wo+j-order+k)
				for j in xrange(order+1)
			])
			for k in xrange(0,order+1)
		]))
		for order in xrange(5)
	])
	maz.check("IntroducingWo", "k'=order-k (%i)",
	[
		wo * exp(-I*z*wo) * (I/2)**order * (
		sum([
			BN(order,(-k+order)) * exp((order-2*(-k+order))*I*z)
			/ prod([ 
				(wo+j-order+(-k+order))
				for j in xrange(order+1)
			])
			for k in xrange(0,order+1)
		]))
		for order in xrange(5)
	])
	maz.check("IntroducingWo", "simplify (%i)",
	[
		wo * exp(-I*z*wo) * (I/2)**order * (
		sum([
			exp((-order+2*k)*I*z)
			* BN(order,k)
			/ prod([ 
				(wo+j-k)
				for j in xrange(order+1)
			])
			for k in xrange(0,order+1)
		]))
		for order in xrange(5)
	])
	maz.check("IntroducingWo", "partial fraction (%i)",
	# See integer.py
	#	wo / prod([ (wo+j-r) for j in xrange(n+1) ])
	#	= sum([ (r-j) * N(-1)**(j) / f(j) / f(n-j) /(wo+j-r) for j in xrange(n+1) ])
	[
		wo * exp(-I*z*wo) * (I/2)**order * (
		sum([
			exp((-order+2*k)*I*z)
			* BN(order,k)
			* sum([ (k-j) * N(-1)**(j) / f(j) / f(order-j) / (wo+j-k) / wo for j in xrange(order+1) ])
			for k in xrange(0,order+1)
		]))
		for order in xrange(5)
	])

	maz.check("IntroducingWo", "common factor wo (%i)",
	[
		exp(-I*z*wo) * (I/2)**order * (
		sum([
			exp((-order+2*k)*I*z)
			* BN(order,k)
			* sum([
				(k-j)
				* N(-1)**(j)
				/ f(j)
				/ f(order-j)
				/ (wo+j-k)
			for j in xrange(order+1) ])
		for k in xrange(order+1) ]))
	for order in xrange(5) ])

	maz.check("IntroducingWo", "using binomial (%i)",
	[
		exp(-I*z*wo) * (I/2)**order * (
		sum([
			exp((2*k-order)*I*z)
			* BN(order,k)
			* sum([
				(k-j)
				* N(-1)**(j)
				* BN(order,j)
				/ f(order)
				/ (wo+j-k)
			for j in xrange(order+1) ])
		for k in xrange(order+1) ]))
	for order in xrange(5) ])

	maz.check("IntroducingWo", "joining sums (%i)",
	[
		exp(-I*z*wo) * (I/2)**order * (
		sum([
		sum([
			exp((2*k-order)*I*z)
			* (k-j)
			* N(-1)**(j)
			* BN(order,j)
			* BN(order,k)
			/ f(order)
			/ (wo+j-k)
		for j in xrange(order+1) ])
		for k in xrange(order+1) ]))
	for order in xrange(5) ])

	def Inwo(z, wo, order) :
		return exp(-I*z*wo) * (I/N(2))**N(order) * (
		sum([
		sum([
			exp(N(2*k-order)*I*z)
			* N(k-j)
			* N(-1)**N(j)
			* BN(order,j)
			* BN(order,k)
			/ f(order)
			/ (wo+j-k)
		for j in xrange(order+1) ])
		for k in xrange(order+1) ])
		)

	maz.check("IntroducingWo", "final function (%i)",
		[ Inwo(z,wo,order) for order in xrange(5)]
	)

	print Inwo(z,wo,0)
	print Inwo(z,wo,10)


	"""
	# Recalling the original integral

	integrate(cos(x)*Pl(sin(x))*exp(-1j*w*to*x))
		= integrate(cos(x)*Pl(sin(x))*exp(-1j*w*to*x))
		= sum ( [ (-1)**i  * In(i+1)(g(x),x) * di(Pl(x)) for i in xrange(l+1) ] )
		x = sin(z) = ((I*exp(-I*z)-I*exp(I*z))/2)
	"""

	maz.check("ShadowedJoint","initial, order 1+%i", [
		sum ([
			N(-1)**N(i)
			* Inwo(z,wo,i+1)
			* Pli(l,i,((I*exp(-I*z)-I*exp(I*z))/2))
			for i in xrange(l+1) ] )
		for l in xrange(5) ])

	maz.check("ShadowedJoint","substitute Inwo, order 1+%i", [
		sum ([
			N(-1)**N(i)
			* exp(-I*z*wo)
			* I**N(i+1)
			/ N(2)**N(i+1)
			* 
				sum([
				sum([
					exp(N(2*k-i-1)*I*z)
					* N(k-j)
					* N(-1)**N(j)
					* BN(i+1,j)
					* BN(i+1,k)
					/ f(i+1)
					/ (wo+j-k)
				for j in xrange(i+1+1) ])
				for k in xrange(i+1+1) ])
			* Pli(l,i,((I*exp(-I*z)-I*exp(I*z))/2))
			for i in xrange(l+1) ] )
		for l in xrange(5) ])

	maz.check("ShadowedJoint","k->m (for later use), order 1+%i", [
		sum ([
			N(-1)**N(i)
			* exp(-I*z*wo)
			* I**N(i+1)
			/ N(2)**N(i+1)
			* 
				sum([
				sum([
					exp(N(2*m-i-1)*I*z)
					* N(m-j)
					* N(-1)**N(j)
					* BN(i+1,j)
					* BN(i+1,m)
					/ f(i+1)
					/ (wo+j-m)
				for j in xrange(i+1+1) ])
				for m in xrange(i+1+1) ])
			* Pli(l,i,((I*exp(-I*z)-I*exp(I*z))/2))
			for i in xrange(l+1) ] )
		for l in xrange(5) ])

	maz.check("ShadowedJoint","isolating exp(-I*z), order 1+%i", [
		sum ([
			N(-1)**N(i)
			* exp(-I*z)**wo
			* I**N(i+1)
			/ N(2)**N(i+1)
			* 
				sum([
				sum([
					exp(-I*z)**N(i+1-2*m)
					* N(m-j)
					* N(-1)**N(j)
					* BN(i+1,j)
					* BN(i+1,m)
					/ f(i+1)
					/ (wo+j-m)
				for j in xrange(i+1+1) ])
				for m in xrange(i+1+1) ])
			* Pli(l,i,((I*exp(-I*z)-I/exp(-I*z))/2))
			for i in xrange(l+1) ] )
		for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp", "original, ShadowedJoint with x=exp(-I*z)", [
		eq.subs(z,I*sympy.ln(x))
		for eq in maz.recall("ShadowedJoint") ])

	maz.check("ShadowJointWithoutExp","inlining, order 1+%i", [
		sum ([
			N(-1)**N(i)
			* exp(wo*sympy.ln(x))
			* I**N(i+1)
			/ N(2)**N(i+1)
			* 
				sum([
				sum([
					x**N(i+1-2*m)
					* N(m-j)
					* N(-1)**N(j)
					* BN(i+1,j)
					* BN(i+1,m)
					/ f(i+1)
					/ (wo+j-m)
				for j in xrange(i+2) ])
				for m in xrange(i+2) ])
			* Pli(l,i,((I*x-I/x)/2))
			for i in xrange(l+1) ] )
		for l in xrange(5) ])

	maz.skip("ShadowJointWithoutExp","sympy doesn't see that x**wo=exp(wo*ln(x)), order 1+%i", [
		sum ([
			N(-1)**N(i)
			* x**wo
			* I**N(i+1)
			/ N(2)**N(i+1)
			* 
				sum([
				sum([
					x**N(i+1-2*m)
					* N(m-j)
					* N(-1)**N(j)
					* BN(i+1,j)
					* BN(i+1,m)
					/ f(i+1)
					/ (wo+j-m)
				for j in xrange(i+2) ])
				for m in xrange(i+2) ])
			* Pli(l,i,((I*x-I/x)/2))
			for i in xrange(l+1) ] )
		for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","substitute Pli, order 1+%i", [
		sum ([
			N(-1)**N(i)
			* x**wo
			* I**N(i+1)
			/ N(2)**N(i+1)
			* 
				sum([
				sum([
					x**N(i+1-2*m)
					* N(m-j)
					* N(-1)**N(j)
					* BN(i+1,j)
					* BN(i+1,m)
					/ f(i+1)
					/ (wo+j-m)
				for j in xrange(i+2) ])
				for m in xrange(i+2) ])
			* sum([
				(-1)**N(k//2)
				/ N(2)**N(l)
				* f(2*l-k)
				/ f(l-k//2)
				/ f(l-k-i)
				/ f(k//2)
				* ((I*x-I/x)/N(2))**N(l-k-i)
				for k in xrange(0,l-i+1,2) ] )
			for i in xrange(l+1) ] )
		for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","extracting 2's and I's, order 1+%i", [
		sum ([
			N(-1)**N(i)
			* x**wo
			* I**N(i+1)
			/ N(2)**N(i+1)
			* 
				sum([
				sum([
					x**N(i+1-2*m)
					* N(m-j)
					* N(-1)**N(j)
					* BN(i+1,j)
					* BN(i+1,m)
					/ f(i+1)
					/ (wo+j-m)
				for j in xrange(i+2) ])
				for m in xrange(i+2) ])
			* sum([
				(-1)**N(k//2)
				/ N(2)**N(l)
				* f(2*l-k)
				/ f(l-k//2)
				/ f(l-k-i)
				/ f(k//2)
				* I**N(l-k-i)
				/ N(2)**N(l-k-i)
				* (x-1/x)**N(l-k-i)
				for k in xrange(0,l-i+1,2) ] )
			for i in xrange(l+1) ] )
		for l in xrange(5) ], skip=[3,4])

	maz.check("ShadowJointWithoutExp","simplifying 2's and I's, order 1+%i", [
		sum ([
			N(-1)**N(i)
			* x**wo
			* I**N(l+1)
			/ N(2)**N(2*l+1)
			* 
				sum([
				sum([
					x**N(i+1-2*m)
					* N(m-j)
					* N(-1)**N(j)
					* BN(i+1,j)
					* BN(i+1,m)
					/ f(i+1)
					/ (wo+j-m)
				for j in xrange(i+2) ])
				for m in xrange(i+2) ])
			*
				sum([
					1
					* N(2)**N(k)
					* f(2*l-k)
					/ f(l-k//2)
					/ f(l-k-i)
					/ f(k//2)
					* (x-1/x)**N(l-k-i)
				for k in xrange(0,l-i+1,2) ])
			for i in xrange(l+1) ] )
		for l in xrange(5) ], skip=[3,4])

	maz.check("ShadowJointWithoutExp","extracting constants, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum ([
			N(-1)**N(i)
			* 
				sum([
				sum([
					x**N(i+1-2*m)
					* N(m-j)
					* N(-1)**N(j)
					* BN(i+1,j)
					* BN(i+1,m)
					/ f(i+1)
					/ (wo+j-m)
				for j in xrange(i+2) ])
				for m in xrange(i+2) ])
			*
				sum([
					1
					* N(2)**N(k)
					* f(2*l-k)
					/ f(l-k//2)
					/ f(l-k-i)
					/ f(k//2)
					* (x-1/x)**N(l-k-i)
				for k in xrange(0,l-i+1,2) ])
			for i in xrange(l+1) ] )
		for l in xrange(5) ], skip=[3,4])

	maz.check("ShadowJointWithoutExp","expand binomial power, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum ([
			N(-1)**N(i)
			* 
				sum([
				sum([
					x**N(i+1-2*m)
					* N(m-j)
					* N(-1)**N(j)
					* BN(i+1,j)
					* BN(i+1,m)
					/ f(i+1)
					/ (wo+j-m)
				for j in xrange(i+2) ])
				for m in xrange(i+2) ])
			*
				sum([
					1
					* N(2)**N(k)
					* f(2*l-k)
					/ f(l-k//2)
					/ f(l-k-i)
					/ f(k//2)
					* sum([
						BN(l-k-i, r)
						* x**N(r-(l-k-i-r))
						* N(-1)**N(l-k-i-r)
						for r in xrange(l-k-i+1) ])
				for k in xrange(0,l-i+1,2) ])
			for i in xrange(l+1) ] )
		for l in xrange(5) ], skip=[3,4])

	maz.check("ShadowJointWithoutExp","reagrouping summatories, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(i)
			* N(2)**N(k)
			* f(2*l-k)
			/ f(l-k//2)
			/ f(l-k-i)
			/ f(k//2)
			* BN(l-k-i, r)
			* x**N(2*r-l+k+i)
			* N(-1)**N(l-k-i-r)
			* x**N(i+1-2*m)
			* N(m-j)
			* N(-1)**N(j)
			* BN(i+1,j)
			* BN(i+1,m)
			/ f(i+1)
			/ (wo+j-m)
		for r in xrange(l-k-i+1) ])
		for k in xrange(0,l-i+1,2) ])
		for j in xrange(i+2) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
	for l in xrange(5) ], skip=[3,4])

	maz.check("ShadowJointWithoutExp","sorting terms, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(i)
			* N(-1)**N(l-k-i-r)
			* N(-1)**N(j)
			* N(m-j)
			* N(2)**N(k)
			* BN(l-k-i, r)
			* BN(i+1,j)
			* BN(i+1,m)
			* f(2*l-k)
			/ f(l-k//2)
			/ f(l-k-i)
			/ f(k//2)
			/ f(i+1)
			* x**N(2*r-l+k+i)
			* x**N(i+1-2*m)
			/ (wo+j-m)
		for r in xrange(l-k-i+1) ])
		for k in xrange(0,l-i+1,2) ])
		for j in xrange(i+2) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","joining terms, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(l-k-r+j)
			* N(m-j)
			* N(2)**N(k)
			* BN(l-k-i, r)
			* BN(i+1,j)
			* BN(i+1,m)
			* f(2*l-k)
			/ f(l-k//2)
			/ f(l-k-i)
			/ f(k//2)
			/ f(i+1)
			* x**N(2*r-l-2*m+k+2*i+1)
			/ (wo+j-m)
		for r in xrange(l-k-i+1) ])
		for k in xrange(0,l-i+1,2) ])
		for j in xrange(i+2) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","introducing 2*p==2*r-2*m+k+2*i+2, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(l-k-r+j)
			* N(m-j)
			* N(2)**N(k)
			* BN(l-k-i, r)
			* BN(i+1,j)
			* BN(i+1,m)
			* f(2*l-k)
			/ f(l-k//2)
			/ f(l-k-i)
			/ f(k//2)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+j-m)
			if 2*p == 2*r-2*m+k+2*i+2
			else 0
		for r in xrange(l-k-i+1) ])
		for k in xrange(0,l-i+1,2) ])
		for j in xrange(i+2) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","r bounds to conditions and relaxing r, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(l-k-r+j)
			* N(m-j)
			* N(2)**N(k)
			* BN(l-k-i, r)
			* BN(i+1,j)
			* BN(i+1,m)
			* f(2*l-k)
			/ f(l-k//2)
			/ f(l-k-i)
			/ f(k//2)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+j-m)
			if True
			and 2*p == 2*r-2*m+k+2*i+2
			and r <= l-k-i
			and r >=0
			else 0
		for k in xrange(0,l-i+1,2) ])
		for r in xrange(l-i+1) ])
		for j in xrange(i+2) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","substituting k, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(l+j+r)
			* N(m-j)
			* N(2)**N(2*(p-r+m-i-1))
			* BN(l+i-2*(p-r+m-1), r)
			* BN(i+1,j)
			* BN(i+1,m)
			* f(2*l-2*(p-r+m-i-1))
			/ f(l-(p-r+m-i-1))
			/ f(l-(2*p-2*r+2*m-2*i-2)-i)
			/ f(p-r+m-i-1)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+j-m)
			if True
			and k == 2*p-2*r+2*m-2*i-2
			and (2*p-2*r+2*m-2*i-2) <= l-r-i
			and r >=0
			else 0
		for k in xrange(0,l-i+1,2) ])
		for r in xrange(l-i+1) ])
		for j in xrange(i+2) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","k loop as conditions and removed, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(l+j+r)
			* N(m-j)
			* N(2)**N(2*(p-r+m-i-1))
			* BN(l+i-2*(p-r+m-1), r)
			* BN(i+1,j)
			* BN(i+1,m)
			* f(2*l-2*(p-r+m-i-1))
			/ f(l-(p-r+m-i-1))
			/ f(l-(2*p-2*r+2*m-2*i-2)-i)
			/ f(p-r+m-i-1)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+j-m)
			if True
			and 2*p-2*r+2*m-2*i-2 >= 0
			and 2*p-2*r+2*m-2*i-2 <= l-i
			and (2*p-2*r+2*m-2*i-2) <= l-r-i
			and r >=0
			else 0
		for r in xrange(l-i+1) ])
		for j in xrange(i+2) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","rewrite conditions, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(l+j+r)
			* N(m-j)
			* N(2)**N(2*(p-r+m-i-1))
			* BN(l+i-2*(p-r+m-1), r)
			* BN(i+1,j)
			* BN(i+1,m)
			* f(2*l-2*(p-r+m-i-1))
			/ f(l-(p-r+m-i-1))
			/ f(l-(2*p-2*r+2*m-2*i-2)-i)
			/ f(p-r+m-i-1)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+j-m)
			if True
			and p-r+m-i-1 >= 0
			and 2*p-2*r+2*m-2*i-2 <= l-i
			and 2*p-2*r+2*m-2*i-2 <= l-i-r
			and r >=0
			else 0
		for r in xrange(l-i+1) ])
		for j in xrange(i+2) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","r is positive so last condition is more restrictive, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(l+j+r)
			* N(m-j)
			* N(2)**N(2*(p-r+m-i-1))
			* BN(l+i-2*(p-r+m-1), r)
			* BN(i+1,j)
			* BN(i+1,m)
			* f(2*l-2*(p-r+m-i-1))
			/ f(l-(p-r+m-i-1))
			/ f(l-(2*p-2*r+2*m-2*i-2)-i)
			/ f(p-r+m-i-1)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+j-m)
			if True
			and p-r+m-i-1 >= 0
			and 2*p+2*m-2 <= l+r+i
			and r >=0
			else 0
		for r in xrange(l-i+1) ])
		for j in xrange(i+2) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","setting a stricter condition for r, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(l+j+r)
			* N(m-j)
			* N(2)**N(2*(p-r+m-i-1))
			* BN(l+i-2*(p-r+m-1), r)
			* BN(i+1,j)
			* BN(i+1,m)
			* f(2*l-2*(p-r+m-i-1))
			/ f(l-(p-r+m-i-1))
			/ f(l-(2*p-2*r+2*m-2*i-2)-i)
			/ f(p-r+m-i-1)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+j-m)
		for r in xrange(max(0,2*p+2*m-2-l-i), p+m-i) ])
		for j in xrange(i+2) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])


	maz.check("ShadowJointWithoutExp","just the (w+(2*p-l-1)) fractions are non zero according to the goal, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(l+j+r)
			* N(m-j)
			* N(2)**N(2*(p-r+m-i-1))
			* BN(l+i-2*(p-r+m-1), r)
			* BN(i+1,j)
			* BN(i+1,m)
			* f(2*l-2*(p-r+m-i-1))
			/ f(l-(p-r+m-i-1))
			/ f(l-(2*p-2*r+2*m-2*i-2)-i)
			/ f(p-r+m-i-1)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+j-m)
			if True
			and +j-m == 2*p-l-1
			else 0
		for r in xrange(max(0,2*p+2*m-2-l-i), p+m-i) ])
		for j in xrange(i+2) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","replacing j=2*p-l-1+m , order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(l+(2*p-l-1+m)+r)
			* N(m-(2*p-l-1+m))
			* N(2)**N(2*(p-r+m-i-1))
			* BN(l+i-2*(p-r+m-1), r)
			* BN(i+1,(2*p-l-1+m))
			* BN(i+1,m)
			* f(2*l-2*(p-r+m-i-1))
			/ f(l-(p-r+m-i-1))
			/ f(l-(2*p-2*r+2*m-2*i-2)-i)
			/ f(p-r+m-i-1)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+(2*p-l-1+m)-m)
			if True
			and +j == 2*p-l-1+m
			else 0
		for r in xrange(max(0,2*p+2*m-2-l-i), p+m-i) ])
		for j in xrange(i+2) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","j bounds to conditions, j summ removed, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(l+(2*p-l-1+m)+r)
			* N(m-(2*p-l-1+m))
			* N(2)**N(2*(p-r+m-i-1))
			* BN(l+i-2*(p-r+m-1), r)
			* BN(i+1,(2*p-l-1+m))
			* BN(i+1,m)
			* f(2*l-2*(p-r+m-i-1))
			/ f(l-(p-r+m-i-1))
			/ f(l-(2*p-2*r+2*m-2*i-2)-i)
			/ f(p-r+m-i-1)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+(2*p-l-1+m)-m)
			if True
			and 2*p-l-1+m >= 0
			and 2*p-l-1+m <= i+1
			else 0
		for r in xrange(max(0,2*p+2*m-2-l-i), p+m-i) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","expression simplification, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(-1+m+r)
			* N(-2*p+l+1)
			* N(2)**N(2*(p-r+m-i-1))
			* BN(l+i-2*(p-r+m-1), r)
			* BN(i+1,(2*p-l-1+m))
			* BN(i+1,m)
			* f(2*(l-p+r-m+i+1))
			/ f(l-p+r-m+i+1)
			/ f(l-2*p+2*r-2*m+i+2)
			/ f(p-r+m-i-1)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+(2*p-l-1))
			if True
			and 2*p-l-1+m >= 0
			and 2*p-l-1+m <= i+1
			else 0
		for r in xrange(max(0,2*p+2*m-2-l-i), p+m-i) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","r'=r-p-m, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(-1+r+p)
			* N(-2*p+l+1)
			* N(2)**N(2*(-r-i-1))
			* BN(l+i+2*r+2, r+p+m)
			* BN(i+1,2*p-l-1+m)
			* BN(i+1,m)
			* f(2*(l+r+i+1))
			/ f(l+r+i+1)
			/ f(l+2*r+i+2)
			/ f(-r-i-1)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+(2*p-l-1))
			if True
			and 2*p-l-1+m >= 0
			and 2*p-l-1+m <= i+1
			and r >= -1 +p+m-1 -l-i
	#		and r <= -i-1
			else 0
		for r in xrange(-p-m, -i) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","r'=r+i, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(-1+r-i+p)
			* N(-2*p+l+1)
			* N(2)**N(2*(-r-1))
			* BN(l+2*r-i+2, r-i+p+m)
			* BN(i+1,2*p-l-1+m)
			* BN(i+1,m)
			* f(2*(l+r+1))
			/ f(l+r+1)
			/ f(l+2*r-i+2)
			/ f(-r-1)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+(2*p-l-1))
			if True
			and 2*p-l-1+m >= 0
			and 2*p-l-1+m <= i+1
			and r >= -1 +p+m-1 -l
	#		and r <= -1
			else 0
		for r in xrange(-p-m-i, 0) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","r'=r+1, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(r-2-i+p)
			* N(-2*p+l+1)
			* N(2)**N(2*(-r))
			* BN(l+2*r-i, r-1-i+p+m)
			* BN(i+1, 2*p-l-1+m)
			* BN(i+1,m)
			* f(2*(l+r))
			/ f(l+r)
			/ f(l+2*r-i)
			/ f(-r)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+(2*p-l-1))
			if True
			and 2*p-l-1+m >= 0
			and 2*p-l-1+m <= i+1
			and r >= +p+m-1 -l
	#		and r <= 0
			else 0
		for r in xrange(-p-m-i+1, 1) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","r'=-r, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(-r-i+p)
			* N(-2*p+l+1)
			* N(2)**N(2*r)
			* BN(l-2*r-i, p+m-r-1-i)
			* BN(i+1, 2*p-l-1+m)
			* BN(i+1,m)
			* f(2*(l-r))
			/ f(l-r)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			* x**N(2*p-l-1)
			/ (wo+(2*p-l-1))
			if True
			and 2*p-l-1+m >= 0
			and 2*p-l-1+m <= i+1
			and -r >= +p+m-1 -l
	#		and -r <= 0
			else 0
		for r in xrange(0,+p+m+i) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(l+2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","2*p'=p, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(-r-i+p//2)
			* N(-p+l+1)
			* N(2)**N(2*r)
			* BN(l-2*r-i, p//2+m-r-1-i)
			* BN(i+1, p-l-1+m)
			* BN(i+1, m)
			* f(2*(l-r))
			/ f(l-r)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			* x**N(p-l-1)
			/ (wo+(p-l-1))
			if True
			and p-l-1+m >= 0
			and p-l-1+m <= i+1
	#		and r <= +p//2 +m+i+1
			and r <= l+1 -p//2 -m
	#		and -r <= 0
			else 0
		for r in xrange(0,+p//2+m+i) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(0,2*l+2+1,2) ])
	for l in xrange(5) ])

	maz.check("ShadowJointWithoutExp","boost order, unchanged formula, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(-r-i+p//2)
			* N(-p+l+1)
			* N(2)**N(2*r)
			* BN(l-2*r-i, p//2+m-r-1-i)
			* BN(i+1, p-l-1+m)
			* BN(i+1, m)
			* f(2*(l-r))
			/ f(l-r)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			* x**N(p-l-1)
			/ (wo+(p-l-1))
			if True
			and p-l-1+m >= 0
			and p-l-1+m <= i+1
	#		and r <= +p//2 +m+i+1
			and r <= l+1 -p//2 -m
	#		and -r <= 0
			else 0
		for r in xrange(0,+p//2+m+i) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		for p in xrange(0,2*l+2+1,2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","moved p,l dependent out, unchanged formula, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(-r-i+p//2)
			* N(2)**N(2*r)
			* BN(l-2*r-i, p//2+m-r-1-i)
			* BN(i+1, p-l-1+m)
			* BN(i+1, m)
			* f(2*(l-r))
			/ f(l-r)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			if True
			and p+2*m-2*r-2-2*i >= 0
			and l-2*r-i >=0
			and p+2*m-2 <= 2*l-2*r
			and p-l-1+m >= 0
			and p-l-1+m <= i+1
	#		and r <= +p//2 +m+i+1
			and r <= l+1 -p//2 -m
	#		and -r <= 0
			else 0
		for r in xrange(0,+p//2+m+i) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		* N(-p+l+1)
		* x**N(p-l-1)
		/ (wo+(p-l-1))
		for p in xrange(0,2*l+2+1,2) ])
	for l in xrange(8) ])


	maz.check("ShadowJointWithoutExp","n=p-l-1, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			N(-1)**N(-r-i+(n+l+1)//2)
			* N(2)**N(2*r)
			* BN(l-2*r-i, (n+l+1)//2+m-r-1-i)
			* BN(i+1, n+m)
			* BN(i+1, m)
			* f(2*(l-r))
			/ f(l-r)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			if True
			and n+l-1+2*m-2*i >= 2*r
			and l-2*r-i >=0
			and 2*r <= l-n-2*m+1
			and n+m >= 0
			and n+m <= i+1
	#		and r <= +(n+l+1)//2 +m+i+1
			and r <= -(n-l-1)//2 -m
	#		and -r <= 0
			else 0
		for r in xrange(0,+(n+l+1)//2+m+i) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		* N(-n)
		* x**N(n)
		/ (wo+(n))
		for n in xrange(-l-1, l+1+1, 2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","obtaining 2*r, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			I**N(-2*r-2*i+n+l+1)
			* N(2)**N(2*r)
			* BN(l-2*r-i, (n+l-1)//2+m-r-i)
			* BN(i+1, n+m)
			* BN(i+1, m)
			* f(2*(l-r))
			/ f(l-r)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			if True
			and n+m >= 0
			and n+m <= i+1
	#		and m >= 0
			and m <= i+1
			and 2*r <= l -i
			and 2*r <= l -n +1 -2*m
			and 2*r <= l +n -1 +2*m -2*i
	#		and r >= 0
			else 0
		for r in xrange(l+1) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		* N(-n)
		* x**N(n)
		/ (wo+(n))
		for n in xrange(-l-1, l+1+1, 2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","spliting n sum by sign, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		(
		sum([
		sum([
		sum([
		sum([
			I**N(-2*r-2*i+n+l+1)
			* N(2)**N(2*r)
			* BN(l-2*r-i, (n+l-1)//2+m-r-i)
			* BN(i+1, n+m)
			* BN(i+1, m)
			* f(2*(l-r))
			/ f(l-r)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			if True
			and n+m >= 0
			and n+m <= i+1
	#		and m >= 0
			and m <= i+1
			and 2*r <= l -i
			and 2*r <= l -n +1 -2*m
			and 2*r <= l +n -1 +2*m -2*i
	#		and r >= 0
			else 0
		for r in xrange(l+1) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		* N(-n)
		* x**N(n)
		/ (wo+(n))
		for n in xrange(l+1, 0, -2) ])
		+
		sum([
		sum([
		sum([
		sum([
			I**N(-2*r-2*i+n+l+1)
			* N(2)**N(2*r)
			* BN(l-2*r-i, (n+l-1)//2+m-r-i)
			* BN(i+1, n+m)
			* BN(i+1, m)
			* f(2*(l-r))
			/ f(l-r)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			if True
			and n+m >= 0
			and n+m <= i+1
	#		and m >= 0
			and m <= i+1
			and 2*r <= l -i
			and 2*r <= l -n +1 -2*m
			and 2*r <= l +n -1 +2*m -2*i
	#		and r >= 0
			else 0
		for r in xrange(l+1) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		* N(-n)
		* x**N(n)
		/ (wo+(n))
		for n in xrange(-l-1, 0, 2) ])
		)
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","n'=-n, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		(
		sum([
		sum([
		sum([
		sum([
			I**N(-2*r-2*i+n+l+1)
			* N(2)**N(2*r)
			* BN(l-2*r-i, (n+l-1)//2+m-r-i)
			* BN(i+1, n+m)
			* BN(i+1, m)
			* f(2*(l-r))
			/ f(l-r)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			if True
	#		and m >= 0
			and m >= 0-n
			and m <= i+1-n
			and m <= i+1
			and 2*r <= l -i
			and 2*r <= l -n +1 -2*m
			and 2*r <= l +n -1 +2*m -2*i
	#		and r >= 0
			else 0
		for r in xrange(l+1) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		* N(-n)
		* x**N(n)
		/ (wo+(n))
		for n in xrange(l+1, 0, -2) ])
		+
		sum([
		sum([
		sum([
		sum([
			I**N(-2*r-2*i-n+l+1)
			* N(2)**N(2*r)
			* BN(l-2*r-i, (-n+l-1)//2+m-r-i)
			* BN(i+1, -n+m)
			* BN(i+1, m)
			* f(2*(l-r))
			/ f(l-r)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			if True
	#		and m >= 0
			and m >= n
			and m <= i+1+n
			and m <= i+1
			and 2*r <= l -i
			and 2*r <= l +n +1 -2*m
			and 2*r <= l -n -1 +2*m -2*i
	#		and r >= 0
			else 0
		for r in xrange(l+1) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		* N(-(-n))
		* x**N((-n))
		/ (wo+((-n)))
		for n in xrange(+l+1, 0, -2) ])
		)
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","Removing loose constraints on m, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		(
		sum([
		sum([
		sum([
		sum([
			I**N(-2*r-2*i+n+l+1)
			* N(2)**N(2*r)
			* BN(l-2*r-i, (n+l-1)//2+m-r-i)
			* BN(i+1, n+m)
			* BN(i+1, m)
			* f(2*(l-r))
			/ f(l-r)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			if True
			and m >= 0
			and m <= i+1-n
			and 2*r <= l -i
			and 2*r <= l -n +1 -2*m
			and 2*r <= l +n -1 +2*m -2*i
	#		and r >= 0
			else 0
		for r in xrange(l+1) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		* N(-n)
		* x**N(n)
		/ (wo+(n))
		for n in xrange(l+1, 0, -2) ])
		+
		sum([
		sum([
		sum([
		sum([
			I**N(-2*r-2*i-n+l+1)
			* N(2)**N(2*r)
			* BN(l-2*r-i, (-n+l-1)//2+m-r-i)
			* BN(i+1, -n+m)
			* BN(i+1, m)
			* f(2*(l-r))
			/ f(l-r)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			if True
			and m >= n
			and m <= i+1
			and 2*r <= l -i
			and 2*r <= l +n +1 -2*m
			and 2*r <= l -n +1 +2*m -2*i -2
	#		and r >= 0
			else 0
		for r in xrange(l+1) ])
		for m in xrange(i+2) ])
		for i in xrange(l+1) ])
		* N(-(-n))
		* x**N((-n))
		/ (wo+((-n)))
		for n in xrange(+l+1, 0, -2) ])
		)
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","m restrictions into sum bounds, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		(
		sum([
		sum([
		sum([
		sum([
			I**N(-2*r-2*i+n+l+1)
			* N(2)**N(2*r)
			* BN(l-2*r-i, (n+l-1)//2+m-r-i)
			* BN(i+1, n+m)
			* BN(i+1, m)
			* f(2*(l-r))
			/ f(l-r)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			if True
			and 2*r <= l -i
			and 2*r <= l +n +1 +2*m -2*n -4*m
			and 2*r <= l +n +1 +2*m -2*i -2
	#		and r >= 0
			else 0
		for r in xrange(l+1) ])
		for m in xrange(i-n+2) ])
		for i in xrange(l+1) ])
		* N(-n)
		* x**N(n)
		/ (wo+(n))
		for n in xrange(l+1, 0, -2) ])
		+
		sum([
		sum([
		sum([
		sum([
			I**N(-2*r-2*i-n+l+1)
			* N(2)**N(2*r)
			* BN(l-2*r-i, (+n+l-1)//2+m-r-i)
			* BN(i+1, m)
			* BN(i+1, m+n)
			* f(2*(l-r))
			/ f(l-r)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			if True
			and 2*r <= l -i
			and 2*r <= l -n +1 -2*m
			and 2*r <= l +n +1 +2*m -2*i -2
	#		and r >= 0
			else 0
		for r in xrange(l+1) ])
		for m in xrange(i-n+2) ])
		for i in xrange(l+1) ])
		* N(-(-n))
		* x**N((-n))
		/ (wo+((-n)))
		for n in xrange(+l+1, 0, -2) ])
		)
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","unify summatories, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			(
				(
					I**N(-2*r-2*i+n+l+1)
					* N(2)**N(2*r)
					* BN(l-2*r-i, (n+l-1)//2+m-r-i)
					* BN(i+1, n+m)
					* BN(i+1, m)
					* f(2*(l-r))
					/ f(l-r)
					/ f(l-2*r-i)
					/ f(r)
					/ f(i+1)
				)
				* N(-n)
				* x**N(n)
				/ (wo+n)
				+
				(
					I**N(-2*r-2*i-n+l+1)
					* N(2)**N(2*r)
					* BN(l-2*r-i, (+n+l-1)//2+m-r-i)
					* BN(i+1, m)
					* BN(i+1, m+n)
					* f(2*(l-r))
					/ f(l-r)
					/ f(l-2*r-i)
					/ f(r)
					/ f(i+1)
				)
				* N(+n)
				* x**N(-n)
				/ (wo-n)
			)
			if True
			and 2*r <= l -i
			and 2*r <= l +n +1 +2*m -2*n -4*m
			and 2*r <= l +n +1 +2*m -2*i -2
	#		and r >= 0
			else 0
		for r in xrange(l+1) ])
		for m in xrange(i-n+2) ])
		for i in xrange(l+1) ])
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","common factors from inner sum, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			I**N(-2*r-2*i+l+1)
			* N(2)**N(2*r)
			* BN(l-2*r-i, (n+l-2*r-1)//2+m-i)
			* BN(i+1, m+n)
			* BN(i+1, m)
			/ f(l-2*r-i)
			/ f(r)
			/ f(i+1)
			/ f(l-r)
			* f(2*(l-r))
			* (
				I**N(+n)
				* N(-n)
				* x**N(+n)
				/ (wo+n)
				+
				I**N(-n)
				* N(+n)
				* x**N(-n)
				/ (wo-n)
			)
			if True
			and 2*r <= l -i
			and 2*r <= l +n +1 +2*m -2*n -4*m
			and 2*r <= l +n +1 +2*m -2*i -2
	#		and r >= 0
			else 0
		for m in xrange(i-n+2) ])
		for r in xrange(l+1) ])
		for i in xrange(l+1) ])
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	def debug(result, *args) :
		if args[4] :
			print " ".join((str(arg) for arg in args))
		return result

	maz.check("ShadowJointWithoutExp","binomials -> factorials, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			I**N(-2*r-2*i+l+1)
			* N(2)**N(2*r)
			* f(l-2*r-i)
			/ f((l-2*r+n-1+2*m-2*i)//2)
			/ f((l-2*r-n+1-2*m)//2)
			* f(i+1)
			/ f(i+1-m)
			/ f(i+1-m-n)
			/ f(m)
			/ f(m+n)
			/ f(l-2*r-i)
			/ f(r)
			/ f(l-r)
			* f(2*(l-r))
			* (
				I**N(+n)
				* N(-n)
				* x**N(+n)
				/ (wo+n)
				+
				I**N(-n)
				* N(+n)
				* x**N(-n)
				/ (wo-n)
			)
			if True
			and 2*r <= l -i
			and 2*r <= l -n +1 -2*m
			and 2*r <= l +n -1 +2*m -2*i
			and 2*r <= l -n +1
	#		and r >= 0
			else 0
		for m in xrange(i-n+2) ])
		for i in xrange(l+1) ])
		for r in xrange((l+1-n)//2+1) ])
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","r'=2*r, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			I**N(-r-2*i+l+1)
			* N(2)**N(r)
			/ f((l-r+n-1+2*m-2*i)//2)
			/ f((l-r-n+1-2*m)//2)
			* f(i+1)
			/ f(i+1-m)
			/ f(i+1-m-n)
			/ f(m)
			/ f(m+n)
			/ f(r//2)
			/ f(l-r//2)
			* f(2*l-r)
			* (
				I**N(+n)
				* N(-n)
				* x**N(+n)
				/ (wo+n)
				+
				I**N(-n)
				* N(+n)
				* x**N(-n)
				/ (wo-n)
			)
			if True
			and r <= l -i
			and r <= l -n +1 -2*m
			and r <= l +n -1 +2*m -2*i
			and r <= l -n +1
	#		and r >= 0
			else 0
		for m in xrange(i-n+2) ])
		for i in xrange(l+1) ])
		for r in xrange(0,(l+1-n)+1,2) ])
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","reorder terms, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			I**N(-r-2*i+l+1)
			* N(2)**N(r)
			* f(i+1)
			/ f(i+1-m)
			/ f(i+1-m-n)
			/ f(m)
			/ f(m+n)
			* f(2*l-r)
			/ f(r//2)
			/ f(l-r//2)
			/ f((l-r+n+1)//2+m-i-1)
			/ f((l-r+n+1)//2-m-n)
			* (
				I**N(+n)
				* N(-n)
				* x**N(+n)
				/ (wo+n)
				+
				I**N(-n)
				* N(+n)
				* x**N(-n)
				/ (wo-n)
			)
			if True
			and r <= l -i
			and r <= l -n +1 -2*m
			and r <= l +n -1 +2*m -2*i
			and r <= l -n +1
	#		and r >= 0
			else 0
		for m in xrange(i-n+2) ])
		for i in xrange(l+1) ])
		for r in xrange(0,(l+1-n)+1,2) ])
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","x powers out of inner sums, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			I**N(-r-2*i+l+1)
			* N(2)**N(r)
			* f(i+1)
			/ f(i+1-m)
			/ f(i+1-m-n)
			/ f(m)
			/ f(m+n)
			* f(2*l-r)
			/ f(r//2)
			/ f(l-r//2)
			/ f((l-r+n+1)//2+m-i-1)
			/ f((l-r+n+1)//2-m-n)
			if True
			and 2*i <= l +n -1 +2*m -r
			and 2*m <= l+1 -n -r
			and i <= l -r
			and r <= l +1 -n
	#		and r >= 0
			else 0
		for m in xrange(i-n+2) ])
		for i in xrange(l+1) ])
		for r in xrange(0,(l+1-n)+1,2) ])
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","i<=l is stricter than i<=l+r, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		sum([
			I**N(-2*i-r+l+1)
			* N(2)**N(r)
			* f(i+1)
			/ f(i+1-m)
			/ f(i+1-m-n)
			/ f(m)
			/ f(m+n)
			* f(2*l-r)
			/ f(r//2)
			/ f(l-r//2)
			/ f((l-r+n+1)//2+m-i-1)
			/ f((l-r+n+1)//2-m-n)
			if True
			and 2*m >= ( 0           if l-1 +n -r > 2*i else 2*i -l -n +1 +r)
			and 2*m <= ( 2*i -2*n +2 if l-1 +n -r > 2*i else      l -n +1 -r )
			else 0
		for m in xrange(l+2) ])
		for i in xrange(0,l-r+1) ])
		for r in xrange(0,(l+1-n)+1,2) ])
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","split r sum range, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		(
			sum([
			sum([
			sum([
				I**N(-2*i-r+l+1)
				* N(2)**N(r)
				* f(i+1)
				/ f(i+1-m)
				/ f(i+1-m-n)
				/ f(m)
				/ f(m+n)
				* f(2*l-r)
				/ f(r//2)
				/ f(l-r//2)
				/ f((l-r+n+1)//2+m-i-1)
				/ f((l-r+n+1)//2-m-n)
				if True
				and l-1 +n -r > 2*i
				and 2*m >= 0           
				and 2*m <= 2*i -2*n +2
				else 0
			for m in xrange(l+2) ])
			for i in xrange(0,l-r+1) ])
			for r in xrange(0,(l+1-n)+1,2) ])
			+
			sum([
			sum([
			sum([
				I**N(-2*i-r+l+1)
				* N(2)**N(r)
				* f(i+1)
				/ f(i+1-m)
				/ f(i+1-m-n)
				/ f(m)
				/ f(m+n)
				* f(2*l-r)
				/ f(r//2)
				/ f(l-r//2)
				/ f((l-r+n+1)//2+m-i-1)
				/ f((l-r+n+1)//2-m-n)
				if True
				and l-1 +n -r <= 2*i
				and 2*m >= 2*i -l -n +1 +r
				and 2*m <=      l -n +1 -r
				else 0
			for m in xrange(l+2) ])
			for i in xrange(0,l-r+1) ])
			for r in xrange(0,(l+1-n)+1,2) ])
		)
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","first term conditions into sum bounds, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		(
			sum([
			sum([
			sum([
				I**N(-2*i-r+l+1)
				* N(2)**N(r)
				* f(i+1)
				/ f(i+1-m)
				/ f(i+1-m-n)
				/ f(m)
				/ f(m+n)
				* f(2*l-r)
				/ f(r//2)
				/ f(l-r//2)
				/ f((l-r+n+1)//2+m-i-1)
				/ f((l-r+n+1)//2-m-n)
			for m in xrange(i-n+2) ])
			for i in xrange(0,(l-r-2+n)//2+1) ])
			for r in xrange(0,(l+1-n)+1,2) ])
			+
			sum([
			sum([
			sum([
				I**N(-2*i-r+l+1)
				* N(2)**N(r)
				* f(i+1)
				/ f(i+1-m)
				/ f(i+1-m-n)
				/ f(m)
				/ f(m+n)
				* f(2*l-r)
				/ f(r//2)
				/ f(l-r//2)
				/ f((l-r+n+1)//2+m-i-1)
				/ f((l-r+n+1)//2-m-n)
				if True
				and 2*m >= 2*i -l -n +1 +r
				and 2*m <=      l -n +1 -r
				and r <= l+1-n
				else 0
			for m in xrange(l+2) ])
			for i in xrange((l+n-r-1)//2,l-r+1) ])
			for r in xrange(0,(l+1-n)+1,2) ])
		)
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","i'=2*i, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		(
			sum([
			sum([
			sum([
				I**N(-i-r+l+1)
				* N(2)**N(r)
				* f((i+2) //2)
				/ f((i-2*m+2) //2)
				/ f((i-2*m-2*n+2) //2)
				/ f(m)
				/ f(m+n)
				* f(2*l-r)
				/ f(r//2)
				/ f((2*l-r) //2)
				/ f((l-r+n+1+2*m-i-2) //2)
				/ f((l-r-n+1-2*m) //2)
			for m in xrange(i//2-n+2) ])
			for i in xrange(0,(l-r-2+n)+1,2) ])
			for r in xrange(0,(l+1-n)+1,2) ])
			+
			sum([
			sum([
			sum([
				I**N(-(i+l+n-r-1)-r+l+1)
				* N(2)**N(r)
				* f((l+i-r+1+n)//2)
				/ f((l+i-r+1+n-2*m)//2)
				/ f((l+i-r+1-n-2*m)//2)
				/ f(m)
				/ f(m+n)
				* f(2*l-r)
				/ f(r//2)
				/ f((2*l-r)//2)
				/ f((2*m-i)//2)
				/ f((l-r-n+1-2*m)//2)
			for m in xrange(i//2,(l-n-r+1)//2+1) ])
			for i in xrange(0,(-n+1+l-r)+1,2) ])
			for r in xrange(0,(l+1-n)+1,2) ])
		)
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","r and i summatories joined, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		(
			sum([
			sum([
				I**N(-i-r+l+1)
				* f((i+2) //2)
				/ f((i-2*m+2) //2)
				/ f((i-2*m-2*n+2) //2)
				/ f(m)
				/ f(m+n)
				/ f((l-r+n+1+2*m-i-2) //2)
				/ f((l-r-n+1-2*m) //2)
			for m in xrange(0,(i-2*n+2)//2+1) ])
			for i in xrange(0,(l-r+n-2)+1,2) ])
			+
			sum([
			sum([
				I**N(-i-n+2)
				* f((l+i-r+1+n)//2)
				/ f((l+i-r+1+n-2*m)//2)
				/ f((l+i-r+1-n-2*m)//2)
				/ f(m)
				/ f(m+n)
				/ f((2*m-i)//2)
				/ f((l-r-n+1-2*m)//2)
			for m in xrange(i//2,(l-n-r+1)//2+1) ])
			for i in xrange(0,(l-r-n+1)+1,2) ])
		)
			* N(2)**N(r)
			* f(2*l-r)
			/ f(r//2)
			/ f((2*l-r)//2)
		for r in xrange(0,(l+1-n)+1,2) ])
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","double m, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		(
			sum([
			sum([
				I**N(-i-r+l+1)
				* f((i+2) //2)
				/ f((i-m+2) //2)
				/ f((i-m-2*n+2) //2)
				/ f((l-r-i-1+n+m) //2)
				/ f((l-r  +1-n-m) //2)
				/ f(m//2)
				/ f(m//2+n)
			for m in xrange(0,(i-2*n+2)+1,2) ])
			for i in xrange(0,(l+n-r-2)+1,2) ])
			+
			sum([
			sum([
				I**N(-i-n+2)
				* f((l-r+i+1+n)//2)
				/ f((l-r+i+1+n-m)//2)
				/ f((l-r+i+1-n-m)//2)
				/ f((l-r  +1-n-m)//2)
				/ f(m//2)
				/ f(m//2+n)
				/ f((m-i)//2)
			for m in xrange(i,(l-n-r+1)+1,2) ])
			for i in xrange(0,(l-n-r+1)+1,2) ])
		)
			* N(2)**N(r)
			* f(2*l-r)
			/ f(r//2)
			/ f((2*l-r)//2)
		for r in xrange(0,(l+1-n)+1,2) ])
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])


	maz.check("ShadowJointWithoutExp","inverting i and m in the first term, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		(
			sum([
			sum([
				I**N(-i-r+l+1)
				* f((i+2) //2)
				/ f((i-m+2) //2)
				/ f((i-m-2*n+2) //2)
				/ f((l-r-i-1+n+m) //2)
				/ f((l-r  +1-n-m) //2)
				/ f(m//2)
				/ f(m//2+n)
			for i in xrange(m+2*n-2,(l+n-r-2)+1,2) ])
			for m in xrange(0,(l-n-r+1)+1,2) ])
			+
			sum([
			sum([
				I**N(-i-n+2)
				* f((l-r+i+1+n)//2)
				/ f((l-r+i+1+n-m)//2)
				/ f((l-r+i+1-n-m)//2)
				/ f((l-r  +1-n-m)//2)
				/ f((m-i)//2)
				/ f(m//2)
				/ f(m//2+n)
			for m in xrange(i,(l-n-r+1)+1,2) ])
			for i in xrange(0,(l-n-r+1)+1,2) ])
		)
			* N(2)**N(r)
			* f(2*l-r)
			/ f(r//2)
			/ f((2*l-r)//2)
		for r in xrange(0,(l+1-n)+1,2) ])
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","extracting the only non symmetric sumand, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		(
			# TODO: not algebraically checked that are antisymetric but, for the orders they are
			sum([
				I**N(-(s+2*n-2)-r+l+1)
				* f((s+2*n) //2)
				/ f((s+2*n-i) //2)
				/ f((l-r-s-n+1+i) //2)
				/ f((l-r  -n+1-i) //2)
				/ f((s-i) //2)
				/ f(i//2)
				/ f(i//2+n)
			for s in xrange(i,(l-n-r)+1,2) ])
			+
			sum([
				I**N(-i-n+2)
				* f((l-r+i+1+n)//2)
				/ f((l-r+i+1+n-m)//2)
				/ f((l-r+i+1-n-m)//2)
				/ f((l-r  +1-n-m)//2)
				/ f((m-i)//2)
				/ f(m//2)
				/ f(m//2+n)
			for m in xrange(i,(l-n-r)+1,2) ])
			+
			I**N(-i-n+2)
			* f((l-r+i+1+n)//2)
			/ f((i+2*n)//2)
			/ f((i)//2)
			/ f((l-n-r+1-i)//2)
			/ f((l-n-r+1)//2)
			/ f((l+n-r+1)//2)
		)
		for i in xrange(0,(l-n-r+1)+1,2) ])
			* N(2)**N(r)
			* f(2*l-r)
			/ f(r//2)
			/ f((2*l-r)//2)
		for r in xrange(0,(l+1-n)+1,2) ])
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","discarded symmetric terms, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		(
			I**N(-i-n+2)
			* f((l+n-r+1+i)//2)
			/ f((l-n-r+1-i)//2)
			/ f((l-n-r+1)//2)
			/ f((l+n-r+1)//2)
			/ f((i+2*n)//2)
			/ f((i)//2)
		)
		for i in xrange(0,(l-n-r+1)+1,2) ])
			* N(2)**N(r)
			* f(2*l-r)
			/ f(r//2)
			/ f((2*l-r)//2)
		for r in xrange(0,(l+1-n)+1,2) ])
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp","extracted i independent, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		(
			I**N(-i-n+2)
			* f((l-n-r+1+i+2*n)//2)
			/ f((l-n-r+1-i)//2)
			/ f((i+2*n)//2)
			/ f((i)//2)
		)
		for i in xrange(0,(l-n-r+1)+1,2) ])
			* N(2)**N(r)
			* f(2*l-r)
			/ f((l-n-r+1)//2)
			/ f((l-n-r+1+2*n)//2)
			/ f(r//2)
			/ f((2*l-r)//2)
		for r in xrange(0,(l+1-n)+1,2) ])
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp", "adding factors to create binomials, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
		sum([
		(
			I**N(-i)
			* f((l-n-r+1+i+2*n)//2)
			/ f((l-n-r+1-i)//2)
			/ f((2*i+2*n)//2)
			* f((2*i+2*n)//2)
			/ f((i+2*n)//2)
			/ f((i)//2)
		)
		for i in xrange(0,(l-n-r+1)+1,2) ])
			* I**N(-n+2)
			* N(2)**N(r)
			* f(2*l-r)
			/ f((l-n-r+1)//2)
			/ f((l-n-r+1+2*n)//2)
			/ f(r//2)
			/ f((2*l-r)//2)
		for r in xrange(0,(l+1-n)+1,2) ])
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp", "halve i, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
			I**(l-n-r+1)
			* I**N(-n+2)
			* N(2)**N(r)
			* f(2*l-r)
			/ f((l-n-r+1)//2)
			/ f((l-n-r+1+2*n)//2)
			/ f(r//2)
			/ f((2*l-r)//2)
		for r in xrange(0,(l+1-n)+1,2) ])
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp", "i summatory is just I**(l-n-r+1), order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
			I**(l+2*n-r-1)
			* N(2)**N(r)
			* f(2*l-r)
			/ f((l-n-r+1)//2)
			/ f((l-n-r+1+2*n)//2)
			/ f(r//2)
			/ f((2*l-r)//2)
		for r in xrange(0,(l+1-n)+1,2) ])
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp", "halve r, order 1+%i", [
		x**wo
		* I**N(l+1)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
			I**(l+2*n-2*r-1)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			+ I**N(-n) * x**N(-n) / (wo-n)
			- I**N(+n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowJointWithoutExp", "redistribute fasors, order 1+%i", [
		x**wo
		* I**N(2*l)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
			N(-1)**N(r)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			+ I**N(+n) * x**N(-n) / (wo-n)
			- I**N(-n) * x**N(+n) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	for formula in maz.recall("ShadowJointWithoutExp") : print formula; print

	maz.check("ShadowedJoint", "retaking old formula with x=exp(-I*z), %s", [
		eq.subs(x,exp(-I*z))
		for eq in maz.recall("ShadowJointWithoutExp") ], skip=[4])

	maz.check("ShadowedJoint", "inlining substituted formula, %s", [
		exp(-I*z*wo)
		* I**N(2*l)
		/ N(2)**N(2*l+1)
		*
		sum([
		sum([
			N(-1)**N(r)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange((l+1-n)//2+1) ])
		* N(n) 
		* (
			+ I**N(+n) * exp(+I*z*N(n)) / (wo-n)
			- I**N(-n) * exp(-I*z*N(n)) / (wo+n)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowedJoint", "dropping a 2 factor to the exponentials, %s", [
		exp(-I*z*wo)
		* I**N(2*l)
		/ N(2)**N(2*l)
		*
		sum([
		sum([
			N(-1)**N(r)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			+ I**N(+n) * exp(+I*z*N(n)) / (wo-n) / N(2)
			- I**N(-n) * exp(-I*z*N(n)) / (wo+n) / N(2)
		)
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowedJoint", "adding exponential fractions, %s", [
		exp(-I*z*wo)
		* I**N(2*l)
		/ N(2)**N(2*l)
		*
		sum([
		sum([
			N(-1)**N(r)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			+ wo*(exp(+I*n*(z+sympy.pi/N(2))) - exp(-I*n*(z+sympy.pi/N(2))))
			+ n *(exp(+I*n*(z+sympy.pi/N(2))) + exp(-I*n*(z+sympy.pi/N(2))))
		) / N(2) / (wo*wo-N(n*n))
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.skip("ShadowedJoint", "exponential into cos and sin, [sympy unable to tell] %s", [
		exp(-I*z*wo)
		* I**N(2*l)
		/ N(2)**N(2*l)
		*
		sum([
		sum([
			N(-1)**N(r)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			+ wo*I*sin(n*(z+sympy.pi/N(2)))
			+ n   *cos(n*(z+sympy.pi/N(2)))
		) / (wo*wo-N(n*n))
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	for formula in maz.recall("ShadowedJoint") : print formula; print

	maz.check("ShadowDefiniteIntegralWo", "taking from ShadowedJoint, %i", [
		eq.subs(z,sympy.pi/N(2)) - eq.subs(z,0)
		for eq in maz.recall("ShadowedJoint") ], skip=[4])

	maz.check("ShadowDefiniteIntegralWo", "inlining ShadowedJoint, %i", [
		exp(-I*sympy.pi*wo/N(2))
		* I**N(2*l)
		/ N(2)**N(2*l)
		*
		sum([
		sum([
			N(-1)**N(r)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			+ wo*I*sin(n*(sympy.pi))
			+ n   *cos(n*(sympy.pi))
		) / (wo*wo-N(n*n))
		for n in xrange(l+1, 0, -2) ])
		-
		I**N(2*l)
		/ N(2)**N(2*l)
		*
		sum([
		sum([
			N(-1)**N(r)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			+ wo*I*sin(n*(sympy.pi/N(2)))
			+ n   *cos(n*(sympy.pi/N(2)))
		) / (wo*wo-N(n*n))
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowDefiniteIntegralWo", "joint parts not being the same, %i", [
		I**N(2*l)
		/ N(2)**N(2*l)
		*
		sum([
		sum([
			N(-1)**N(r)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			+ exp(-I*sympy.pi/N(2)*wo) * wo*I*sin(n*(sympy.pi))
			+ exp(-I*sympy.pi/N(2)*wo) * n   *cos(n*(sympy.pi))
		) / (wo*wo-N(n*n))
		for n in xrange(l+1, 0, -2) ])
		-
		I**N(2*l)
		/ N(2)**N(2*l)
		*
		sum([
		sum([
			N(-1)**N(r)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			+ wo*I*sin(n*(sympy.pi/N(2)))
			+ n   *cos(n*(sympy.pi/N(2)))
		) / (wo*wo-N(n*n))
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowDefiniteIntegralWo", "add the two sides, %i", [
		I**N(2*l)
		/ N(2)**N(2*l)
		*
		sum([
		sum([
			N(-1)**N(r)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			+ exp(-I*sympy.pi/N(2)*wo) * wo*I*sin(n*(sympy.pi))
			+ exp(-I*sympy.pi/N(2)*wo) * n   *cos(n*(sympy.pi))
			- wo*I*sin(n*(sympy.pi/N(2)))
			- n   *cos(n*(sympy.pi/N(2)))
		) / (wo*wo-N(n*n))
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowDefiniteIntegralWo", "sin n*pi always zero, %i", [
		I**N(2*l)
		/ N(2)**N(2*l)
		*
		sum([
		sum([
			N(-1)**N(r)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			+ n   *cos(n*sympy.pi)     *exp(-I*sympy.pi/N(2)*wo)
			- wo*I*sin(n*sympy.pi/N(2))
			- n   *cos(n*sympy.pi/N(2))
		) / (wo*wo-N(n*n))
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	for formula in maz.recall("ShadowDefiniteIntegralWo") : print formula; print


# TODO: HERE! FIX THIS!!
	maz.check("ShadowDefiniteIntegral", "taking from ShadowDefiniteIntegral wo=-I*w*to, %i", [
		eq.subs(wo,w*to)
		for eq in maz.recall("ShadowDefiniteIntegralWo") ], skip=[4])

	maz.check("ShadowDefiniteIntegral", "inline substitution, %i", [
		I**N(2*l)
		/ N(2)**N(2*l)
		*
		sum([
		sum([
			N(-1)**N(r)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			+ n   *cos(n*(sympy.pi)) * exp(-I*sympy.pi/N(2)*(w*to))
			- (w*to)*I*sin(n*(sympy.pi/N(2)))
			- n   *cos(n*(sympy.pi/N(2)))
		) / ((w*to)*(w*to)-N(n*n))
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowDefiniteIntegral", "simplify, %i", [
		I**N(2*l)
		/ N(2)**N(2*l)
		*
		sum([
		sum([
			N(-1)**N(r)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			+ n   *cos(n*(sympy.pi)) * exp(-sympy.pi/N(2)*I*w*to)
			- I*w*to*sin(n*(sympy.pi/N(2)))
			- n   *cos(n*(sympy.pi/N(2)))
		) / (-(I*w*I*w*to*to)-N(n*n))
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowDefiniteIntegral", "sign arrangement, %i", [
		I**N(2*l)
		/ N(2)**N(2*l)
		*
		sum([
		sum([
			N(-1)**N(r)
			* N(2)**N(2*r)
			* f(2*l-2*r)
			/ f((l-n+1)//2   -r)
			/ f((l-n+1)//2 +n-r)
			/ f(r)
			/ f(l-r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			- n   *cos(n*(sympy.pi)) * exp(-sympy.pi/N(2)*I*w*to)
			+ n   *cos(n*(sympy.pi/N(2)))
			+ I*w*to*sin(n*(sympy.pi/N(2)))
		) / ((I*w*I*w*to*to)+N(n*n))
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowDefiniteIntegral", "double factorial, %i", [
		I**N(2*l)
		/ N(2)**N(2*l)
		*
		sum([
		sum([
			I**N(2*r)
			* N(2)**(2*l+1)
			/ R(l-n+1-2*r)
			/ R(l+n+1-2*r)
			* R(2*l-2*r-1)
			/ R(2*r)
		for r in xrange(0,(l+1-n)//2+1) ])
		* N(n) 
		* (
			- n   *cos(n*(sympy.pi)) * exp(-sympy.pi/N(2)*I*w*to)
			+ n   *cos(n*(sympy.pi/N(2)))
			+ I*w*to*sin(n*(sympy.pi/N(2)))
		) / ((I*w*I*w*to*to)+N(n*n))
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowDefiniteIntegral", "r'=2*r, %i", [
		I**N(2*l)
		/ N(2)**N(2*l)
		*
		sum([
		sum([
			I**N(r)
			* N(2)**(2*l+1)
			/ R(l-n+1-r)
			/ R(l+n+1-r)
			* R(2*l-r-1)
			/ R(r)
		for r in xrange(0,(l+1-n)+1,2) ])
		* N(n) 
		* (
			- n   *cos(n*(sympy.pi)) * exp(-sympy.pi/N(2)*I*w*to)
			+ n   *cos(n*(sympy.pi/N(2)))
			+ I*w*to*sin(n*(sympy.pi/N(2)))
		) / ((I*w*I*w*to*to)+N(n*n))
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowDefiniteIntegral", "two step -> phasor, %i", [
		I**N(2*l)
		*
		sum([
		sum([
			I**N(r)
			* (1+N(-1)**N(r))
			* R(2*l-r-1)
			/ R(l-n+1-r)
			/ R(l+n+1-r)
			/ R(r)
		for r in xrange(0,(l+1-n)+1,2) ])
		* N(n) 
		* (
			- n   *cos(n*(sympy.pi)) * exp(-sympy.pi/N(2)*I*w*to)
			+ n   *cos(n*(sympy.pi/N(2)))
			+ I*w*to*sin(n*(sympy.pi/N(2)))
		) / ((I*w*I*w*to*to)+N(n*n))
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowDefiniteIntegral", "r'=-r, %i", [
		I**N(2*l)
		*
		sum([
		sum([
			debug(1,l,n,r, "+-", l+r-2, "/", r-n, r+n, l+1-r) *
			I**N(-r+l+1)
			* (1+N(-1)**N(l+1-r))
			* R(l+r-2)
			/ R(r-n)
			/ R(r+n)
			/ R(l-r+1)
		for r in xrange(n,l+2) ])
		* N(n) 
		* (
			- n   *cos(n*(sympy.pi)) * exp(-sympy.pi/N(2)*I*w*to)
			+ n   *cos(n*(sympy.pi/N(2)))
			+ I*w*to*sin(n*(sympy.pi/N(2)))
		) / ((I*w*I*w*to*to)+N(n*n))
		for n in xrange(l+1, 0, -2) ])
	for l in xrange(8) ])

	maz.check("ShadowDefiniteIntegral", "for all n, %i", [
		I**N(2*l)
		/ N(2)
		*
		sum([
		sum([
			I**N(-r+l+1)
			* (1+N(-1)**N(l+1-r))
			* (1+N(-1)**N(n-l-1))
			* R(l+r-2)
			/ R(r-n)
			/ R(r+n)
			/ R(l-r+1)
		for r in xrange(n,l+2) ])
		* N(n) 
		* (
			- n   *cos(n*(sympy.pi)) * exp(-sympy.pi/N(2)*I*w*to)
			+ n   *cos(n*(sympy.pi/N(2)))
			+ I*w*to*sin(n*(sympy.pi/N(2)))
		) / ((I*w*I*w*to*to)+N(n*n))
		for n in xrange(1,l+1+1) ])
	for l in xrange(8) ])

	maz.check("ShadowDefiniteIntegral", "trig to signs, %i", [
		I**N(2*l)
		/ N(2)
		*
		sum([
		sum([
			I**N(-r+l+1)
			* (1+N(-1)**N(l+1-r))
			* (1+N(-1)**N(n-l-1))
			* R(l+r-2)
			/ R(r-n)
			/ R(r+n)
			/ R(l-r+1)
		for r in xrange(n,l+2) ])
		* N(n) 
		* (
			- n   *(-1)**N(n) * exp(-sympy.pi/N(2)*I*w*to)
			+ n   *(I**N(n)+I**N(-n)) / N(2)
			+ w*to*(I**N(n)-I**N(-n)) / N(2)
		) / ((I*w*I*w*to*to)+N(n*n))
		for n in xrange(0,l+1+1) ])
	for l in xrange(8) ])

	maz.check("ShadowDefiniteIntegral", "join phases, %i", [
		I**N(2*l)
		/ N(2)
		*
		sum([
		sum([
			I**N(-r+l+1)
			* (1+N(-1)**N(l+1-r)+N(-1)**N(n-l-1)+N(-1)**N(n-r))
			* R(l+r-2)
			/ R(r-n)
			/ R(r+n)
			/ R(l-r+1)
		for r in xrange(n,l+2) ])
		* N(n) 
		* (
			- n   *(-1)**N(n) * exp(-sympy.pi/N(2)*I*w*to)
			+ n   *(I**N(n)+I**N(-n)) / N(2)
			- I*I*w*to*(I**N(n)-I**N(-n)) / N(2)
		) / ((I*w*I*w*to*to)+N(n*n))
		for n in xrange(0,l+1+1) ])
	for l in xrange(8) ])

	maz.check("ShadowDefiniteIntegral", "r'=l-r+1 and n starts at 1, %i", [
		I**N(2*l)
		*
		sum([
		sum([
			0 if r&1 or (l+1+n)&1 else
			I**N(r)
			* R(2*l-1-r)
			/ R(l+1+n-r)
			/ R(l+1-n-r)
			/ R(r)
		for r in xrange(0, l+1-n+1) ])
		* N(n) 
		* (
			- N(2)*n*(-1)**N(n) * exp(-sympy.pi/N(2)*I*w*to)
			+ n     *(I**N(n)+I**N(-n))
			- I*I*w*to*(I**N(n)-I**N(-n))
		) / ((I*w*I*w*to*to)+N(n*n))
		for n in xrange(1,l+1+1) ])
	for l in xrange(8) ])

	maz.check("ShadowDefiniteIntegral", "l phasor on -1, %i", [
		N(-1)**N(l)
		* N(2)
		*
		sum([
		sum([
			0 if r&1 or (l+1+n)&1 else
			I**N(r)
			* N(n) 
			* R(2*l-1-r)
			/ R(l+1+n-r)
			/ R(l+1-n-r)
			/ R(r)
		for r in xrange(0, l+1-n+1) ])
		* (
			- n     *(-1)**N(n) * exp(-sympy.pi/N(2)*I*w*to)
			+ n     *(I**N(n)+I**N(-n)) / N(2)
			- I*I*w*to*(I**N(n)-I**N(-n)) / N(2)
		) / ((I*w*I*w*to*to)+N(n*n))
		for n in xrange(1,l+1+1) ])
	for l in xrange(8) ])

def ShadowDefiniteIntegral(l,w) :
	return (
		N(-1)**N(l)
		* N(2)
		*
		sum([
		sum([
			0 if r&1 or (l+1+n)&1 else
			I**N(r)
			* N(n) 
			* R(2*l-1-r)
			/ R(l+1+n-r)
			/ R(l+1-n-r)
			/ R(r)
		for r in xrange(0, l+1-n+1) ])
		* (
			- n     *(-1)**N(n) * exp(-sympy.pi/N(2)*I*w*to)
			+ n     *(I**N(n)+I**N(-n)) / N(2)
			- I*I*w*to*(I**N(n)-I**N(-n)) / N(2)
		) / ((I*w*I*w*to*to)+N(n*n))
		for n in xrange(1,l+1+1) ])
	)

if __name__ == "__main__" :
	maz.check("ShadowDefiniteIntegral", "Function, %i", [
		ShadowDefiniteIntegral(l,w)
	for l in xrange(8) ])

	for formula in maz.recall("ShadowDefiniteIntegral") : print formula; print

