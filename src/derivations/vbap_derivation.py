import maz
import sympy as sp

x1, y1, z1 = sp.symbols("x1 y1 z1")
x2, y2, z2 = sp.symbols("x2 y2 z2")
x3, y3, z3 = sp.symbols("x3 y3 z3")
xs, ys, zs = sp.symbols("xs ys zs")

M = sp.Matrix( [
	[x1, y1, z1],
	[x2, y2, z2],
	[x3, y3, z3],
	] )

P = sp.Matrix( [
	[xs, ys, zs],
	] )

G = P.multiply(M.inv())


maz.check("InverseDeterminant", "original",
	M.det()
)
maz.check("InverseDeterminant", "expanded",
	x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1
)
maz.check("InverseDeterminant", "divided by x1",
	(
		+ (x1*z3 - x3*z1)*(x1*y2 - x2*y1)
		- (x1*z2 - x2*z1)*(x1*y3 - x3*y1)
	) / x1
)
maz.check("InverseDeterminant", "divided by y1",
	(
		+ (y1*x3 - y3*x1)*(y1*z2 - y2*z1)
		- (y1*x2 - y2*x1)*(y1*z3 - y3*z1)
	) / y1
)
maz.check("InverseDeterminant", "divided by z1",
	(
		+ (z1*y3 - z3*y1)*(z1*x2 - z2*x1)
		- (z1*y2 - z2*y1)*(z1*x3 - z3*x1)
	) / z1
)
maz.check("InverseDeterminant", "As dot of v1 of cross product of v2 to v3",
	M[0,:].dot(M[1,:].cross(M[2,:]))
)

maz.check("InverseDeterminant", "As dot of v2 of cross product of v3 to v1",
	M[1,:].dot(M[2,:].cross(M[0,:]))
)
maz.check("InverseDeterminant", "As dot of v3 of cross product of v1 to v2",
	M[2,:].dot(M[0,:].cross(M[1,:]))
)
maz.check("InverseDeterminant", "Inverting order inverts sign",
	-M[2,:].dot(M[1,:].cross(M[0,:]))
)


D = maz.recall("InverseDeterminant",2)

"""
maz.check( "VbapGains", "Original Formula, Gain %i",
	[G[0], G[1], G[2]]
	)
maz.check( "VbapGains", "Result expanded by sympy, Gain %i",
	[
		+ xs*(
			- 1
				* (+ x1*x2*y3 - y1*x2*x3 - x1*y2*x3 + y1*x2*x3)
				* (+ y1*z1*x2 + z1*x1*y2 - y1*x1*z2 - z1*x2*y1)
			+ (x2*y1 + x1*y2 - x2*y1)
				* (
						+ (x1*z3 - x3*z1)*(x1*y2 - x2*y1)
						- (x1*y3 - x3*y1)*(x1*z2 - x2*z1)
				)
			)
			/ x1
			/ (x1*y2 - x2*y1)
			/ (
					+ (x1*z3 - x3*z1)*(x1*y2 - x2*y1)
					- (x1*y3 - x3*y1)*(x1*z2 - x2*z1)
			)
		+ ys*(
			-(z2 - x2*z1/x1)*(x2*(y3 - x3*y1/x1)/(x1*(y2 - x2*y1/x1)) - x3/x1)/((y2 - x2*y1/x1)*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)) - x2/(x1*(y2 - x2*y1/x1)))
		+ zs*(
			x2*(y3 - x3*y1/x1)/(x1*(y2 - x2*y1/x1)) - x3/x1)/(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1),
		xs*((y3 - x3*y1/x1)*(-y1*(z2 - x2*z1/x1)/(x1*(y2 - x2*y1/x1)) + z1/x1)/((y2 - x2*y1/x1)*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)) - y1/(x1*(y2 - x2*y1/x1))) + ys*(1/(y2 - x2*y1/x1) + (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/((y2 - x2*y1/x1)**2*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1))) - zs*(y3 - x3*y1/x1)/((y2 - x2*y1/x1)*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)),
		xs*(y1*(z2 - x2*z1/x1)/(x1*(y2 - x2*y1/x1)) - z1/x1)/(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1) - ys*(z2 - x2*z1/x1)/((y2 - x2*y1/x1)*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)) + zs/(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)
	])
maz.check( "VbapGains", "Expand inner multiplication, Gain %i",
	[
		+ xs*x1
			* (y2*z3 - y3*z2)
			/ (
					+ (x1*z3 - x3*z1)*(x1*y2 - x2*y1)
					- (x1*z2 - x2*z1)*(x1*y3 - x3*y1)
			)
		+ ys*(
			-(z2 - x2*z1/x1)*(x2*(y3 - x3*y1/x1)/(x1*(y2 - x2*y1/x1)) - x3/x1)/((y2 - x2*y1/x1)*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)) - x2/(x1*(y2 - x2*y1/x1)))
		+ zs*(
			x2*(y3 - x3*y1/x1)/(x1*(y2 - x2*y1/x1)) - x3/x1)/(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1),
		xs*((y3 - x3*y1/x1)*(-y1*(z2 - x2*z1/x1)/(x1*(y2 - x2*y1/x1)) + z1/x1)/((y2 - x2*y1/x1)*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)) - y1/(x1*(y2 - x2*y1/x1))) + ys*(1/(y2 - x2*y1/x1) + (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/((y2 - x2*y1/x1)**2*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1))) - zs*(y3 - x3*y1/x1)/((y2 - x2*y1/x1)*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)),
		xs*(y1*(z2 - x2*z1/x1)/(x1*(y2 - x2*y1/x1)) - z1/x1)/(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1) - ys*(z2 - x2*z1/x1)/((y2 - x2*y1/x1)*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)) + zs/(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)
	])

maz.check( "VbapGains", "generalizing other terms by rotating x,y,z, Gain %i",
	[
		+ xs*x1
			* (y2*z3 - y3*z2)
			/ (
					+ (x1*z3 - x3*z1)*(x1*y2 - x2*y1)
					- (x1*z2 - x2*z1)*(x1*y3 - x3*y1)
			)
		+ ys*y1
			* (z2*x3 - z3*x2)
			/ (
					+ (y1*x3 - y3*x1)*(y1*z2 - y2*z1)
					- (y1*x2 - y2*x1)*(y1*z3 - y3*z1)
			)
		+ zs*z1
			* (x2*y3 - x3*y2)
			/ (
					+ (z1*y3 - z3*y1)*(z1*x2 - z2*x1)
					- (z1*y2 - z2*y1)*(z1*x3 - z3*x1)
			)
		,

		xs*((y3 - x3*y1/x1)*(-y1*(z2 - x2*z1/x1)/(x1*(y2 - x2*y1/x1)) + z1/x1)/((y2 - x2*y1/x1)*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)) - y1/(x1*(y2 - x2*y1/x1))) + ys*(1/(y2 - x2*y1/x1) + (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/((y2 - x2*y1/x1)**2*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1))) - zs*(y3 - x3*y1/x1)/((y2 - x2*y1/x1)*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)),
		xs*(y1*(z2 - x2*z1/x1)/(x1*(y2 - x2*y1/x1)) - z1/x1)/(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1) - ys*(z2 - x2*z1/x1)/((y2 - x2*y1/x1)*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)) + zs/(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)
	])

maz.check( "VbapGains", "generalizing other terms by rotating 1,2,3, Gain %i",
	[
		+ xs*x1
			* (y2*z3 - y3*z2)
			/ (
					+ (x1*z3 - x3*z1)*(x1*y2 - x2*y1)
					- (x1*z2 - x2*z1)*(x1*y3 - x3*y1)
			)
		+ ys*y1
			* (z2*x3 - z3*x2)
			/ (
					+ (y1*x3 - y3*x1)*(y1*z2 - y2*z1)
					- (y1*x2 - y2*x1)*(y1*z3 - y3*z1)
			)
		+ zs*z1
			* (x2*y3 - x3*y2)
			/ (
					+ (z1*y3 - z3*y1)*(z1*x2 - z2*x1)
					- (z1*y2 - z2*y1)*(z1*x3 - z3*x1)
			)
		,
		+ xs*x1
			* (y3*z1 - y1*z3)
			/ (
					+ (z3*x1 - x3*z1) * (x1*y2 - x2*y1)
					- (y3*x1 - x3*y1) * (x1*z2 - x2*z1)
				)
		+ ys*x1
			* (x1*z3 - x3*z1)
			/ (
					+ (z3*x1 - x3*z1) * (x1*y2 - x2*y1)
					- (y3*x1 - x3*y1) * (x1*z2 - x2*z1)
				)
		+ zs*x1
			* (x3*y1 - x1*y3)
			/ (
					+ (z3*x1 - x3*z1) * (x1*y2 - x2*y1)
					- (y3*x1 - x3*y1) * (x1*z2 - x2*z1)
			),

		xs*(y1*(z2 - x2*z1/x1)/(x1*(y2 - x2*y1/x1)) - z1/x1)/(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1) - ys*(z2 - x2*z1/x1)/((y2 - x2*y1/x1)*(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)) + zs/(z3 - (y3 - x3*y1/x1)*(z2 - x2*z1/x1)/(y2 - x2*y1/x1) - x3*z1/x1)
	])
"""
maz.check( "VbapGains", "Using determinant, Gain %i",
	[
		+ xs * (y2*z3 - y3*z2) / D
		+ ys * (z2*x3 - z3*x2) / D
		+ zs * (x2*y3 - x3*y2) / D
		,
		+ xs * (y3*z1 - y1*z3) / D
		+ ys * (z3*x1 - z1*x3) / D
		+ zs * (x3*y1 - x1*y3) / D
		,
		+ xs * (y1*z2 - y2*z1) / D
		+ ys * (z1*x2 - z2*x1) / D
		+ zs * (x1*y2 - x2*y1) / D
		,
	])

maz.check("Denormalized", "Taking previous series, Gain %i",
	[
		D*f
		for f in maz.recall("VbapGains")
	])
maz.check( "Denormalized", "Without D, Gain %i",
	[
		+ xs * (y2*z3 - y3*z2)
		+ ys * (z2*x3 - z3*x2)
		+ zs * (x2*y3 - x3*y2)
		,
		+ xs * (y3*z1 - y1*z3)
		+ ys * (z3*x1 - z1*x3)
		+ zs * (x3*y1 - x1*y3)
		,
		+ xs * (y1*z2 - y2*z1)
		+ ys * (z1*x2 - z2*x1)
		+ zs * (x1*y2 - x2*y1)
		,
	])


ca1, sa1, ce1, se1 = sp.symbols("ca1 sa1 ce1 se1")
ca2, sa2, ce2, se2 = sp.symbols("ca2 sa2 ce2 se2")
ca3, sa3, ce3, se3 = sp.symbols("ca3 sa3 ce3 se3")

maz.check("VbapGainsTrigonometric", "Taking previous series, Gain %i",
	[
		(D*f).subs(dict(
			x1=ce1*ca1, y1=ce1*sa1, z1=se1,
			x2=ce2*ca2, y2=ce2*sa2, z2=se2,
			x3=ce3*ca3, y3=ce3*sa3, z3=se3,
			))
		for f in maz.recall("VbapGains")
	])

maz.check("VbapGainsTrigonometric", "Substitution, Gain %i",
	[
		+ xs * ( ce2 * sa2 * se3 - ce3 * sa3 * se2 )
		+ ys * ( se2 * ce3 * ca3 - se3 * ce2 * ca2 )
		+ zs * ( ce2 * ce3 * ca2 * sa3 - ce2 * ce3 * ca3 * sa2 )
		,
		+ xs * ( ce3 * sa3 * se1 - ce1 * sa1 * se3 )
		+ ys * ( se3 * ce1 * ca1 - se1 * ce3 * ca3 )
		+ zs * ( ce3 * ce1 * ca3 * sa1 - ce3 * ce1 * ca1 * sa3 )
		,
		+ xs * ( ce1 * sa1 * se2 - ce2 * sa2 * se1 )
		+ ys * ( se1 * ce2 * ca2 - se2 * ce1 * ca1 )
		+ zs * ( ce1 * ce2 * ca1 * sa2 - ce1 * ce2 * ca2 * sa1 )
		,
	])

maz.check("VbapGainsTrigonometric", "Associative, Gain %i",
	[
	+ ce3 * se2 * (ca3 * ys - sa3 * xs)
	- ce2 * (
			+ se3 * (ca2 * ys - sa2 * xs)
	 		+ ce3 * (ca3 * sa2 - ca2 * sa3 ) * zs
			)
	,
	+ ce1 * se3 * (ca1 * ys - sa1 * xs)
	- ce3 * (
			+ se1 * (ca3 * ys - sa3 * xs)
	 		+ ce1 * (ca1 * sa3 - ca3 * sa1 ) * zs
			)
	,
	+ ce2 * se1 * (ca2 * ys - sa2 * xs)
	- ce1 * (
			+ se2 * (ca1 * ys - sa1 * xs)
	 		+ ce2 * (ca2 * sa1 - ca1 * sa2 ) * zs
			)
	])

# g3 =

# + xs ( ce1 sa1 se2 - ce2 sa2 se1 )
# + ys ( se1 ce2 ca2 - se2 ce1 ca1 )
# + zs ( ce1 ce2 ca1 sa2 - ce1 ce2 ca2 sa1 )

# + xs ( (se2 ce1) sa1 - (se1 ce2) sa2 )
# + ys ( (se1 ce2) ca2 - (se2 ce1) ca1 )
# + zs ( ce1 ce2 (ca1 sa2 - ca2 sa1) )

# + (se1 ce2) ca2 ys - (se2 ce1) ca1 ys
# + (se2 ce1) sa1 xs - (se1 ce2) sa2 xs
# + zs ce1 ce2 (ca1 sa2 - ca2 sa1)

# + (se1 ce2) (ca2 ys - sa2 xs)
# + (se2 ce1) (sa1 xs - ca1 ys)
# + zs ce1 ce2 (ca1 sa2 - ca2 sa1)

# + (se1 ce2) (ca2 ys - sa2 xs)
# + (se2 ce1) (sa1 xs - ca1 ys)
# + zs ce1 ce2 (ca1 sa2 - ca2 sa1)

# + (ce2 se1) (ca2 ys - sa2 xs)
# - (ce1 se2) (ca1 ys - sa1 xs)
# + ce1 zs ce2 (ca1 sa2 - ca2 sa1)

# + ce2 se1 (ca2 ys - sa2 xs)
# - ce1 (
#		+ se2 (ca1 ys - sa1 xs)
# 		+ ce2 (ca2 sa1 - ca1 sa2 ) zs
#		)


(
	+ 2*ys*zs * (
		+ x1*x2*y1*z2
		+ x1*x2*y2*z1
		- x1**2*y2*z2
		- x2**2*y1*z1

		+ x1*x3*y1*z3
		+ x1*x3*y3*z1
		- x1**2*y3*z3
		- x3**2*y1*z1

		+ x2*x3*y2*z3
		+ x2*x3*y3*z2
		- x2**2*y3*z3
		- x3**2*y2*z2
	)
	+ 2*ys*xs * (
		+ x2*y1*z1*z2
		- x1*y1*z2**2
		- x1*y1*z3**2
		+ x1*y2*z1*z2
		+ x1*y3*z1*z3
		- x2*y2*z1**2
		- x2*y2*z3**2
		+ x2*y3*z2*z3
		+ x3*y1*z1*z3
		+ x3*y2*z2*z3
		- x3*y3*z1**2
		- x3*y3*z2**2
	)
	+ 2*xs*zs * (
		+ x1*y1*y2*z2
		+ x1*y1*y3*z3
		- x1*y2**2*z1
		- x1*y3**2*z1
		- x2*y1**2*z2
		+ x2*y1*y2*z1
		+ x2*y2*y3*z3
		- x3*y1**2*z3
		- x2*y3**2*z2
		+ x3*y1*y3*z1
		- x3*y2**2*z3
		+ x3*y2*y3*z2
	)
	+ xs**2 * (
		+ (y2*z1-y1*z2)**2
		+ (y2*z3-y3*z2)**2
		+ (y1*z3-y3*z1)**2
	)
	+ ys**2 * (
		+ (x2*z1-x1*z2)**2
		+ (x2*z3-x3*z2)**2
		+ (x1*z3-x3*z1)**2
	)
	+ zs**2 * (
		+ (x2*y1-x1*y2)**2
		+ (x2*y3-x3*y2)**2
		+ (x1*y3-x3*y1)**2
	)
)







