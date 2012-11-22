import maz
import sympy as sp

x1, y1, z1 = sp.symbols("x1 y1 z1")
x2, y2, z2 = sp.symbols("x2 y2 z2")
x3, y3, z3 = sp.symbols("x3 y3 z3")
xs, ys, zs = sp.symbols("xs ys zs")


maz.check( "Denormalized", "Without D, Gain %i",
	[
		(sp.Matrix([
			[xs,ys],
			]) *
		sp.Matrix([
			[x1,y1],
			[x2,y2],
			]).inv())[0],
		(sp.Matrix([
			[xs,ys],
			]) *
		sp.Matrix([
			[x1,y1],
			[x2,y2],
			]).inv())[1],
	])

maz.check( "Denormalized", "Without D, Gain %i",
	[
	+ xs*y2/(x1*y2 - x2*y1)
	- ys*x2/(x1*y2 - x2*y1)
	,
	+ ys*x1/(x1*y2 - x2*y1)
	- xs*y1/(x1*y2 - x2*y1)
	])

maz.check( "Denormalized", "Without D, Gain %i",
	[
	(xs*y2 - ys*x2)/(x1*y2 - x2*y1)
	,
	(ys*x1 - xs*y1)/(x1*y2 - x2*y1)
	])

