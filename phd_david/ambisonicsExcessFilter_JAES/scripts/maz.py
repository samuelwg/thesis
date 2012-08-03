
from collections import defaultdict
from sympy import simplify
import sympy

_db = defaultdict(lambda : [])
_cases = defaultdict(lambda : None)

def recall(name, step=-1) :
	return _db[name][step]

def pprint(name, step=-1, orderoffset=0) :
	for i, form in enumerate(recall(name,step)) :
		print "==order %i==="%(i+orderoffset)
		sympy.pprint(form)

def cases(name, cases) :
	_cases[name] = cases

def check(name, stepDescription, expression, skip=None) :
	thread = _db[name]
	fullDescription = name+(", step %i"%len(thread))
	if thread :
		if stepDescription:
			fullDescription+=": "+stepDescription
		assertEquivalent(thread[-1],expression, fullDescription, skip, cases=_cases[name])
	else :
		print "\033[34mStart:\033[0m %s"%fullDescription
	thread.append(expression)

def skip(name, stepDescription, expression, skip=None) : # skip param ignored to match 'check' signature
	check(name, stepDescription, expression, skip=True)

def assertEquivalent(expr1, expr2, description, skip=None, cases=None ) :
	if skip is True:
		print "\033[34;1mSkipped:\033[0m %s"%description
		return
	if type(expr1) is list :
		if cases is None : cases = xrange(len(expr1))
		if not skip : skip = []
		for i, (case,e1,e2) in enumerate(zip(cases,expr1,expr2)) :
			assertEquivalent(e1,e2,description%(case,), i in skip)
		return
	if simplify(expr1-expr2) == 0 :
		print "\033[32mPassed:\033[0m", description
		return
	print "\033[31mFailed %s!!!!!\033[0m"%description
	print "\033[33;1mExpected:\033[0m\n", expr1
	print "\033[33mResult:\033[0m\n", expr2


if __name__ == "__main__" :
	import sympy
	x,y=sympy.var("x y", real=True)
	n=sympy.var("n", integer=True)
	i=sympy.var("i", integer=True)
	check("first", "initial",
		(x+y)*(x+y)**(n-1)
		)
	check("first", "join term",
		(x+y)**n
		)
	print recall("first")
	skip("first", "expand binomial",
		sympy.Sum(sympy.binomial(n,i) * x**i * y**(n-i), (i,0,n))
		)
	print recall("first")
	print recall("first", 0)
	check("first", "bad step",
		sympy.simplify(sympy.Sum(x**i, (i,0,n)) - sympy.Sum(x**i-y**n,(i,0,n)))
		)


