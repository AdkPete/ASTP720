import rootfinder
import interpolation

def f(x):
	## A function to find the roots of.
	## This will find the cube root of 129
	return 129 - x ** 3
	
def fprime(x):
	## Derivative of f
	return -3 * x ** 2

##This will use bisection to find the root of f.
root = rootfinder.bisection(f , 3 , 10)

##To use the secant method, we have:
root = rootfinder.Secant(f , 3 , 10)

##Using the Newton method, we have:

root = rootfinder.Newton(f , fprime , 3)

print (root , root ** 3)

###To test out our interpolation function, we need an array of data points

x = []
y = []
i = -2
while i < 2:
	x.append(i)
	y.append(f(i))
	i += .5
	
### To produce a function to find linear interpolations, use:
lin_interp = interpolation.piecewise_linear(x , y)

print (lin_interp(-1.75))

###To use a cubic spline interpolation:

cubic_interp = interpolation.natural_cubic(x , y)

print (cubic_interp(-1.75))
