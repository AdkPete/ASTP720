
import matrix

def lin_fit(y , t):
	###We are going to solve a matrix problem to find m and b, where y = m * t + b
	n = len(y)
	NY = []
	for i in y:
		NY.append([i])
	Y = matrix.Matrix(( n , 1 ) , NY)
	setup = []
	for i in range(n):
		setup.append([1 , t[i]])
	X = matrix.Matrix( (n , 2) , setup)
	Xdagger = X.transpose()
	XdagX = Xdagger * X
	theta = XdagX.inverse() * Xdagger
	theta *= Y
	
	return theta.elements[0][0] , theta.elements[1][0]
	

