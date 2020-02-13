import unittest
import matplotlib.pyplot as plt
import astropy.units as u

def test_ode(t , y):
	
	return [2 * t]

class solve_ode:
	
	def __init__(self , f , t0 , y0 , h = 1e-3):
		self.f = f
		self.t0 = t0
		self.y0 = y0
		self.h = h
		## f is a function of y and t -- f(t , y)
	
	def Forward_Euler(self, t_end):
		
		
		t = []
		y = []
		
		t.append(self.t0)
		y.append(self.y0)
		
		while t[-1] < t_end:
			
			fi = self.f(t[-1] , y[-1])
			t.append(t[-1] + self.h)
			ny = []
			for i in range(len(y[-1])):
				ny.append( y[-1][i] + fi[i] * self.h)
			y.append(ny)
			
		return t , y

	def Heun(self , t_end):
		
		t = []
		y = []
		
		t.append(self.t0)
		y.append(self.y0)
		
		while t[-1] < t_end:
			
			fi = self.f(t[-1] , y[-1])
			t.append(t[-1] + self.h)
			nyt = []
			for i in range(len(y[-1])):
				nyt.append( y[-1][i] + fi[i] * self.h)
			
			###nyt is the result of the predictor step
			
			fi1 = self.f(t[-1] , nyt)
			ny = []
			for i in range(len(y[-1])):
				ny.append( y[-1][i] + 0.5 * self.h * (fi1[i] + fi[i]))
			y.append(ny)
			
		return t , y
		
def test():
	A = solve_ode(test_ode , 0 , [0] )
	t , y = A.Heun(10)
	
	py = []
	for i in y:
		py.append(i[0])
		
	plt.plot(t , py)
	plt.show()
	
	
test()
