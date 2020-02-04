import unittest


class TestCalc(unittest.TestCase):
	
	def test_centered_diff(self):
		self.assertAlmostEqual( 5 , centered_diff(test_func_1 , 5 , 1e-4))
		
	def test_rectangle(self):
		###note the low accuracy restriction due to the inefficient nature of this method
		self.assertEqual( round(right_rect_integrate(test_func_1 , 0 , 5 , 1e-5) , 3) , 112.5)
	
		
		
def test_func_1(x):
	return 5 * x + 10


def centered_diff(f , x , h):
	
	'''
	this function will return the derivative of the function f at a point x
	f should be a function
	x and h both are numbers
	returns a value for the derivative
	'''
	
	return (f(x + h) - f(x - h)) / (2 * h)
	
def right_rect_integrate(f , a , b , h):

	'''
	this function will return the integral of f from a to b
	f should be a function
	a and b are the limits of integration
	h is the width of the regions
	'''
	
	
	x = a
	integral = 0
	while x < b - h:
		integral += f(x) * h
		x += h
	
	###handles any remaining bit of the region
	if x < b:
		nh = b - x
		integral += f(x) * nh

	return integral
	
if __name__ == "__main__":
	unittest.main()
