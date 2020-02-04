import unittest


class TestCalc(unittest.TestCase):
	
	def test_centered_diff(self):
		self.assertAlmostEqual( 5 , centered_diff(test_func_1 , 5 , 1e-4))
	
		
		
def test_func_1(x):
	return 5 * x + 10


def centered_diff(f , x , h):
	
	'''
	this function will return the derivative of the function f at a point x
	f should be a function
	x nad h both are numbers
	returns a value for the derivative
	'''
	
	return (f(x + h) - f(x - h)) / (2 * h)
	
if __name__ == "__main__":
	unittest.main()
