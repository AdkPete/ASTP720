
import unittest
from least_squares import *

class TestLS(unittest.TestCase):

	def test_least_squares(self):
		def f(t):
			return 2 * t + 5 + np.random.rand() / 10.0
		times = []
		y = []
		for i in range(100):
			times.append(i)
			y.append(f(i))
			
		b , m = lin_fit(y , times)
		self.assertEqual(m , 2)
		self.assertEqual(b , 5)
		
unittest.main()
