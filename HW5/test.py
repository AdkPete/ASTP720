
import unittest
from least_squares import *
import numpy as np

np.random.seed(42)
class TestLS(unittest.TestCase):

	def test_least_squares(self):
		def f(t):
			return 2 * t + 5
		times = []
		y = []
		for i in range(100):
			times.append(i)
			y.append(f(i))
			
		b , m = lin_fit(y , times)
		self.assertAlmostEqual(m , 2)
		self.assertAlmostEqual(b , 5)
		
	def test_with_err(self):
		def f(t):
			return 2 * t + 5 + np.random.rand() / 100.0
		times = []
		y = []
		for i in range(10):
			t = float(i) * 10
			times.append(t)
			y.append(f(t))
		
		b , m = lin_fit(y , times)

		self.assertEqual(round(m , 0) , 2)
		self.assertEqual(round(b , 0) , 5)
		
		
unittest.main()
