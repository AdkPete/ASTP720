
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
		
	def test_2_param(self):
		def f(t1 , t2):
			return 5 + 2 * t1 + 3 * t2
		ti1 = []
		ti2 = []
		y = []
		for i in range(10):
			for k in range(10):
				ti1.append(i)
				ti2.append(k)
				y.append(f(i , k))
		a , b , g = more_general_fit(y , ti1 , ti2)
		mgplot(y , ti1 , ti2 , a , b , g)
		self.assertEqual(round(a , 0) , 5)
		self.assertEqual(round(b , 0) , 2)
		self.assertEqual(round(g , 0) , 3)

		
unittest.main()
