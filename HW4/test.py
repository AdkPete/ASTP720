import unittest
from nbody import *
import numpy as np


class TestNbody(unittest.TestCase):
	
	def test_vect_add(self):
		V1 = Vector([ 1 , 2 , 3])
		V2 = Vector([0 , 3 , 5])
		A = Vector([1 , 5 , 8])
		self.assertEqual(A , V1 + V2)
		
	def test_vect_sub(self):
		V1 = Vector([ 1 , 2 , 3])
		V2 = Vector([ 2 , 3 , 4 ])
		A = Vector([ -1 , -1 , -1 ])
		self.assertEqual(A , V1 - V2)
		
	def test_vect_magnitude(self):
		V = Vector([ 5 , 5 , 5])
		A = np.sqrt(3 * 25)
		self.assertEqual(A , V.mag())
unittest.main()
