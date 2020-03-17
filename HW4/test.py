import unittest
from nbody import *
import numpy as np
import astropy.units as u
import astropy.constants as const

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
		
	def test_vect_mult(self):
		V = Vector( [5 , 5 , 5])
		A = Vector( [ 10 , 10 , 10])
		self.assertEqual(A , V * 2.0)
		
	
	def test_vect_magnitude(self):
		V = Vector([ 5 , 5 , 5])
		A = np.sqrt(3 * 25)
		self.assertEqual(A , V.mag())
		
	def test_acc_nosoft(self):
		P1 = Particle(1 * u.pc , 0 * u.pc , 0 * u.pc , 1 * u.M_sun)
		P2 = Particle( 0 * u.pc , 0 * u.pc , 0 * u.pc , 1 * u.M_sun)
		a = P1.accel(P2 , 0)
		corr = Vector( [ -1 * const.G * 1 * u.M_sun / (1 * u.pc ** 2) , 0 * u.m / (u.s ** 2) , 0 * u.m / (u.s ** 2)])
		self.assertEqual(a , corr)
unittest.main()
