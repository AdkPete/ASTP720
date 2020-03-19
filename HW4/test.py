import unittest
from nbody import *
import numpy as np
import astropy.units as u
import astropy.constants as const

def dsum_acc(IC , t  , Particle): ##test function to compute true acceleraation
	A = Vector( [ 0 , 0 , 0 ] )
	for i in IC:
		if i == Particle:
		
			continue
		
		A += Particle.soft_acc(i , t , 1 * u.pc)
	
	return A
	
	
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
	
	
	def test_soft(self):
		P1 = Particle(1e4 * u.pc , 0 * u.pc , 0 * u.pc , 1 * u.M_sun)
		P2 = Particle( 0 * u.pc , 0 * u.pc , 0 * u.pc , 1 * u.M_sun)
		a = P1.accel(P2 , 0)
		corr = P1.soft_acc(P2 , 0 , 0.0001 * u.pc) ##for very small epsilon, we should get the same asnwer
		self.assertEqual(a , corr)
		
	def test_tree_builder(self):
		P1 = Particle( 1 * u.pc , 1 * u.pc , 1 * u.pc , 1 * u.Msun)
		P2 = Particle( 9 * u.pc , 9 * u.pc , 9 * u.pc , 1 * u.Msun)
		P3 = Particle(3 * u.pc , 1 * u.pc , 3 * u.pc , 1 * u.Msun)
		IC = [P1 , P2 , P3]
		
		Tree = Find_Tree(IC , 0)
		
		self.assertEqual(len(Tree.children) , 8)
		self.assertEqual(len(Tree.children[0].children) , 8)
		self.assertEqual(len(Tree.children[0].particles) , 2)

	
	def test_barnes_hut_acc(self):
		P1 = Particle( 1 * u.pc , 1 * u.pc , 1 * u.pc , 1 * u.Msun)
		P2 = Particle( 9 * u.pc , 9 * u.pc , 9 * u.pc , 1 * u.Msun)
		P3 = Particle(3 * u.pc , 1 * u.pc , 3 * u.pc , 1 * u.Msun)
		P4 = Particle( 8.5 * u.pc , 9 * u.pc , 9 * u.pc , 1 * u.Msun)
		P5 = Particle( 8.5 * u.pc , 9.1 * u.pc , 9 * u.pc , 1 * u.Msun)
		IC = [P1 , P2 , P3]# , P4 , P5]
		Tree = Find_Tree(IC , 0)
		
		self.assertEqual(BH_Acceleration(IC , 0 ,Tree , P1) , dsum_acc(IC , 0  , P1))
	
unittest.main()
