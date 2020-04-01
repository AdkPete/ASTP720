import unittest
from nbody import *
test_mode()
set_params("params.txt")
import numpy as np
import astropy.units as u
import astropy.constants as const
import time

def gen_ic(numb):
	IC = []
	for i in range(numb):
		P = Particle(np.random.rand() * 5 * u.Mpc , np.random.rand() * 5 * u.Mpc , np.random.rand() * 5 * u.Mpc , 1e12 * u.M_sun , i)
		P.r.append(Vector([np.random.rand() * 5 * u.Mpc , np.random.rand() * 5 * u.Mpc , np.random.rand() * 5 * u.Mpc]))
		IC.append(P)
	return IC

def dsum_acc(IC , Particle): ##test function to compute true acceleraation
	A = Vector( [ 0 , 0 , 0 ] )
	for i in IC:
		if i == Particle:
		
			continue
		
		A += Particle.soft_acc(i)
	
	return A
	
def timing(numb):
	IC = gen_ic(numb)
	s = time.time()
	N = Barnes_Hut(IC , 1000 * u.yr , t_end = 100 * u.yr)
	print ("Time required for Barnes Hut is {}".format(str(time.time() - s)))
	
	s = time.time()
	N = direct_summation(100 * u.yr , 1000 * u.yr , IC)
	print ("Time required for Direct Summation is {}".format(str(time.time() - s)))
	
	
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
		P1 = Particle(1 * u.pc , 0 * u.pc , 0 * u.pc , 1 * u.M_sun , 1)
		P2 = Particle( 0 * u.pc , 0 * u.pc , 0 * u.pc , 1 * u.M_sun , 1)
		a = P1.accel(P2)
		corr = Vector( [ -1 * const.G * 1 * u.M_sun / (1 * u.pc ** 2) , 0 * u.m / (u.s ** 2) , 0 * u.m / (u.s ** 2)])
		self.assertEqual(a , corr)
	
	
		
	def test_tree_builder(self):
		P1 = Particle( 1 * u.pc , 1 * u.pc , 1 * u.pc , 1 * u.Msun , 1)
		P2 = Particle( 9 * u.pc , 9 * u.pc , 9 * u.pc , 1 * u.Msun , 2)
		P3 = Particle(3 * u.pc , 1 * u.pc , 3 * u.pc , 1 * u.Msun , 3)
		IC = Sim([P1 , P2 , P3])
		
		Tree = Find_Tree(IC)
		
		self.assertEqual(len(Tree.children) , 8)
		self.assertEqual(len(Tree.children[0].children) , 8)
		self.assertEqual(len(Tree.children[0].particles) , 2)

	
	def test_barnes_hut_acc(self):
		P1 = Particle( 1 * u.pc , 1 * u.pc , 1 * u.pc , 1 * u.Msun , 1)
		P2 = Particle( 9 * u.pc , 9 * u.pc , 9 * u.pc , 1 * u.Msun , 2)
		P3 = Particle(3 * u.pc , 1 * u.pc , 3 * u.pc , 1 * u.Msun , 3)
		P4 = Particle( 8.5 * u.pc , 9 * u.pc , 9 * u.pc , 1 * u.Msun , 4)
		P5 = Particle( 8.5 * u.pc , 9.1 * u.pc , 9 * u.pc , 1 * u.Msun , 5)
		IC = Sim([P1 , P2 , P3])
		Tree = Find_Tree(IC)
		
		self.assertEqual(BH_Acceleration(IC ,Tree , P1) , dsum_acc(IC , P1))
		
	def test_bhacc(self):
	
		P1 = Particle( 1 * u.pc , 1 * u.pc , 1 * u.pc , 1 * u.Msun , 1)
		P2 = Particle( 9 * u.pc , 9 * u.pc , 9 * u.pc , 1 * u.Msun , 2)
		P3 = Particle(3 * u.pc , 1 * u.pc , 3 * u.pc , 1 * u.Msun , 3)
		P4 = Particle( 8.5 * u.pc , 9 * u.pc , 9 * u.pc , 1 * u.Msun , 4)
		P5 = Particle( 8.5 * u.pc , 9.1 * u.pc , 9 * u.pc , 1 * u.Msun , 5)
		IC = Sim([P1 , P2 , P3 , P4 , P5])
		##This set of particles has two particles close together, so the Barnes - Hut algorithm
		##Will not give the same answer as the direct sum
		Tree = Find_Tree(IC)
		
		self.assertNotEqual(BH_Acceleration(IC ,Tree , P1) , dsum_acc(IC  , P1))
	
unittest.main()
#timing(908)
