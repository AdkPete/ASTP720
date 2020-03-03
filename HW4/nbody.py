import numpy as np
import unittest
import astropy.constants as const
import astropy.units as u
import unittest


class Particle:
	
	def __init__(self , x , y , z , M):
		self.r =  [ [x , y , z] ]
		self.M = M
	
	def __eq__(self , other):
		if self.r[0] == other.r[0]:
			return True
		return False
		
	def accel(self , other , tstep):
		"""
		Calculates the acceleration on this particle due to another particle
		Takes in self and another particle object
		Does not include force softening
		returns an acceleration vector
		"""
		r = []
		
		for i in range(3):
			r.append(other.r[tstep][i] - self.r[tstep][i])
		
		###Magnitude of r12
		mag = np.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)
		
		###r12^ (unit vector)
		unit_r = []
		unit_r.append(r[0] / mag)
		unit_r.append(r[1] / mag)
		unit_r.append(r[2] / mag)
		a_mag = other.M * const.G / mag ** 2
		
		acc = []
		for i in range(len(unit_r)):
			acc.append(a_mag * unit_r[i])
		return acc
		
	def update_r(self , acc , h , t_step):
		##Updates r given an acceleration and step size and t_step
		nr = [ 0 , 0 , 0 ]
		for i in range(3):
			nr[i] = 2 * self.r[t_step][i] - self.r[t_step - 1][i] + h ** 2 * acc[i]
		self.r.append(nr)
		
def direct_summation(t_end , h , IC):
	'''
	This is an nbody solver that uses direct summation
	Not useful for our problem, but it is good for testing purposes
	takes in an end time and step size, t_end and h
	takes in a list of particles, IC
	returns a new list of particles
	'''
	
	t = 0 * t_end ##Preserves units
	t_step = 1
	while t < t_end:
		for P1 in IC:
			total_acc = [0 , 0 , 0]
			for P2 in IC:
				if P1 == P2:
					continue
				acc = P1.accel(P2 , t_step)
				for i in range(3):
					total_acc[i] += acc[i]
			P1.update_r(total_acc , h , t_step)
		t_step += 1
		t += h
	return IC
	
