import numpy as np
import unittest
import astropy.constants as const
import astropy.units as u
import unittest

class Particle:
	
	def __init__(self , x , y , z , vx , vy , vz , M):
		self.r = [x , y , z]
		self.v = [vx , vy , vz]
		self.M = M
		
	def accel(self , other):
		"""
		Calculates the acceleration on this particle due to another particle
		Takes in self and another particle object
		Does not include force softening
		returns an acceleration vector
		"""
		r = []
		
		for i in range(3):
			r.append(other.r[i] - self.r[i])
		
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
		

P1 = Particle( 0 * u.m , 0 * u.m , 0 * u.m , 5 * u.m / u.s, 5 * u.m , 5 * u.m , 1e12 * u.kg)
P2 = Particle( 1 * u.m , 1 * u.m , 0 * u.m , 5 * u.m / u.s, 5 * u.m , 5 * u.m , 1e12 * u.kg)

print (P1.accel(P2))		
	
