import numpy as np
import unittest
import astropy.constants as const
import astropy.units as u
import unittest
import copy




def S(r , eps):
	return 1 / np.sqrt(r * r + eps ** 2)
	

class Node:
	def __init__(self , x , y , z , slength):
		self.x = x
		self.y = y
		self.z = z
		self.L = slength
		self.children = []
		self.particles = []
	
	def in_node(self , particle , t):
		
		px = particle.r[t][0]
		py = particle.r[t][1]
		pz = particle.r[t][2]
		
		if px < self.x + self.L and px >= self.x:
			if py < self.y + self.L and py >= self.y:
				if pz < self.z + self.L and pz >= self.z:
					return True
					
		return False

class Vector:		

	def __init__(self , L):
		self.elements = L
		self.m = None
		
	def __getitem__(self , index):
		return self.elements[index]
		
	def __setitem__(self , index , val):
		self.elements[index] = val
	
	def __len__(self):
		return len(self.elements)
	
	def __add__(self , other):
	
		N = Vector([0] * len(self))
		for i in range(len(self)):
			N[i] = self[i] + other[i]
		return N
		
	def __sub__(self , other):
	
		N = Vector([0] * len(self))
		for i in range(len(self)):
			N[i] = self[i] - other[i]
			
		return N
		
	def __mul__(self , other):
		if type(other) == type(42) or type(other) == type(42.42) or type(other) == type(42.01 * u.pc):
			N = Vector( [0] * len(self))
			for i in range(len(N)):
				N[i] = self[i] * other
				
		elif type(other) == type(self):
			N = 0
			
			for i in range(len(self)):
				N += self[i] * other[i]
			
		return N
		
	def __eq__(self , other):
	
		if self.elements == other.elements:
			return True
		return False
	
	def print(self):
		print (self.elements)
		
		
	def unit(self):
		if self.m == None:
			self.mag()
		N = Vector([0] * len(self))
		for i in range(len(N)):
			N[i] = self[i] / self.m
		self.uv = N
		return N
		
		
	def mag(self):
		M = 0
		for i in self.elements:
			M += i ** 2
		self.m = np.sqrt(M)
		return self.m
		
	
class Particle:
	
	def __init__(self , x , y , z , M , id = None):
		self.r =  [ Vector([x , y , z])]
		self.M = M
		self.id = id
	def __eq__(self , other):
		if self.id == other.id and self.id != None and other.id != None:
			return True
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
		
		
		
		rel = other.r[tstep] - self.r[tstep]
		
		###Magnitude of r12
		rel.mag()
		
		###r12^ (unit vector)
		rel.unit()
		a_mag = other.M * const.G / rel.m ** 2
		
		acc = rel.uv * a_mag
		return acc
		
	def soft_acc(self , other , tstep , epsilon):
		
		
		rel = other.r[tstep] - self.r[tstep]
		
		rel.mag()
		
		rel.unit()
		
		a_mag = other.M * const.G * S(rel , epsilon) / rel.m
		
		acc = rel.uv * a_mag
		return acc
		
	def update_r(self , acc , h , t_step):
		##Updates r given an acceleration and step size and t_step
		nr = copy.deepcopy(self.r[t_step])
		nr *= 2
		nr -= self.r[t_step - 1]
		nr += nr + acc * h ** 2

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
			total_acc = Vector([0 , 0 , 0])
			for P2 in IC:
				if P1 == P2:
					continue
				acc = P1.accel(P2 , t_step)
				total_acc += acc
			P1.update_r(total_acc , h , t_step)
		t_step += 1
		t += h
	return IC
	
	
def Barnes_Hut(IC):
	###start by finding parent node
	x = []
	
	y = []
	
	z = []
	for i in IC:
	

		x.append(i.r[0][0])
		y.append(i.r[0][1])
		z.append(i.r[0][2])
	
	lunit = i.r[0][0].unit
	
	sx = min(x) - 1 * lunit
	sy = min(y) - 1 * lunit
	sz = min(z) - 1 * lunit
	
	L = max([ max(x) - sx + 1 * lunit , max(y) - sy+ 1 * lunit , max(z) - sz + 1 * lunit])
	Pnode = Node(sx , sy , sz , L)
	
	tp = Particle( 5 * u.pc , 8 * u.pc , 3 * u.pc , 1 * u.Msun)
	print (Pnode.in_node(tp , 0))
	
def test():
	P1 = Particle( 1 * u.pc , 1 * u.pc , 1 * u.pc , 1 * u.Msun)
	P2 = Particle( 9 * u.pc , 9 * u.pc , 9 * u.pc , 1 * u.Msun)
	
	IC = [P1 , P2]
	Barnes_Hut(IC)
test()
