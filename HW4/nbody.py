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
		
	def calc_all_com(self , t): ###rewrite later. This goes in the wrong direction (top down instead of bottom up)
	
								###Less efficient than it should be
	
		TM = 0
		com = Vector( [ 0 , 0 , 0 ] )
		for i in self.particles:
			TM += i.M
			com += i.r[t] * i.M
		if TM == 0:
			self.com = com
			self.TM = 0
		else:
			self.com = com * (1 / TM)
			self.TM = TM
		if len(self.children) > 0:
			for i in self.children:
				i.calc_all_com(t)
				
				
		
	def reproduce(self , t):
	
		NL = self.L / 2.0
		LLL = Node(self.x , self.y , self.z , NL)
		LLU = Node(self.x , self.y , self.z + NL , NL)
		LUL = Node(self.x , self.y + NL , self.z , NL)
		ULL = Node(self.x + NL , self.y , self.z , NL)
		LUU = Node(self.x , self.y + NL , self.z + NL , NL)
		ULU = Node(self.x + NL , self.y , self.z + NL , NL)
		UUL = Node(self.x + NL , self.y + NL , self.z , NL)
		UUU = Node(self.x + NL , self.y + NL , self.z + NL , NL)
		###Now we need to place the existing particles into one of the child nodes.
		
		self.children = [LLL , LLU , LUL , ULL , LUU , ULU , UUL , UUU]
		for i in self.children:
			if i.in_node(self.particles[0] , t):
				i.add_particle(self.particles[0] , t)
		
		
	def add_particle(self , particle , t):
		if len(self.particles) >= 1:
			self.particles.append(particle)
			if len(self.children) == 0:
				self.reproduce(t)
			
			
			for i in self.children:
			
				if i.in_node(particle , t):
				
					i.add_particle(particle , t)
					break
					
		
		elif len(self.particles) == 0:
		
			self.particles.append(particle)
		
				

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
		
		
	def transform(self , unit):
	
		for i in range(len(self.elements)):
			self.elements[i] = self.elements[i].to(unit)
	
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
	'''
	Class used to store information for a single particle
	takes in the x , y , z positions. Defaults to units of pc, unless inputs include astropy units
	Also takes in a mass, assumed to be in solar masses
	Can optionally take in a particle id
	'''
	
	def __init__(self , x , y , z , M , id = None):
		ct = type(5 * u.m)
		warn = False ##Triggers a unit warning if set to True
		if type(x) != ct:
		
			x *= u.pc
			warn = True
			
		if type(y) != ct:
			y *= u.pc
			warn = True
			
			
		if type(z) != ct:
			z *= u.pc
			warn = True
		
		if type(M) != ct:
			M *= u.Msun
			warn = True
			
			
		if warn:
			print ("Warning: At least on quantity lacked units, so the default was assumed")
		self.r =  [ Vector([x , y , z])]
		self.M = M
		self.id = id
		
		
	def __eq__(self , other):
		'''
		defines equality of particles. If id numbers are the same, then particles are the same
		if either particle lacks an id, then we use the starting position
		'''
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
		
		"""
		Calculates the acceleration on this particle due to another particle
		Takes in self and another particle object
		Includes force softening
		returns an acceleration vector
		"""
		
		rel = other.r[tstep] - self.r[tstep]
		
		rel.mag()
		rel.unit()
		
		a_mag = other.M * const.G * S(rel , epsilon) / rel.m
		
		acc = rel.uv * a_mag
		
		return acc
		
	def update_r(self , acc , h , t_step):
		'''
		Takes in an acceleration vector acc, step size h, and time step t_step
		updates the position of the particle based on this acceleration
		'''
		
		
		##Updates r given an acceleration and step size and t_step
		nr = copy.deepcopy(self.r[t_step])
		nr *= 2
		nr -= self.r[t_step - 1]
		nr += nr + acc * h ** 2

		self.r.append(nr)
		
	def distance(self , other , t):
		#Complutes the distance between two particles at a time t
		r = self.r[t] - other.r[t]
		return r.mag()
		
	
		
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
	
	
def Find_Tree(IC , t):
	'''
	This takes in a list of particle objects, IC, and a time t.
	t should be an integer describing which time step we are interested in
	This function generates and returns a tree
	returns a node object
	'''
	x = []
	
	y = []
	
	z = []
	for i in IC:
	

		x.append(i.r[t][0])
		y.append(i.r[t][1])
		z.append(i.r[t][2])
	
	lunit = i.r[t][0].unit
	
	sx = min(x) - 1 * lunit
	sy = min(y) - 1 * lunit
	sz = min(z) - 1 * lunit
	
	L = max([ max(x) - sx + 1 * lunit , max(y) - sy+ 1 * lunit , max(z) - sz + 1 * lunit])
	Pnode = Node(sx , sy , sz , L)
	
	
	for i in IC:
	
		Pnode.add_particle( i , t )
	
	
	Pnode.calc_all_com(t)
	return Pnode
	
def BH_Acceleration(IC , t , Tree , Part):

	'''
	Calculates an acceleration using the Barnes - Hut tree algorithm
	Takes in a list of particles IC and time step t
	Takes in a tree (Node object) for the particle list
	Part is a particle object
	returns the acceleration vector for Part
	'''
	
	epsilon = 1 * u.pc
	theta_limit = 1
	Acc = Vector ( [ 0 , 0 , 0 ] )
	L = Tree.L
	D = Part.distance(Particle(Tree.com[0] , Tree.com[1] , Tree.com[2] , 0) , 0)
	T = L / D
	if len(Tree.children) == 0 and len(Tree.particles) == 1:
		if Tree.particles[0] == Part:
			return Acc
			
		else:
			
			Acc = Part.soft_acc(Tree.particles[0] , t , epsilon)
			
	elif T < theta_limit:
		Tree_part = Particle( Tree.com[0] , Tree.com[1] , Tree.com[2] , Tree.TM)
		Acc += Part.soft_acc(Tree_part , t , epsilon)
		
	elif len(Tree.children) > 0:
		for i in Tree.children:
			Acc += BH_Acceleration(IC , t , i , Part)
	return Acc
	

		
def Barnes_Hut(IC):


	'''
	Given a list of particles, will run an n-body simulation using the Barnes - Hut tree algorith
	IC is a list of particles
	will return the list of particles complete with updated positions
	
	'''
	t = 0
	Tree = Find_Tree(IC , t)
	theta = 1
	for i in Tree.particles:
		acc = BH_Acceleration(IC , t , Tree , i)
		acc.transform(u.m / (u.s ** 2))
		acc.print()
	return 0
	

